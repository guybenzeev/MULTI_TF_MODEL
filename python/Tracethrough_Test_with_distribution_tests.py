from pathlib import Path
import re


import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.stats import chisquare, kstest
except ImportError:  # Keep the script usable even if scipy is not installed.
    chisquare = None
    kstest = None


REPO_ROOT = Path(__file__).resolve().parents[1]
LOG_PATH = REPO_ROOT / "data" / "multitf_cout.log"
HOLDING_TIME_OUTPUT_PATH = REPO_ROOT / "data" / "scaled_holding_times_hist.png"
THRESHOLD_OUTPUT_PATH = REPO_ROOT / "data" / "scaled_thresholds_hist.png"

EDIT_LINE_RE = re.compile(
    r"^Edit type: (?P<edit>.*), rate: (?P<rate>[-+0-9.eE]+), "
    r"cumulativeRate: (?P<cumulative>[-+0-9.eE]+), "
    r"Threshold: (?P<threshold>[-+0-9.eE]+)$"
)


def iter_trace_lines(log_path=LOG_PATH):
    with log_path.open("r", encoding="utf-8") as log_file:
        for line_number, line in enumerate(log_file, start=1):
            yield line_number, line.rstrip("\n")


def iter_steps(log_path=LOG_PATH):
    current_step = None
    current_lines = []

    for line_number, line in iter_trace_lines(log_path):
        if line.startswith("STEP "):
            if current_step is not None:
                yield current_step, current_lines

            _, step_number = line.split(maxsplit=1)
            current_step = int(step_number)
            current_lines = [(line_number, line)]
        elif current_step is not None:
            current_lines.append((line_number, line))

    if current_step is not None:
        yield current_step, current_lines


def parse_value_after_equals(line):
    return float(line.split("=", maxsplit=1)[1])


def parse_float_list(line):
    values = line.split("=", maxsplit=1)[1]
    return [float(value) for value in values.split(",") if value]


def parse_key_values(line):
    _, payload = line.split(maxsplit=1)
    result = {}

    for part in payload.split(","):
        key, value = part.strip().split("=", maxsplit=1)
        result[key] = value

    return result


def parse_step(step_number, lines):
    step = {
        "step": step_number,
        "edit_candidates": [],
    }

    for _, line in lines:
        if line.startswith("RATES total_exit_rate="):
            step["total_exit_rate"] = parse_value_after_equals(line)
        elif line.startswith("RATES subchain_exit_rates="):
            step["subchain_exit_rates"] = parse_float_list(line)
        elif line.startswith("TIME holding_time="):
            step["holding_time"] = parse_value_after_equals(line)
        elif line.startswith("THRESHOLD raw="):
            step["raw_threshold"] = parse_value_after_equals(line)
        elif line.startswith("THRESHOLD local_after_subchain_selection="):
            step["local_threshold"] = parse_value_after_equals(line)
        elif line.startswith("SELECTED subchain_index="):
            values = parse_key_values(line)
            step["selected_subchain_index"] = int(values["subchain_index"])
            step["selected_subchain_exit_rate"] = float(values["subchain_exit_rate"])
        elif line.startswith("EDIT selected="):
            selected, protein = line[len("EDIT selected="):].rsplit(", protein=", maxsplit=1)
            step["selected_edit"] = selected
            step["selected_edit_protein"] = int(protein)
        else:
            edit_match = EDIT_LINE_RE.match(line)
            if edit_match:
                step["edit_candidates"].append({
                    "edit": edit_match.group("edit"),
                    "rate": float(edit_match.group("rate")),
                    "cumulative": float(edit_match.group("cumulative")),
                    "threshold": float(edit_match.group("threshold")),
                })

    return step


def iter_parsed_steps(log_path=LOG_PATH):
    for step_number, lines in iter_steps(log_path):
        yield parse_step(step_number, lines)


def extract_scaled_holding_times(log_path=LOG_PATH):
    scaled_holding_times = []

    for step in iter_parsed_steps(log_path):
        total_exit_rate = step.get("total_exit_rate")
        holding_time = step.get("holding_time")

        if total_exit_rate is None or holding_time is None:
            continue

        scaled_holding_times.append((step["step"], holding_time * total_exit_rate))

    return scaled_holding_times


def extract_scaled_thresholds(log_path=LOG_PATH):
    scaled_thresholds = []

    for step in iter_parsed_steps(log_path):
        total_exit_rate = step.get("total_exit_rate")
        raw_threshold = step.get("raw_threshold")

        if total_exit_rate is None or raw_threshold is None:
            continue

        scaled_thresholds.append((step["step"], raw_threshold / total_exit_rate))

    return scaled_thresholds


def plot_scaled_holding_times(scaled_holding_times, save_path=HOLDING_TIME_OUTPUT_PATH):
    if not scaled_holding_times:
        raise ValueError("No complete steps with holding time and total exit rate were found.")

    values = np.array([value for _, value in scaled_holding_times])
    x_max = max(values.max(), 6.0)
    x = np.linspace(0.0, x_max, 500)
    exp_pdf = np.exp(-x)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(
        values,
        bins=50,
        density=True,
        alpha=0.65,
        edgecolor="black",
        label="holding_time * total_exit_rate",
    )
    ax.plot(x, exp_pdf, color="red", linewidth=2, label="Exp(lambda=1)")

    ax.set_title("Scaled Holding Times vs Exponential(1)")
    ax.set_xlabel("holding time * total exit rate")
    ax.set_ylabel("density")
    ax.legend()
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    if plt.get_backend().lower() != "agg":
        plt.show()


def plot_scaled_thresholds(scaled_thresholds, save_path=THRESHOLD_OUTPUT_PATH):
    if not scaled_thresholds:
        raise ValueError("No complete steps with raw threshold and total exit rate were found.")

    values = np.array([value for _, value in scaled_thresholds])
    x = np.linspace(0.0, 1.0, 2)
    uniform_pdf = np.ones_like(x)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(
        values,
        bins=50,
        range=(0.0, 1.0),
        density=True,
        alpha=0.65,
        edgecolor="black",
        label="raw_threshold / total_exit_rate",
    )
    ax.plot(x, uniform_pdf, color="red", linewidth=2, label="U(0,1)")

    ax.set_title("Scaled Thresholds vs Uniform(0,1)")
    ax.set_xlabel("raw threshold / total exit rate")
    ax.set_ylabel("density")
    ax.set_xlim(0.0, 1.0)
    ax.legend()
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    if plt.get_backend().lower() != "agg":
        plt.show()




def describe_samples(values):
    """Return basic summary statistics for a 1D numeric sample."""
    values = np.asarray(values, dtype=float)
    return {
        "n": int(values.size),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
        "variance": float(np.var(values, ddof=1)) if values.size > 1 else float("nan"),
    }


def print_sample_summary(name, values, expected_mean, expected_variance):
    summary = describe_samples(values)
    print(f"\n{name} summary:")
    print(f"  n:        {summary['n']}")
    print(f"  min:      {summary['min']:.10g}")
    print(f"  max:      {summary['max']:.10g}")
    print(f"  mean:     {summary['mean']:.10g}  expected about {expected_mean:.10g}")
    print(f"  variance: {summary['variance']:.10g}  expected about {expected_variance:.10g}")


def empirical_ks_statistic(values, cdf):
    """Compute the one-sample Kolmogorov-Smirnov D statistic."""
    values = np.sort(np.asarray(values, dtype=float))
    n = values.size
    if n == 0:
        raise ValueError("Cannot compute a KS statistic with no samples.")

    cdf_values = np.asarray([cdf(value) for value in values], dtype=float)
    empirical_upper = np.arange(1, n + 1) / n
    empirical_lower = np.arange(0, n) / n
    return float(np.max(np.maximum(empirical_upper - cdf_values, cdf_values - empirical_lower)))


def approximate_ks_pvalue(ks_statistic, sample_size):
    """
    Approximate the one-sample KS p-value.

    This is a fallback used only when scipy is unavailable. It is good enough for
    a diagnostic script, but scipy.stats.kstest is preferred when available.
    """
    if sample_size <= 0:
        return float("nan")

    adjusted = (np.sqrt(sample_size) + 0.12 + 0.11 / np.sqrt(sample_size)) * ks_statistic
    terms = [(-1) ** (k - 1) * np.exp(-2 * (k ** 2) * (adjusted ** 2)) for k in range(1, 101)]
    return float(np.clip(2 * np.sum(terms), 0.0, 1.0))


def run_ks_test(values, distribution_name):
    """
    Run a KS test against U(0,1) or Exp(1).

    Returns (statistic, pvalue, used_scipy).
    """
    values = np.asarray(values, dtype=float)

    if kstest is not None:
        statistic, pvalue = kstest(values, distribution_name)
        return float(statistic), float(pvalue), True

    if distribution_name == "uniform":
        statistic = empirical_ks_statistic(values, lambda x: min(max(x, 0.0), 1.0))
    elif distribution_name == "expon":
        statistic = empirical_ks_statistic(values, lambda x: 0.0 if x < 0.0 else 1.0 - np.exp(-x))
    else:
        raise ValueError(f"Unsupported distribution for fallback KS test: {distribution_name}")

    return statistic, approximate_ks_pvalue(statistic, values.size), False


def run_chi_square_uniform_test(values, num_bins=10):
    values = np.asarray(values, dtype=float)
    counts, _ = np.histogram(values, bins=num_bins, range=(0.0, 1.0))
    expected = np.full(num_bins, values.size / num_bins)

    if chisquare is not None:
        statistic, pvalue = chisquare(counts, expected)
        return counts, float(statistic), float(pvalue), True

    statistic = float(np.sum((counts - expected) ** 2 / expected))
    return counts, statistic, None, False


def run_chi_square_exponential_test(values, num_bins=10):
    """
    Chi-square test for Exp(1) using equal-probability bins.

    Equal-probability bins make the expected count the same in each bin. The last
    bin is open-ended, so very large exponential samples are handled correctly.
    """
    values = np.asarray(values, dtype=float)
    probabilities = np.linspace(0.0, 1.0, num_bins + 1)
    edges = -np.log1p(-probabilities[:-1])
    edges = np.append(edges, np.inf)
    counts, _ = np.histogram(values, bins=edges)
    expected = np.full(num_bins, values.size / num_bins)

    if chisquare is not None:
        statistic, pvalue = chisquare(counts, expected)
        return counts, edges, float(statistic), float(pvalue), True

    statistic = float(np.sum((counts - expected) ** 2 / expected))
    return counts, edges, statistic, None, False


def print_ks_result(label, statistic, pvalue, used_scipy, alpha=0.05):
    source = "scipy.stats.kstest" if used_scipy else "fallback KS approximation"
    print(f"  KS test ({source}):")
    print(f"    statistic: {statistic:.10g}")
    print(f"    p-value:   {pvalue:.10g}")
    if pvalue < alpha:
        print(f"    result:    FAIL to match expected distribution at alpha={alpha}")
    else:
        print(f"    result:    pass; no strong evidence against expected distribution at alpha={alpha}")


def print_chi_square_result(counts, statistic, pvalue, used_scipy, alpha=0.05):
    source = "scipy.stats.chisquare" if used_scipy else "manual statistic only; install scipy for p-value"
    print(f"  Chi-square bin test ({source}):")
    print(f"    bin counts: {counts.tolist()}")
    print(f"    statistic:  {statistic:.10g}")
    if pvalue is None:
        print("    p-value:    unavailable without scipy")
        return
    print(f"    p-value:    {pvalue:.10g}")
    if pvalue < alpha:
        print(f"    result:     FAIL to match expected bin counts at alpha={alpha}")
    else:
        print(f"    result:     pass; no strong evidence against expected bin counts at alpha={alpha}")


def test_scaled_holding_times(scaled_holding_times, alpha=0.05, num_bins=10):
    if not scaled_holding_times:
        raise ValueError("No complete steps with holding time and total exit rate were found.")

    values = np.array([value for _, value in scaled_holding_times], dtype=float)
    print_sample_summary("Scaled holding times: holding_time * total_exit_rate", values, 1.0, 1.0)

    if np.any(values < 0.0):
        print("  Bounds check: FAIL; exponential samples should never be negative.")
    else:
        print("  Bounds check: pass; all scaled holding times are nonnegative.")

    ks_statistic, ks_pvalue, used_scipy_ks = run_ks_test(values, "expon")
    print_ks_result("Exp(1)", ks_statistic, ks_pvalue, used_scipy_ks, alpha=alpha)

    counts, edges, chi_statistic, chi_pvalue, used_scipy_chi = run_chi_square_exponential_test(
        values,
        num_bins=num_bins,
    )
    print_chi_square_result(counts, chi_statistic, chi_pvalue, used_scipy_chi, alpha=alpha)


def test_scaled_thresholds(scaled_thresholds, alpha=0.05, num_bins=10):
    if not scaled_thresholds:
        raise ValueError("No complete steps with raw threshold and total exit rate were found.")

    values = np.array([value for _, value in scaled_thresholds], dtype=float)
    print_sample_summary("Scaled thresholds: raw_threshold / total_exit_rate", values, 0.5, 1.0 / 12.0)

    if np.any(values < 0.0) or np.any(values > 1.0):
        bad_count = int(np.sum((values < 0.0) | (values > 1.0)))
        print(f"  Bounds check: FAIL; {bad_count} scaled thresholds are outside [0, 1].")
    else:
        print("  Bounds check: pass; all scaled thresholds are in [0, 1].")

    ks_statistic, ks_pvalue, used_scipy_ks = run_ks_test(values, "uniform")
    print_ks_result("U(0,1)", ks_statistic, ks_pvalue, used_scipy_ks, alpha=alpha)

    counts, chi_statistic, chi_pvalue, used_scipy_chi = run_chi_square_uniform_test(
        values,
        num_bins=num_bins,
    )
    print_chi_square_result(counts, chi_statistic, chi_pvalue, used_scipy_chi, alpha=alpha)


def select_subchain_from_threshold(raw_threshold, exit_rates):
    cumulative = 0.0

    for index, rate in enumerate(exit_rates):
        previous = cumulative
        cumulative += rate
        if raw_threshold < cumulative:
            return index, raw_threshold - previous

    return None, None


def select_edit_from_threshold(local_threshold, edit_candidates):
    for candidate in edit_candidates:
        if candidate["cumulative"] >= local_threshold:
            return candidate["edit"]

    return None


def close_enough(left, right, tolerance=1e-3):
    return abs(left - right) <= tolerance


def validate_trace(log_path=LOG_PATH):
    issues = []
    checked_subchains = 0
    checked_edits = 0

    for step in iter_parsed_steps(log_path):
        required_subchain_fields = {
            "total_exit_rate",
            "subchain_exit_rates",
            "raw_threshold",
            "local_threshold",
            "selected_subchain_index",
            "selected_subchain_exit_rate",
        }

        if required_subchain_fields.issubset(step):
            expected_index, expected_local_threshold = select_subchain_from_threshold(
                step["raw_threshold"],
                step["subchain_exit_rates"],
            )
            checked_subchains += 1

            if expected_index != step["selected_subchain_index"]:
                issues.append(
                    f"STEP {step['step']}: selected subchain {step['selected_subchain_index']} "
                    f"but expected {expected_index}"
                )

            if expected_local_threshold is not None and not close_enough(
                expected_local_threshold,
                step["local_threshold"],
            ):
                issues.append(
                    f"STEP {step['step']}: local threshold {step['local_threshold']} "
                    f"but expected {expected_local_threshold}"
                )

            selected_rate = step["subchain_exit_rates"][step["selected_subchain_index"]]
            if not close_enough(selected_rate, step["selected_subchain_exit_rate"]):
                issues.append(
                    f"STEP {step['step']}: selected subchain rate "
                    f"{step['selected_subchain_exit_rate']} but expected {selected_rate}"
                )

        if "local_threshold" in step and "selected_edit" in step and step["edit_candidates"]:
            expected_edit = select_edit_from_threshold(
                step["local_threshold"],
                step["edit_candidates"],
            )
            checked_edits += 1

            if expected_edit != step["selected_edit"]:
                issues.append(
                    f"STEP {step['step']}: selected edit {step['selected_edit']} "
                    f"but expected {expected_edit}"
                )

    return issues, checked_subchains, checked_edits


def main():
    if not LOG_PATH.exists():
        raise FileNotFoundError(f"Trace log not found: {LOG_PATH}")

    scaled_holding_times = extract_scaled_holding_times()
    scaled_thresholds = extract_scaled_thresholds()
    issues, checked_subchains, checked_edits = validate_trace()

    print(f"Parsed {len(scaled_holding_times)} complete holding-time samples.")
    print(f"Parsed {len(scaled_thresholds)} complete threshold samples.")
    print(f"Saving holding-time histogram to: {HOLDING_TIME_OUTPUT_PATH}")
    plot_scaled_holding_times(scaled_holding_times)
    print(f"Saving threshold histogram to: {THRESHOLD_OUTPUT_PATH}")
    plot_scaled_thresholds(scaled_thresholds)

    test_scaled_holding_times(scaled_holding_times)
    test_scaled_thresholds(scaled_thresholds)

    print(f"\nChecked {checked_subchains} subchain selections.")
    print(f"Checked {checked_edits} edit selections.")
    if issues:
        print("Trace validation issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("Trace validation passed.")


if __name__ == "__main__":
    main()
