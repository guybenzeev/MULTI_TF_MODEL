from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


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

    print(f"Checked {checked_subchains} subchain selections.")
    print(f"Checked {checked_edits} edit selections.")
    if issues:
        print("Trace validation issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("Trace validation passed.")


if __name__ == "__main__":
    main()
