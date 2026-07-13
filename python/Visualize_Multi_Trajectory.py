from __future__ import annotations

import os
from collections.abc import Mapping

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle
import numpy as np


class Visualize_Multi_Trajectory:
    """Visualize trajectories that use the C++ multi-protein state IDs."""

    def __init__(self, trajectory):
        self.trajectory = trajectory

    def _num_states(self):
        """Return FREE plus the NS/S states for every protein and side."""
        return (
            1
            + 2
            * self.trajectory.get_num_proteins()
            * self.trajectory.get_num_sides()
        )

    def _state_label(self, state_id):
        if state_id == 0:
            return "Free"

        num_sides = self.trajectory.get_num_sides()
        num_proteins = self.trajectory.get_num_proteins()
        states_per_binding_type = num_sides * num_proteins

        state_offset = state_id - 1
        if state_offset < states_per_binding_type:
            binding_type = "NS"
        else:
            binding_type = "S"
            state_offset -= states_per_binding_type

        protein = state_offset // num_sides + 1
        side = state_offset % num_sides + 1
        return f"Protein {protein}, {binding_type}, side {side}"

    def _default_palette(self):
        """Assign a distinct color to every state in the trajectory schema."""
        num_states = self._num_states()
        palette = {0: "white"}

        # Sample a continuous map so the palette also works when there are more
        # protein/side combinations than the categorical maps contain.
        color_map = plt.get_cmap("turbo")
        for state_id in range(1, num_states):
            fraction = (state_id - 1) / max(num_states - 2, 1)
            palette[state_id] = color_map(fraction)

        return palette

    def _complete_palette(self, palette):
        if palette is None:
            return self._default_palette()

        if not isinstance(palette, Mapping):
            raise TypeError("palette must be a mapping from state IDs to colors")

        # A partial custom palette remains supported. Unspecified valid states
        # are gray instead of being omitted from the pcolormesh normalization.
        return {
            state_id: palette.get(state_id, "gray")
            for state_id in range(self._num_states())
        }

    def _add_state_legend(self, ax, palette):
        handles = [
            Patch(
                facecolor=palette[state_id],
                edgecolor="black",
                linewidth=0.3,
                label=self._state_label(state_id),
            )
            for state_id in range(self._num_states())
        ]
        ax.legend(
            handles=handles,
            title="State",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0,
        )

    @staticmethod
    def _event_time(line):
        if len(line) < 5:
            raise ValueError(
                "Expected CSV row with five columns: "
                "time,index,edit,protein,side"
            )
        return float(line[0])

    def visualize_blockwise(self, palette=None, show_legend=True):
        length = self.trajectory.get_length()
        run_time = self.trajectory.get_runtime()
        palette = self._complete_palette(palette)
        state = [0] * length

        fig, ax = plt.subplots(figsize=(12, 6))
        current_time = 0.0
        self.trajectory.init_iterator()

        while self.trajectory.has_next():
            line = self.trajectory.next_line()
            time = self._event_time(line)
            holding_time = time - current_time

            # Draw the state occupied during [current_time, time), then apply
            # the event. applyEdit decodes its protein and side fields.
            for index, state_id in enumerate(state):
                ax.add_patch(
                    Rectangle(
                        (current_time, index),
                        holding_time,
                        1,
                        facecolor=palette[state_id],
                        edgecolor="black",
                        linewidth=0.2,
                    )
                )

            self.trajectory.applyEdit(line, state)
            current_time = time

        if current_time < run_time:
            for index, state_id in enumerate(state):
                ax.add_patch(
                    Rectangle(
                        (current_time, index),
                        run_time - current_time,
                        1,
                        facecolor=palette[state_id],
                        edgecolor="black",
                        linewidth=0.2,
                    )
                )

        ax.set_xlim(0, run_time)
        ax.set_ylim(length, 0)
        ax.set_xlabel("Time")
        ax.set_ylabel("Index")
        ax.set_title("Blockwise Multi-Protein Trajectory")
        if show_legend:
            self._add_state_legend(ax, palette)

        fig.tight_layout()
        plt.show()
        return fig, ax

    def visualize_blockwise_fast(
        self,
        palette=None,
        filename="trajectory.png",
        show_legend=True,
    ):
        length = self.trajectory.get_length()
        run_time = self.trajectory.get_runtime()
        palette = self._complete_palette(palette)

        self.trajectory.init_iterator()
        lines = []
        while self.trajectory.has_next():
            lines.append(self.trajectory.next_line())

        matrix = np.zeros((length, len(lines) + 1), dtype=int)
        state = [0] * length
        times = [0.0]

        for interval_index, line in enumerate(lines):
            matrix[:, interval_index] = state
            self.trajectory.applyEdit(line, state)
            times.append(self._event_time(line))

        matrix[:, -1] = state
        times.append(run_time)

        num_states = self._num_states()
        colors = [palette[state_id] for state_id in range(num_states)]
        color_map = ListedColormap(colors)
        norm = BoundaryNorm(np.arange(num_states + 1) - 0.5, color_map.N)

        fig, ax = plt.subplots(figsize=(12, 6))
        x_edges = np.asarray(times)
        y_edges = np.arange(length + 1)
        ax.pcolormesh(
            x_edges,
            y_edges,
            matrix,
            cmap=color_map,
            norm=norm,
            shading="flat",
            edgecolors="none",
        )

        for time in x_edges:
            ax.axvline(time, color="black", linewidth=0.2, alpha=0.3)

        ax.set_xlim(0, run_time)
        ax.set_ylim(length, 0)
        ax.set_xlabel("Time")
        ax.set_ylabel("Index")
        ax.set_title("Blockwise Multi-Protein Trajectory")
        if show_legend:
            self._add_state_legend(ax, palette)

        fig.tight_layout()
        os.makedirs("figures", exist_ok=True)
        save_path = os.path.join("figures", filename)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved to {save_path}")
        plt.show()
        return fig, ax

    def holding_time_visualization(
        self,
        stateIndex,
        bins=20,
        figsize=(5, 4),
        rate=None,
        log_scale=False,
        filename="holding_time.png",
    ):
        """Plot visits with exactly one protein in the requested state ID."""
        if stateIndex < 0 or stateIndex >= self._num_states():
            raise ValueError(
                f"stateIndex must be between 0 and {self._num_states() - 1}"
            )

        holding_times = []
        length = self.trajectory.get_length()
        run_time = self.trajectory.get_runtime()
        state = [0] * length
        current_time = 0.0

        self.trajectory.init_iterator()
        while self.trajectory.has_next():
            line = self.trajectory.next_line()
            time = self._event_time(line)
            occupied_indices = np.flatnonzero(state)

            if stateIndex == 0 and occupied_indices.size == 0:
                holding_times.append(time - current_time)
            elif occupied_indices.size == 1:
                protein_index = int(occupied_indices[0])
                if (
                    state[protein_index] == stateIndex
                    and protein_index not in {0, length - 1}
                ):
                    holding_times.append(time - current_time)

            self.trajectory.applyEdit(line, state)
            current_time = time

        # Account for the final state, which persists until the declared runtime.
        if current_time < run_time:
            occupied_indices = np.flatnonzero(state)
            if stateIndex == 0 and occupied_indices.size == 0:
                holding_times.append(run_time - current_time)
            elif occupied_indices.size == 1:
                protein_index = int(occupied_indices[0])
                if (
                    state[protein_index] == stateIndex
                    and protein_index not in {0, length - 1}
                ):
                    holding_times.append(run_time - current_time)

        label = self._state_label(stateIndex)
        print(f"Average holding time for {label}: {self._average(holding_times)}")
        self._plot_holding_times(
            holding_times,
            title=f"Holding Time Distribution for {label}",
            bins=bins,
            figsize=figsize,
            rate=rate,
            log_scale=log_scale,
            filename=filename or f"holding_time_state_{stateIndex}.png",
        )
        return holding_times

    def holding_time_for_exact_state(
        self,
        target_state,
        bins=20,
        figsize=(5, 4),
        rate=None,
        log_scale=False,
        filename="holding_time_exact_state.png",
    ):
        """Plot holding times for an exact list of multi-protein state IDs."""
        length = self.trajectory.get_length()
        run_time = self.trajectory.get_runtime()
        target_state = list(target_state)

        if len(target_state) != length:
            raise ValueError("target_state must match trajectory length")

        invalid_ids = [
            state_id
            for state_id in target_state
            if state_id < 0 or state_id >= self._num_states()
        ]
        if invalid_ids:
            raise ValueError(
                "target_state contains invalid state IDs: "
                f"{sorted(set(invalid_ids))}"
            )

        holding_times = []
        state = [0] * length
        current_time = 0.0
        self.trajectory.init_iterator()

        while self.trajectory.has_next():
            line = self.trajectory.next_line()
            time = self._event_time(line)
            if state == target_state:
                holding_times.append(time - current_time)

            self.trajectory.applyEdit(line, state)
            current_time = time

        if current_time < run_time and state == target_state:
            holding_times.append(run_time - current_time)

        print(f"Target state: {target_state}")
        print(f"Average holding time: {self._average(holding_times)}")
        self._plot_holding_times(
            holding_times,
            title=f"Holding Time Distribution for State {target_state}",
            bins=bins,
            figsize=figsize,
            rate=rate,
            log_scale=log_scale,
            filename=filename,
        )
        return holding_times

    @staticmethod
    def _average(values):
        return float(np.mean(values)) if values else 0.0

    @staticmethod
    def _plot_holding_times(
        holding_times,
        title,
        bins,
        figsize,
        rate,
        log_scale,
        filename,
    ):
        fig, ax = plt.subplots(figsize=figsize)
        if holding_times:
            ax.hist(
                holding_times,
                bins=bins,
                density=True,
                alpha=0.6,
                edgecolor="black",
                label="Simulated holding times",
            )

        if log_scale:
            ax.set_yscale("log")

        if rate is not None and holding_times:
            x = np.linspace(0, max(holding_times), 300)
            ax.plot(
                x,
                rate * np.exp(-rate * x),
                linewidth=2,
                label=f"Exp(λ={rate})",
            )

        ax.set_xlabel("Holding time")
        ax.set_ylabel("Density")
        ax.set_title(title)
        if holding_times:
            ax.legend()

        fig.tight_layout()
        os.makedirs("figures", exist_ok=True)
        save_path = os.path.join("figures", filename)
        fig.savefig(save_path, dpi=300)
        print(f"Figure saved to {save_path}")
        plt.show()

    def test(self):
        """Print and apply the first ten multi-protein events."""
        length = self.trajectory.get_length()
        state = [0] * length
        current_time = 0.0
        self.trajectory.init_iterator()

        for _ in range(10):
            if not self.trajectory.has_next():
                break
            line = self.trajectory.next_line()
            time = self._event_time(line)
            print(f"Processing line: {line}")
            print(f"Current time: {current_time}, State: {state}")
            self.trajectory.applyEdit(line, state)
            current_time = time
