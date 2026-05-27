import os

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import os

class Visualize_Multi_Trajectory:
    def __init__(self, trajectory):
        self.trajectory = trajectory

    def visualize_blockwise(self, palette=None):
        length = self.trajectory.get_length()
        num_sides = self.trajectory.get_num_sides()
        run_time = self.trajectory.get_runtime()

        if palette is None:
            palette = {
                0: "white",
                1: "skyblue",
                2: "orange",
                3: "green",
            }

        state = [0] * length

        fig, ax = plt.subplots(figsize=(12, 6))

        current_time = 0.0
        self.trajectory.init_iterator()

        while self.trajectory.has_next():
            line = self.trajectory.next_line()
            time, index, edit_type, side = line

            time = float(time)
            holding_time = time - current_time

            # Draw current state during interval [current_time, time]
            for i in range(length):
                block_state = state[i]
                color = palette.get(block_state, "gray")

                rect = Rectangle(
                    (current_time, i),      # bottom-left corner
                    holding_time,           # width in time units
                    1,                      # height of one block
                    facecolor=color,
                    edgecolor="black",
                    linewidth=0.2,
                )
                ax.add_patch(rect)

            # Apply edit AFTER drawing interval
            self.trajectory.applyEdit(line, state)

            current_time = time

        # Draw final state from last event to run_time
        if current_time < run_time:
            holding_time = run_time - current_time

            for i in range(length):
                block_state = state[i]
                color = palette.get(block_state, "gray")

                rect = Rectangle(
                    (current_time, i),
                    holding_time,
                    1,
                    facecolor=color,
                    edgecolor="black",
                    linewidth=0.2,
                )
                ax.add_patch(rect)

        ax.set_xlim(0, run_time)
        ax.set_ylim(0, length)

        ax.set_xlabel("Time")
        ax.set_ylabel("Index")
        ax.set_title("Blockwise Trajectory Visualization")

        ax.invert_yaxis()  # optional: index 0 at top
        plt.tight_layout()
        plt.show()

    def visualize_blockwise_fast(self, palette=None, filename="trajectory.png"):
        length = self.trajectory.get_length()
        run_time = self.trajectory.get_runtime()

        if palette is None:
            palette = {
                0: "white",
                1: "skyblue",
                2: "orange",
                3: "green",
            }

        self.trajectory.init_iterator()
        lines = []

        while self.trajectory.has_next():
            lines.append(self.trajectory.next_line())

        num_intervals = len(lines) + 1
        matrix = np.zeros((length, num_intervals), dtype=int)

        state = [0] * length
        times = [0.0]

        self.trajectory.init_iterator()

        j = 0
        while self.trajectory.has_next():
            line = self.trajectory.next_line()
            time, index, edit_type, side = line
            time = float(time)

            matrix[:, j] = state
            self.trajectory.applyEdit(line, state)

            times.append(time)
            j += 1

        matrix[:, j] = state
        times.append(run_time)

        y_edges = np.arange(length + 1)
        x_edges = np.array(times)

        max_state = max(palette.keys())
        colors = [palette.get(i, "gray") for i in range(max_state + 1)]

        cmap = ListedColormap(colors)
        norm = BoundaryNorm(np.arange(max_state + 2) - 0.5, cmap.N)

        fig, ax = plt.subplots(figsize=(12, 6))

        ax.pcolormesh(
            x_edges,
            y_edges,
            matrix,
            cmap=cmap,
            norm=norm,
            shading="flat",
            edgecolors="none",
        )

        for t in x_edges:
            ax.axvline(t, color='black', linewidth=0.2, alpha=0.3)

        ax.set_xlim(0, run_time)
        ax.set_ylim(length, 0)

        ax.set_xlabel("Time")
        ax.set_ylabel("Index")
        ax.set_title("Blockwise Trajectory Visualization")

        plt.tight_layout()

        os.makedirs("figures", exist_ok=True)  # create folder if it doesn't exist
        save_path = os.path.join("figures", filename)
        plt.savefig(save_path, dpi=300)

        print(f"Figure saved to {save_path}")

        plt.show()

    def holding_time_visualization(self, stateIndex, bins=20, figsize=(5, 4), rate=None, log_scale=False, filename="holding_time.png"):
        holding_times = []
        length = self.trajectory.get_length()
        run_time = self.trajectory.get_runtime()

        state = [0] * length
        current_time = 0.0
        current_state = 0

        self.trajectory.init_iterator()

        while self.trajectory.has_next():
            count = np.count_nonzero(state)
            #print(f"num_tfs={count}, state={state}")
            line = self.trajectory.next_line()
            time, index, edit_type, side = line
            #print(f"edit: time={time}, index={index}, type={edit_type}, side={side}")

            time = float(time)
            print(f"time at check 1: {time}")
            holding_time = time - current_time
            
            if(count == 0):
                current_state = 0
            elif(count == 1):
                current_state = state[np.nonzero(state)[0][0]]
                tf_index = np.nonzero(state)[0][0]
            elif(count > 1):
                current_state = -1

            #print(f"current_state: {current_state}")
            if count == 0 and stateIndex == 0:
                print(f"adding to holding time for state {current_state} with num_tfs={count}")
                holding_times.append(holding_time)
            elif count == 1 and stateIndex == current_state:
                if state[0] != 1 and state[length-1] != 1:
                    #if holding_time > 0.7:
                        #print(
                        #"BIG STATE 1 HT:",
                      #  "holding_time=", holding_time,
                      #  "state before=", state.copy(),
                      #  "tf_index=", tf_index,
                      #  "next line=", line
                     #   )
                    holding_times.append(holding_time)
            
            current_time = time
            self.trajectory.applyEdit(line, state)

        avg = np.mean(holding_times) if holding_times else 0.0
        print(f"Average holding time: {avg}")

        # Plot histogram
        fig, ax = plt.subplots(figsize=figsize)

        ax.hist(
            holding_times,
            bins=bins,
            density=True,
            alpha=0.6,
            edgecolor="black",
            label="Simulated holding times"
        )
        if log_scale:
            ax.set_yscale('log')

        if rate is not None and len(holding_times) > 0:
            x = np.linspace(0, max(holding_times), 300)
            y = rate * np.exp(-rate * x)
            ax.plot(x, y, linewidth=2, label=f"Exp(λ={rate})")

        ax.set_xlabel("Holding time")
        ax.set_ylabel("Density")
        ax.set_title(f"Holding Time Distribution for State {stateIndex}")
        ax.legend()
        plt.tight_layout()

        os.makedirs("figures", exist_ok=True)

        if filename is None:
            filename = f"holding_time_state_{stateIndex}.png"

        save_path = os.path.join("figures", filename)
        plt.savefig(save_path, dpi=300)

        print(f"Figure saved to {save_path}")

    
        plt.show()

        return holding_times

    def holding_time_for_exact_state(
        self,
        target_state,
        bins=20,
        figsize=(5, 4),
        rate=None,
        log_scale=False,
        filename="holding_time_exact_state.png"
    ):
        """
        target_state: list of length L (e.g., [0,1,0,2,0])
        """

        holding_times = []
        length = self.trajectory.get_length()
        run_time = self.trajectory.get_runtime()

        if len(target_state) != length:
            raise ValueError("target_state must match trajectory length")

        state = [0] * length
        print(f"Initial state: {state}")
        current_time = 0.0

        self.trajectory.init_iterator()

        while self.trajectory.has_next():
            line = self.trajectory.next_line()
            #print(line)
            time, index, edit_type, side = line
            time = float(time)

            holding_time = time - current_time

            if state == target_state:
                #print(f"Adding holding time {holding_time} for state {state} at time {time}")
                holding_times.append(holding_time)

            # advance state
            self.trajectory.applyEdit(line, state)
            current_time = time

        avg = np.mean(holding_times) if holding_times else 0.0
        print(f"Target state: {target_state}")
        print(f"Average holding time: {avg}")

        # Plot
        fig, ax = plt.subplots(figsize=figsize)

        ax.hist(
            holding_times,
            bins=bins,
            density=True,
            alpha=0.6,
            edgecolor="black",
            label="Simulated holding times"
        )

        if log_scale:
            ax.set_yscale('log')

        if rate is not None and len(holding_times) > 0:
            x = np.linspace(0, max(holding_times), 300)
            y = rate * np.exp(-rate * x)
            ax.plot(x, y, linewidth=2, label=f"Exp(λ={rate})")

        ax.set_xlabel("Holding time")
        ax.set_ylabel("Density")
        ax.set_title(f"Holding Time Distribution for State {target_state}")
        ax.legend()
        plt.tight_layout()
        os.makedirs("figures", exist_ok=True)

        save_path = os.path.join("figures", filename)
        plt.savefig(save_path, dpi=300)

        print(f"Figure saved to {save_path}")

        plt.show()

        return holding_times
    


    def test(self):
        length = self.trajectory.get_length()
        num_sides = self.trajectory.get_num_sides()
        run_time = self.trajectory.get_runtime()
        state = [0] * length

        fig, ax = plt.subplots(figsize=(12, 6))

        current_time = 0.0

        for i in range(10):  # Just process the first 10 lines for testing
            line = self.trajectory.next_line()
            time, index, edit_type, side = line
            print(f"Processing line: {line}")
            print(f"Current time: {current_time}, State: {state}")

            time = float(time)
            holding_time = time - current_time
            self.trajectory.applyEdit(line, state)

            current_time = time