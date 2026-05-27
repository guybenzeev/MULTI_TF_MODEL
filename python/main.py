from Multi_Trajectory import Multi_Trajectory
from Visualize_Multi_Trajectory import Visualize_Multi_Trajectory


def main() -> None:

    events_trajectory = Multi_Trajectory("events.csv")
    events_visualizer = Visualize_Multi_Trajectory(events_trajectory)
    events_visualizer.visualize_blockwise_fast(filename="events_blockwise.png")

    josh_trajectory = Multi_Trajectory("josh's_events.csv")
    josh_visualizer = Visualize_Multi_Trajectory(josh_trajectory)
    josh_visualizer.visualize_blockwise_fast(filename="josh_blockwise.png")

    reduced_trajectory = Multi_Trajectory("reduced_events.csv")
    reduced_visualizer = Visualize_Multi_Trajectory(reduced_trajectory)
    reduced_visualizer.holding_time_visualization(2, bins=20, rate = 0.054, log_scale=False, filename="reduced_holding_time_SB.png")

    #reduced_visualizer.holding_time_visualization(0, bins=20, rate = 0.005, log_scale=True)
    #reduced_visualizer.holding_time_visualization(1, bins=20, rate = 4.604, log_scale=True)

    mini_trajectory = Multi_Trajectory("mini_events.csv")
    mini_visualizer = Visualize_Multi_Trajectory(mini_trajectory)
    mini_visualizer.holding_time_for_exact_state([0, 1, 1], bins=20, rate=4.2, log_scale=False)



if __name__ == "__main__":
    main()
