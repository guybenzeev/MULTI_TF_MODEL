# Multi-TF Model

Multi-TF Model simulates multiple protein types binding to, unbinding from, and
sliding along a DNA strand. The C++ simulator uses a Gillespie-style
continuous-time Markov chain. Local rate recalculation and a Fenwick tree keep
event selection efficient as the strand grows.

The repository also contains Python utilities for reconstructing trajectories
from the simulator's CSV output and creating state and holding-time plots.

## Requirements

- A C++17 compiler
- Doxygen (optional, for API documentation)
- Python 3 with NumPy and Matplotlib (optional, for visualization)

## Build

From the repository root:

```sh
mkdir -p build
c++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
  -Iinclude src/*.cpp -o build/multitf-model
```

## Run

```sh
./build/multitf-model [length] [run_time] [output.csv]
```

For example:

```sh
./build/multitf-model 10 1000 data/events.csv
```

Arguments are optional and positional:

1. `length` is the number of sites in the strand.
2. `run_time` is the simulated time limit.
3. `output.csv` is the trajectory output path.

The CSV header records the strand length, number of sides, number of protein
types, and requested run time. Each remaining row contains the event time,
zero-based head-site index, edit type, one-based protein ID, and one-based
strand side.

The executable currently constructs a single-sided topology with `ProteinA`
(width 3) and `ProteinB` (width 5). Add another `Protein` subclass and pass an
instance to `SingleStrandTopo` to define another kinetic model.

## API documentation

The public headers in `include/` use Doxygen conventions. Generate the HTML
reference with:

```sh
cd docs
doxygen Doxyfile
```

Open `docs/html/index.html` after generation.

`open docs/html/index.html`

## Python visualization

Install the optional dependencies:

```sh
python3 -m pip install numpy matplotlib
```

`python/Multi_Trajectory.py` parses simulator output, and
`python/Visualize_Multi_Trajectory.py` provides plotting helpers. Update the
input filename and desired visualization in `python/main.py`, then run:

```sh
cd python
python3 main.py
```

Generated figures are written under `python/figures/` by default.
