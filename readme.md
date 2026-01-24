# High-Performance Trajectory Design and Optimization in Julia

**Course:** Computer Science I (Programming)  
**Program:** PhD in Statistics and Computer Science (Bocconi University)

---

## Overview

This project implements a **high-performance scientific computing framework in Julia** for interplanetary trajectory design. The primary goal is to write **efficient, allocation-aware, and scalable Julia code**, in line with the objectives of the course.

The framework supports:

- Earth–Mars transfer optimization via porkchop plots
- Multi-leg missions with gravity assists (patched conics)
- Voyager-2–like trajectories as a compact benchmark

The emphasis is not on mission-grade astrodynamics, but on **explicit numerical implementations, performance engineering, and clean software structure**.

---

## Modeling Assumptions

- Spacecraft motion is modeled using the **two-body approximation** with the Sun as the central body.
- Interplanetary transfers are computed by solving **Lambert boundary-value problems**.
- Gravity assists are modeled using a **patched-conic approximation**.
- Planetary ephemerides are obtained from **NAIF SPICE kernels (DE440)** via `SPICE.jl`.

Importantly, SPICE is used **only to precompute planetary states**. All performance-critical loops operate exclusively on precomputed data.

---

## Code Structure

The project is organized as a modular Julia codebase:

- `PorkchopSolver.jl` – main module definition and public API
- `types.jl` – core data structures and shared type definitions
- `stumpff.jl`, `kepler.jl` – low-level numerical primitives
- `lambert.jl` – Lambert solver implementation
- `ephemeris.jl` – ephemeris backends (SPICE and analytic)
- `porkchop.jl` – high-performance porkchop solvers
- `flyby.jl`, `gravity_assist.jl` – patched-conic gravity-assist modeling
- `viz.jl`, `viz_traj.jl` – trajectory and solution visualization utilities
- `zoom_analysis.jl` – local refinement and high-resolution analysis tools
- `demo*.jl` – reproducible numerical experiments and benchmarks
- `make_plots.jl` – automated plot generation
- `run_project.jl` – top-level driver script that executes all experiments end-to-end

Numerical kernels, data handling, and visualization are clearly separated to enable focused optimization.

## Performance-Oriented Design

The implementation explicitly targets performance:

- **StaticArrays** (`SVector{3,Float64}`) are used for all 3D vectors to avoid heap allocations.
- All large arrays are **preallocated** and reused via dedicated workspace structures.
- Core numerical kernels are **allocation-free**; remaining allocations originate from setup, I/O, or plotting.
- The outer grid loops are parallelized using `Threads.@threads`, achieving ~80–85% CPU utilization.

---

## Numerical Experiments and Results

### Earth–Mars Porkchop (2026 Opportunity)

- Time resolution: **3 hours**
- Grid size: **5840 × 8296** ≈ **48.4 million** departure/arrival combinations
- Lambert solves executed: **≈ 41.1 million** (short- and long-way)

**Measured performance (Julia, multi-threaded):**

- Short-way Lambert scan: **6.03 s**
- Long-way Lambert scan: **9.47 s**
- Total allocations: **≈ 56k**
- Peak allocated memory: **≈ 3.0 MB**

This corresponds to **millions of Lambert solves per second** with negligible memory pressure, confirming a strongly compute-bound implementation.

**Best trajectory found:**

- Departure: **2026-10-31 06:00 UTC**
- Arrival: **2027-09-07 09:00 UTC**
- Time of flight: **311.12 days** (7467 hours)
- Total Δv: **5.59 km/s**

Representative porkchop, Δv, C3, and $v_\infty$ plots are automatically generated and stored in `plots/`.

---

### Flyby Mission (Earth–Jupiter–Saturn, 2020 Window)

- Scan resolution: **5 days**
- Runtime: **2.59 s**
- Allocations: **≈ 2.8k**
- Allocated memory: **≈ 2.7 MB**

Despite involving conditional logic, geometry checks, and patched-conic transitions, the solver remains highly efficient and allocation-light.

Garbage collection remains below **1%**, indicating a compute-bound execution even at higher resolutions.

A high-resolution zoom of the optimal region is also generated for post-analysis.

---

### Voyager-2–Like Trajectory (1977 Benchmark)

This scenario reproduces a simplified version of the historical Voyager-2 mission using Lambert arcs and patched-conic gravity assists.

**Measured performance:**

- Runtime: **0.04 s**
- Allocations: **726**
- Memory usage: **160 kB**

**Computed optimal trajectory:**

- Earth departure: **1977-09-03**
- Jupiter flyby: **1979-07-29**
- Saturn arrival: **1981-09-30**

**Mission duration:**

- Earth → Jupiter: **694 days**
- Jupiter → Saturn: **794 days**
- Total duration: **4.07 years**

**Δv budget (km/s):**

- Launch Δv (C3 proxy): **6.76**
- Flyby Δv: **0.0**
- Total mission Δv: **6.76**

---

## Comparison with the Real Voyager 2 Mission

| Quantity                | Real Voyager 2  | This Project    |
| ----------------------- | --------------- | --------------- |
| Launch date             | 1977-08-20      | 1977-09-03      |
| Jupiter flyby           | 1979-07-09      | 1979-07-29      |
| Saturn flyby            | 1981-08-26      | 1981-09-30      |
| Earth → Saturn duration | ~4.0 years      | **4.07 years**  |
| Launch C3               | ~6.5–7.0 km²/s² | **6.76 km²/s²** |
| Flyby modeling          | Full N-body     | Patched conics  |

Despite using a **two-body Sun-centered model** and **patched-conic gravity assists**, the computed trajectory:

- matches mission dates within **weeks**
- matches total mission duration within **~2%**
- matches launch energy within historical uncertainty

This validates both the **numerical correctness** of the solvers and the **physical plausibility** of the simplified modeling assumptions, while keeping the computational cost extremely low.

---

## GPU Acceleration (Exploratory / Future Work)

The computational structure of the porkchop solvers, consisting of large, embarrassingly parallel grids of independent Lambert problems, makes the framework naturally amenable to GPU acceleration.

An initial GPU-oriented refactor was explored, targeting a CUDA-based implementation using `CUDA.jl`, with the goal of accelerating wider and higher-resolution grids beyond what is practical on a CPU alone.

However, GPU execution could not be fully validated due to hardware and environment constraints:

- No local CUDA-capable GPU was available during development.
- Attempts to use Google Colab GPUs were unsuccessful, as Julia–CUDA integration could not be reliably established in that environment.

As a result, all reported benchmarks correspond to optimized multi-threaded CPU execution.
Nevertheless, the solver architecture (pure kernels, minimal branching, allocation-free inner loops) was designed with GPU portability in mind, and a GPU backend could be integrated with limited structural changes given appropriate hardware access.

---

## How to Run

Clone the repository, then:

```bash
julia --project
pkg> instantiate
```

To run all experiments and generate plots:

```bash
julia -t auto --project=. run_project.jl
```

This command:

- runs all demos (`demo.jl`, `demo_flyby.jl`, `demo_voyager2*.jl`)
- saves numerical results to `.jld2` files
- generates plots in the `plots/` directory

---

## Reproducibility Notes

The project includes both `Project.toml` and `Manifest.toml`. It has been tested by running:

```julia
deleteat!(LOAD_PATH, 2)
```

to ensure no reliance on globally installed packages.

NAIF kernels and `.jld2` output files are not tracked in the repository and are downloaded/generated automatically when needed.

---

This project was developed as a final project for **Computer Science I (Programming)** and is primarily intended to demonstrate **efficient Julia programming for scientific computing**.
