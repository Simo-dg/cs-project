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

- `stumpff.jl`, `kepler.jl` – low-level numerical primitives
- `lambert.jl` – Lambert solver implementation
- `ephemeris.jl` – ephemeris backends (SPICE and analytic)
- `porkchop.jl` – high-performance porkchop solvers
- `flyby.jl`, `gravity_assist.jl` – gravity-assist modeling
- `demo*.jl` – reproducible numerical experiments
- `make_plots.jl` – automated plot generation

Numerical kernels, data handling, and visualization are clearly separated to enable focused optimization.

---

## Performance-Oriented Design

The implementation explicitly targets performance:

- **StaticArrays** (`SVector{3,Float64}`) are used for all 3D vectors to avoid heap allocations.
- All large arrays are **preallocated** and reused via dedicated workspace structures.
- Core numerical kernels are **allocation-free**; remaining allocations originate from setup, I/O, or plotting.
- The outer grid loops are parallelized using `Threads.@threads`, achieving ~80–85% CPU utilization.

---

## Numerical Experiments and Results

### Earth–Mars Porkchop

- Time resolution: **3 hours**
- Grid size: ~48.4 million departure/arrival combinations
- Lambert solves: ~41.1 million (short- and long-way)

**Best solution found:**

- Departure: 2026-10-31 06:00 UTC
- Arrival: 2027-09-07 09:00 UTC
- Time of flight: 311.12 days
- Total Δv: **5.59 km/s**

Representative plots are stored in `plots/`.

---

### Flyby Mission (Earth–Jupiter–Saturn)

| Resolution | Runtime | Allocations |   Memory |
| ---------- | ------: | ----------: | -------: |
| 5 days     |    ~3 s |       ~2.7k |  ~2.7 MB |
| 2 day      |   ~42 s |       ~6.7k | ~15.6 MB |

Garbage collection remains below 1% even at the highest resolution, indicating a compute-bound implementation.

---

### Voyager-2 Scenario

- Runtime: ~0.04 s
- Allocations: ~700
- Memory: ~160 kB
- Best mission Δv: 6.76 km/s

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

- runs all demos (`demo.jl`, `demo_flyby.jl`, `demo_voyager*.jl`)
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
