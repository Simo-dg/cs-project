# run_project.jl
using Dates

function run_script(path)
    file_name = basename(path)
    println("\n" * "="^60)
    println("🚀 RUNNING: $file_name")
    println("="^60)
    
    try
        # ADDED "-t auto" HERE 👇
        run(`julia -t auto --project=. $path`)
        println("✅ SUCCESS: $file_name completed.")
    catch e
        println("❌ FAILED: $file_name encountered an error.")
        println(e)
    end
end

println("Starting Project Execution at $(now())\n")

# 1. Physics Solvers (Heavy Computation)
run_script(joinpath("src", "demo.jl"))           # Earth -> Mars
run_script(joinpath("src", "demo_flyby.jl"))     # Flyby 2020
run_script(joinpath("src", "demo_voyager2.jl"))  # Voyager 2

# 2. Plotter (Visualization)
# Plotting doesn't strictly need threads, but it doesn't hurt.
run_script(joinpath("src", "make_plots.jl"))

println("\n" * "="^60)
println("🎉 PROJECT COMPLETE.")
println("="^60)