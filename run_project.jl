using Dates

function run_script(path)
    file_name = basename(path)
    println("\n" * "="^60)
    println("RUNNING: $file_name")
    println("="^60)
    
    try
        run(`julia -t auto --project=. $path`)
        println("SUCCESS: $file_name completed.")
    catch e
        println("FAILED: $file_name encountered an error.")
        println(e)
    end
end

println("Starting Project Execution at $(now())\n")

# 1. Solvers
run_script(joinpath("src", "demo.jl"))           
run_script(joinpath("src", "demo_flyby.jl"))     
run_script(joinpath("src", "demo_voyager2.jl"))  

# 2. Plotter 
run_script(joinpath("src", "make_plots.jl"))

println("\n" * "="^60)
println("PROJECT COMPLETE.")
println("="^60)