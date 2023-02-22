using DrWatson
@quickactivate "TumorSim"
using TumorSim


using BenchmarkTools
fitness=Dict([0,0,0]=>0.027, 
            [1,0,0]=>0.035,
            [0,1,0]=>0.032,
            [1,1,0]=>0.040,
            [1,1,1]=>0.032)

scenario_3D = create_scenario((100,100,100),10,"center")
adaptive_therapy = create_treatment(3000, 0.65, 0.5, 3, 0.75) 

println("Starting Benchmark:")
@benchmark simulate(parameters,3000) setup=(parameters = Dict("cost_of_resistance" => 0.2, "dr" => 0.55, "mutation_rate" => 0.01, "scenario" => scenario_3D, "fitness" => fitness4, "treatment" => adaptive_therapy, "seed" => abs(rand(Int64)))) seconds=120

