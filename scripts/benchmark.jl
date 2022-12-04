using DrWatson
@quickactivate "TumorSim"
using TumorSim


using BenchmarkTools
fitness4=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.2)

scenario_3D = create_scenario((100,100,100),10,"center")
adaptive_therapy = create_treatment(3000, 2000, 1000, 3, 0.75) 


@benchmark simulate(parameters,3000) setup=(parameters = Dict("pr" => 0.027, "dr" => 0.55, "mr" => 0.01, "scenario" => scenario_3D, "fitness" => fitness4, "treatment" => adaptive_therapy, "seed" => abs(rand(Int64)))) seconds=1