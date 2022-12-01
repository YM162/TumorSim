using DrWatson
@quickactivate "TumorSim"

using TumorSim
using Distributed
using ProgressMeter
# Here you may include files from the source directory
fitness=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.2)



scenario_0D = create_scenario((1000000),10)
scenario_1D = create_scenario((1000000,),10,"center")
scenario_2D = create_scenario((1000,1000),10,"center")
scenario_3D = create_scenario((100,100,100),10,"center")

adaptive_therapy = Treatment(3000, 2000, 1000, 3, 0.75, false, false) 
continuous_therapy = Treatment(3000, 2000, 0, 3, 0.75, false, false) 

scenario = create_scenario((1,1,1),10,"center")
treatment = Treatment(3000, 2000, 0, 3, 0.75, false, false)

parameters = Dict(
    "pr" => 0.027,
    "dr" => 0.015,
    "mr" => [0.01,0.01],   
    "scenario" => scenario, 
    "fitness" => fitness,
    "treatment" => treatment,
    "seed" => 0
)

parameter_combinations = dict_list(parameters)

steps=10
results = @showprogress pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)))

for (i, d) in enumerate(parameter_combinations)
    safesave(datadir("simulations", savename(d, "jld2")), results[i])
end