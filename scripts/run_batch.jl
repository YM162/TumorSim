using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
#We create the fitness, scenarios and therapies we want to use
fitness=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.2)

scenario_0D = create_scenario((1000000),10)
scenario_1D = create_scenario((1000000,),10,"center")
scenario_2D = create_scenario((1000,1000),10,"center")
scenario_3D = create_scenario((100,100,100),10,"center")

adaptive_therapy = create_treatment(3000, 2000, 1000, 3, 0.75) 
continuous_therapy = create_treatment(3000, 2000, 0, 3, 0.75) 

#This would be cool to do, but we need the cluster
parameters = Dict(
    "pr" => [0.005,0.01,0.0015,0.02,0.025],
    "dr" => [0,0.25,0.5,0.75],
    "mr" => [0.001,0.005,0.01,0.015,0.03,0.05],   
    "scenario" => scenario_3D, 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy,continuous_therapy],
    "seed" => map(abs,rand(Int64,100))
)
#we can do this instead
parameters = Dict(
    "pr" => 0.027,
    "dr" => 0.55,
    "mr" => 0.01,   
    "scenario" => scenario_3D, 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy,continuous_therapy],
    "seed" => map(abs,rand(Int64,100))
)

parameter_combinations = dict_list(parameters)

steps=3000

results = @showprogress pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)))

for (i, d) in enumerate(parameter_combinations)
    safesave(datadir("simulations", savename(d, "jld2")), results[i])
end