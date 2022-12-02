#We load the project
using DrWatson
@quickactivate "TumorSim"
#Import dependencies
using TumorSim


#We create the fitness, scenarios and therapies we want to use
fitness=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.2)

scenario = create_scenario((100,100,100),10,"center")
treatment = create_treatment(3000, 2000, 1000, 3, 0.75)

parameters = Dict(
    "pr" => 0.027,
    "dr" => 0.015,
    "mr" => 0.01,   
    "scenario" => scenario, 
    "fitness" => fitness,
    "treatment" => treatment,
    "seed" => 0
)

launch_interactive_simulation(parameters)