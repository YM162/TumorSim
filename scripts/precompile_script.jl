using DrWatson
@quickactivate "TumorSim"

using TumorSim
using Distributed
#We run a simple simulation to preload the functions needed.
fitness=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.2)

scenario = create_scenario((100,100,100),10,"center")
treatment = create_treatment(3000, 2000, 1000, 3, 0.75)

params = Dict(
    "pr" => 0.027,
    "dr" => 0.5,
    "mr" => [0.01,0.1],   
    "scenario" => scenario, 
    "fitness" => fitness,
    "treatment" => treatment,
    "seed" => 0
)

dicts = dict_list(params)

steps=100
results = pmap(simulate,dicts,fill(steps,length(dicts)))