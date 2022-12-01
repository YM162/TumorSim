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

scenario = create_scenario((1,1,1),10,"center")
treatment = Treatment(3000, 2000, 0, 3, 0.75, false, false)

params = Dict(
    "pr" => 0.027,
    "dr" => 0.015,
    "mr" => [0.01,0.01],   
    "scenario" => scenario, 
    "fitness" => fitness,
    "treatment" => treatment,
    "seed" => 0
)

dicts = dict_list(params)

steps=10
results = pmap(simulate,dicts,fill(steps,length(dicts)))