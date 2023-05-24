using DrWatson
@quickactivate "TumorSim"
using TumorSim


using BenchmarkTools
using ProfileView

fitness=Dict([0,0,0]=>0.027, 
            [1,0,0]=>0.030,
            [0,1,0]=>0.033,
            [1,1,0]=>0.036,
            [1,1,1]=>0.032)

adaptive_therapy = create_treatment(3000, 1, 0.5, 3, 0.75) 

parameters = Dict(
    "death_rate" => [0.0],
    "mutation_rate" => 0.007,
    "scenario" => [create_scenario((100,100),1,"center",false)], 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy],
    "cost_of_resistance" => [0.0],
    "migration_rate" => [0.3],
    "interaction_rule" => [:hierarchical_voter],
    "seed" => 0
)
parameter_combinations = dict_list(parameters)
println("Starting Benchmark:")

#We test the performance of the model.
@benchmark simulate(parameter_combinations[1],3000) seconds=60

#@profview simulate(parameter_combinations[1],3000)
#1.6s