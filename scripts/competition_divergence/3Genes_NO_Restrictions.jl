using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using BSON
using DataFrames
using Dates

#We test adaptive and continuous therapy

restrictions= []

ngenes = 3
base_pr = 0.027
#Base pr multiplicative, except for the last one, which is the resistant gene.
mult_pr = [1.16,1.35]
cost_of_resistance = 0.15

fitness = build_fitness_table(restrictions,base_pr,mult_pr,cost_of_resistance,ngenes)

adaptive_therapy = create_treatment(3000, 1, 0.5, 3, 0.75) 
continuous_therapy = create_treatment(3000, 1, 0.0, 3, 0.75) 

parameters = Dict(
    "death_rate" => [0.3],
    "mutation_rate" => 0.0033,
    "scenario" => [create_scenario((100,100),100,"center",false)], 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy],
    "migration_rate" => [0.1],
    "interaction_rule" => [:contact],
    "seed" => map(abs,rand(Int64,100))
)

parameter_combinations = dict_list(parameters)
println("Number of simulations: ",length(parameter_combinations))
steps=5000

filename = "3Genes_NO_Restrictions_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")

p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

println("Saving simulations...")
df = DataFrame(results)
filepath = datadir("simulations","competition_divergence",filename*".bson")

bson(filepath,Dict("df" => df))
