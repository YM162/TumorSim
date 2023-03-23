using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using BSON
using DataFrames
using Dates

#We test adaptive and continuous therapy

restrictions= [[1,2,1],
                [1,3,1],
                [2,4,1],
                [2,5,2],
                [3,5,2],
                [3,6,1],
                [4,8,1],
                [5,8,1],
                [6,7,1]]

ngenes = 7
base_pr = 0.027
#Base pr multiplicative, except for the last one, which is the resistant gene.
mult_pr = [1.06,1.11,1.10,1.04,1.09,1.07]
cost_of_resistance = 0.3

fitness = build_fitness_table(restrictions,base_pr,mult_pr,cost_of_resistance,ngenes)

adaptive_therapy = create_treatment(5000, 1, 0.5, 7, 0.75) 
continuous_therapy = create_treatment(5000, 1, 0.0, 7, 0.75) 

parameters = Dict(
    "death_rate" => [0.3],
    "mutation_rate" => 0.02,
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

filename = "7Genes_Example2_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")

p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

println("Saving simulations...")
df = DataFrame(results)
filepath = datadir("simulations","competition_divergence",filename*".bson")

bson(filepath,Dict("df" => df))


println(df[!,"TTP"])
println(df[!,"Resistant_on_detection"])
stack(df[!,"Genotypes"][n],names(df[!,"Genotypes"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:zero}, color="variable:n")

