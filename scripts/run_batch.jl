#a = run(Cmd(`julia -q --sysimage build/TumorSim.so -p auto scripts/run_batch.jl`,detach=true))

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
                [3,4,1]]

ngenes = 3
base_pr = 0.027
#Base pr multiplicative, except for the last one, which is the resistant gene.
mult_pr = [1.16,1.35]
cost_of_resistance = 0.15

fitness = build_fitness_table(restrictions,base_pr,mult_pr,cost_of_resistance,ngenes)
println(fitness)

adaptive_therapy = create_treatment(3000, 1, 0.5, 3, 0.75) 
continuous_therapy = create_treatment(3000, 1, 0.0, 3, 0.75) 

parameters = Dict(
    "death_rate" => [0.3],
    "mutation_rate" => 0.01,
    "scenario" => [create_scenario((100,100),100,"center",false)], 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy,continuous_therapy],
    "migration_rate" => [0.1],
    "interaction_rule" => [:contact],
    "seed" => map(abs,rand(Int64,1))
)

parameter_combinations = dict_list(parameters)
println("Number of simulations: ",length(parameter_combinations))
steps=5000

filename = "simulations_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")

p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

println("Saving simulations...")
df = DataFrame(results)
filepath = datadir("simulations",filename*".bson")

bson(filepath,Dict("df" => df))

println(df[!,"TTP"])
println(df[!,"Resistant_on_detection"])


n=2
using VegaLite
#stack(df[!,"Genotypes"][n],names(df[!,"Genotypes"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")
#stack(df[!,"Genotypes"][n],names(df[!,"Genotypes"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:zero}, color="variable:n")

#stack(df[!,"Resistant_inhibited_by"][n],names(df[!,"Resistant_inhibited_by"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")

#df[!,"Divergence"][n] |> @vlplot(:line, x=:step, y=:jenshen_shannon, color=:variable)













#new_inhibited = DataFrame()
#fr = min(length(df[!,"Resistant_inhibited_by"][n][!,"step"]),length(df[!,"Resistant_inhibited_by"][n-1][!,"step"]))
#for i in names(df[!,"Resistant_inhibited_by"][n])[2:end]
#    resta = df[!,"Resistant_inhibited_by"][n-1][!,i][1:fr]-df[!,"Resistant_inhibited_by"][n][!,i][1:fr]
#    insertcols!(new_inhibited,length(names(new_inhibited))+1,i => resta)
#end
#insertcols!(new_inhibited,1,"step" => df[!,"Resistant_inhibited_by"][n][1:fr,"step"])
#stack(new_inhibited,names(df[!,"Resistant_inhibited_by"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:zero}, color="variable:n")
