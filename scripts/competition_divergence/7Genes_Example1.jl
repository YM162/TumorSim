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
                [2,3,1],
                [2,4,1],
                [2,5,1],
                [3,6,2],
                [4,6,2],
                [4,7,1],
                [5,7,1],
                [6,8,2],
                [7,8,2]]

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
    "seed" => map(abs,rand(Int64,10000))
)

parameter_combinations = dict_list(parameters)
println("Number of simulations: ",length(parameter_combinations))
steps=5000

filename = "7Genes_Example1_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")

p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

println("Saving simulations...")
df = DataFrame(results)
filepath = datadir("simulations","competition_divergence",filename*".bson")

for i in 1:length(eachrow(results))
    
end

#bson(filepath,Dict("df" => df))
using StatsBase

newdf=DataFrame(step=Int64[])
for i in df[!,"Divergence"]
    global newdf = outerjoin(newdf,i[:,1:2],on=:step,makeunique=true)
end
sort!(newdf,:step)

enddf = DataFrame(step=Int64[],jenshen_shannon_mean=Float64[],jenshen_shannon_sd=Float64[])
for i in eachrow(newdf)
    push!(enddf,Dict(:step=>i[1],:jenshen_shannon_mean=>mean(skipmissing(i[2:end])),:jenshen_shannon_sd=>std(skipmissing(i[2:end]))))
end

#Remove every row where any of the columns has an undefined value
filter(row -> all(x -> !(x isa Number && isnan(x)), row), enddf)
finaldf = filter(row -> all(x -> !(x isa Number && isnan(x)), row), enddf)

bson(datadir("simulations","competition_divergence","cleanup",filename),Dict("divergence" => finaldf, "divergence_raw" => Matrix(newdf), "TTP" => df[!,"TTP"], "detecting_time" => [sim[!,"step"][findfirst(sim[!,"status"])] for sim in df[!,"Treatment_status"]]))

