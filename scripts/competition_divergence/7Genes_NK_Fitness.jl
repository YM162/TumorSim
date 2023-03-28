using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using BSON
using DataFrames
using Dates
using CSV

#We test adaptive and continuous therapy

#R Generation of landscape using NK model
#rnk <- rfitness(6, K = 3, model = "NK", scale = c(0,0.055,0.027))
#write.csv2(rnk,"C://Users/yomis/TFG/TumorSim/scripts/competition_divergence/test_NK.csv")


csv_reader = CSV.File(projectdir("scripts","competition_divergence","test_NK.csv"))
fitness = Dict([[x[i] for i in 2:length(x)-1]=>parse(Float64,replace(x[end],","=>".")) for x in csv_reader])
cost_of_resistance = 0.5
fitness = Dict(vcat([vcat(x[1],0)=>x[2] for x in fitness],[vcat(x[1],1)=>x[2]*cost_of_resistance for x in fitness]))



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

filename = "7Genes_NO_Restrictions_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")

p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

println("Saving simulations...")
df = DataFrame(results)
filepath = datadir("simulations","competition_divergence",filename*".bson")

bson(filepath,Dict("df" => df))
