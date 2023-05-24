using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using BSON
using DataFrames
using Dates
using CSV

#We test NK fitness landscapes with 7 genes and 3 interactions per gene under the NK model

#R Generation of landscape using NK model
#rnk <- rfitness(6, K = 3, model = "NK", scale = c(0,0.055,0.027))
#write.csv2(rnk,"C://Users/yomis/TFG/TumorSim/scripts/AT_in_cancer_the_role_of_restrictions_2023/test_NK.csv")


csv_reader = CSV.File(projectdir("scripts","AT_in_cancer_the_role_of_restrictions_2023","NK_Fitness_2.csv"))
fitness = Dict([[x[i] for i in 2:length(x)-1]=>parse(Float64,replace(x[end],","=>".")) for x in csv_reader if parse(Float64,replace(x[end],","=>"."))>0.0001])
#cost_of_resistance = 0.5
#fitness = Dict(vcat([vcat(x[1],0)=>x[2] for x in fitness],[vcat(x[1],1)=>x[2]*cost_of_resistance for x in fitness]))



adaptive_therapy = create_treatment(5000, 1, 0.5, 7, 0.75) 
continuous_therapy = create_treatment(5000, 1, 0.0, 7, 0.75) 

parameters = Dict(
    "death_rate" => [0.3],
    "mutation_rate" => 0.01,
    "scenario" => [create_scenario((100,100),100,"center",false)], 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy],
    "migration_rate" => [0.1],
    "interaction_rule" => [:contact],
    "seed" => map(abs,rand(Int64,1000))
)

parameter_combinations = dict_list(parameters)
println("Number of simulations: ",length(parameter_combinations))
steps=5000

filename = "7Genes_NK_Fitness_2_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")

p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

println("Saving simulations...")
df = DataFrame(results)
filepath = datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023",filename*".bson")


#bson(filepath,Dict("df" => df))
using StatsBase

newdf=DataFrame(step=Int64[])
for i in df[!,"Divergence"]
    global newdf = outerjoin(newdf,i[:,1:2],on=:step,makeunique=true)
end
sort!(newdf,:step)

enddf = DataFrame(step=Int64[],jensen_shannon_mean=Float64[],jensen_shannon_sd=Float64[])
for i in eachrow(newdf)
    push!(enddf,Dict(:step=>i[1],:jensen_shannon_mean=>mean(skipmissing(i[2:end])),:jensen_shannon_sd=>std(skipmissing(i[2:end]))))
end

#Remove every row where any of the columns has an undefined value
filter(row -> all(x -> !(x isa Number && isnan(x)), row), enddf)
finaldf = filter(row -> all(x -> !(x isa Number && isnan(x)), row), enddf)

#bson(datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023","cleanup",filename),Dict("divergence" => finaldf, "divergence_raw" => Matrix(newdf), "TTP" => df[!,"TTP"], "detecting_time" => [sim[!,"step"][findfirst(sim[!,"status"])] for sim in df[!,"Treatment_status"]]))
println("Calculating stats of f0 for",filename)
println("Mean: ",mean(skipmissing(df[!,"Resistant_on_detection"])))
println("Standard deviation: ",std(skipmissing(df[!,"Resistant_on_detection"])))