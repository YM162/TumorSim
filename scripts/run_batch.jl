#a = run(Cmd(`julia -q --sysimage build/TumorSim.so -p auto scripts/run_batch.jl`,detach=true))

using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using BSON
using DataFrames
using Dates
#With the same fitness we can change the cost of resistance


#We test adaptive and continuous therapy
fitness=Dict([0,0,0]=>0.027, 
            [1,0,0]=>0.03,
            [0,1,0]=>0.033,
            [1,1,0]=>0.035,
            [1,1,1]=>0.035)

scenario = create_scenario((100,100),10)

adaptive_therapy = create_treatment(3000, 0.8, 0.5, 3, 0.75) 
continuous_therapy = create_treatment(3000, 0.8, 0.0, 3, 0.75) 

parameters = Dict(
    "dr" => 0.0,
    "mr" => 0.01,
    "scenario" => scenario, 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy,continuous_therapy],
    "cr" => 0.3,
    "seed" => map(abs,rand(Int64,10))
)

parameter_combinations = dict_list(parameters)
println("Number of simulations: ",length(parameter_combinations))
steps=3000

filename = "simulations_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")

p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

println("Saving simulations...")
df = DataFrame(results)
filepath = datadir("simulations",filename*".bson")

bson(filepath,Dict("df" => df))

println(df[!,"TTP"])



