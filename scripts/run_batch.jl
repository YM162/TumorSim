#a = run(Cmd(`julia -q --sysimage build/TumorSim.so -p auto scripts/run_batch.jl`,detach=true))

using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using JLD2
using DataFrames
using Dates
#With the same fitness we can change the cost of resistance
fitness1=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>0.6)

fitness2=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>0.8)

fitness3=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1)

fitness4=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.5)

fitness5=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.5)
                
#We test the dimensionality.
scenario_0D = create_scenario((1000000),10)
scenario_1D = create_scenario((1000000,),10,"center")
scenario_2D = create_scenario((1000,1000),10,"center")
scenario_3D = create_scenario((100,100,100),10,"center")

#We test adaptive and continuous therapy
adaptive_therapy = create_treatment(3000, 0.65, 0.5, 3, 0.75) 
continuous_therapy = create_treatment(3000, 0.65, 0.0, 3, 0.75) 

parameters = Dict(
    "pr" => [0.01,0.015,0.02,0.025,0.03,0.04,0.1],
    "dr" => [0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7],
    "mr" => [0.01,0.1,0.2],  
    "scenario" => [scenario_3D], 
    "fitness" => [fitness4],
    "cr" => [0.0,0.1,0.3,0.5,0.7],
    "treatment" => [adaptive_therapy,continuous_therapy],
    "seed" => map(abs,rand(Int64,10))
)

parameter_combinations = dict_list(parameters)
println("Number of simulations: ",length(parameter_combinations))
steps=3000

filename = "simulations_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")
open("logs/progress/"*filename*".log", "w") do io
    p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
    #Hay que cambiar Progress para que solamente de el numero y el tiempo, para sacarlo facilmente del log
    results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

    println("Saving simulations...")
    df = DataFrame(results)
    filepath = datadir("simulations",filename*".jld2")

    jldopen(filepath, "w") do file
        file["df"] = df
    end;
    println(df[!,"TTP"])
end



