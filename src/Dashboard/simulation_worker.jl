#julia -p 8 src/Dashboard/simulation_worker.jl 0.027 0.001 0.027 0.5 0.05 0.5 0.01 0.01 0.01 1000000 3 10 3000 100 3000 0.65 0.1 0.65 0.5 0.1 0.5 0.75 0.2 0.1 0.2 1

using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using BSON
using DataFrames
using Dates

#This is not ideal, we should use JSON for this in the future, making it not order-dependent and allowing for better compatibility.

pr_low = parse(Float64,ARGS[1])
pr_step = parse(Float64,ARGS[2])
pr_high = parse(Float64,ARGS[3])

dr_low = parse(Float64,ARGS[4])
dr_step = parse(Float64,ARGS[5])
dr_high = parse(Float64,ARGS[6])

mr_low = parse(Float64,ARGS[7])
mr_step = parse(Float64,ARGS[8])
mr_high = parse(Float64,ARGS[9])

s_size = parse(Int64,ARGS[10])
s_dim = parse(Int64,ARGS[11])
s_initial_cells = parse(Int64,ARGS[12])

t_detecting_size_low = parse(Int64,ARGS[13])
t_detecting_size_step = parse(Int64,ARGS[14]) #Absolute value
t_detecting_size_high = parse(Int64,ARGS[15])

t_starting_size_low = parse(Float64,ARGS[16])
t_starting_size_step = parse(Float64,ARGS[17]) #Relative value (%) of detecting size
t_starting_size_high = parse(Float64,ARGS[18])

t_pausing_size_low = parse(Float64,ARGS[19])
t_pausing_size_step = parse(Float64,ARGS[20]) #Relative value (%) of starting size
t_pausing_size_high = parse(Float64,ARGS[21])

t_kill_rate = parse(Float64,ARGS[22])

cr_low = parse(Float64,ARGS[23]) #% de penalización por tener la mutación de resistencia
cr_step = parse(Float64,ARGS[24]) 
cr_high = parse(Float64,ARGS[25]) 

repetitions = parse(Int64,ARGS[26])


#Hay que pensar en como pasar el fitness landscape


if s_dim == 0
    scenario=create_scenario(s_size,s_initial_cells)
elseif s_dim == 1
    scenario=create_scenario((s_size,),s_initial_cells)
elseif s_dim == 2
    scenario=create_scenario((round(Int64,s_size^(1/2)),round(Int64,s_size^(1/2))),s_initial_cells)
elseif s_dim == 3
    scenario=create_scenario((round(Int64,s_size^(1/3)),round(Int64,s_size^(1/3)),round(Int64,s_size^(1/3))),s_initial_cells)
end

fitness=Dict([0,0,0]=>1, 
            [1,0,0]=>1.3,
            [0,1,0]=>1.2,
            [1,1,0]=>1.5,
            [1,1,1]=>1.5)

adaptive_therapy = [create_treatment(t_detecting_size, t_starting_size, t_pausing_size, 3, t_kill_rate) for t_detecting_size in t_detecting_size_low:t_detecting_size_step:t_detecting_size_high
                                                                            for t_starting_size in t_starting_size_low:t_starting_size_step:t_starting_size_high
                                                                            for t_pausing_size in t_pausing_size_low:t_pausing_size_step:t_pausing_size_high]

continuous_therapy = [create_treatment(t_detecting_size, t_starting_size, 0.0, 3, t_kill_rate) for t_detecting_size in t_detecting_size_low:t_detecting_size_step:t_detecting_size_high
                                                                            for t_starting_size in t_starting_size_low:t_starting_size_step:t_starting_size_high]

parameters = Dict(
    "pr" => collect(pr_low:pr_step:pr_high),
    "dr" => collect(dr_low:dr_step:dr_high),
    "mr" => collect(mr_low:mr_step:mr_high), 
    "scenario" => scenario,
    "fitness" => fitness,
    "cr" => collect(cr_low:cr_step:cr_high),
    "treatment" => append!(adaptive_therapy,continuous_therapy),
    "seed" => map(abs,rand(Int64,repetitions))
)

parameter_combinations = dict_list(parameters)

steps=3000
println("Starting simulations...")

filename = "simulations_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")

println("Number of simulations: ",length(parameter_combinations))

open(projectdir("logs","progress",filename*".log"), "w") do io
    p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=io,desc="",barlen=0)
    results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

    println("Saving simulations...")

    df = DataFrame(results)

    filepath = datadir("simulations",filename*".bson")


    bson(filepath,Dict("df" => df))
    
    println(df[!,"TTP"])
    #Meter un mensaje al log para señalizar que ya está terminado.
end

