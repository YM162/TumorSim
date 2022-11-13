#We import everything we need
cd(@__DIR__) #src
using Agents, Random
using Agents.DataFrames, Agents.Graphs
using Distributions: Poisson, DiscreteNonParametric
using DrWatson: @dict
using CairoMakie
CairoMakie.activate!()


using InteractiveDynamics
using StatsBase
using ColorSchemes
using DataStructures
using VegaLite
using DataFrames

using RCall

using Images
using FileIO


#We create the Cell agent
@agent Cell GridAgent{2} begin
        time_alive::Int  # Time the cell has been alive
        near_cells::Int # Number of cells in the neighborhood
        genotype::BitArray # Genotype of the cell
end

#We initialize the model according to some parameters.
function  model_init(;pr,dr,mr,seed,h,w,ngenes,fitness,cell_pos,wall_pos,treatment)
    
    rng = MersenneTwister(seed)

    space = GridSpace((h, w))
    #we create the walls visualization matrix
    wall_matrix=zeros(Int8, h, w)
    for t in wall_pos
        i,j=t
        wall_matrix[i,j]=1
    end

    properties=@dict(ngenes,pr,dr,mr,fitness,wall_pos,wall_matrix,treatment)
    model = ABM(Cell, space;properties, rng) 
    #we create each cell
    for cell in cell_pos
        add_agent!((cell[1],cell[2]),model,0,0,BitArray([false for x in 1:ngenes])) # With this one we use the scenario
    end


    return model
end

#Function to get the number of cells "near" each cell. As i dont know what should count as near (only the space or a mean of the 8 surrounding spaces? for example) i prefer to define it in a single place and then change it.
function get_near!(agent,model)
    length(ids_in_position(agent, model))
end

#Step evey agent, updating its parameters and then reproducing, moving and dying.
function agent_step!(agent, model)
    if agent.time_alive == 0 #We mutate the cell if it has just been born
        mutate!(agent,model)
    end
    agent.time_alive += 1
    agent.near_cells = get_near!(agent,model)
    
    if reproduce!(agent, model) #We want to stop doing things if the cell has died.
        return
    end
    move!(agent, model)
    if die!(agent, model)
        return
    end
end

#We use the model step to evaluate the treatment
function model_step!(model)
    current_size = length(model.agents)
    if current_size < model.treatment.pausing_size
        model.treatment.active = false
    end
    if current_size > model.treatment.starting_size
        model.treatment.active = true
    end
end

#If the cell is susceptible to the treatment, and treatment is active, it dies. Returns true if the cell has dies
function treat!(agent,model)
    if model.treatment.active && agent.genotype[model.treatment.resistance_gene]!=1
        kill_agent!(agent,model)
        return true
    end
    return false
end

#with a probability p choose a random non mutated gene and mutate it.
function mutate!(agent,model)
    genes=findall(agent.genotype .!=1)
    if genes!=[] && rand(model.rng) < model.mr
        agent.genotype[rand(model.rng,genes)]=true
    end
end

#Reproduce, creating a new cell in the same space with a probability that decreases with how many cells are already in its space.
#With a probability (the kill rate of the treatment), the cell is subjected to a treatment check.
#Returns true if the cell has died.
function reproduce!(agent,model)
    pr = model.pr*model.fitness[bit_2_int(agent.genotype)]
    pid = agent.pos
    newgenom = copy(agent.genotype)
    if rand(model.rng) < pr/(get_near!(agent,model)^2)
        if rand(model.rng) < model.treatment.kill_rate
            if treat!(agent,model)
                return true
            end
        end
        add_agent!(pid,model,0,0,newgenom)
    end
    return false
end

#Move every cell to a random nearby space ONLY if your space is "crowded", crowded for example is more than 1 cell in your space 
function move!(agent, model)
    pos = agent.pos
    nearby = [x for x in nearby_positions(agent,model,1)]
    setdiff!(nearby,model.wall_pos)
    if length(nearby)==0
        return
    end
    newpos = rand(model.rng, nearby)
    if length(ids_in_position(agent, model)) > 1
        move_agent!(agent,newpos, model)
    end
end

#die, with a probability that increases with the number of cells that are in its space. returns true if the cell has died.
function die!(agent, model)
    if rand(model.rng) < model.dr*(get_near!(agent,model)^2)
        kill_agent!(agent, model)
        return true
    end
    return false
end

#A generator that returns a list of functions that each get the number of cells of each genotype given a number of genes.
function genotype_fraction_function_generator(ngenes)
    functions = []
    for i in 0:((2^ngenes)-1)
        compare = reverse(digits(i, base=2, pad=ngenes))
        func = function get_perc(x)
            len = length(findall([string(y)[5:end] for y in x] .== string(compare)))
            return len
        end
        push!(functions,func)
    end
    return functions
end

#Function to create a random fitness landscape using the OncoSimulR library.
#How can i feed it a seed??
function OncoSimulR_rfitness(;g,c,sd)
    R"library(OncoSimulR)"
    fitness = R"rfitness(g=$g ,c=$c ,sd=$sd )"
    rows=2^g
    dictionary=Dict()
    for i in 1:rows
        genotype=BitArray([])
        for j in 1:g
            push!(genotype,fitness[i,j])
        end
        push!(dictionary,(bit_2_int(genotype)=>Float64(fitness[i,g+1])))
    end
    return dictionary
end

#Function to go from BitArray to Int. Taken from https://discourse.julialang.org/t/parse-an-array-of-bits-bitarray-to-an-integer/42361/5
function bit_2_int(arr)
    arr = reverse(arr)
    sum(((i, x),) -> Int(x) << ((i-1) * sizeof(x)), enumerate(arr.chunks))
end

#Function to read a scenario from a file. Reads cells and walls and defines the size of the grid.
#IDEA: associate each color with a genotype, to allow for more flexible scenarios.
function get_scenario(filepath)
    scenario = FileIO.load(filepath)
    dims = size(scenario)
    cell_pos=Tuple[]
    wall_pos=Tuple[]
    for i in 1:dims[1]
        for j in 1:dims[2]
            if scenario[i,j]==RGB{N0f8}(0.133,0.694,0.298) #Initial cells are painted with the default deep green color from MS Paint.
                push!(cell_pos,(j,(dims[1]+1)-i))
            elseif scenario[i,j]==RGB{N0f8}(0,0,0) #Walls are painted black
                push!(wall_pos,(j,(dims[1]+1)-i))
            end
        end
    end

    return dims[2],dims[1],cell_pos,wall_pos #We return height, width, the positions of the cells and the positions of the walls
end

#We define what a treatment is
mutable struct Treatment
    pausing_size::Int 
    starting_size::Int
    resistance_gene::Int
    kill_rate::Float16
    active::Bool
end