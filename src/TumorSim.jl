module TumorSim
export create_scenario, model_init, agent_step!, model_step!,genotype_fraction_function_generator, bit_2_int, launch_simulation, plot_genotypes, Treatment
#We import everything we need
cd(@__DIR__) #src

 using Agents, Random
 using Agents.DataFrames, Agents.Graphs
 using Distributions
 using DrWatson: @dict
 using GLMakie
GLMakie.activate!()


 using InteractiveDynamics
 using StatsBase
 using ColorSchemes
 using DataStructures
 using VegaLite
 using DataFrames

 using RCall

 using Images
 using FileIO

 using CSV

#We create the Cell agent
 @agent Cell GridAgent{3} begin
        time_alive::Int  # Time the cell has been alive
        near_cells::Number # Number of cells in the neighborhood
        genotype::BitArray # Genotype of the cell
end

#We initialize the model according to some parameters.
 function  model_init(;seed,pr,dr,mr,fitness,scenario,treatment)
    #We need to do this to reuse the treatment in paramscan
    treatment = Treatment(treatment.detecting_size,
                            treatment.starting_size,
                            treatment.pausing_size,
                            treatment.resistance_gene,
                            treatment.kill_rate,
                            treatment.active,
                            treatment.detected)

    x = scenario.x
    y = scenario.y
    z = scenario.z
    cell_pos = scenario.cell_pos
    wall_pos = scenario.wall_pos
    current_size=length(cell_pos)

    ngenes=length(collect(keys(fitness))[1])
    
    fitness = Dict(zip([BitArray(i) for i in keys(fitness)],[fitness[i] for i in keys(fitness)]))
    fitness=DefaultDict(0,fitness)

    rng = MersenneTwister(seed)

    if y!=0 && z!=0
        space = GridSpace((x, y, z),periodic=false)
        #we create the walls visualization matrix
        wall_matrix=zeros(Int8, x, y, z)
        for t in wall_pos
            i,j,k=t
            wall_matrix[i,j,k]=1
        end

        properties=@dict(pr,dr,mr,fitness,wall_pos,wall_matrix,treatment,scenario,current_size)
        model = ABM(Cell, space;properties, rng) 
        #we create each cell
        for cell in cell_pos
            add_agent!((cell[1],cell[2],cell[3]),model,0,0,BitArray([false for x in 1:ngenes])) # With this one we use the scenario
        end
    else
        space = GridSpace((1,1,1),periodic=false)
        wall_matrix=zeros(Int8, 1)
        properties=@dict(pr,dr,mr,fitness,wall_pos,wall_matrix,treatment,scenario,current_size)
        model = ABM(Cell, space;properties, rng) 
        for cell in cell_pos
            add_agent!(model,0,0,BitArray([false for x in 1:ngenes])) # With this one we use the scenario
        end
    end

    return model
end

#Function to get the number of cells "near" each cell.
 function get_near!(agent,model)
    if model.scenario.z==0 #If we are in 0D

        bin = Binomial(length(model.agents),1/model.scenario.x)
        return (length(model.agents)/model.scenario.x)/(1-pdf(bin,0)) #We calculate the mean number of cells in each cell´s space using a binomial distribution.

    end
    return length(ids_in_position(agent, model))
end

#Step evey agent, updating its parameters and then reproducing, moving and dying.
 function agent_step!(agent, model)
    agent.time_alive += 1
    agent.near_cells = get_near!(agent,model)
    
    if reproduce!(agent, model) #We want to stop doing things if the cell has died.
        return
    end
    if model.scenario.z!=0 #We dont move if we are in 0D
        move!(agent, model)
    end
    if die!(agent, model)
        return
    end
end

#We use the model step to evaluate the treatment
 function model_step!(model)
    current_size = length(model.agents)
    model.current_size = current_size
    if model.treatment.detected
        if current_size < model.treatment.pausing_size
            model.treatment.active = false
        end
        if current_size > model.treatment.starting_size
            model.treatment.active = true
        end
    else
        if current_size > model.treatment.detecting_size
            model.treatment.detected = true
        end
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
    pr = model.pr*model.fitness[agent.genotype]
    pid = agent.pos
    newgenom = copy(agent.genotype)
    prob = pr/(get_near!(agent,model)^2)
    if rand(model.rng) < prob/(1+prob)
        if rand(model.rng) < model.treatment.kill_rate
            if treat!(agent,model)
                return true
            end
        end
        newagent = add_agent!(pid,model,0,0,newgenom)
        mutate!(newagent,model)
        kill_non_viable!(newagent, model)

        mutate!(agent,model)
        if kill_non_viable!(agent, model)
            return true
        end
    end
    return false
end

#Move every cell to a random nearby space ONLY if your space is "crowded", crowded for example is more than 1 cell in your space 
 function move!(agent, model)
    pos = agent.pos
    nearby = [x for x in nearby_positions(agent,model,1)]

    setdiff!(nearby,model.wall_pos)
    nearby_density = [1/(length(ids_in_position(x,model))+1) for x in nearby]
    
    push!(nearby,pos)
    push!(nearby_density,1/length(ids_in_position(agent,model)))

    newpos = sample(model.rng, nearby, Weights(nearby_density))
    if length(ids_in_position(agent, model)) > 1
        move_agent!(agent,newpos, model)
    end
end

#die, with a probability that increases with the number of cells that are in its space. returns true if the cell has died.
 function die!(agent, model)
    prob = model.dr*(get_near!(agent,model)^2)
    if rand(model.rng) < prob/(1+prob)
        kill_agent!(agent, model)
        return true
    end
    return false
end

#we kill all non viable agents instantly to make our data cleaner
 function kill_non_viable!(agent, model)
    if !(agent.genotype in keys(model.fitness))
        kill_agent!(agent,model)
        return true
    end
    return false
end
#A generator that returns a list of functions that each get the number of cells of each genotype given a number of genes.
 function genotype_fraction_function_generator(fitness)
    functions = []
    for i in [x for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
        compare = i
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
        genotype::Vector{Int64}=[]
        for j in 1:g
            push!(genotype,Int(fitness[i,j]))
        end
        push!(dictionary,(genotype=>Real(fitness[i,g+1])))
    end
    return dictionary
end

#Function to read a scenario from a file. Reads cells and walls and defines the size of the grid.
#IDEA: associate each color with a genotype, to allow for more flexible scenarios.
 function get_2D_scenario_from_bmp(filepath)
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
    return Scenario(dims[2],dims[1],1,cell_pos,wall_pos)
end

#Function to go from BitArray to Int. Taken from https://discourse.julialang.org/t/parse-an-array-of-bits-bitarray-to-an-integer/42361/5
 function bit_2_int(arr)
    arr = reverse(arr)
    sum(((i, x),) -> Int(x) << ((i-1) * sizeof(x)), enumerate(arr.chunks))
end

#We define what a treatment is
 mutable struct Treatment
    detecting_size::Int
    starting_size::Int
    pausing_size::Int 
    resistance_gene::Int
    kill_rate::Float16
    active::Bool
    detected::Bool
end

#We define the scenario and the functions to create it using multiple dispatch.
 struct Scenario
    x::Int
    y::Int
    z::Int
    cell_pos::Vector
    wall_pos::Vector
end

#0D
 function create_scenario(size::Int64,ncells::Int)   
    return Scenario(size,0,0,[(1,1,1) for i in 1:ncells],[])
end

#1D
 function create_scenario(size::Tuple{Int64},ncells::Int,cell_pos::String="center",wall_pos::Vector=[])
    if cell_pos=="random"
        cell_pos = [((rand(1:size[1])),1,1) for i in 1:ncells]
    elseif cell_pos=="center"
        cell_pos = [(Int(floor(size[1]/2))+1,1,1) for i in 1:ncells]
    else
        print("Error, cell_pos option not found.")
        return
    end
        return Scenario(size[1],1,1,cell_pos,wall_pos)
end

 function create_scenario(size::Tuple{Int64},cell_pos::Vector{Tuple{Int64}},wall_pos=[])
        return Scenario(size[1],1,1,[(x[1],1,1) for x in cell_pos],wall_pos)
end

#2D
 function create_scenario(size::Tuple{Int64, Int64},ncells::Int,cell_pos::String="center",wall_pos::Vector=[])
    if cell_pos=="random"
        cell_pos = [((rand(1:size[1])),(rand(1:size[2])),1) for i in 1:ncells]
    elseif cell_pos=="center"
        cell_pos = [(Int(floor(size[1]/2))+1,Int(floor(size[2]/2))+1,1) for i in 1:ncells]
    else
        print("Error, cell_pos option not found.")
        return
    end
        return Scenario(size[1],size[2],1,cell_pos,wall_pos)
end

 function create_scenario(size::Tuple{Int64,Int64},cell_pos::Vector{Tuple{Int64, Int64}},wall_pos::Vector=[])
        return Scenario(size[1],size[2],1,[(x[1],x[2],1) for x in cell_pos],wall_pos)
end

#3D
 function create_scenario(size::Tuple{Int64, Int64, Int64},ncells::Int,cell_pos::String="center",wall_pos::Vector=[])
    if cell_pos=="random"
        cell_pos = [((rand(1:size[1])),(rand(1:size[2])),(rand(1:size[3]))) for i in 1:ncells]
    elseif cell_pos=="center"
        cell_pos = [(Int(floor(size[1]/2))+1,Int(floor(size[2]/2))+1,Int(floor(size[3]/2))+1) for i in 1:ncells]
    else
        print("Error, cell_pos option not found.")
        return
    end
        return Scenario(size[1],size[2],size[3],cell_pos,wall_pos)
end

 function create_scenario(size::Tuple{Int64,Int64,Int64},cell_pos::Vector{Tuple{Int64, Int64,Int64}},wall_pos::Vector=[])
        return Scenario(size[1],size[2],size[3],[(x[1],x[2],x[3]) for x in cell_pos],wall_pos)
end

function launch_simulation(model;fitness,agent_step!,model_step!)
    #We plot the genotype of each cell with a different color.
    genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
    genotype_bits = [x for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
    order = Dict(zip(genotype_bits,1:length(genotype_bits)))

    #Functions to get a different color for each genotype
    genotypecolor(a) = get(colorschemes[:hsv], order[a.genotype], (1,length(genotypes)+1))
    genotypecolor_legend(a) = get(colorschemes[:hsv], a, (1,length(genotypes)+1))

    #We make a static plot
    #figure, _ = abmplot(model;ac = genotypecolor,as=8,am='■',heatarray,heatkwargs)

    #We make a dynamic plot
    figure, _ = abmplot(model;agent_step! = agent_step!,model_step! = model_step!,ac = genotypecolor,as=0.5)


    #We create a legend for the genotypes
    Legend(figure[1, 2],
        [MarkerElement(color = genotypecolor_legend(a), marker = '■', markersize = 15, strokecolor = :black) for a in 1:length(genotypes)],
        genotypes,
        patchsize = (20, 20), rowgap = 1)


    #We display the figure
    display(figure)
end

function plot_genotypes(adata,fitness,mode="absolute")
    #And lastly we can make plots of both the total number of cells of each genotype
    genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
    stacked = stack(adata,genotypes)
    if mode=="absolute"
        stacked |>
        @vlplot(:area, x=:step, y={:value, stack=:zero}, color="variable:n")
    elseif mode=="relative"
        #And the relative number of cells of each genotype
        stacked |>
        @vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")
    else
        println("Mode not supported")
    end
end

end
