module TumorModel
    export model_init, create_stop_function, agent_step!, model_step!

    using DrWatson
    using Agents, Random
    using Agents.DataFrames, Agents.Graphs
    using Distributions
    using InteractiveDynamics
    using DataStructures
    using StatsBase

    using TumorSim.Scenario
    using TumorSim.Treatment

    #We create the Cell agent
    @agent Cell GridAgent{3} begin
            time_alive::Int  # Time the cell has been alive
            genotype::BitArray # Genotype of the cell
            phylogeny::Array{Int} # Phylogeny of the cell
            inhibited_by::Array{BitArray} # List of cells inhibiting this cell
    end

    #We initialize the model according to some parameters.
    function model_init(;seed::Int64,death_rate::Float64,mutation_rate::Float64,migration_rate::Float64,fitness::Dict{Vector{Int64},Float64},scenario::ScenarioObject,treatment::TreatmentObject,interaction_rule::Symbol)
        #We need to do this to reuse the treatment in paramscan
        treatment = TreatmentObject(treatment.detecting_size,
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

        current_size::Int=length(cell_pos)

        ngenes::Int=length(collect(keys(fitness))[1])
        
        fitness = copy(fitness)
        fitness = DefaultDict(0,Dict(zip([BitArray(i) for i in keys(fitness)],[fitness[i] for i in keys(fitness)])))

        rng = MersenneTwister(seed)

        space = GridSpace((x, y, z),periodic=false,metric=:euclidean) 
        
        abs_death_rate = fitness[BitArray([false for x in 1:ngenes])]*death_rate

        status = copy(treatment.active)

        properties=@dict(death_rate,mutation_rate,fitness,treatment,scenario,current_size,ngenes,migration_rate,abs_death_rate,interaction_rule,status)

        scheduler = Schedulers.Randomly()

        model = ABM(Cell, space;properties, rng, scheduler) 
        #we create each cell
        for cell in cell_pos
            add_agent!((cell[1],cell[2],cell[3]),model,0,BitArray([false for x in 1:ngenes]),[],[]) # With this one we use the scenario
        end

        return model
    end

    #Function to get the number of cells "near" each cell.

    #Step evey agent, updating its parameters and then reproducing, moving and dying.
    function agent_step!(agent, model)
        agent.time_alive += 1
        
        if die!(agent, model)
            return
        end

        if reproduce!(agent, model) #We want to stop doing things if the cell has died.
            return
        end
        
        move!(agent, model)

    end

    #We use the model step to evaluate the treatment and randomize cell positions
    function model_step!(model)
        current_size::Int = nagents(model)
        model.current_size = current_size
        if model.treatment.detected
            if current_size < (model.treatment.starting_size * model.treatment.detecting_size * model.treatment.pausing_size)
                model.treatment.active = false
            end
            if current_size > (model.treatment.starting_size * model.treatment.detecting_size)
                model.treatment.active = true
            end
        else
            if current_size > model.treatment.detecting_size
                model.treatment.detected = true
            end
        end
        #We randomize positions if mix is active
        if model.scenario.mix
            agents = allagents(model)
            usedpositions = [agent.pos for agent in agents]
            shuffle!(model.rng,usedpositions)
            for (agent,pos) in zip(agents,usedpositions)
                move_agent!(agent,pos,model)
            end
        end
        #We record the status of the treatment for this step
        model.status = copy(model.treatment.active)
    end

    #We stop if any of this conditions are met.
    function create_stop_function(steps::Int,stop_size::Int)
        function step(model,s)
            if nagents(model)==0
                return true
            end
            if nagents(model)>=stop_size && model.treatment.detected
                return true
            end
            if s==steps
                return true
            end
                return false
        end
        return step
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
        if genes!=[] 
            for gene in genes
                if rand(model.rng) < model.mutation_rate
                    agent.genotype[gene]=true
                    push!(agent.phylogeny,gene)
                end
            end
        end
    end

    #Reproduce, creating a new cell in the same space with a probability that decreases with how many cells are already in its space.
    #With a probability (the kill rate of the treatment), the cell is subjected to a treatment check.
    #Returns true if the cell has died.
    function reproduce!(agent,model)
        if rand(model.rng) < model.fitness[agent.genotype]
            npos = nearby_positions(agent,model,1)
            if model.interaction_rule==:contact
                aviable_pos = [pos for pos in npos if isempty(pos,model)]
            elseif model.interaction_rule==:hierarchical_voter
                aviable_pos = [pos for pos in npos if isempty(pos,model) || model.fitness[collect(agents_in_position(pos,model))[1].genotype]<model.fitness[agent.genotype]]
            end
            
            if aviable_pos!=[]
                agent.inhibited_by = []
                if rand(model.rng) < model.treatment.kill_rate
                    if treat!(agent,model)
                        return true
                    end
                end
                newgenom::BitArray = copy(agent.genotype)
                newphylo::Array{Int64} = copy(agent.phylogeny)
                newpos = sample(model.rng,aviable_pos)
                if model.interaction_rule==:hierarchical_voter
                    if !isempty(newpos,model)
                        kill_agent!(collect(agents_in_position(newpos,model))[1],model)
                    end
                end
                newagent = add_agent!(newpos,model,0,newgenom,newphylo,[])
                mutate!(newagent,model)
                kill_non_viable!(newagent, model)
                mutate!(agent,model)
                if kill_non_viable!(agent, model)
                    return true
                end
            else
                agent.inhibited_by = [collect(agents_in_position(pos,model))[1].genotype for pos in npos]
            end
        end
        return false
    end

    #Move every cell to a random nearby space
    function move!(agent, model)
        if rand(model.rng) < model.migration_rate
            npos = nearby_positions(agent,model,1)
            empty_pos = [pos for pos in npos if isempty(pos,model)]
            if empty_pos!=[]
                newpos = sample(model.rng,empty_pos)
                move_agent!(agent,newpos,model)
            end
        end
    end

    #die, with a probability that increases with the number of cells that are in its space. returns true if the cell has died.
    function die!(agent, model)

        #Base apoptosis rate (Turnover)
        if rand(model.rng) < model.abs_death_rate
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

end
using .TumorModel