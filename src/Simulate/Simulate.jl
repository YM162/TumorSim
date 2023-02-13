module Simulate
    export simulate

    using DataFrames
    using DrWatson
    using Agents, Random
    using Agents.DataFrames, Agents.Graphs
    using StatsBase

    using TumorSim.Fitness
    using TumorSim.TumorModel
    using TumorSim.Analysis

    function simulate(d::Dict,max_steps::Int)
        @unpack seed, dr, mr, fitness, cr, scenario, treatment = d
        
        fulld::Dict = copy(d)
        agent_collect::Array = [(:genotype, f) for f in genotype_fraction_function_generator(fitness)]
        model = model_init(dr=dr, mr=mr, scenario=scenario, fitness=fitness,treatment=treatment,cr = cr, seed=seed)
        #We stop (not a typo, stop != step) early if a size of max or 0 is reached
        step = create_stop_function(max_steps,Int(floor(treatment.detecting_size*1.5)))

        adata::DataFrame, _ = run!(model, agent_step!, model_step!, step; adata = agent_collect)
        
        phylogeny = countmap([x.phylogeny for x in allagents(model)])
        fulld["Phylogeny"] = phylogeny

        fulld["TTP"] = get_TTP(adata,Int(floor(treatment.detecting_size*1.2)))
        fulld["Diversity"] = get_diversity(adata)

        genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
        pushfirst!(genotypes,"step")
        rename!(adata,genotypes)
        fulld["Genotypes"] = adata

        if scenario.x * scenario.y * scenario.z == 0
            fulld["s_dim"] = 0
        else
            fulld["s_dim"] = Int(scenario.x!=1) + Int(scenario.y!=1) + Int(scenario.z!=1)
        end
        
        
        fulld["s_size"] = string((scenario.x,scenario.y,scenario.z))
        fulld["s_initial_cells"] = length(scenario.cell_pos)

        fulld["t_detecting_size"] = treatment.detecting_size
        fulld["t_starting_size"] = treatment.starting_size
        fulld["t_pausing_size"] = treatment.pausing_size
        fulld["t_resistance_gene"] = treatment.resistance_gene
        fulld["t_kill_rate"] = treatment.kill_rate

        return fulld
    end
end
using .Simulate