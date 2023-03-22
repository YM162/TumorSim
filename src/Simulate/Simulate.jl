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
        @unpack seed, death_rate, mutation_rate, fitness, scenario, treatment, migration_rate, interaction_rule = d
        
        fulld::Dict = copy(d)
        f1 = genotype_fraction_function_generator(fitness)
        f2 = inhibited_by_function_generator(fitness)
        agent_collect::Array = [(:genotype, f1),(:inhibited_by, f2,x->x.genotype[treatment.resistance_gene]==1)]
        

        model_collect = [:status]
        
        model = model_init(death_rate=death_rate, mutation_rate=mutation_rate, scenario=scenario, fitness=fitness,interaction_rule=interaction_rule,treatment=treatment,migration_rate=migration_rate, seed=seed)
        #We stop (not a typo, stop != step) early if a size of max or 0 is reached
        step = create_stop_function(max_steps,Int(floor(treatment.detecting_size*1.5)))

        adata::DataFrame, mdata::DataFrame = run!(model, agent_step!, model_step!, step; adata = agent_collect, mdata = model_collect)
        
        fulld["Treatment_status"] = mdata

        gen_df=eachcol(adata)[2]
        ngen_df=DataFrame()
        for row in gen_df
            push!(ngen_df,row)
        end
        insertcols!(ngen_df,1,"step" => eachcol(adata)[1])

        fulld["Genotypes"] = ngen_df

        gen_df2=eachcol(adata)[3]
        ngen_df2=DataFrame()
        for row in gen_df2
            push!(ngen_df2,row)
        end
        insertcols!(ngen_df2,1,"step" => eachcol(adata)[1])
        fulld["Resistant_inhibited_by"] = ngen_df2

        phylogeny = countmap([x.phylogeny for x in allagents(model)])
        fulld["Phylogeny"] = phylogeny
        

        

        fulld["TTP"] = get_TTP(ngen_df,Int(floor(treatment.detecting_size*1.2)))
        fulld["Diversity"] = get_diversity(ngen_df)
        fulld["Divergence"] = get_divergence(ngen_df,ngen_df2)
        fulld["Resistant_on_detection"] = get_resistant_fraction_on_detection(ngen_df,treatment.detecting_size,treatment.resistance_gene)

        fulld["s_dim"] = Int(scenario.x!=1) + Int(scenario.y!=1) + Int(scenario.z!=1)
        
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