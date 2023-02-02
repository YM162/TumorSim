function simulate(d::Dict,max_steps::Int)
    @unpack seed, pr, dr, mr, fitness, scenario, treatment = d
    
    fulld::Dict = copy(d)
    agent_collect::Array = [(:genotype, f) for f in genotype_fraction_function_generator(fitness)]
    model = model_init(pr=pr, dr=dr, mr=mr, scenario=scenario, fitness=fitness,treatment=treatment, seed=seed)
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

    fulld["s_dims"] = string((scenario.x,scenario.y,scenario.z))
    fulld["s_initial_cells"] = length(scenario.cell_pos)

    fulld["t_detecting_size"] = treatment.detecting_size
    fulld["t_starting_size"] = treatment.starting_size
    fulld["t_pausing_size"] = treatment.pausing_size
    fulld["t_resistance_size"] = treatment.resistance_gene
    fulld["t_kill_rate"] = treatment.kill_rate
    return fulld
end