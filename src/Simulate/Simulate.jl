using GLMakie
GLMakie.activate!()
using StatsBase

function launch_interactive_simulation(d::Dict)
    @unpack seed, pr, dr, mr, fitness, scenario, treatment = d

    model = model_init(pr=0.027, dr=0.015, mr=0.01, scenario=scenario, fitness=fitness,treatment=treatment, seed=0)

    #We plot the genotype of each cell with a different color.
    genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
    genotype_bits = [x for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
    order = Dict(zip(genotype_bits,1:length(genotype_bits)))

    #Functions to get a different color for each genotype
    genotypecolor(a) = get(colorschemes[:hsv], order[a.genotype], (1,length(genotypes)+1))
    genotypecolor_legend(a) = get(colorschemes[:hsv], a, (1,length(genotypes)+1))

    #We make a dynamic plot
    figure, _ = abmplot(model;agent_step! = agent_step!,model_step! = model_step!,ac = genotypecolor,as=0.5)

    #We create a legend for the genotypes
    Legend(figure[1, 2],
        [MarkerElement(color = genotypecolor_legend(a), marker = '■', markersize = 15, strokecolor = :black) for a in 1:length(genotypes)],
        genotypes, patchsize = (20, 20), rowgap = 1)

    #We display the figure
    display(figure)
end

function simulate(d::Dict,steps)
    @unpack seed, pr, dr, mr, fitness, scenario, treatment = d
    
    fulld = copy(d)
    agent_collect = [(:genotype, f) for f in genotype_fraction_function_generator(fitness)]
    model = model_init(pr=pr, dr=dr, mr=mr, scenario=scenario, fitness=fitness,treatment=treatment, seed=seed)
    #We stop (not a typo, stop != step) early if a size of max or 0 is reached
    step = create_stop_function(steps,treatment.detecting_size*1.5)

    adata, _ = run!(model, agent_step!, model_step!, step; adata = agent_collect)
    
    phylogeny = countmap([x.phylogeny for x in allagents(model)])
    fulld["Phylogeny"] = phylogeny

    fulld["TTP"] = get_TTP(adata,treatment.detecting_size*1.2)
    fulld["Diversity"] = get_diversity(adata)

    genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
    pushfirst!(genotypes,"step")
    rename!(adata,genotypes)
    fulld["Genotypes"]=adata

    

    return fulld
end