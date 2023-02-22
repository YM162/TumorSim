using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using BSON
using DataFrames
using Dates
using Agents
using GLMakie
using DrWatson
using Agents, Random
using Agents.DataFrames, Agents.Graphs
using Distributions
using InteractiveDynamics
using DataStructures
using StatsBase
using ColorSchemes
GLMakie.activate!()

fitness=Dict([0,0,0]=>0.027, 
            [1,0,0]=>0.031,
            [0,1,0]=>0.035,
            [1,1,0]=>0.040,
            [1,1,1]=>0.040)

scenario = create_scenario((100,100),10,"center",false)

adaptive_therapy = create_treatment(3000, 1, 0.5, 3, 0.75) 
continuous_therapy = create_treatment(3000, 1, 0.0, 3, 0.75) 


model = TumorSim.TumorModel.model_init(cost_of_resistance=0.1, death_rate=0.1, mutation_rate=0.01, migration_rate=0.05,
                scenario=scenario, fitness=fitness,treatment=adaptive_therapy, 
                seed=0)

#We plot the genotype of each cell with a different color.
genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> TumorSim.Fitness.bit_2_int(BitArray(x)))]
genotype_bits = [x for x in sort!([x for x in keys(fitness)],by=x -> TumorSim.Fitness.bit_2_int(BitArray(x)))]
order = Dict(zip(genotype_bits,1:length(genotype_bits)))

#Functions to get a different color for each genotype
genotypecolor(a) = get(colorschemes[:hsv], order[a.genotype], (1,length(genotypes)+1))
genotypecolor_legend(a) = get(colorschemes[:hsv], a, (1,length(genotypes)+1))

#We make a dynamic plot
figure, _ = abmplot(model;agent_step! = TumorSim.TumorModel.agent_step!,model_step! = TumorSim.TumorModel.model_step!,ac = genotypecolor,as=0.5)

#We create a legend for the genotypes
Legend(figure[1, 2],
    [MarkerElement(color = genotypecolor_legend(a), marker = '■', markersize = 15, strokecolor = :black) for a in 1:length(genotypes)],
    genotypes, patchsize = (20, 20), rowgap = 1)

#We display the figure
display(figure)
