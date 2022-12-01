using Agents
using TumorSim
using DataFrames

fitness=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.2)

scenario_3D = create_scenario((1,1,1),10,"center")

continuous_therapy = Treatment(3000, 2000, 0, 3, 0.75, false, false)

agent_collect = [(:genotype, f) for f in genotype_fraction_function_generator(fitness)]

parameters = Dict(
        :pr => 0.027,
        :dr => 0.015,
        :mr => 0.01,
        :fitness => fitness,
        :scenario => scenario_3D,
        :treatment => continuous_therapy,
        :seed => 3
)

adata2, _ = paramscan(parameters, model_init; adata = agent_collect, agent_step! = agent_step!, model_step! = model_step!, n = 1, parallel = false, showprogress = false,include_constants=true)

model = model_init(pr=0.027, #Proliferation rate
                    dr=0.015, #Death rate
                    mr=0.01, #Mutation rate
                    scenario=scenario_3D, #The scenario
                    fitness=fitness, #The fitness of each genotype
                    treatment=continuous_therapy, #The treatment we are going to use for the simulation
                    seed=0) #Seed to get reproducible results

steps=1
adata, _ = run!(model, agent_step!, model_step!, steps; adata = agent_collect)

genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
pushfirst!(genotypes,"step")
rename!(adata,genotypes)

#We need to add the other things like plotting, but not before we include them in the package as functions.