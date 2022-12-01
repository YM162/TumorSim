#We import the common functions and model.
using Distributed
cd(@__DIR__)
@everywhere begin
    using TumorSim
    using Agents
    using DataFrames
    using CSV
    #We create a fitness landscape
    fitness=Dict([0,0,0]=>1, 
                [1,0,0]=>1.3,
                [0,1,0]=>1.2,
                [1,1,0]=>1.5,
                [1,1,1]=>1.2)

    #Scenarios for every dimension with 1.000.000 spaces each.
    scenario_0D = create_scenario((1000000),10)
    scenario_1D = create_scenario((1000000,),10,"center")
    scenario_2D = create_scenario((1000,1000),10,"center")
    scenario_3D = create_scenario((100,100,100),10,"center")

    #We create the therapies with detecting size 3000, starting size 2000 and pausing size 1000/0. Gene of resistance is 3. kill rate is 0.75
    adaptive_therapy = Treatment(3000, 2000, 1000, 3, 0.75, false, false) 

    continuous_therapy = Treatment(3000, 2000, 0, 3, 0.75, false, false) 

    #We define what we want to collect
    agent_collect = [(:genotype, f) for f in genotype_fraction_function_generator(fitness)]

    #We define the parameters, if it is a vector, it will make a simulation with each of it's values.
    parameters = Dict(
        :pr => 0.027,
        :dr => 0.015,
        :mr => 0.01,
        :fitness => fitness,
        :scenario => [scenario_0D, scenario_1D, scenario_2D, scenario_3D],
        :treatment => [adaptive_therapy, continuous_therapy],
        :seed => 3
    )

end

#We run the simulation
adata, _ = paramscan(parameters, model_init; adata = agent_collect, agent_step! = agent_step!, model_step! = model_step!, n = 1000, parallel = true, showprogress = true,include_constants=true)

#we rename the columns to get a clean "data" DataFrame.
genotypes = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
pushfirst!(genotypes,"step")
#We need to push the parameters we scan to the names
push!(genotypes,["proliferation_rate","death_rate","mutation_rate","fitness","scenario","treatment","seed"]...)

rename!(adata,genotypes)

#We save the data
CSV.write("./results/paramscan_result_1.csv",adata)