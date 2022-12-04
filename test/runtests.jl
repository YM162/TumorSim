using DrWatson, Test
@quickactivate "TumorSim"

using Distributed
using Agents, Random
using Agents.DataFrames, Agents.Graphs
using Distributions
using InteractiveDynamics
using StatsBase
using ColorSchemes
using DataStructures
using DataFrames

# Run test suite
println("Starting tests")
ti = time()

#We include everything we need
include(srcdir("Fitness/Fitness.jl"))
include(srcdir("Scenario/Scenario.jl"))
include(srcdir("Treatment/Treatment.jl"))
include(srcdir("TumorModel/TumorModel.jl"))
include(srcdir("Simulate/Simulate.jl"))
include(srcdir("Analysis/Analysis.jl"))

@testset "Fitness tests" begin
    #genotype_fraction_function_generator tests
    fitness=Dict([0,0,0]=>1, [1,0,0]=>1.3)
    functions = genotype_fraction_function_generator(fitness)
    @test length(functions) == 2
    @test functions[1]([BitArray([0,0,0]),BitArray([1,0,0]),BitArray([0,0,0])]) == 2
    @test functions[2]([BitArray([0,0,0]),BitArray([1,0,0]),BitArray([0,0,0])]) == 1

    #bit_2_int tests
    @test bit_2_int(BitArray([1,0,0,0,0])) == 16
    @test bit_2_int(BitArray([0,0,0,0,0])) == 0
    @test bit_2_int(BitArray([1,1,1])) == 7
end

@testset "Scenario tests" begin
    #0D
    scenario = create_scenario(10,3)
    @test scenario.x == 10
    @test scenario.y == 0
    @test scenario.z == 0
    @test length(scenario.cell_pos) == 3
    #1D
    scenario = create_scenario((40,),4)
    @test scenario.x == 40
    @test scenario.y == 1
    @test scenario.z == 1
    @test length(scenario.cell_pos) == 4

    scenario = create_scenario((40,),[(1,)])
    @test scenario.x == 40
    @test scenario.y == 1
    @test scenario.z == 1
    @test length(scenario.cell_pos) == 1
    #2D
    scenario = create_scenario((15,20),5)
    @test scenario.x == 15
    @test scenario.y == 20
    @test scenario.z == 1
    @test length(scenario.cell_pos) == 5

    scenario = create_scenario((15,20),[(1,1),(2,2)])
    @test scenario.x == 15
    @test scenario.y == 20
    @test scenario.z == 1
    @test length(scenario.cell_pos) == 2
    #3D
    scenario = create_scenario((16,17,18),19)
    @test scenario.x == 16
    @test scenario.y == 17
    @test scenario.z == 18
    @test length(scenario.cell_pos) == 19

    scenario = create_scenario((16,17,18),[(1,1,1),(2,2,2),(3,3,3)])
    @test scenario.x == 16
    @test scenario.y == 17
    @test scenario.z == 18
    @test length(scenario.cell_pos) == 3
end

@testset "Treatment tests" begin
    treatment = create_treatment(3000,2000,1000,3,0.75)
    @test treatment.detecting_size == 3000
    @test treatment.starting_size == 2000
    @test treatment.pausing_size == 1000
    @test treatment.resistance_gene == 3
    @test treatment.kill_rate == 0.75
end

@testset "TumorModel tests" begin
    #0D init
    #3D
    treatment = create_treatment(3000,2000,1000,3,0.75)
    scenario = create_scenario(10,5)
    fitness=Dict([0,0,0]=>1, [1,0,0]=>1.3)
    model = model_init(pr=0.027, dr=0.55, mr=0.01, scenario=scenario, fitness=fitness,treatment=treatment, seed=0)
    agent = collect(allagents(model))[1]
    @test nagents(model) == 5
    @test agent.time_alive == 0
    @test agent.genotype == BitArray([0,0,0])
    @test agent.phylogeny == []
    @test get_near!(agent,model) ≈ 1.220971404849699
    #3D
    scenario = create_scenario((10,10,10),5)
    model = model_init(pr=0.027, dr=0.55, mr=0.01, scenario=scenario, fitness=fitness,treatment=treatment, seed=0)
    agent = collect(allagents(model))[1]
    @test nagents(model) == 5
    @test agent.time_alive == 0
    @test agent.genotype == BitArray([0,0,0])
    @test agent.phylogeny == []
    @test get_near!(agent,model) == 5
end
#Prepare a simulation
treatment = create_treatment(3000,2000,1000,3,0.75)
scenario = create_scenario((100,100,100),5)
fitness=Dict([0,0,0]=>1, 
            [1,0,0]=>1.3,
            [0,1,0]=>1.2,
            [1,1,0]=>1.5,
            [1,1,1]=>1.2)

params = Dict(
    "pr" => 0.027,
    "dr" => [0.55,2],
    "mr" => 0.1,
    "scenario" => scenario, 
    "fitness" => fitness,
    "treatment" => treatment,
    "seed" => 0
)

dicts = dict_list(params)
steps=3000
results = pmap(simulate,dicts,fill(steps,length(dicts)))
adata1 = results[1]["Genotypes"]
adata2 = results[2]["Genotypes"]

@testset "Simulate tests" begin
    @test length(results) == 2
    @test length(results[1]) == 18
    @test length(eachcol(adata1)) == 6
    @test length(eachrow(adata1)) != 0
    @test length(eachrow(adata1)) > length(eachrow(adata2))
end

@testset "Analysis tests" begin
    TTP1 = get_TTP(adata1,3200)
    TTP2 = get_TTP(adata2,3200)
    @test TTP1 != -1
    @test TTP2 == -1

    diversity = get_diversity(adata1)
    @test length(eachcol(diversity)) == 3
    @test length(eachrow(diversity)) != 0
    @test diversity[!,"species_richness"][1] == 1
    @test diversity[!,"shannon_index"][1] == 0
    @test diversity[!,"evenness"][1] == 0

    @test diversity[!,"species_richness"][500] > 1
    @test diversity[!,"shannon_index"][500] > 0
    @test diversity[!,"evenness"][500] > 0
    @test diversity[!,"evenness"][500] < 1

end

using TumorSim

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
