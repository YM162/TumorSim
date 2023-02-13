using DrWatson, Test
@quickactivate "TumorSim"

using TumorSim

using Distributed
using Agents, Random
using Agents.DataFrames, Agents.Graphs
using Distributions
using InteractiveDynamics
using StatsBase
using ColorSchemes
using DataStructures
using DataFrames

using Genie

# Run test suite
println("Starting tests")
ti = time()

@testset "Fitness tests" begin
    #genotype_fraction_function_generator tests
    fitness=Dict([0,0,0]=>0.027, [1,0,0]=>0.035)
    functions = TumorSim.genotype_fraction_function_generator(fitness)
    @test length(functions) == 2
    @test functions[1]([BitArray([0,0,0]),BitArray([1,0,0]),BitArray([0,0,0])]) == 2
    @test functions[2]([BitArray([0,0,0]),BitArray([1,0,0]),BitArray([0,0,0])]) == 1

    #bit_2_int tests
    @test TumorSim.bit_2_int(BitArray([1,0,0,0,0])) == 16
    @test TumorSim.bit_2_int(BitArray([0,0,0,0,0])) == 0
    @test TumorSim.bit_2_int(BitArray([1,1,1])) == 7
end

@testset "Scenario tests" begin
    #0D
    scenario = create_scenario(10,3)
    @test scenario.x == 10
    @test scenario.y == 0
    @test scenario.z == 0
    @test length(scenario.cell_pos) == 3
    #1D
    scenario = create_scenario((40,),4,"random",[(10,)])
    scenario = create_scenario((40,),4,"center",[(10,)])
    @test scenario.x == 40
    @test scenario.y == 1
    @test scenario.z == 1
    @test length(scenario.cell_pos) == 4

    scenario = create_scenario((40,),[(1,)],[(5,)])
    @test scenario.x == 40
    @test scenario.y == 1
    @test scenario.z == 1
    @test length(scenario.cell_pos) == 1
    #2D
    scenario = create_scenario((15,20),5,"random",[(1,1)])
    scenario = create_scenario((15,20),5,"center",[(1,1)])
    @test scenario.x == 15
    @test scenario.y == 20
    @test scenario.z == 1
    @test length(scenario.cell_pos) == 5

    scenario = create_scenario((15,20),[(1,1),(2,2)],[(5,5)])
    @test scenario.x == 15
    @test scenario.y == 20
    @test scenario.z == 1
    @test length(scenario.cell_pos) == 2
    #3D
    scenario = create_scenario((16,17,18),19,"random",[(1,1,1)])
    scenario = create_scenario((16,17,18),19,"center",[(1,1,1)])
    @test scenario.x == 16
    @test scenario.y == 17
    @test scenario.z == 18
    @test length(scenario.cell_pos) == 19

    scenario = create_scenario((16,17,18),[(1,1,1),(2,2,2),(3,3,3)],[(5,5,5)])
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
    treatment = create_treatment(3000,2000,1000,3,0.75)
    scenario = create_scenario(10,5)
    fitness=Dict([0,0,0]=>0.027, [1,0,0]=>0.035)
    model = TumorSim.model_init(dr=0.55, mr=0.01, scenario=scenario, fitness=fitness,treatment=treatment, cr = 0.2, seed=0)
    agent = collect(allagents(model))[1]
    @test nagents(model) == 5
    @test agent.time_alive == 0
    @test agent.genotype == BitArray([0,0,0])
    @test agent.phylogeny == []
    #3D
    scenario = create_scenario((10,10,10),5,"center",[(1,1,1)])
    model = TumorSim.model_init(dr=0.55, mr=0.01, scenario=scenario, fitness=fitness,treatment=treatment, cr = 0.2, seed=0)
    agent = collect(allagents(model))[1]
    @test nagents(model) == 5
    @test agent.time_alive == 0
    @test agent.genotype == BitArray([0,0,0])
    @test agent.phylogeny == []
end
#Prepare a simulation
treatment = create_treatment(3000,0.65,0.5,3,0.75)
scenario = create_scenario((100,100,100),5,"center",[(1,1,1)])
fitness=Dict([0,0,0]=>0.027, 
            [1,0,0]=>0.035,
            [0,1,0]=>0.032,
            [1,1,0]=>0.040,
            [1,1,1]=>0.032)

params = Dict(
    "dr" => [0.55,2],
    "mr" => 0.1,
    "scenario" => scenario, 
    "fitness" => fitness,
    "treatment" => treatment,
    "cr" => 0.2,
    "seed" => 0
)

dicts = dict_list(params)
steps=3000
results = pmap(simulate,dicts,fill(steps,length(dicts)))
adata1 = results[1]["Genotypes"]
adata2 = results[2]["Genotypes"]

@testset "Simulate tests" begin
    @test length(results) == 2
    @test length(results[1]) == 19
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

@testset "Dashboard tests" begin
    @test launch_dashboard().webserver._isexception == false
    #Maybe we can test the buttons somehow?
    sleep(5)
    @test kill_dashboard() == Genie.Server.ServersCollection[]
end

@testset "Worker tests" begin
    rm(projectdir("logs","progress","test"),force=true)
    rm(projectdir("data","simulations","test.bson"),force=true)
    
    worker_path = srcdir("Dashboard","simulation_worker.jl")

    worker = `julia --check-bounds=yes --project=$(projectdir()) $worker_path 0.5 0.05 0.5 0.01 0.01 0.01 1000000 3 10 3000 100 3000 0.65 0.1 0.65 0.5 0.1 0.5 0.75 0.05 0.75 0.2 0.1 0.2 1 test`
    run(worker)

    @test isfile(projectdir("logs","progress","test")) == true
    @test isfile(projectdir("data","simulations","test.bson")) == true
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
