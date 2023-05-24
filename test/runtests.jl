using DrWatson, Test
@quickactivate "TumorSim"

using TumorSim

using Distributed
using Agents, Random
using Agents.DataFrames, Agents.Graphs
using Distributions
using StatsBase
using DataStructures
using DataFrames


# Run test suite
println("Starting tests")
ti = time()

@testset "Fitness tests" begin
    #bit_2_int tests
    @test TumorSim.bit_2_int(BitArray([1,0,0,0,0])) == 16
    @test TumorSim.bit_2_int(BitArray([0,0,0,0,0])) == 0
    @test TumorSim.bit_2_int(BitArray([1,1,1])) == 7
    #Generate fitness from restrictions
    restrictions_AND= [[1,2,1],[1,3,1],[2,4,1],[3,4,1]]
    mitosis_probabilities = TumorSim.Fitness.build_fitness_table(restrictions_AND,0.027,[1.16,1.35],0.15,3)
    @test round(mitosis_probabilities[[0,0,0]],digits=10)==round(0.027,digits=10)
    @test round(mitosis_probabilities[[1,0,0]],digits=10)==round(0.027*1.16,digits=10)
    @test round(mitosis_probabilities[[1,1,1]],digits=10)==round(0.027*1.16*1.35*(1-0.15),digits=10)
    restrictions_OR= [[1,2,1],[1,3,1],[2,4,2],[3,4,2]]
    mitosis_probabilities = TumorSim.Fitness.build_fitness_table(restrictions_OR,0.027,[1.16,1.35],0.15,3)
    @test round(mitosis_probabilities[[0,0,0]],digits=10)==round(0.027,digits=10)
    @test round(mitosis_probabilities[[1,0,0]],digits=10)==round(0.027*1.16,digits=10)
    @test round(mitosis_probabilities[[1,0,1]],digits=10)==round(0.027*1.16*(1-0.15),digits=10)
end

@testset "Scenario tests" begin
    #1D
    scenario = create_scenario((40,),4,"random")
    scenario = create_scenario((40,),4,"center")
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
    scenario = create_scenario((15,20),5,"random")
    scenario = create_scenario((15,20),5,"center",true)
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
    scenario = create_scenario((16,17,18),19,"random")
    scenario = create_scenario((16,17,18),19,"center",false)
    @test scenario.x == 16
    @test scenario.y == 17
    @test scenario.z == 18
    @test length(scenario.cell_pos) == 19

    scenario = create_scenario((16,17,18),[(1,1,1),(2,2,2),(3,3,3)])
    @test scenario.x == 16
    @test scenario.y == 17
    @test scenario.z == 18
    @test length(scenario.cell_pos) == 3
    #Errors
    @test create_scenario((10,10),1000) === nothing
    @test create_scenario((10,),1000,"test") === nothing
    @test create_scenario((10,10),1000,"test") === nothing
    @test create_scenario((5,5,5),1000,"test") === nothing

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
    #2D init
    treatment = create_treatment(3000,0.65,0.5,3,0.75)
    scenario = create_scenario((10,10),5,"center",false)
    fitness=Dict([0,0,0]=>0.027, [1,0,0]=>0.035)
    model = TumorSim.model_init(interaction_rule=:contact,death_rate=0.55, mutation_rate=0.01, migration_rate = 0.05, scenario=scenario, fitness=fitness,treatment=treatment, seed=0)
    agent = collect(allagents(model))[1]
    @test nagents(model) == 5
    @test agent.time_alive == 0
    @test agent.genotype == BitArray([0,0,0])
    @test agent.phylogeny == []
    #3D
    scenario = create_scenario((10,10,10),5,"center",true)
    model = TumorSim.model_init(interaction_rule=:hierarchical_voter,death_rate=0.55, mutation_rate=0.01, migration_rate = 0.05, scenario=scenario, fitness=fitness,treatment=treatment, seed=0)
    agent = collect(allagents(model))[1]
    @test nagents(model) == 5
    @test agent.time_alive == 0
    @test agent.genotype == BitArray([0,0,0])
    @test agent.phylogeny == []
end
#Prepare a simulation
treatment = create_treatment(3000,0.65,0.5,3,0.75)
scenario = create_scenario((100,100,100),5,"center")
scenario2 = create_scenario((100,100,100),5,"center",true)
fitness=Dict([0,0,0]=>0.027, 
            [1,0,0]=>0.035,
            [0,1,0]=>0.032,
            [1,1,0]=>0.040,
            [1,1,1]=>0.032)

params = Dict(
    "death_rate" => [0.55,2],
    "interaction_rule"=>[:contact,:hierarchical_voter],
    "mutation_rate" => 0.1,
    "scenario" => [scenario,scenario2], 
    "fitness" => fitness,
    "treatment" => treatment,
    "migration_rate" => 0.05,
    "seed" => 0
)

dicts = dict_list(params)
steps=100
results_short = pmap(simulate,dicts,fill(steps,length(dicts)))
steps=3000
results = pmap(simulate,dicts,fill(steps,length(dicts)))
adata1 = results[1]["Genotypes"]
adata2 = results[2]["Genotypes"]
adata_short = results_short[1]["Genotypes"]
resistant_inhibited1 = results[1]["Resistant_inhibited_by"]

@testset "Simulate tests" begin
    @test length(results) == 8
    @test length(results[1]) == 24
    @test length(eachcol(adata1)) == 6
    @test length(eachrow(adata1)) != 0
    @test length(eachrow(adata1)) > length(eachrow(adata2))
end

@testset "Analysis tests" begin
    TTP1 = get_TTP(adata1,3200)
    TTP2 = get_TTP(adata2,3200)
    TTP_short = get_TTP(adata_short,100)

    @test TTP1 != -1
    @test TTP2 == -2

    diversity = TumorSim.Analysis.get_diversity(adata1)
    @test length(eachcol(diversity)) == 3
    @test length(eachrow(diversity)) != 0
    @test diversity[!,"species_richness"][1] == 1
    @test diversity[!,"shannon_index"][1] == 0
    @test diversity[!,"evenness"][1] == 0

    @test diversity[!,"species_richness"][400] > 1
    @test diversity[!,"shannon_index"][400] > 0
    @test diversity[!,"evenness"][400] > 0
    @test diversity[!,"evenness"][400] < 1

    #Test for resistance fraction
    resistant_fraction1 = TumorSim.Analysis.get_resistant_fraction_on_detection(adata1,3200,3)
    resistant_fraction2 = TumorSim.Analysis.get_resistant_fraction_on_detection(adata2,3200,3)
    @test resistant_fraction1 != -1
    @test resistant_fraction2 == -1
    @test resistant_fraction1 > 0
    @test resistant_fraction1 < 1
    
    #Test for divergence
    divergence = TumorSim.Analysis.get_divergence(adata1,resistant_inhibited1)
    @test length(eachcol(divergence)) == 4
    @test length(eachrow(divergence)) != 0
    @test minimum(divergence[!,"jensen_shannon"]) >= 0
    @test maximum(divergence[!,"jensen_shannon"]) <= 1

    #Test for jensen shannon divergence
    @test TumorSim.Analysis.jensen_shannon([1.0,0.0],[0.0,1.0]) == 1.0
    @test TumorSim.Analysis.jensen_shannon([1.0,0.0],[1.0,0.0]) == 0.0
    @test TumorSim.Analysis.jensen_shannon([0.5,0.5],[0.5,0.5]) == 0.0
    @test round(TumorSim.Analysis.jensen_shannon([0.5,0.5],[0.0,1.0]),digits=4) == 0.3113
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
