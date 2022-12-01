using DrWatson, Test
@quickactivate "TumorSim"

using Distributed


# Run test suite
println("Starting tests")
ti = time()

#Fitness tests
include(srcdir("Fitness/Fitness.jl"))

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

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
