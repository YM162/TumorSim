#A generator that returns a list of functions that each get the number of cells of each genotype given a number of genes.
#Now fitness is birth rate for each genotype.

module Fitness
    export genotype_fraction_function_generator, bit_2_int, inhibited_by_function_generator
    
    using DataStructures
    function genotype_fraction_function_generator(fitness::Dict{Vector{Int64},Float64},)

        genotypes = [BitArray(x) for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
        genotypes_names = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
        func = function get_count(x::Base.Generator)
            c = counter(x)
            return (; zip([Symbol(x) for x in genotypes_names],[c[x] for x in genotypes])...)
        end
        return func

    end

    function inhibited_by_function_generator(fitness::Dict{Vector{Int64},Float64},)
        genotypes = [BitArray(x) for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
        genotypes_names = [replace(string(x)," "=>"") for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
        func = function get_count(x::Base.Generator)
            inhibited = collect(Iterators.flatten(x))
            c = counter(inhibited)
            return (; zip([Symbol(x) for x in genotypes_names],[c[x] for x in genotypes])...)
        end
        return func
    end


    #Function to go from BitArray to Int. Taken from https://discourse.julialang.org/t/parse-an-array-of-bits-bitarray-to-an-integer/42361/5
    function bit_2_int(arr::BitArray)
        arr = reverse(arr)
        sum(((i, x),) -> Int(x) << ((i-1) * sizeof(x)), enumerate(arr.chunks))
    end
end
using .Fitness