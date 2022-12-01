#We include the OncoSimulR interface
include("OncoSimulR.jl")

#A generator that returns a list of functions that each get the number of cells of each genotype given a number of genes.
function genotype_fraction_function_generator(fitness)
    functions = []
    for i in [x for x in sort!([x for x in keys(fitness)],by=x -> bit_2_int(BitArray(x)))]
        compare = i
        func = function get_perc(x)
            len = length(findall([string(y)[5:end] for y in x] .== string(compare)))
            return len
        end
        push!(functions,func)
    end
    return functions
end

#Function to go from BitArray to Int. Taken from https://discourse.julialang.org/t/parse-an-array-of-bits-bitarray-to-an-integer/42361/5
function bit_2_int(arr)
    arr = reverse(arr)
    sum(((i, x),) -> Int(x) << ((i-1) * sizeof(x)), enumerate(arr.chunks))
end

