#A generator that returns a list of functions that each get the number of cells of each genotype given a number of genes.
#Now fitness is birth rate for each genotype.

module Fitness
    export genotype_fraction_function_generator, bit_2_int, inhibited_by_function_generator, build_fitness_table
    
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

    function does_it_follow_rules(Genotype,DAG)
    
        for i in 1:length(Genotype)
            if Genotype[i] == 0
                continue
            else
                #Gene is mutated, but is it valid?
                #If we find anything wrong, we return false. If we dont, we just continue
                rules = DAG[i+1,2:end]
                if all(rules .== 0) && DAG[i+1:1] == 0
    
                    return false
                end
                #We check for simple dependencies.
                for j in 1:length(rules)
                    if rules[j] == 1
                        if Genotype[j] == 0
    
                            return false
                        end
                    end
                end
                if any(rules .== 2) && !(any(Genotype[findall(rules .== 2)].==1))
    
                    return false
                end
            end
        end
        return true
    end

    function build_fitness_table(restrictions,base_pr,mult_pr,cost_of_resistance,ngenes)
    
        DAG = zeros(Int,ngenes+1,ngenes+1)
        for (i,j,value) in restrictions
            DAG[j,i] = value
        end
    
        Genotypes = [digits(x, base=2, pad=ngenes) for x in 0:2^ngenes-1]
    
        fitness::Dict{Vector{Int64},Float64} = Dict()
        for genotype in Genotypes
            if does_it_follow_rules(genotype,DAG)
                fitness[genotype] = base_pr*prod(mult_pr[findall(genotype[1:end-1].==1)])*(1-cost_of_resistance)^(genotype[end]==1)
            end
        end
        return fitness
    end
end
using .Fitness