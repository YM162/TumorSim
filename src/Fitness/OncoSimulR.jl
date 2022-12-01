using RCall

#Function to create a random fitness landscape using the OncoSimulR library.
function OncoSimulR_rfitness(;g,c,sd)
    R"library(OncoSimulR)"
    fitness = R"rfitness(g=$g ,c=$c ,sd=$sd )"
    rows=2^g
    dictionary=Dict()
    for i in 1:rows
        genotype::Vector{Int64}=[]
        for j in 1:g
            push!(genotype,Int(fitness[i,j]))
        end
        push!(dictionary,(genotype=>Real(fitness[i,g+1])))
    end
    return dictionary
end