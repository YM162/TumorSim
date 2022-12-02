using StatsBase
using DataFrames
function get_TTP(adata,TTP_size)
    for row in eachrow(adata)
        tumor_size = sum(row[2:end])
        if tumor_size>TTP_size
            return row["step"]
        end
    end
    return -1
end

function get_diversity(adata)
    diversity = DataFrame(species_richness=Float64[],shannon_index=Float64[],evenness=Float64[])
    for row in eachrow(adata)
        countmap(row[2:end])
        species_richness = count([x!=0 for x in row[2:end]])
        relative_frequencies = collect(row[2:end])./collect(sum(row[2:end]))
        shannon_index = -sum([x*log(x) for x in relative_frequencies])
        evenness = shannon_index/log(species_richness)
        if shannon_index === NaN
            shannon_index=0
            evenness=0
        end
        push!(diversity,[species_richness,shannon_index,evenness])
    end
    return(diversity)
end
