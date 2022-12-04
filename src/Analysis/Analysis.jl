function get_TTP(adata::DataFrame,TTP_size::Int64)
    for row in eachrow(adata)
        tumor_size::Int = sum(row[2:end])
        if tumor_size>TTP_size
            return row["step"]
        end
    end
    return -1
end

function get_diversity(adata::DataFrame)
    diversity = DataFrame(species_richness=Int64[],shannon_index=Float64[],evenness=Float64[])
    for row in eachrow(adata)
        countmap(row[2:end])
        species_richness::Int64 = count([x!=0 for x in row[2:end]])
        relative_frequencies::Vector{Float64} = collect(row[2:end])./collect(sum(row[2:end]))
        shannon_index::Float64 = -sum([x*log(x) for x in relative_frequencies])
        evenness::Float64 = shannon_index/log(species_richness)
        if shannon_index === NaN
            shannon_index=0
            evenness=0
        end
        push!(diversity,[species_richness,shannon_index,evenness])
    end
    return(diversity)
end
