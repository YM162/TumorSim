module Analysis
    export get_TTP, get_diversity, get_resistant_fraction_on_detection

    using DataFrames
    using StatsBase

    function get_resistant_fraction_on_detection(adata::DataFrame, detection_size::Int64, resistance_gene::Int64)
        for row in eachrow(adata)
            tumor_size::Int = sum(row[2:end])
            if tumor_size>detection_size
                resistant = 0
                for i in 2:length(row)
                    if eval(Meta.parse(names(row)[i]))[resistance_gene] == 1
                        resistant += row[i]
                    end
                end
                return resistant/tumor_size
            end
        end
        return -1
    end

    function get_TTP(adata::DataFrame,TTP_size::Int64)
        for row in eachrow(adata)
            tumor_size::Int = sum(row[2:end])
            if tumor_size>TTP_size
                return row["step"]
            end
        end
        if sum(eachrow(adata)[end][2:end]) == 0
            return -2 #Tumor is gone
        else
            return -1 #Tumor is growing and max steps is reached
        end
    end

    function get_diversity(adata::DataFrame)
        diversity = DataFrame(species_richness=Int64[],shannon_index=Float64[],evenness=Float64[])
        for row in eachrow(adata)
            countmap(row[2:end])
            species_richness::Int64 = count([x!=0 for x in row[2:end]])
            relative_frequencies::Vector{Float64} = collect(row[2:end])./collect(sum(row[2:end]))
            
            filter!((x) -> x != 0, relative_frequencies)

            shannon_index::Float64 = -sum([x*log(x) for x in relative_frequencies])
            evenness::Float64 = shannon_index/log(species_richness)
            if shannon_index === NaN || shannon_index == -0.0
                shannon_index=0
                evenness=0
            end
            push!(diversity,[species_richness,shannon_index,evenness])
        end
        return(diversity)
    end
end

using .Analysis