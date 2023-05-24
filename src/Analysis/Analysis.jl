module Analysis
    export get_TTP, get_diversity, get_resistant_fraction_on_detection, get_divergence

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

    function kl_divergence(p::Vector{Float64}, q::Vector{Float64})
        non_zero_p = findall(x->x!=0,p)
        return sum([p[i]*log2(p[i]/q[i]) for i in non_zero_p])
    end

    function jensen_shannon(p::Vector{Float64}, q::Vector{Float64})
        m = 0.5*(p+q)
        return 0.5*(kl_divergence(p,m)+kl_divergence(q,m))
    end

    function get_divergence(adata::DataFrame,resistant_inhibited_by::DataFrame)
        divergence = DataFrame(step=Int64[],jensen_shannon=Float64[],Kullback_Leibler_relfreq_relinhib=Float64[],Kullback_Leibler_relinhib_relfreq=Float64[])
        adata_rows = eachrow(adata)
        resistant_inhibited_by_rows = eachrow(resistant_inhibited_by)
        for (i,row) in enumerate(eachrow(adata))

            if sum(adata_rows[i][2:end]) == 0 || sum(resistant_inhibited_by_rows[i][2:end]) == 0
                continue
            end

            relative_frequencies::Vector{Float64} = collect(adata_rows[i][2:end])./collect(sum(adata_rows[i][2:end]))
            relative_inhibition::Vector{Float64} = collect(resistant_inhibited_by_rows[i][2:end])./collect(sum(resistant_inhibited_by_rows[i][2:end]))

            push!(divergence,[row[1],jensen_shannon(relative_frequencies,relative_inhibition),kl_divergence(relative_frequencies,relative_inhibition),kl_divergence(relative_inhibition,relative_frequencies)])
        end
        return divergence
    end
end

using .Analysis