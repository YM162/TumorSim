using DrWatson
@quickactivate "TumorSim"

using TumorSim

using BSON
using DataFrames
using Dates
using CSV
using StatsBase
using Plots

function get_clean_divergence(filename)
    df = BSON.load(datadir("simulations","competition_divergence",filename))["df"]
    newdf=DataFrame(step=Int64[])
    for i in df[!,"Divergence"]
        newdf = outerjoin(newdf,i[:,1:2],on=:step,makeunique=true)
    end
    sort!(newdf,:step)

    enddf = DataFrame(step=Int64[],jenshen_shannon_mean=Float64[],jenshen_shannon_sd=Float64[])
    for i in eachrow(newdf)
        push!(enddf,Dict(:step=>i[1],:jenshen_shannon_mean=>mean(skipmissing(i[2:end])),:jenshen_shannon_sd=>std(skipmissing(i[2:end]))))
    end

    #Remove every row where any of the columns has an undefined value
    filter(row -> all(x -> !(x isa Number && isnan(x)), row), enddf)
    finaldf = filter(row -> all(x -> !(x isa Number && isnan(x)), row), enddf)
    bson(datadir("simulations","competition_divergence","cleanup",filename),Dict("divergence" => finaldf))
end

function plot_divergence(filename)
    finaldf = BSON.load(datadir("simulations","competition_divergence","cleanup",filename))["divergence"]
    plot!(finaldf[!,"step"],finaldf[!,"jenshen_shannon_mean"],grid=false,ribbon=finaldf[!,"jenshen_shannon_sd"],fillalpha=.5,label=filename)
end

plot(xlabel = "Steps",ylabel="Jenshen Shannon",title="Deviation of real competition from perfect mixing model.")

#plot_divergence("3Genes_AND_Restrictions_28.3.2023.15.56.38.582.bson")
#plot_divergence("3Genes_OR_Restrictions_28.3.2023.16.46.47.769.bson")
#plot_divergence("3Genes_NO_Restrictions_28.3.2023.16.30.28.822.bson")

plot_divergence("7Genes_NO_Restrictions_30.3.2023.15.0.46.206.bson")
plot_divergence("7Genes_Example1_30.3.2023.15.52.29.549.bson")
plot_divergence("7Genes_Example2_30.3.2023.16.42.51.333.bson")
plot_divergence("7Genes_NK_Fitness_1_2.4.2023.16.20.40.75.bson")
plot_divergence("7Genes_NK_Fitness_2_2.4.2023.17.45.15.901.bson")




