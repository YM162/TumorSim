using DrWatson
@quickactivate "TumorSim"

using DataFrames
using TumorSim
using Statistics

using HypothesisTests
using VegaLite
function plot_genotypes(adata::DataFrame,mode::String="absolute")
    #And lastly we can make plots of both the total number of cells of each genotype
    genotypes = names(adata)[2:end]
    stacked = stack(adata,genotypes)
    if mode=="absolute"
        stacked |>
        @vlplot(:area, x=:step, y={:value, stack=:zero}, color="variable:n")
    elseif mode=="relative"
        #And the relative number of cells of each genotype
        stacked |>
        @vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")
    else
        println("Mode not supported")
    end
end

df = collect_results(datadir("simulations"))
println(df)

println("Adaptive therapy:")
adaptive = filter("t_pausing_size" => n -> n == 1000, df)
adaptive_TTP = filter(n -> n !=-1,adaptive[!,"TTP"])
println("Mean ",mean(adaptive_TTP))
println("Continuous therapy:")
continuous = filter("t_pausing_size" => n -> n == 0, df)
continuous_TTP = filter(n -> n !=-1,adaptive[!,"TTP"])
println("Mean ",mean(continuous_TTP))

println("t-test p-value:")
pvalue(EqualVarianceTTest([x for x in adaptive_TTP], [x for x in continuous_TTP]))
