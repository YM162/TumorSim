using DrWatson
@quickactivate "TumorSim"

using DataFrames
using TumorSim
using Statistics

using HypothesisTests

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
