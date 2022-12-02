using DrWatson
@quickactivate "TumorSim"

using DataFrames
using TumorSim

df = collect_results(datadir("simulations"))
print(df)

#for (i,sim) in enumerate(eachrow(df))
#    simdict = Dict(names(sim) .=> values(sim))
#    plot_genotypes(sim["Genotypes"]) |> save(plotsdir(savename(string(i),simdict,"svg",ignores="path")))   
#end