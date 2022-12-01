using DrWatson
@quickactivate "TumorSim"

using DataFrames
using TumorSim: Scenario, Treatment
df = collect_results(datadir("simulations"))
print(df)