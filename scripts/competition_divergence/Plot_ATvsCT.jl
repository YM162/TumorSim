using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "TumorSim"

@everywhere using TumorSim

using ProgressMeter
using BSON
using DataFrames
using Dates

#We test adaptive and continuous therapy



fitness = Dict([0,0,0]=>0.027,
                [1,0,0]=>0.030,
                [0,1,0]=>0.033,
                [1,1,0]=>0.046,
                [1,1,1]=>0.025)

adaptive_therapy = create_treatment(10000, 1, 0.5, 3, 0.75) 
continuous_therapy = create_treatment(10000, 1, 0.0, 3, 0.75) 

parameters = Dict(
    "death_rate" => [0.3],
    "mutation_rate" => 0.005,
    "scenario" => [create_scenario((200,200),100,"center",false)], 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy,continuous_therapy],
    "migration_rate" => [0.1],
    "interaction_rule" => [:hierarchical_voter],
    "seed" => map(abs,rand(Int64,1))
)

parameter_combinations = dict_list(parameters)
println("Number of simulations: ",length(parameter_combinations))
steps=5000

filename = "3Genes_AND_Restrictions_"*Dates.format(now(),"d.m.yyyy.H.M.S.s")
println("Starting simulations...")

p = Progress(length(parameter_combinations), barglyphs=BarGlyphs("[=> ]"),output=stdout,barlen=50)
results = progress_pmap(simulate,parameter_combinations,fill(steps,length(parameter_combinations)),progress=p)

println("Saving simulations...")
df = DataFrame(results)
filepath = datadir("simulations","competition_divergence",filename*".bson")

#bson(filepath,Dict("df" => df))
using StatsBase


println(df[!,"TTP"])
println(df[!,"TTP"][1]-df[!,"TTP"][2])
println(df[!,"Resistant_on_detection"])


n=1
using VegaLite
#stack(df[!,"Genotypes"][n],names(df[!,"Genotypes"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")
#stack(df[!,"Genotypes"][n],names(df[!,"Genotypes"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:zero}, color="variable:n")

#stack(df[!,"Resistant_inhibited_by"][n],names(df[!,"Resistant_inhibited_by"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")

#df[!,"Divergence"][n] |> @vlplot(:line, x=:step, y=:jenshen_shannon, color=:variable)

using Plots

using Plots.PlotMeasures
using LaTeXStrings

pgfplotsx()

f10 = Plots.font("Computer Modern", 10)
f12 = Plots.font("Computer Modern", 12)

TTP_AT = df[!,"TTP"][1]
TTP_CT = df[!,"TTP"][2]
Tstart = [sim[!,"step"][findfirst(sim[!,"status"])] for sim in df[!,"Treatment_status"]][1]

plot(topmargin = 5mm, 
        bottommargin = 5mm,
        leftmargin = 5mm,
        rightmargin = 5mm,
        legend = :topright,
        legend_font_halign= :left,
        size = (800,500), 
        ylims = (0,15000),
        xtickfont=f10, ytickfont=f10, legendfont=f10, guidefont=f10, titlefont=f12,
        xticks = ([Tstart,TTP_AT,TTP_CT],[L"\mathrm{\ Treatment\ start}",L"\mathrm{TTP\\ (AT)}",L"\mathrm{TTP\\ (CT)}"]),
        yticks = ([5000,10000,12000],[L"\mathrm{\tfrac{Acceptable\ burden}{2}}",L"\mathrm{Acceptable\ burden}",L"\mathrm{Tumor\ progression}"]),
        xlabel = L"\mathbf{Time\ (Arbitrary\ units)}",
        ylabel=L"\mathbf{Tumor\ size}",
        tex_output_standalone = true
    )

    plot!(df[!,"Genotypes"][1][!,"step"],[sum(x) for x in eachrow(df[!,"Genotypes"][1][:,2:end])],
    grid=false,
    linewidth = 1.5,
    color = :coral2,
    label=L"\mathrm{Adaptive\ Therapy\ (AT)}",
    linestyle=:solid)

    plot!(df[!,"Genotypes"][2][!,"step"],[sum(x) for x in eachrow(df[!,"Genotypes"][2][:,2:end])],
    grid=false,
    linewidth = 1.5,
    color = :royalblue1,
    label=L"\mathrm{Continuous\ Therapy\ (CT)}",
    linestyle=:dash)

    plot!([5000], seriestype = :hline, color = :gray, linewidth = 2, linestyle=:dash, label=nothing)
    plot!([10000], seriestype = :hline, color = :gray, linewidth = 2, linestyle=:dash, label=nothing)
    plot!([12000], seriestype = :hline, color = :gray, linewidth = 2, linestyle=:dash, label=nothing)

    plot!([Tstart], seriestype = :vline, color = :green, linewidth = 2, label=nothing)
    plot!([TTP_AT], seriestype = :vline, color = :red, linewidth = 2, label=nothing)
    plot!([TTP_CT], seriestype = :vline, color = :skyblue, linewidth = 2, label=nothing)

#[:darkolivegreen4 :skyblue :coral2 :darkorange :royalblue1]

#savefig(datadir("figures","competition_divergence","Plot_Explanation.tex"))
