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
                [1,0,0]=>0.032,
                [0,1,0]=>0.037,
                [1,1,0]=>0.042,
                [1,1,1]=>0.036)

adaptive_therapy = create_treatment(3000, 1, 0.5, 3, 0.75) 
continuous_therapy = create_treatment(3000, 1, 0.0, 3, 0.75) 

parameters = Dict(
    "death_rate" => [0.3],
    "mutation_rate" => 0.01,
    "scenario" => [create_scenario((100,100),100,"center",false)], 
    "fitness" => fitness,
    "treatment" => [adaptive_therapy],
    "migration_rate" => [0.0],
    "interaction_rule" => [:contact],
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
filepath = datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023",filename*".bson")

#bson(filepath,Dict("df" => df))
using StatsBase


println(df[!,"TTP"])
println(df[!,"Resistant_on_detection"])


n=1
using VegaLite
#stack(df[!,"Genotypes"][n],names(df[!,"Genotypes"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")
#stack(df[!,"Genotypes"][n],names(df[!,"Genotypes"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:zero}, color="variable:n")

#stack(df[!,"Resistant_inhibited_by"][n],names(df[!,"Resistant_inhibited_by"][n])[2:end]) |> @vlplot(:area, x=:step, y={:value, stack=:normalize}, color="variable:n")

#df[!,"Divergence"][n] |> @vlplot(:line, x=:step, y=:jensen_shannon, color=:variable)

using Plots

using Plots.PlotMeasures
using LaTeXStrings

pgfplotsx()

f10 = Plots.font("Computer Modern", 10)
f12 = Plots.font("Computer Modern", 12)

TTP = df[!,"TTP"][1]
Tstart = [sim[!,"step"][findfirst(sim[!,"status"])] for sim in df[!,"Treatment_status"]][1]

p1 = areaplot(df[!,"Genotypes"][1][!,"step"], [df[!,"Genotypes"][1][!,"[1,1,1]"] df[!,"Genotypes"][1][!,"[1,1,0]"] df[!,"Genotypes"][1][!,"[0,1,0]"] df[!,"Genotypes"][1][!,"[1,0,0]"] df[!,"Genotypes"][1][!,"[0,0,0]"]], labels=reshape(reverse(names(df[!,"Genotypes"][1])[2:end]), (1,5)), seriescolor = [:darkolivegreen4 :skyblue :coral2 :darkorange :royalblue1],
        topmargin = 0mm, 
        bottommargin = 0mm,
        leftmargin = 5mm,
        rightmargin = 0mm,
        legend_font_halign= :left,
        ylims = (0,5000),
        xlims = (0,df[!,"Genotypes"][1][!,"step"][end]),
        xtickfont=f10, ytickfont=f10, legendfont=f10, guidefont=f10, titlefont=f12,
        xlabel = nothing,
        xticks = nothing,
        ylabel=L"\parbox{5em}{\textbf{Cell\ count\\ (Absolute)}}")

p2 = areaplot(df[!,"Genotypes"][1][!,"step"], [df[!,"Genotypes"][1][!,"[1,1,1]"]./[sum(x) for x in eachrow(df[!,"Genotypes"][1][:,2:end])] df[!,"Genotypes"][1][!,"[1,1,0]"]./[sum(x) for x in eachrow(df[!,"Genotypes"][1][:,2:end])] df[!,"Genotypes"][1][!,"[0,1,0]"]./[sum(x) for x in eachrow(df[!,"Genotypes"][1][:,2:end])] df[!,"Genotypes"][1][!,"[1,0,0]"]./[sum(x) for x in eachrow(df[!,"Genotypes"][1][:,2:end])] df[!,"Genotypes"][1][!,"[0,0,0]"]./[sum(x) for x in eachrow(df[!,"Genotypes"][1][:,2:end])]], seriescolor = [:darkolivegreen4 :skyblue :coral2 :darkorange :royalblue1],
        labels = nothing,        
        topmargin = 0mm, 
        bottommargin = 0mm,
        leftmargin = 5mm,
        rightmargin = 0mm,
        legend_font_halign= :left,
        ylims = (0,1),
        xlims = (0,df[!,"Genotypes"][1][!,"step"][end]),
        xtickfont=f10, ytickfont=f10, legendfont=f10, guidefont=f10, titlefont=f12,
        xlabel = nothing,
        xticks = nothing,
        ylabel=L"\parbox{5em}{\textbf{Cell\ count\\ (Relative)}}")

p3 = areaplot(df[!,"Resistant_inhibited_by"][1][!,"step"], [df[!,"Resistant_inhibited_by"][1][!,"[1,1,1]"]./[sum(x) for x in eachrow(df[!,"Resistant_inhibited_by"][1][:,2:end])] df[!,"Resistant_inhibited_by"][1][!,"[1,1,0]"]./[sum(x) for x in eachrow(df[!,"Resistant_inhibited_by"][1][:,2:end])] df[!,"Resistant_inhibited_by"][1][!,"[0,1,0]"]./[sum(x) for x in eachrow(df[!,"Resistant_inhibited_by"][1][:,2:end])] df[!,"Resistant_inhibited_by"][1][!,"[1,0,0]"]./[sum(x) for x in eachrow(df[!,"Resistant_inhibited_by"][1][:,2:end])] df[!,"Resistant_inhibited_by"][1][!,"[0,0,0]"]./[sum(x) for x in eachrow(df[!,"Resistant_inhibited_by"][1][:,2:end])]], seriescolor = [:darkolivegreen4 :skyblue :coral2 :darkorange :royalblue1],
        labels = nothing,        
        topmargin = 0mm, 
        bottommargin = 0mm,
        leftmargin = 5mm,
        rightmargin = 0mm,
        legend_font_halign= :left,
        ylims = (0,1),
        xlims = (0,df[!,"Genotypes"][1][!,"step"][end]),
        xtickfont=f10, ytickfont=f10, legendfont=f10, guidefont=f10, titlefont=f12,
        xticks = nothing,
        ylabel=L"\parbox{5em}{\textbf{Relative\ resistant\\cell\ inhibition}}")

p4 = plot(df[!,"Divergence"][1][!,"step"], df[!,"Divergence"][1][!,"jensen_shannon"], seriescolor = [:royalblue1],
        grid=false,
        label = nothing, 
        linewidth = 2,       
        topmargin = 0mm, 
        bottommargin = 5mm,
        leftmargin = 5mm,
        rightmargin = 0mm,
        legend_font_halign= :left,
        ylims = (-0.05,1),
        xlims = (0,df[!,"Genotypes"][1][!,"step"][end]),
        xtickfont=f10, ytickfont=f10, legendfont=f10, guidefont=f10, titlefont=f12,
        xlabel = L"\mathbf{Time\ (Arbitrary\ units)}",
        xticks = ([Tstart,TTP],[L"\mathrm{\ Treatment\ start}",L"\mathrm{\ Tumor\ progression}",]),
        ylabel=L"\parbox{5em}{\textbf{Jensen-Shannon\ divergence}}")

p = plot(p1, p2, p3, p4, layout=(4,1),bottom_margin = -12mm, link = :y)



vline!(p[1], [Tstart], color = :green, linewidth = 2, label=nothing)
vline!(p[1], [TTP], color = :red, linewidth = 2, label=nothing)

vline!(p[2], [Tstart], color = :green, linewidth = 2, label=nothing)
vline!(p[2], [TTP], color = :red, linewidth = 2, label=nothing)

vline!(p[3], [Tstart], color = :green, linewidth = 2, label=nothing)
vline!(p[3], [TTP], color = :red, linewidth = 2, label=nothing)

vline!(p[4], [Tstart], color = :green, linewidth = 2, label=nothing)
vline!(p[4], [TTP], color = :red, linewidth = 2, label=nothing,
tex_output_standalone = true)

savefig(datadir("figures","AT_in_cancer_the_role_of_restrictions_2023","Plot_Explanation.tex"))
