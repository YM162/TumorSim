using DrWatson
@quickactivate "TumorSim"

using TumorSim

using BSON
using DataFrames
using Dates
using CSV
using StatsBase
using Plots
using Plots.PlotMeasures
using LaTeXStrings

pgfplotsx()


function plot_divergence(finaldf,legendname,linestyle)
    plot!(finaldf[!,"step"],finaldf[!,"jensen_shannon_mean"],grid=false,ribbon=finaldf[!,"jensen_shannon_sd"],fillalpha=.5,label=legendname,linestyle=linestyle)
end

function load_scaled_data()
    _7Genes_NK1 = BSON.load(datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023","cleanup","7Genes_NK_Fitness_1_3.4.2023.6.35.13.179.bson"))
    _7Genes_NK1_df = _7Genes_NK1["divergence"]
    _7Genes_NK1_df[!,"step"] = (_7Genes_NK1_df[!,"step"] .- mean(_7Genes_NK1["detecting_time"])) ./ (mean(_7Genes_NK1["TTP"]) .- mean(_7Genes_NK1["detecting_time"]))

    _7Genes_NK2 = BSON.load(datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023","cleanup","7Genes_NK_Fitness_2_7.4.2023.21.29.47.666.bson"))
    _7Genes_NK2_df = _7Genes_NK2["divergence"]
    _7Genes_NK2_df[!,"step"] = (_7Genes_NK2_df[!,"step"] .- mean(_7Genes_NK2["detecting_time"])) ./ (mean(_7Genes_NK2["TTP"]) .- mean(_7Genes_NK2["detecting_time"]))

    _7Genes_NO = BSON.load(datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023","cleanup","7Genes_NO_Restrictions_12.4.2023.18.21.50.966.bson"))
    _7Genes_NO_df = _7Genes_NO["divergence"]
    _7Genes_NO_df[!,"step"] = (_7Genes_NO_df[!,"step"] .- mean(_7Genes_NO["detecting_time"])) ./ (mean(_7Genes_NO["TTP"]) .- mean(_7Genes_NO["detecting_time"]))

    return (_7Genes_NK1_df, _7Genes_NK2_df, _7Genes_NO_df)
end

_7Genes_NK1_df, _7Genes_NK2_df, _7Genes_NO_df = load_scaled_data()

function plot_scaled_figure(df_NK1,df_NK2,df_NO)
    f10 = Plots.font("Computer Modern", 10)
    f12 = Plots.font("Computer Modern", 12)
    plot(topmargin = 5mm, 
        bottommargin = 5mm,
        leftmargin = 5mm,
        rightmargin = 5mm,
        legend = :topright,
        legend_font_halign= :left,
        size = (800,500), 
        ylims = (-0.05,0.7),
        xtickfont=f10, ytickfont=f10, legendfont=f10, guidefont=f10, titlefont=f12,
        xticks = ([0,1],[L"\mathrm{\ Treatment\ start}",L"\mathrm{\ Tumor\ progression}",]),
        xlabel = L"\mathbf{Time\ (Arbitrary\ units)}",
        ylabel=L"\mathbf{Jensen-Shannon\ divergence}",
        title=L"\mathbf{Deviation\ of\ real\ competition\ from\ perfect\ mixing\ model}", 
        tex_output_standalone = true
    )

    plot_divergence(df_NK1, L"\mathrm{7\ Genes\ w/\ NK\ Fitness\ Landscape\ 1}",:solid)
    plot_divergence(df_NK2, L"\mathrm{7\ Genes\ w/\ NK\ Fitness\ Landscape\ 2}",:dash)
    plot_divergence(df_NO, L"\mathrm{7\ Genes\ w/\ NO\ Restrictions}",:dashdot)

    plot!([0], seriestype = :vline, color = :green, linewidth = 2, label=nothing,ribbon=0.5)
    plot!([1], seriestype = :vline, color = :red, linewidth = 2, label=nothing)
    
end

plot_scaled_figure(_7Genes_NK1_df,_7Genes_NK2_df,_7Genes_NO_df)

savefig(datadir("figures","AT_in_cancer_the_role_of_restrictions_2023","7Genes_NK_Divergence_Scaled.tex"))