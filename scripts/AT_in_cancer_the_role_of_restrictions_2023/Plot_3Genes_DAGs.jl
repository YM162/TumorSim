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
    _3Genes_AND = BSON.load(datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023","cleanup","3Genes_AND_Restrictions_3.4.2023.2.49.14.497.bson"))
    _3Genes_AND_df = _3Genes_AND["divergence"]
    _3Genes_AND_df[!,"step"] = (_3Genes_AND_df[!,"step"] .- mean(_3Genes_AND["detecting_time"])) ./ (mean(_3Genes_AND["TTP"]) .- mean(_3Genes_AND["detecting_time"]))

    _3Genes_OR = BSON.load(datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023","cleanup","3Genes_OR_Restrictions_3.4.2023.3.15.28.889.bson"))
    _3Genes_OR_df = _3Genes_OR["divergence"]
    _3Genes_OR_df[!,"step"] = (_3Genes_OR_df[!,"step"] .- mean(_3Genes_OR["detecting_time"])) ./ ( mean(_3Genes_OR["TTP"]) .- mean(_3Genes_OR["detecting_time"]))

    _3Genes_NO = BSON.load(datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023","cleanup","3Genes_NO_Restrictions_3.4.2023.2.27.15.692.bson"))
    _3Genes_NO_df = _3Genes_NO["divergence"]
    _3Genes_NO_df[!,"step"] = (_3Genes_NO_df[!,"step"] .- mean(_3Genes_NO["detecting_time"])) ./ (mean(_3Genes_NO["TTP"]) .- mean(_3Genes_NO["detecting_time"]))

    return (_3Genes_AND_df, _3Genes_OR_df, _3Genes_NO_df)
end

#_3Genes_AND_df, _3Genes_OR_df, _3Genes_NO_df = load_scaled_data()

function plot_scaled_figure(df_AND,df_OR,df_NO)
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

    plot_divergence(df_AND, L"\mathrm{3\ Genes\ w/\ AND\ Restrictions}",:solid)
    plot_divergence(df_OR, L"\mathrm{3\ Genes\ w/\ OR\ Restrictions}",:dash)
    plot_divergence(df_NO, L"\mathrm{3\ Genes\ w/\ NO\ Restrictions}",:dashdot)

    plot!([0], seriestype = :vline, color = :green, linewidth = 2, label=nothing,ribbon=0.5)
    plot!([1], seriestype = :vline, color = :red, linewidth = 2, label=nothing)
    
end

plot_scaled_figure(_3Genes_AND_df,_3Genes_OR_df,_3Genes_NO_df)

savefig(datadir("figures","AT_in_cancer_the_role_of_restrictions_2023","3Genes_Divergence_Scaled.tex"))



function print_all_data()

    filename = "7Genes_NO_Restrictions_12.4.2023.18.21.50.966.bson"
    thisdf = BSON.load(datadir("simulations","AT_in_cancer_the_role_of_restrictions_2023","cleanup",filename))
    println("Simulation ",filename," has:")
    println("TTP: ",mean(thisdf["TTP"]))
    println("Detecting time: ",mean(thisdf["detecting_time"]))

end
