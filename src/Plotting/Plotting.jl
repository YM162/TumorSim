using VegaLite

function plot_genotypes(adata,mode="absolute")
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