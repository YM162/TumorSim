function get_TTP(adata,TTP_size)
    for i in eachrow(adata)
        tumor_size = sum(i[2:end])
        if tumor_size>TTP_size
            return i["step"]
        end
    end
    return -1
end

function get_diversity(adata)
    return 1
end