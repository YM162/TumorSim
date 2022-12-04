#We define what a treatment is
mutable struct Treatment
    detecting_size::Int64
    starting_size::Int64
    pausing_size::Int64
    resistance_gene::Int64
    kill_rate::Float64
    active::Bool
    detected::Bool
end

function create_treatment(detecting_size,starting_size,pausing_size,resistance_gene,kill_rate)
    return Treatment(detecting_size,starting_size,pausing_size,resistance_gene,kill_rate,false,false)
end