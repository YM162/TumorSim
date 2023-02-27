#We define what a treatment is
module Treatment
export TreatmentObject, create_treatment
    mutable struct TreatmentObject
        detecting_size::Int64
        starting_size::Float64
        pausing_size::Float64
        resistance_gene::Int64
        kill_rate::Float64
        active::Bool
        detected::Bool
    end

    function create_treatment(detecting_size,starting_size,pausing_size,resistance_gene,kill_rate)
        return TreatmentObject(detecting_size,starting_size,pausing_size,resistance_gene,kill_rate,false,false)
    end
end
using .Treatment