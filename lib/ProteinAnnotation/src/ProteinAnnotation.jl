module ProteinAnnotation

using CSV, DataFrames
using Downloads: download
using TaxMax

include("Antibody.jl")
include("SCOP.jl")
include("import_pss.jl")

export merge_annotations

function merge_annotations(dataframes::DataFrame...)
    outerjoin(dataframes...; on = [:pdb_id, :chain_id])
end

include("shared.jl")
include("V2.jl")

end
