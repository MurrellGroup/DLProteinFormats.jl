module Flatom

using BioStructures

include("utils.jl")
include("renumber.jl")
include("flatoms.jl")

export flatoms
export MMCIFDict
export PDBFormat, MMCIFFormat

function flatoms_pdb(pdb_id::AbstractString, ba_number::Int=0)
    mktempdir() do dir
        structure = retrievepdb(pdb_id; dir, ba_number)
        return flatoms(structure)
    end
end

macro flatoms_str(string::String, ba_number::Int=0)
    flatoms_pdb(string, ba_number)
end

export flatoms_pdb, @flatoms_str

end