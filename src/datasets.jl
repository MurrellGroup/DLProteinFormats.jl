using HuggingFaceApi
using ProteinChains: ProteinStructureStore
using Serialization: deserialize
using MergedArrays

export PSSDataset, SerializedDataset

abstract type AbstractProteinDataset end


struct PSSDataset <: AbstractProteinDataset
    hfurl::HuggingFaceURL
end

PSSDataset(args...) = PSSDataset(HuggingFaceURL(args...; repo_type="datasets"))

load(dataset::PSSDataset) = ProteinStructureStore(cached_download(dataset.hfurl), "r")


struct SerializedDataset <: AbstractProteinDataset
    hfurl::HuggingFaceURL
end

SerializedDataset(args...) = SerializedDataset(HuggingFaceURL(args...; repo_type="datasets"))

load(dataset::SerializedDataset) = deserialize(cached_download(dataset.hfurl))


export PDBStore, PDBStoreSubset, PDBSimpleFlat, PDBSimpleFlatSubset
export SwissProtStore, SwissProtStoreSubset, SwissProtSimpleFlat, SwissProtSimpleFlatSubset

const PDBStore = PSSDataset("MurrellLab/ProteinChains", "pdb.pss")
const PDBStoreSubset = PSSDataset("MurrellLab/ProteinChains", "pdb-500.pss")
const PDBSimpleFlat = SerializedDataset("MurrellLab/ProteinChains", "pdb-simple-flat.jls")
const PDBSimpleFlatSubset = SerializedDataset("MurrellLab/ProteinChains", "pdb-simple-flat-500.jls")

const SwissProtStore = PSSDataset("MurrellLab/ProteinChains", "swissprot-split.pss")
const SwissProtStoreSubset = PSSDataset("MurrellLab/ProteinChains", "swissprot-split-500.pss")
const SwissProtSimpleFlat = SerializedDataset("MurrellLab/ProteinChains", "swissprot-simple-flat.jls")
const SwissProtSimpleFlatSubset = SerializedDataset("MurrellLab/ProteinChains", "swissprot-simple-flat-500.jls")

remove_cache(dataset::AbstractProteinDataset) = HuggingFaceApi.remove_cache(dataset.hfurl)
