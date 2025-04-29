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


export PDBStore, PDBStore500, PDBSimpleFlat, PDBSimpleFlat500
export SwissProtSplitStore, SwissProtSplitStore500, SwissProtSplitSimpleFlat, SwissProtSplitSimpleFlat500

const PDBStore = PSSDataset("MurrellLab/ProteinChains", "pdb.pss")
const PDBStore500 = PSSDataset("MurrellLab/ProteinChains", "pdb-500.pss")
const PDBSimpleFlat = SerializedDataset("MurrellLab/ProteinChains", "pdb-simple-flat.jls")
const PDBSimpleFlat500 = SerializedDataset("MurrellLab/ProteinChains", "pdb-simple-flat-500.jls")

const SwissProtSplitStore = PSSDataset("MurrellLab/ProteinChains", "swissprot-split.pss")
const SwissProtSplitStore500 = PSSDataset("MurrellLab/ProteinChains", "swissprot-split-500.pss")
const SwissProtSplitSimpleFlat = SerializedDataset("MurrellLab/ProteinChains", "swissprot-split-simple-flat.jls")
const SwissProtSplitSimpleFlat500 = SerializedDataset("MurrellLab/ProteinChains", "swissprot-split-simple-flat-500.jls")

remove_cache(dataset::AbstractProteinDataset) = HuggingFaceApi.remove_cache(dataset.hfurl; now=true)
remove_cache() = HuggingFaceApi.remove_cache("MurrellLab/ProteinChains")
