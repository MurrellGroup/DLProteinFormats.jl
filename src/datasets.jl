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


using CSV, DataFrames

struct CSVDataset <: AbstractProteinDataset
    hfurl::HuggingFaceURL
end

CSVDataset(args...) = CSVDataset(HuggingFaceURL(args...; repo_type="datasets"))

load(dataset::CSVDataset, sink=DataFrame) = CSV.read(cached_download(dataset.hfurl), sink)


using JLD2

struct JLD2Dataset <: AbstractProteinDataset
    hfurl::HuggingFaceURL
end

JLD2Dataset(args...) = JLD2Dataset(HuggingFaceURL(args...; repo_type="datasets"))

load(dataset::JLD2Dataset) = JLD2.load(cached_download(dataset.hfurl))["data"]


export PDBStore, PDBStore500, PDBSimpleFlat, PDBSimpleFlat500, PDBSimpleFlatV2, PDBSimpleFlatV2_500
export PDBAtom14, PDBAtom14_500
export PDBClusters, PDBMethods, PDBTable
export PDBFlatom169K
export SwissProtSplitStore, SwissProtSplitStore500, SwissProtSplitSimpleFlat, SwissProtSplitSimpleFlat500
export QM9Crude1000

const PDBStore = PSSDataset("MurrellLab/ProteinChains", "pdb.pss")
const PDBStore500 = PSSDataset("MurrellLab/ProteinChains", "pdb-500.pss")
const PDBSimpleFlat = SerializedDataset("MurrellLab/ProteinChains", "pdb-simple-flat.jls")
const PDBSimpleFlat500 = SerializedDataset("MurrellLab/ProteinChains", "pdb-simple-flat-500.jls")
const PDBSimpleFlatV2 = JLD2Dataset("MurrellLab/ProteinChains", "flat-v2.jld2")
const PDBSimpleFlatV2_500 = JLD2Dataset("MurrellLab/ProteinChains", "flat-v2-500.jld2")

const PDBAtom14 = JLD2Dataset("MurrellLab/ProteinChains", "pdb-atom14.jld2")
const PDBAtom14_500 = JLD2Dataset("MurrellLab/ProteinChains", "pdb-atom14-500.jld2")

const PDBClusters = SerializedDataset("MurrellLab/ProteinChains", "pdb-clusters.jls")
const PDBMethods = CSVDataset("MurrellLab/ProteinChains", "pdb-methods.csv")
const PDBTable = CSVDataset("MurrellLab/ProteinChains", "pdb-table.csv")

using StaticArrays, StaticStrings
const PDBFlatom169K = SerializedDataset("MurrellLab/ProteinChains", "flatom-169k.jls")

const SwissProtSplitStore = PSSDataset("MurrellLab/ProteinChains", "swissprot-split.pss")
const SwissProtSplitStore500 = PSSDataset("MurrellLab/ProteinChains", "swissprot-split-500.pss")
const SwissProtSplitSimpleFlat = SerializedDataset("MurrellLab/ProteinChains", "swissprot-split-simple-flat.jls")
const SwissProtSplitSimpleFlat500 = SerializedDataset("MurrellLab/ProteinChains", "swissprot-split-simple-flat-500.jls")

const QM9Crude1000 = SerializedDataset("MurrellLab/ProteinChains", "qm9-crude-1000.jls")

remove_cache(dataset::AbstractProteinDataset) = HuggingFaceApi.remove_cache(dataset.hfurl; now=true)
remove_cache() = HuggingFaceApi.remove_cache("MurrellLab/ProteinChains")
