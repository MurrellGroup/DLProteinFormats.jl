using Statistics: mean
import DLProteinFormats

struct AnnotatedV2 end

const host_class_cats = Categorizer(Dict("Gammaproteobacteria" => 1, "Mammalia" => 2, "Insecta" => 3), 3)
const gene_class_cats = Categorizer(Dict("Gammaproteobacteria" => 1, "Mammalia" => 2, "Insecta" => 3), 3)
const gene_superkingdom_cats = Categorizer(Dict("Eukaryota" => 1, "Bacteria" => 2, "Viruses" => 3, "Archaea" => 4), 4)
const therm_cats = Categorizer([50.0, 65.0, 75.0, 85.0])
const ss_cats = Categorizer([0.2, 0.4, 0.6, 0.8])
const scop_type_cats = Categorizer(
    Dict("globular" => 1, "membrane" => 2, "unstructured" => 3, "fibrous" => 4), 5)

function features(::AnnotatedV2, row)
    ss = onehotcats([row.helix_count, row.sheet_count] ./ row.length, ss_cats)
    therm = onehotcats(row.avg_temp, therm_cats)
    gene_class = onehotcats(row.gene_class, gene_class_cats)
    host_class = onehotcats(row.host_class, host_class_cats)
    gene_superkingdom = onehotcats(row.gene_superkingdom, gene_superkingdom_cats)
    (; ss, therm, gene_class, host_class, gene_superkingdom)
end

flat_features(::AnnotatedV2, features) = mapreduce(vec, vcat, features)

function index_table(df)
    T = NamedTuple{Tuple(propertynames(df)), Tuple{eltype.(eachcol(df))...}}
    dict = Dict{String,Vector{T}}()
    for (i, pdb_id) in enumerate(df.pdb_id)
        push!(get!(dict, pdb_id, T[]), NamedTuple(df[i,:]))
    end
    return dict
end

function DLProteinFormats.flatten(pt::AnnotatedV2, structure::ProteinStructure; T=Float32)
    len = sum(length, structure)
    chainids = reduce(vcat, [fill(i, l) for (i, l) in enumerate(length.(structure))], init=Int[])
    locs = zeros(T, 3, 1, len)
    rots = zeros(T, 3, 3, len)
    resinds = zeros(Int, len)
    AAs = zeros(Int, len)
    pos = 1
    chain_labels = String[]
    for (i, chain) in enumerate(structure)
        push!(chain_labels, chain.id)
        f = ProteinChains.Frames(chain)
        locs[:, 1, pos:pos+length(chain)-1] .= f.translations
        rots[:, :, pos:pos+length(chain)-1] .= f.rotations
        resinds[pos:pos+length(chain)-1] .= chain.numbering
        AAs[pos:pos+length(chain)-1] .= DLProteinFormats.aa_to_ints(chain.sequence)
        pos += length(chain)
    end
    @assert pos-1 == len == length(AAs)
    return (; record_type = :AnnotatedV2, structure.name, chain_labels, len,
        locs = (locs .- mean(locs, dims = 3)) ./ DLProteinFormats.unit_scaling, rots, AAs,
        resinds, chainids)
end
