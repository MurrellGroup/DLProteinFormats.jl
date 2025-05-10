struct AnnotatedV1 end
AnnotatedV1NT = NamedTuple{<:Any,<:Tuple{AnnotatedV1,Vararg{Any}}}
struct BatchedAnnotatedV1 end
BatchedAnnotatedV1NT = NamedTuple{<:Any,<:Tuple{BatchedAnnotatedV1,Vararg{Any}}}

#Add protein repeats to 
#########Chain level:
const host_class_cats = Categorizer(Dict("Gammaproteobacteria" => 1, "Mammalia" => 2, "Insecta" => 3), 3)
const gene_class_cats = Categorizer(Dict("Gammaproteobacteria" => 1, "Mammalia" => 2, "Insecta" => 3), 3)
const gene_superkingdom_cats = Categorizer(Dict("Eukaryota" => 1, "Bacteria" => 2, "Viruses" => 3, "Archaea" => 4), 4)
const therm_cats = Categorizer([50.0, 65.0, 75.0, 85.0])
const ss_cats = Categorizer([0.2, 0.4, 0.6, 0.8])
#These should be chain level, but they're record level - will be fixed.
#const scop_type_cats = Categorizer(Dict("0" => 1, "1" => 2, "2" => 3, "3" => 4), 4)
#const scop_type_cluster_cats = Categorizer(Dict("0" => 1, "1" => 2, "2" => 3, "3" => 4, "4" => 5), 5)
##########Record level:
const soluble_cats = Categorizer(Dict(false => 1, true => 2), 2) #From the list of soluble PDBs for solubleMPNN


function featurize_rec(::AnnotatedV1, rec::ProteinStructure)
    sol = onehotcats(rec.soluble, soluble_cats)
    return sol
end

function featurize_chain(::AnnotatedV1, ch::ProteinChain)
    ss = onehotcats(mean(onehotbatch(unwrap(ch.secondary_structure), 1:3), dims = 2), ss_cats)[:] #Because there were 3 values
    therm = onehotcats(ch.temp["max"], therm_cats)
    gene_class = onehotcats(ch.gene_lineage["class"], gene_class_cats)
    host_class = onehotcats(ch.host_lineage["class"], host_class_cats)
    gene_superkingdom = onehotcats(ch.gene_lineage["superkingdom"], gene_superkingdom_cats)
    vcat(ss[:], therm, gene_class, host_class, gene_superkingdom)
end

#Current plan: since there aren't that many chains, broadcast the rec features over all chains.
function featurize_chains(pt::AnnotatedV1, rec::ProteinStructure)
    vcat(stack([featurize_rec(pt, rec) for _ in 1:length(rec)]), stack(featurize_chain.((pt,), rec)))
end

"""
    flatten(rec::ProteinStructure; T = Float32, extras::NamedTuple = (;))

Takes a ProteinStructure and returns a tuple of the translations, rotations, residue indices, and features for each chain.
`extras` must be a `NamedTuple`, so if you have a single extra value, use eg:
`flatten(rec, extras = (dataset = "PDB",))`
"""
function flatten(pt::AnnotatedV1, rec::ProteinStructure; T = Float32, cluster_offset = 0, dataset = "")
    L = sum(length.(values(rec)))
    chainids = vcat([fill(i, l) for (i, l) in enumerate(length.(rec))]...)
    broadcasted_feats = featurize_chains(pt, rec)[:,chainids] 
    translations = zeros(T, 3, 1, L)
    rotations = zeros(T, 3, 3, L)
    resinds = zeros(Int, L)
    aas = zeros(Int, L)
    pos = 1
    for (i, ch) in enumerate(values(rec))
        translations[:, 1, pos:pos+length(ch)-1] .= Frames(ch).translations
        rotations[:, :, pos:pos+length(ch)-1] .= Frames(ch).rotations
        resinds[pos:pos+length(ch)-1] .= ch.numbering
        aas[pos:pos+length(ch)-1] .= aa_to_ints(ch.sequence)
        pos += length(ch)
    end
    @assert pos-1 == L == length(aas)
    #Make named tuple
    return (;record_type = AnnotatedV1(), rec.name, dataset = dataset, len = L, cluster = rec.cluster + cluster_offset,
            locs = (translations .- mean(translations, dims = 3)) ./ unit_scaling, rots = rotations, AAs = aas,
            resinds = resinds, chainids = chainids, chainfeats = broadcasted_feats)
end

function batch_flatrecs(recs::AbstractVector{<:AnnotatedV1NT})
    maxlen = maximum([r.len for r in recs])
    b = length(recs)
    locs = zeros(Float32, 3, 1, maxlen, b)
    rots = zeros(Float32, 3, 3, maxlen, b)
    rots[1,1,:,:] .= 1
    rots[2,2,:,:] .= 1
    rots[3,3,:,:] .= 1
    AAs = zeros(Int, maxlen, b) .+ 21
    resinds = zeros(Int, maxlen, b)
    chainids = zeros(Int, maxlen, b)
    aa_cmask = ones(Bool, maxlen, b) #<- Poor sep of concerns
    res_cmask = ones(Bool, maxlen, b) #<- Poor sep of concerns
    padmask = zeros(Bool, maxlen, b) #<- If mask bug, come back to this and make it an Int again. But then you'll have to change this in the batched mask sampler: inds = b.padmask[:,i]
    chainfeats = similar(recs[1].chainfeats, size(recs[1].chainfeats,1), maxlen, b) .= 0
    for i in 1:b
        flat = recs[i]
        locs[:,:,1:flat.len,i] .= flat.locs
        rots[:,:,1:flat.len,i] .= flat.rots
        AAs[1:flat.len,i] .= flat.AAs
        resinds[1:flat.len,i] .= flat.resinds
        chainids[1:flat.len,i] .= flat.chainids
        chainfeats[:,1:flat.len,i] .= flat.chainfeats
        padmask[1:flat.len,i] .= 1
        if haskey(flat, :aa_cmask)
            aa_cmask[1:flat.len,i] .= flat.aa_cmask
        end
        if haskey(flat, :res_cmask)
            res_cmask[1:flat.len,i] .= flat.res_cmask
        end
    end
    return (;record_type = BatchedAnnotatedV1(), locs, rots, AAs, resinds, chainids, padmask, chainfeats, aa_cmask, res_cmask)
end


"""
Takes in a flat record, and drops out some blocks of features for each chain, and some residues entirely.
"""
function feature_and_residue_dropout(flatrec; block_drop_p = 0.2, block_size_dist = Poisson(5), feature_drop_p = 0.3, all_feature_drop_p = 0.2)
    L = length(flatrec.chainids)
    drop_mask = trues(L)
    if rand() < block_drop_p
        l = rand(block_size_dist)
        pos = rand((1-l):L) #We allow the start of the mask to fall before the start of the chain, and the end to extend beyong the end of the chain
        st,en = (max(1,pos),min(pos+l,L))
        drop_mask[st:en] .= false
    end
    if sum(drop_mask) < 15
        drop_mask .= true
    end
    chain_feats = copy(flatrec.chainfeats)
    chain_ids = unique(flatrec.chainids)
    for c in chain_ids
        drop_rows = rand(size(chain_feats,1)) .< feature_drop_p
        if rand() < all_feature_drop_p
            drop_rows .= true
        end
        chain_feats[drop_rows, flatrec.chainids .== c] .= 0
    end
    return (;len = sum(drop_mask), cluster = flatrec.cluster, locs = flatrec.locs[:,:,drop_mask], rots = flatrec.rots[:,:,drop_mask], AAs = flatrec.AAs[drop_mask],
    resinds = flatrec.resinds[drop_mask], chainids = flatrec.chainids[drop_mask], chainfeats = chain_feats[:,drop_mask])
end

