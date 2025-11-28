const unit_scaling = 10
const AAs = collect("ACDEFGHIKLMNPQRSTVWY")
const AA2ind = Dict(zip(AAs,1:length(AAs)))
aa_to_int(aa::Char) = get(AA2ind, aa, 21)
int_to_aa(i::Int) = i == 21 ? 'X' : AAs[i]
aa_to_ints(aa::AbstractString) = [aa_to_int(a) for a in aa]
ints_to_aa(ints::AbstractVector) = join(int_to_aa.(ints))

"""
    flatten(rec::ProteinStructure; T = Float32)

Takes a ProteinStructure and returns a tuple of the translations, rotations, residue indices, and features for each chain.
"""
function flatten(rec::ProteinStructure; T = Float32)
    L = sum(length.(values(rec)))
    chainids = vcat([fill(i, l) for (i, l) in enumerate(length.(rec))]...)
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
    return (;len = L, cluster = rec.cluster, locs = (translations .- mean(translations, dims = 3)) ./ unit_scaling, rots = rotations, AAs = aas, resinds = resinds, chainids = chainids)
end

function batch_flatrecs(recs::AbstractVector)
    maxlen = maximum([r.len for r in recs])
    b = length(recs)
    locs = zeros(Float32, 3, 1, maxlen, b)
    rots = zeros(Float32, 3, 3, maxlen, b)
    rots[1,1,:,:] .= 1
    rots[2,2,:,:] .= 1
    rots[3,3,:,:] .= 1
    aas = zeros(Int, maxlen, b) .+ 21
    resinds = zeros(Int, maxlen, b)
    chainids = zeros(Int, maxlen, b)
    padmask = zeros(Int, maxlen, b)
    for i in 1:b
        flat = recs[i]
        locs[:,:,1:flat.len,i] .= flat.locs
        rots[:,:,1:flat.len,i] .= flat.rots
        aas[1:flat.len,i] .= flat.AAs
        resinds[1:flat.len,i] .= flat.resinds
        chainids[1:flat.len,i] .= flat.chainids
        padmask[1:flat.len,i] .= 1
    end
    return (;locs, rots, aas, resinds, chainids, padmask)
end

one_ind_per_cluster(clusters) = [rand(findall(clusters .== c)) for c in unique(clusters)]
length2batch(length_max, fac) = L -> floor(Int, length_max^fac/L^fac)

"""
    sample_batched_inds(flatrecs; l2b = length2batch(1000, 1.9))

Takes a vector of (flattened) protein structures, and returns a vector of indices into the original array, with each batch containing a random sample of one protein from each cluster.
"""
function sample_batched_inds(lens, clusters; l2b = length2batch(1000, 1.9))
    sampled_inds = filter(ind -> l2b(lens[ind]) > 0, one_ind_per_cluster(clusters))
    indices_lengths_jitter = [(ind, lens[ind], lens[ind] + 2randn()) for ind in sampled_inds]
    sort!(indices_lengths_jitter, by = x -> x[3])
    batch_inds = Vector{Int}[]
    current_batch = Int[]
    current_max_len = 0
    for (sampled_idx, original_len, _) in indices_lengths_jitter
        potential_max_len = max(current_max_len, original_len)
        if isempty(current_batch) || (length(current_batch) + 1 <= l2b(potential_max_len))
            push!(current_batch, sampled_idx)
            current_max_len = potential_max_len
        else
            push!(batch_inds, current_batch)
            current_batch = [sampled_idx]
            current_max_len = original_len
        end
    end
    if !isempty(current_batch)
        push!(batch_inds, current_batch)
    end
    return shuffle(batch_inds)
end

sample_batched_inds(flatrecs::MergedVector; l2b = length2batch(1000, 1.9)) = sample_batched_inds(flatrecs.len, flatrecs.cluster, l2b = l2b)

"""
    unflatten(locs, rots, seqints, chainids, resnums)
    unflatten(locs, rots, seqhots, chainids, resnums)  
    unflatten(locs, rots, seq, chainids, resnums)

Converts flattened protein structure data back into ProteinChain objects.

# Arguments
- `locs`: Array of translations/locations (3×1×L or 3×1×L×B for batched)
- `rots`: Array of rotations (3×3×L or 3×3×L×B for batched) 
- `seqints`/`seqhots`/`seq`: Sequence data as integers, one-hot encoding, or generic sequence
- `chainids`: Chain identifiers for each residue
- `resnums`: Residue numbers for each position

# Returns
- Vector of `ProteinChain` objects (or vector of vectors for batched input)

The function reconstructs protein chains from flattened representations, applying unit scaling
to locations and converting sequence integers back to amino acid strings.
"""
function unflatten(locs::AbstractArray{T,3}, rots::AbstractArray{T,3}, seqints::AbstractVector, chainids, resnums) where T
    seqstrs = ints_to_aa(seqints)
    return [ProteinChain(string(i),
                get_atoms(Frames(rots[:,:,chainids .== i], reshape(locs,3,:)[:,chainids .== i] .* unit_scaling)(ProteinChains.STANDARD_RESIDUE)),
                seqstrs[findall(chainids .== i)], resnums[findall(chainids .== i)]) 
            for i in unique(chainids)]
end

function unflatten(locs::AbstractArray{T,3}, rots::AbstractArray{T,3}, seqhots::AbstractArray, chainids, resnums) where T
    seqints = [argmax(x[:]) for x in eachslice(seqhots, dims = 2)]
    return unflatten(locs, rots, seqints, chainids, resnums)
end

function unflatten(locs::AbstractArray{T,4}, rots::AbstractArray{T,4}, seq, chainids, resnums) where T
    return [unflatten(locs[:,:,:,i], rots[:,:,:,i], selectdim(seq, ndims(seq), i), chainids[:,i], resnums[:,i]) for i in 1:size(locs, 4)]
end
