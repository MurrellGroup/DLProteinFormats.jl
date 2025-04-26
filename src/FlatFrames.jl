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
function sample_batched_inds(flatrecs; l2b = length2batch(1000, 1.9))
    sampled_inds = filter(ind -> l2b(flatrecs.len[ind]) > 0, one_ind_per_cluster([r.cluster for r in flatrecs]))
    indices_lengths_jitter = [(ind, flatrecs.len[ind], flatrecs.len[ind] + 2randn()) for ind in sampled_inds]
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

#Add function to convert a flattened protein into a ProteinChains struct for plotting etc