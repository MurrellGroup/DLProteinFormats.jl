########Conventions########
const unit_scaling = 10
const AAs = collect("ACDEFGHIKLMNPQRSTVWY")
const AA2ind = Dict(zip(AAs,1:length(AAs)))
aa_to_int(aa::Char) = get(AA2ind, aa, 21)
int_to_aa(i::Int) = i == 21 ? 'X' : AAs[i]
aa_to_ints(aa::AbstractString) = [aa_to_int(a) for a in aa]
ints_to_aa(ints::AbstractVector) = join(int_to_aa.(ints))

########QC########
"""
    qc_filter(rec::ProteinStructure)
    qc_filter(chain::ProteinChain)

Determines whether or not to discard a protein structure or chain.
"""
function qc_filter(rec::ProteinStructure)
    flag = true
    #flag &= !rec.truncated #For now we'll allow these
    for ch in values(rec)
        flag &= qc_filter(ch)
    end
    return flag
end

function qc_filter(chain::ProteinChain)
    flag = true
    flag &= chain.complete_renumber
    flag &= chain.is_sorted
    return flag
end

########Categorizing properties########
struct Categorizer{T}
    cats::T
    max::Int
    default::Int
end

Categorizer(c,m) = Categorizer(c,m,0)
Categorizer(a::AbstractVector{<:Real}) = Categorizer(v -> bin(v, a), length(a)+1)

"""
    bin(v, a::AbstractVector)

Returns the index of the bin that `v` falls into in `a`. If v is below the lowest value, it returns 1. If it is above the highest value, it returns `length(a) + 1`.
"""
bin(v, a::AbstractVector) = searchsortedfirst(vcat(a, [Inf]), v)

"""
    categorize(v::Nothing, f::Function) = 0    
    categorize(v, d::Dict) = get(d, v, 0)
    categorize(v, f::Function) = f(v)
"""
categorize(v::Nothing, d::Categorizer) = d.default
categorize(v::Union{Number, String, Char, Symbol, Bool}, d::Categorizer{<:Dict}) = get(d.cats, v, d.default)
categorize(v::Union{Number, String, Char, Symbol, Bool}, d::Categorizer{<:Function}) = d.cats(v)
categorize(v::AbstractArray, d::Categorizer{<:Function}) = categorize.(v, (d,))
categorize(v::AbstractArray, d::Categorizer{<:Dict}) = categorize.(v, (d,))

"""
    onehotcats(v, d::Categorizer)
    
Convert a value or array of values into a one-hot BitMatrix, fitting the categorization scheme.
"""
onehotcats(v, d::Categorizer) = onehotbatch(categorize(v, d), 0:d.max) .* true

########Batch sampling########
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

########Handling repeats########
export get_repeats_graph, get_residue_to_master, get_reduced_chainids, expand_mask

minimal_unflatten(fr) = ProteinStructure("",
    map(findall(==(i), fr.chainids) for i in unique(fr.chainids)) do is
        ProteinChain("", [Atom{Float32}[] for _ in is], join(AAs[fr.AAs[is]]), fr.resinds[is])
    end
)

get_repeats_graph(fr; kwargs...) = detect_repeats(minimal_unflatten(fr); connect_fully=true, kwargs...)

function get_master_map(chainids::Vector{Int}, g::Graphs.SimpleGraph)
    master_map = Dict{Int,Int}()
    for chainid in Random.shuffle(unique(chainids))
        chainid in keys(master_map) && continue
        matches = findall(==(chainid), chainids)
        copies = chainids[Graphs.neighbors(g, rand(matches))]
        for copy in copies
            master_map[copy] = chainid
        end
    end
    return master_map
end

function get_residue_to_master(chainids::Vector{Int}, g::Graphs.SimpleGraph)
    masters = Set(values(get_master_map(chainids, g)))
    residue_to_master = fill(-1, length(chainids))  # -1 indicates "no master residue found"
    for idx in eachindex(chainids)
        if chainids[idx] in masters
            residue_to_master[idx] = idx  # master residue maps to itself
        else
            neighbor_idxs = Graphs.neighbors(g, idx)
            master_neighbor_idx = findfirst(n -> chainids[n] in masters, neighbor_idxs)
            if master_neighbor_idx !== nothing
                residue_to_master[idx] = neighbor_idxs[master_neighbor_idx]
            end
        end
    end
    return residue_to_master
end

function get_reduced_chainids(chainids::Vector{Int}, residue_to_master::Vector{Int})
    master_residue_idxs = [idx for idx in eachindex(residue_to_master) if residue_to_master[idx] == idx]
    return chainids[master_residue_idxs]
end

function expand_mask(master_mask::Vector{Bool}, residue_to_master::Vector{Int}, chainids::Vector{Int})
    n = length(residue_to_master)
    expanded_mask = falses(n)
    mask_broadcasted = falses(n)

    master_residues = findall(i -> residue_to_master[i] == i, 1:n)
    master_idx_to_mask = Dict(master_residues .=> master_mask)

    for idx in 1:n
        master_idx = residue_to_master[idx]
        if master_idx == idx
            expanded_mask[idx] = master_idx_to_mask[idx]
            mask_broadcasted[idx] = true
        elseif master_idx > 0
            expanded_mask[idx] = master_idx_to_mask[master_idx]
            mask_broadcasted[idx] = true
        end
    end

    for chainid in unique(chainids)
        chain_idxs = findall(==(chainid), chainids)
        start_idx, end_idx = first(chain_idxs), last(chain_idxs)

        segment_start = nothing
        idx = start_idx
        while idx <= end_idx
            if !mask_broadcasted[idx]
                segment_start = idx
                while idx <= end_idx && !mask_broadcasted[idx]
                    idx += 1
                end
                segment_end = idx - 1

                left_boundary = segment_start - 1
                right_boundary = segment_end + 1

                left_valid = left_boundary ≥ start_idx && mask_broadcasted[left_boundary]
                right_valid = right_boundary ≤ end_idx && mask_broadcasted[right_boundary]

                chosen_mask = false
                if left_valid && right_valid
                    if expanded_mask[left_boundary] == expanded_mask[right_boundary]
                        chosen_mask = expanded_mask[left_boundary]
                    else
                        chosen_mask = rand(Bool) ? expanded_mask[left_boundary] : expanded_mask[right_boundary]
                    end
                elseif left_valid
                    chosen_mask = expanded_mask[left_boundary]
                elseif right_valid
                    chosen_mask = expanded_mask[right_boundary]
                else
                    chosen_mask = false
                end
                expanded_mask[segment_start:segment_end] .= chosen_mask
            else
                idx += 1
            end
        end
    end

    return expanded_mask
end

### </repeat masking> ###
#=
function sample_masks(flat_record, sample_chain_masks) #<- sample_chain_masks(l) takes a length, and returns a tuple of bool vectors, one for pos and one for aa
    g = get_repeats_graph(flat_record)
    chainids = flat_record.chainids
    residue_to_master = get_residue_to_master(chainids, g)
    reduced_chainids = get_reduced_chainids(chainids, residue_to_master)
    reduced_posmask, reduced_aamask = ones(Bool, length(reduced_chainids)), ones(Bool, length(reduced_chainids))
    for c in unique(reduced_chainids)
        chainposmask, chainaamask = sample_chain_masks(sum(reduced_chainids .== c))
        reduced_posmask[reduced_chainids .== c] .= chainposmask
        reduced_aamask[reduced_chainids .== c] .= chainaamask
    end
    if sum(reduced_posmask) + sum(reduced_aamask) == 0 #If we wound up with nothing flowing in all chains, we make everything flow
        reduced_posmask .= true
        reduced_aamask .= true
    end
    posmask = expand_mask(reduced_posmask, residue_to_master, chainids)
    aamask = expand_mask(reduced_aamask, residue_to_master, chainids)
    return posmask, aamask
end
=#

function sample_repeat_aware_masks(flat_record, num_masks, sample_chain_masks) #<- sample_chain_masks(l) takes a length, and returns a tuple of bool vectors, one for pos and one for aa
    g = get_repeats_graph(flat_record)
    chainids = flat_record.chainids
    residue_to_master = get_residue_to_master(chainids, g)
    reduced_chainids = get_reduced_chainids(chainids, residue_to_master)
    reduced_masks = [ones(Bool, length(reduced_chainids)) for _ in 1:num_masks]
    for c in unique(reduced_chainids)
        master_masks = sample_chain_masks(sum(reduced_chainids .== c))
        for i in 1:num_masks
            reduced_masks[i][reduced_chainids .== c] .= master_masks[i]
        end
    end
    if sum(sum.(reduced_masks)) == 0 #If we wound up with nothing flowing in all chains, we make everything flow
        for i in 1:num_masks
            reduced_masks[i] .= true
        end
    end
    expanded_masks = [expand_mask(reduced_masks[i], residue_to_master, chainids) for i in 1:num_masks]
    return expanded_masks
end

function batched_sample_repeat_aware_masks(b, num_masks, sample_chain_masks)
    batched_masks = [similar(b.padmask) .= 0 for _ in 1:num_masks]
    for i in 1:size(b.padmask,2)
        inds = b.padmask[:,i]
        minimal_fr = (; chainids = b.chainids[inds,i], AAs = b.AAs[inds,i], resinds = b.resinds[inds,i])
        chain_masks = sample_repeat_aware_masks(minimal_fr, num_masks, sample_chain_masks)
        for j in 1:num_masks
            batched_masks[j][inds,i] .= chain_masks[j]
        end
    end
    return batched_masks
end



#Getting from flat locs, rots, etc to ProteinStructure
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
