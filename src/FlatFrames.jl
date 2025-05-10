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
