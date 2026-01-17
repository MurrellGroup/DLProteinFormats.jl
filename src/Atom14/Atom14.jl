struct Atom14 end

using ProteinChains
using ProteinChains: threeletter_resname, oneletter_resname

const ref_atoms = Dict{String, Vector{String}}(
    "PAD" => [],
    "UNK" => ["N", "CA", "C", "O", "CB"],
    "-" => [],
    "GLY" => ["N", "CA", "C", "O"],  # 0
    "ALA" => ["N", "CA", "C", "O", "CB"], # 1
    "CYS" => ["N", "CA", "C", "O", "CB", "SG"], # 2
    "SER" => ["N", "CA", "C", "O", "CB", "OG"], # 2
    "PRO" => ["N", "CA", "C", "O", "CB", "CG", "CD"],# 3
    "THR" => ["N", "CA", "C", "O", "CB", "OG1", "CG2"],# 3
    "VAL" => ["N", "CA", "C", "O", "CB", "CG1", "CG2"],# 3
    "ASN" => ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],# 4
    "ASP" => ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],# 4
    "ILE" => ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],# 4
    "LEU" => ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],# 4
    "MET" => ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],# 4
    "GLN" => ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],# 5
    "GLU" => ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],# 5
    "LYS" => ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],# 5
    "HIS" => ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],# 6
    "PHE" => ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],# 7
    "ARG" => ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],# 7
    "TYR" => ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"], # 8
    "TRP" => ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],  # 10 noqa: E501
)

atom_name_to_element(name::String) = name[1:1]

const protein_backbone_atom_names = ["N", "CA", "C", "O"]
const fake_atom_placements = Dict{String, Vector{String}}(
    "UNK" => [".", ".", ".", ".", ".", "N", "N", "N", "N", "N", "N", "N", "N", "N"], # 0
    "GLY" => [".", ".", ".", ".", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O"], # 0
    "ALA" => [".", ".", ".", ".", ".", "O", "O", "O", "O", "O", "O", "O", "O", "O"], # 1
    "CYS" => [".", ".", ".", ".", ".", ".", "O", "O", "O", "O", "O", "O", "O", "O"], # 2
    "SER" => [".", ".", ".", ".", ".", ".", "N", "N", "N", "N", "N", "N", "N", "N"], # 2
    "PRO" => [".", ".", ".", ".", ".", ".", ".", "O", "O", "O", "O", "O", "O", "O"], # 3
    "THR" => [".", ".", ".", ".", ".", ".", ".", "N", "N", "N", "O", "O", "O", "O"], # 3
    "VAL" => [".", ".", ".", ".", ".", ".", ".", "N", "N", "N", "N", "N", "N", "N"], # 3
    "ILE" => [".", ".", ".", ".", ".", ".", ".", ".", "O", "O", "O", "O", "O", "O"], # 4
    "ASN" => [".", ".", ".", ".", ".", ".", ".", ".", "N", "O", "O", "O", "O", "O"], # 4
    "ASP" => [".", ".", ".", ".", ".", ".", ".", ".", "N", "N", "O", "O", "O", "O"], # 4
    "LEU" => [".", ".", ".", ".", ".", ".", ".", ".", "N", "N", "N", "N", "O", "O"], # 4
    "MET" => [".", ".", ".", ".", ".", ".", ".", ".", "N", "N", "N", "N", "N", "N"], # 4
    "GLN" => [".", ".", ".", ".", ".", ".", ".", ".", ".", "O", "O", "O", "O", "O"], # 5
    "GLU" => [".", ".", ".", ".", ".", ".", ".", ".", ".", "N", "N", "O", "O", "O"], # 5
    "LYS" => [".", ".", ".", ".", ".", ".", ".", ".", ".", "N", "N", "N", "N", "N"], # 5
    "HIS" => [".", ".", ".", ".", ".", ".", ".", ".", ".", ".", "O", "O", "O", "O"], # 6
    "PHE" => [".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", "O", "O", "O"], # 7
    "ARG" => [".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", "N", "N", "N"], # 7
    "TYR" => [".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", "O", "O"], # 8
    "TRP" => [".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", "."], # 10
)

_get(x, i, fallback="UNK") = get(x, i, x[fallback])

# Dict("UNK" => [9, 0, 0, 0], "GLY" => [0, 0, 0, 10], "ALA" => [0, 0, 0, 9], "CYS" => [0, 0, 0, 8], ...)
const token_to_placement_count = Dict{String, Vector{Int}}(token_type => [count(==(backbone_atom_name), placement) for backbone_atom_name in protein_backbone_atom_names] for (token_type, placement) in pairs(fake_atom_placements))
const placement_count_to_token = Dict{NTuple{4, Int}, String}(Tuple(v) => k for (k, v) in pairs(token_to_placement_count))

const expected_atom_count = Dict{String, Int}(k => length(v) for (k, v) in ref_atoms)

struct MissingAtomError <: Exception
    atom::String
    residue::String
end

function add_atom14!(fake_atom_coords::AbstractMatrix{<:Real}, atom_coords, atom_names::Vector{<:AbstractString}, residue_type_3letter::String)
    ref_atom_names = _get(ref_atoms, residue_type_3letter)
    placements = _get(fake_atom_placements, residue_type_3letter)
    for (i, placement) in enumerate(placements)
        target_atom = placement == "." ? ref_atom_names[i] : placement
        j = findfirst(==(target_atom), atom_names)
        j === nothing && throw(MissingAtomError(target_atom, residue_type_3letter))
        fake_atom_coords[:, i] .= atom_coords[j]
    end
    return fake_atom_coords
end

function add_atom14_property!(chain::ProteinChain; sequence=collect(chain.sequence))
    L = length(chain)
    atom14_coords = zeros(Float32, 3, 14, L)
    has_complete_side_chain = Vector{Bool}(undef, L)
    @inbounds for i in 1:L
        atom_coords = map(ProteinChains.atom_coords, chain.atoms[i])
        residue_type = sequence[i]
        residue_type_3letter = threeletter_resname(residue_type)
        atom_names = map(strip ∘ ProteinChains.atom_name, chain.atoms[i])
        try
            add_atom14!(view(atom14_coords, :, :, i), atom_coords, atom_names, residue_type_3letter)
            has_complete_side_chain[i] = true
        catch e
            e isa MissingAtomError || rethrow()
            add_atom14!(view(atom14_coords, :, :, i), atom_coords, atom_names, "GLY")
            has_complete_side_chain[i] = false
        end
    end
    chain.atom14_coords = atom14_coords
    return has_complete_side_chain
end

function atom14_to_atoms(atom14_coords::AbstractArray{T,3}, sequence::String) where T
    atoms = Vector{Vector{Atom{T}}}(undef, size(atom14_coords, 3))
    @assert size(atom14_coords, 2) == 14
    @assert size(atom14_coords, 3) == length(sequence)
    for (i, res) in enumerate(sequence)
        ref_atom_names = _get(ref_atoms, threeletter_resname(res), "GLY")
        num_atoms = length(ref_atom_names)
        @assert num_atoms >= 4
        atoms[i] = map(1:num_atoms) do j
            Atom(ref_atom_names[j], atom_name_to_element(ref_atom_names[j]), atom14_coords[:, j, i])
        end
    end
    return atoms
end

function flatten(::Atom14, structure::ProteinStructure; T=Float32)
    len = sum(length, structure)
    chainids = reduce(vcat, [fill(i, l) for (i, l) in enumerate(length.(structure))], init=Int[])
    locs = zeros(T, 3, 1, len)
    rots = zeros(T, 3, 3, len)
    resinds = zeros(Int, len)
    AAs = zeros(Int, len)
    atom14_coords = zeros(T, 3, 14, len)
    complete_side_chain = Vector{Bool}(undef, len)
    pos = 1
    chain_labels = String[]
    for (i, chain) in enumerate(structure)
        push!(chain_labels, chain.id)
        _complete_side_chain = add_atom14_property!(chain)
        chain_len = length(chain)
        chain_range = pos:pos+chain_len-1
        f = ProteinChains.Frames(chain)
        locs[:, 1, chain_range] .= f.translations ./ unit_scaling
        rots[:, :, chain_range] .= f.rotations
        resinds[chain_range] .= chain.numbering
        AAs[chain_range] .= DLProteinFormats.aa_to_ints(chain.sequence)
        atom14_coords[:, :, chain_range] .= chain.atom14_coords ./ unit_scaling
        complete_side_chain[chain_range] .= _complete_side_chain
        pos += chain_len
    end
    @assert pos-1 == len == length(AAs)
    return (; record_type = :Atom14, structure.name, chain_labels, len,
        locs = locs, rots, AAs,
        resinds, chainids, atom14_coords, complete_side_chain)
end

function unflatten(::Atom14, seqints, chainids, resnums, atom14_coords)
    sequence = ints_to_aa(seqints)

    chains = map(unique(chainids)) do i
        chain_range = findall(chainids .== i)
        seq = sequence[chain_range]
        chain = ProteinChain(string(i), atom14_to_atoms(view(atom14_coords,:,:,chain_range) .* unit_scaling, seq), seq, Int.(resnums[chain_range]))
        chain.atom14_coords = view(atom14_coords,:,:,chain_range) .* unit_scaling
        return chain
    end

    return ProteinStructure("structure", chains)
end

function atom14_to_aacodes(atom14_coords::AbstractArray{<:Real,3}; radius::Real=0.02) #Tolerance assuming nm
    r2 = radius^2
    L = size(atom14_coords, 3)
    seq = Vector{Char}(undef, L)

    @inbounds for i in 1:L
        n1, n2, n3 = atom14_coords[1,1,i], atom14_coords[2,1,i], atom14_coords[3,1,i]
        a1, a2, a3 = atom14_coords[1,2,i], atom14_coords[2,2,i], atom14_coords[3,2,i] # CA
        c1, c2, c3 = atom14_coords[1,3,i], atom14_coords[2,3,i], atom14_coords[3,3,i] # C
        o1, o2, o3 = atom14_coords[1,4,i], atom14_coords[2,4,i], atom14_coords[3,4,i] # O

        counts = Int[0, 0, 0, 0]  # N, CA, C, O (virtual placements among atoms 5:14)

        for j in 5:14
            x, y, z = atom14_coords[1,j,i], atom14_coords[2,j,i], atom14_coords[3,j,i]

            d2n  = (x-n1)^2 + (y-n2)^2 + (z-n3)^2
            d2ca = (x-a1)^2 + (y-a2)^2 + (z-a3)^2
            d2c  = (x-c1)^2 + (y-c2)^2 + (z-c3)^2
            d2o  = (x-o1)^2 + (y-o2)^2 + (z-o3)^2

            bestk, bestd2 = 1, d2n
            if d2ca < bestd2; bestk, bestd2 = 2, d2ca; end
            if d2c  < bestd2; bestk, bestd2 = 3, d2c;  end
            if d2o  < bestd2; bestk, bestd2 = 4, d2o;  end

            if bestd2 <= r2
                counts[bestk] += 1
            end
        end

        token = get(placement_count_to_token, (counts[1], counts[2], counts[3], counts[4]), "UNK")
        if token == "UNK" || token == "PAD" || token == "-"
            seq[i] = 'X'
        else
            ol = oneletter_resname(token)
            seq[i] = ol isa Char ? ol : ol[1]
        end
    end

    return String(seq)
end

int_seq_from_at14(at14) = aa_to_ints(atom14_to_aacodes(at14))
