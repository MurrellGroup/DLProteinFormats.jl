tryparse_fallback(T, str, fallback) = something(tryparse(T, str), fallback)

function renumber!(chain::Chain, mmcif_dict::BioStructures.MMCIFDict)
    id = split(chain.id, '-')[1]

    label_seq_ids = mmcif_dict["_atom_site.label_seq_id"]
    pdbx_PDB_ins_codes = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]
    auth_seq_ids = mmcif_dict["_atom_site.auth_seq_id"]
    auth_asym_ids = mmcif_dict["_atom_site.auth_asym_id"]

    label_seq_id_dict = Dict{Tuple{Int,Char},Int}()

    current_residue = ("", ' ')
    for (i, auth_asym_id) in enumerate(auth_asym_ids)
        auth_asym_id == id || continue
        ins_code = pdbx_PDB_ins_codes[i] == "?" ? " " : pdbx_PDB_ins_codes[i]

        auth_num = tryparse_fallback(Int, auth_seq_ids[i], -1)
        residue = (auth_num, first(ins_code))
        current_residue == residue && continue
        current_residue = residue

        label_seq_id_dict[residue] = tryparse_fallback(Int, label_seq_ids[i], auth_num)
    end

    for residue in collectresidues(chain)
        residue isa DisorderedResidue && (residue = defaultresidue(residue))
        residue.number = get(label_seq_id_dict, (resnumber(residue), inscode(residue)), residue.number)
    end

    return chain
end