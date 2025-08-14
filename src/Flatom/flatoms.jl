function get_category(residue)
    return backbone_residue_selector(residue) ? 1 : !ishetero(residue) ? 2 : 3
end

function flatoms(residue::AbstractResidue, chainid=1)
    chainid = chainid % Int16
    residue isa DisorderedResidue && (residue = defaultresidue(residue))
    resnum = resnumber(residue) % Int32
    resname = StaticString{3}(lpad(residue.name, 3))
    category = get_category(residue) % Int8
    map(collectatoms(residue)) do atom
        atom isa DisorderedAtom && (atom = defaultatom(atom))
        atomname = encode_atom_name(atom.name, atom.element)
        element = element_symbol_to_number(atom.element) % Int8
        return (;
            element,
            category,
            chainid,
            resnum,
            resname,
            atomname,
            coords = coords(atom)
        )
    end
end

function flatoms(chain::Chain; mmcif_dict=nothing, chainid=1)
    mmcif_dict isa MMCIFDict && renumber!(chain, mmcif_dict)
    return reduce(vcat, [flatoms(residue, chainid) for residue in collectresidues(chain)])
end

function flatoms(structure::MolecularStructure, mmcif_dict=nothing)
    atoms = reduce(vcat, [flatoms(chain; mmcif_dict, chainid) for (chainid, chain) in enumerate(collectchains(structure))])
    return atoms
end

function flatoms(filename::AbstractString, ::Type{MMCIFFormat})
    mmcif_dict = MMCIFDict(filename)
    atoms = flatoms(MolecularStructure(mmcif_dict; structure_name=basename(filename)), mmcif_dict)
    return atoms
end

flatoms(filename::AbstractString, format) = flatoms(read(filename, format)::MolecularStructure)

function _atom_record(flatom, i)
    AtomRecord(
        flatom.category == 3,
        i,
        flatom.atomname,
        ' ',
        flatom.resname,
        string(Char(64 + flatom.chainid)),
        flatom.resnum,
        ' ',
        flatom.coords,
        1.0,
        0.0,
        NUMBER_TO_ELEMENT_SYMBOL[flatom.element],
        "  ",
    )
end

function write_to_pdb(io::IO, atoms)
    for (i, atom) in enumerate(atoms)
        write(io, pdbline(_atom_record(atom, i)), '\n')
    end
end

write_to_pdb(path::AbstractString, atoms) = open(io -> write_to_pdb(io, atoms), path, "w")

function write_to_cif(path::AbstractString, atoms)
    mktempdir() do dir
        temp_path = joinpath(dir, "temp.pdb")
        write_to_pdb(temp_path, atoms)
        writemmcif(path, read(temp_path, PDBFormat))
    end
end