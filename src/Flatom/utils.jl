using PeriodicTable: elements
using StaticStrings
using StaticArrays

const ELEMENT_SYMBOL_TO_NUMBER = Dict(uppercase(elements[number].symbol) => number for number in 1:118)
const NUMBER_TO_ELEMENT_SYMBOL = Dict(n => s for (s, n) in ELEMENT_SYMBOL_TO_NUMBER)

element_symbol_to_number(element_symbol::AbstractString) = get(ELEMENT_SYMBOL_TO_NUMBER, uppercase(strip(element_symbol)), 0)
number_to_element_symbol(number::Integer) = get(NUMBER_TO_ELEMENT_SYMBOL, number, "X")

function pad_atom_name(name::AbstractString, element_symbol::AbstractString)
    length(name) == 4 && return name
    rpad(" "^(2-length(strip(element_symbol)))*strip(name), 4)
end

encode_atom_name(name::AbstractString, element_symbol::AbstractString) = StaticString{4}(pad_atom_name(name, element_symbol))

coords(atom::BioStructures.AbstractAtom) = SVector{3,Float32}(Float32(BioStructures.x(atom)), Float32(BioStructures.y(atom)), Float32(BioStructures.z(atom)))

const BACKBONE_ATOM_NAMES = ("N", "CA", "C")
backbone_atom_selector(atom::BioStructures.AbstractAtom) = BioStructures.atomnameselector(atom, BACKBONE_ATOM_NAMES)

function backbone_residue_selector(residue::Residue)
    standardselector(residue) && countatoms(residue, backbone_atom_selector) == 3
end

#OUR STANDARDIZED ENCODING SCHEME that tracks backbone atom names, and for everything else, the top 26 elements. 31 is other.
#=
element_counts = countmap(vcat([map(x -> x.element, dat[i].atoms) for i in 1:length(dat)]...));
common_elements = sort([e for e in keys(element_counts) if element_counts[e] > 1450])
=#
const element_coding_dict = Dict(zip(Int8[0, 1, 6, 7, 8, 9, 11, 12, 15, 16, 17, 19, 20, 25, 26, 27, 28, 29, 30, 33, 34, 35, 48, 53, 74, 80], 5:30))
const backbone_coding_dict = Dict(DLProteinFormats.static" N  "4 => 1,
                            DLProteinFormats.static" CA "4 => 2, 
                            DLProteinFormats.static" C  "4 => 3, 
                            DLProteinFormats.static" O  "4 => 4)

#119th element means "too rare, was coded as other".
const element_from_atom_code_dict = Dict(vcat(collect(zip(values(element_coding_dict),keys(element_coding_dict))),[(1,7),(2,6),(3,6),(4,8),(31,119)]))
#1,2,3,4 are for NA, CA, C, O. The rest are the top 26 elements. Then 31 is other. 32 will be reserved for "masked".
function atom_code(atom)
    if atom.category == 1 && haskey(backbone_coding_dict, atom.atomname)
        return backbone_coding_dict[atom.atomname]
    else
        return get(element_coding_dict,atom.element,31)
    end
end
element_from_atom_code(code) = element_from_atom_code_dict[code]
