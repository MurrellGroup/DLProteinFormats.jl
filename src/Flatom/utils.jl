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