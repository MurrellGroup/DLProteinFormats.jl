# DLProteinFormats

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/DLProteinFormats.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/DLProteinFormats.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/DLProteinFormats.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/DLProteinFormats.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/DLProteinFormats.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/DLProteinFormats.jl)

## Installation

```julia
using Pkg
Pkg.Registry.add("https://github.com/MurrellGroup/MurrellGroupRegistry")
Pkg.add("DLProteinFormats")
```

## Quickstart

```julia
using DLProteinFormats

data = DLProteinFormats.load(PDBSimpleFlat500);

flat_chains = data[1];

chains = DLProteinFormats.unflatten(
    flat_chains.locs,
    flat_chains.rots,
    flat_chains.AAs,
    flat_chains.chainids,
    flat_chains.resinds)

DLProteinFormats.writepdb("chains-1.pdb", chains) # view in e.g. chimerax or vscode protein viewer extension
```

## Flat atom PDB dataset

Flatom is a dataset that stores biomolecular structures in a minimal flat atom primitive format,
with each atom represented as a 28-byte NamedTuple with the following fields:

- `element::Int8`: element number
- `category::Int8`: structure category
  - 1: protein residue
  - 2: nucleic residue
  - 3: other (e.g. hetero)
- `chainid::Int16`: chain identifier
- `resnum::Int32`: residue number (can be optionally renumbered using the MMCIF file)
- `resname::StaticStrings.StaticString{3}`: 3-character alphanumeric residue name
- `atomname::StaticStrings.StaticString{4}`: 4-character alphanumeric atom name
- `coords::StaticArrays.SVector{3,Float32}`: 3D coordinates

```julia
using DLProteinFormats

# load 169k PDB structures at flat vectors of atoms (~500 million atoms total, 28 bytes per atom)
# once downloaded, takes ~6 seconds to load
structures = DLProteinFormats.load(PDBFlatom169K);

# atoms of the first structure
atoms = structures[1].atoms

# use the `stack(::Function, ...)` method to get a struct of array vibe
elements = stack(atom -> atom.element, atoms)
coords = stack(atom -> atom.coords, atoms)
# etc. etc.

# split by chain
chainids = map(a -> a.chainid, atoms)
chains = [atoms[findall(==(id), chainids)] for id in unique(chainids)]

using DLProteinFormats.Flatom

# write to PDB
write_to_pdb("chains-1.pdb", atoms)

# write to MMCIF
write_to_cif("chains-1.cif", atoms)
```
