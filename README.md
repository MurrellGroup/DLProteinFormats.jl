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
