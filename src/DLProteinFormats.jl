module DLProteinFormats

using ProteinChains, StatsBase, Random, MergedArrays

include("Flatom/Flatom.jl")

include("FlatFrames.jl")
export flatten, unflatten

include("datasets.jl")

include("shared.jl")

include("features.jl")

include("Atom14/Atom14.jl")
export Atom14


import ProteinChains: writepdb, writecif, pdbentry, @pdb_str
export writepdb, writecif, pdbentry, @pdb_str

"""
    writepdb(path, chains::AbstractVector{<:ProteinChains.ProteinChain})

## Examples

```julia
using DLProteinFormats

data = DLProteinFormats.load(PDBSimpleFlat500);

flat_chains = data[1];

chains = DLProteinFormats.unflatten(
    flat_chains.locs,
    flat_chains.rots,
    flat_chains.AAs,
    flat_chains.chainids,
    flat_chains.resinds) # unflatten the flat data

writepdb("chains-1.pdb", chains) # view in e.g. chimerax or vscode protein viewer extension
```
"""
writepdb

end
