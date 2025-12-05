module TaxMax

using Reexport
@reexport using Taxonomy

using Scratch, Downloads
using Tar, TranscodingStreams, CodecZlib

const TAXDB_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
const SCRATCH_NAME = "taxdump_scratch"
const taxonomy_scratch = Ref{String}()

function get_database(path=taxonomy_scratch[]) 
    nodes_path = joinpath(path, "nodes.dmp")
    names_path = joinpath(path, "names.dmp")
    if !all(isfile, (nodes_path, names_path))
        throw(ArgumentError("Taxonomy database not found at $path"))
    end
    # calling this function automatically caches
    # it as the current database in Taxonomy.jl
    return Taxonomy.DB(nodes_path, names_path)
end

function init_database()
    taxonomy_scratch[] = @get_scratch!(SCRATCH_NAME)
    try
        get_database()
    catch
        @info "Initializing Taxonomy database..."
        delete_scratch!(SCRATCH_NAME)
        taxonomy_scratch[] = @get_scratch!(SCRATCH_NAME)
        download_path = Downloads.download(TAXDB_URL)
        open(download_path, "r") do f
            io = TranscodingStream(CodecZlib.GzipDecompressor(), f)
            Tar.extract(io, taxonomy_scratch[])
        end
        get_database()
    end
end

__init__() = init_database()

export get_ranks, safe_taxon, safe_name

function get_ranks(taxon, ranks::Vector{Symbol})
    l = Lineage(taxon)
    return Dict{Symbol,Union{String,Missing}}(
        rank => get(l, rank, missing) |>
        t -> t isa Taxon ? name(t) : t
        for rank in ranks)
end

function safe_taxon(id::Int)
    try
        Taxonomy.Taxon(id)
    catch e
        #@warn "could not find taxon $id in taxonomy database. using root (1) instead."
        Taxonomy.Taxon(1)
    end
end

safe_taxon(::Missing) = missing

function safe_name(id::Int, rank::Symbol)
    try
        Taxonomy.name(Taxonomy.Lineage(Taxonomy.Taxon(id))[rank])
    catch e
        #@warn "$e"
        missing
    end
end

safe_name(taxon::Taxon, rank::Symbol) = safe_name(taxid(taxon), rank)
safe_name(::Missing, rank::Symbol) = missing

end
