using DataFrames

export get_scop_data

function parse_scopcla(scopcla_str::AbstractString)::Dict{String, String}
    SCOPCLA = Dict{String,String}()
    for (key, value) in split.(split(scopcla_str, ","), "=", limit=2)
        SCOPCLA[key] = value
    end
    return SCOPCLA
end

const PROTEIN_TYPE = Dict("1" => "globular", "2" => "membrane", "3" => "fibrous", "4" => "unstructured")

function parse_scop_line(line::AbstractString)
    fields = split(strip(line))
    length(fields) == 11 || throw(ArgumentError("Line needs to have 11 fields: $line"))
    cla = parse_scopcla(fields[11])
    return (
        FA_DOMID=fields[1],
        FA_PDBID=fields[2],
        FA_PDBREG=fields[3],
        FA_UNIID=fields[4],
        FA_UNIREG=fields[5],
        SF_DOMID=fields[6],
        SF_PDBID=fields[7],
        SF_PDBREG=fields[8],
        SF_UNIID=fields[9],
        SF_UNIREG=fields[10],
        SCOPCLA_TP=PROTEIN_TYPE[cla["TP"]],
        SCOPCLA_SF=cla["CL"],
        SCOPCLA_FA=cla["CF"],
        SCOPCLA_CL=cla["SF"],
        SCOPCLA_CF=cla["FA"]
    )
end

function parse_scop_file(filename::AbstractString)
    entries = NamedTuple[]
    open(filename, "r") do io
        for (line_number, line) in enumerate(eachline(io))
            (startswith(strip(line), "#") || isempty(strip(line))) && continue
            try
                entry = parse_scop_line(line)
                push!(entries, entry)
            catch e
                @warn "Failed to parse line $line_number: $line" exception=e
            end
        end
    end
    return DataFrame(entries)
end

using DataFrames

const TYPE_PRIORITY = Dict("membrane" => 4, "fibrous" => 3, "globular" => 2, "unstructured" => 1)

function resolve_conflicts(types_list)
    unique_types = unique(types_list)
    length(unique_types) == 1 && return unique_types[1]
    best_type = ""
    max_score = -1
    for t in unique_types
        score = get(TYPE_PRIORITY, t, 0)
        if score > max_score
            max_score = score
            best_type = t
        end
    end
    return best_type
end

get_raw_scop_data() = parse_scop_file(download("https://www.ebi.ac.uk/pdbe/scop/files/scop-cla-latest.txt"))

function get_scop_data()
    df = get_raw_scop_data()
    pdb_id = lowercase.(df.FA_PDBID)
    chain_id = first.(split.(df.FA_PDBREG, ":"))
    reduced = DataFrame(; pdb_id, chain_id, protein_type = df.SCOPCLA_TP)
    return combine(groupby(reduced, [:pdb_id, :chain_id])) do gdf
        (; protein_type = resolve_conflicts(gdf.protein_type))
    end
end
