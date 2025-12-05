const URL_SABDAB_ALL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/"
const URL_SABDAB_NANO = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all_nano/"

export get_antibody_data

function get_antibody_data()
    df_std_raw = CSV.read(download(URL_SABDAB_ALL), DataFrame, missingstring=["NA", "None", "nan", "-"], silencewarnings=true)
    df_nano_raw = CSV.read(download(URL_SABDAB_NANO), DataFrame, missingstring=["NA", "None", "nan", "-"], silencewarnings=true)

    df_std  = ingest_rows(df_std_raw, is_nano=false)
    df_nano = ingest_rows(df_nano_raw, is_nano=true)

    claimed = Set(zip(df_nano.pdb_id, df_nano.chain_id))
    df_std_unique = filter(row -> (row.pdb_id, row.chain_id) ∉ claimed, df_std)
    
    return vcat(df_nano, df_std_unique)
end

"""
Converts 'IGHV3' -> 'VH3' (or 'VHH3' if from Nanobody file).
"""
function parse_shorthand(raw_sub, is_nano::Bool)
    ismissing(raw_sub) && return (; locus=missing, subclass=missing)
    s = uppercase(raw_sub)
    # Regex: Capture Letter (H/K/L) and Family Number
    m = match(r"IG([HKL])V(\d+)", s)
    if m !== nothing
        letter, num = m.captures[1], m.captures[2]
        if letter == "H"
            # File source determines the locus
            loc = is_nano ? "VHH" : "Heavy"
            sub = is_nano ? "VHH"*num : "VH"*num
            return (; locus=loc, subclass=sub)
        elseif letter == "K"
            return (; locus="Kappa", subclass="VK"*num)
        elseif letter == "L"
            return (; locus="Lambda", subclass="VL"*num)
        end
    else
        if s in ("ANTIGEN", "PROTEIN", "PEPTIDE")
            return (; locus="Antigen", subclass="Antigen")
        end
    end
    return (; locus=missing, subclass=missing)
end

function ingest_rows(raw_df; is_nano::Bool)
    long_df = DataFrame(
        pdb_id = String[], 
        chain_id = String[], 
        locus = Vector{Union{String, Missing}}(), 
        subclass = Vector{Union{String, Missing}}()
    )

    for row in eachrow(raw_df)
        pdb = lowercase(row.pdb)
        function add(cid, raw_sub)
            if !ismissing(cid) && cid != "NA" && cid != "None"
                meta = parse_shorthand(raw_sub, is_nano)
                push!(long_df, (;
                    pdb_id = pdb, 
                    chain_id = String(cid), 
                    locus = meta.locus, 
                    subclass = meta.subclass
                ))
            end
        end

        add(row.Hchain, row.heavy_subclass)
        add(row.Lchain, row.light_subclass)

        # Split Antigens (A | B)
        if !ismissing(row.antigen_chain) && row.antigen_chain != "NA"
            for ag in split(row.antigen_chain, r"\s*\|\s*")
                add(ag, "ANTIGEN")
            end
        end
    end

    return unique(long_df)
end
