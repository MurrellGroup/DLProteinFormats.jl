using ProteinChains
using ProgressMeter

export import_from_pss

function import_from_pss(store::ProteinChains.ProteinStructureStore)
    data = []
    @showprogress for (name, structure) in store
        for chain in structure
            contains('-')(chain.id) && continue
            push!(data, (;
                pdb_id = lowercase(name),
                chain_id = chain.id,
                length = length(chain) |> Int32,
                coil_count = count(==(1), chain.secondary_structure |> unwrap) |> Int32,
                helix_count = count(==(2), chain.secondary_structure |> unwrap) |> Int32,
                sheet_count = count(==(3), chain.secondary_structure |> unwrap) |> Int32,
                gene_taxid = chain.gene_taxid |> x -> isone(x) ? missing : x,
                host_taxid = chain.host_taxid |> x -> isone(x) ? missing : x,
                max_temp = something(chain.temp["max"], missing),
                min_temp = something(chain.temp["min"], missing),
                avg_temp = something(chain.temp["avg"], missing)
            ))
        end
    end
    return DataFrame(data)
end
