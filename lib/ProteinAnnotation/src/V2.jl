struct AnnotatedV2 end

const host_class_cats = Categorizer(Dict("Gammaproteobacteria" => 1, "Mammalia" => 2, "Insecta" => 3), 3)
const gene_class_cats = Categorizer(Dict("Gammaproteobacteria" => 1, "Mammalia" => 2, "Insecta" => 3), 3)
const gene_superkingdom_cats = Categorizer(Dict("Eukaryota" => 1, "Bacteria" => 2, "Viruses" => 3, "Archaea" => 4), 4)
const therm_cats = Categorizer([50.0, 65.0, 75.0, 85.0])
const ss_cats = Categorizer([0.2, 0.4, 0.6, 0.8])
const scop_type_cats = Categorizer(
    Dict("globular" => 1, "membrane" => 2, "unstructured" => 3, "fibrous" => 4), 5)

function features(::AnnotatedV2, row)
    ss = onehotcats([row.helix_count, row.sheet_count] ./ row.length, ss_cats)
    therm = onehotcats(row.avg_temp, therm_cats)
    gene_class = onehotcats(row.gene_class, gene_class_cats)
    host_class = onehotcats(row.host_class, host_class_cats)
    gene_superkingdom = onehotcats(row.gene_superkingdom, gene_superkingdom_cats)
    (; ss, therm, gene_class, host_class, gene_superkingdom)
end
