#=
#V1 feature set. Intervals are derived by quantiles of *unique* values.
bin_feature_names = ["contact_order_over_sqrt_length",
                          "radius_of_gyration_over_cube_root_length", 
                          "interchain_atom_prop", 
                          "helix_proportion", 
                          "sheet_proportion",
                          "cys_prop",
                          "pairedcys_ratio",
                          "avg_temp",
                          "scTM",
                          "interchain_residue_prop",
                          "length"]

for f in bin_feature_names
    print(f," = ")
    print("(col = \"$(f)\", d = ")
    print(quantile(union(skipmissing(chainfeatures[!,f])), [0.2, 0.4, 0.6, 0.8]))
    println("),")
end
=#
const CHAIN_FEATS_V1 = (
    ig =        (col = "locus", d = Dict("Heavy" => 1, "Lambda" => 2, "Kappa" => 2, "VHH" => 3)),
    ag =        (col = "locus", d = Dict("Antigen" => 1)),
    genetax =   (col = "gene_superkingdom", d = Dict("Eukaryota" => 1, "Bacteria" => 2, "Archaea" => 3)),
    exprtax =   (col = "host_superkingdom", d = Dict("Eukaryota" => 1, "Bacteria" => 2, "Archaea" => 3)),
    prottype =  (col = "protein_type", d = Dict("globular" => 1, "membrane" => 2)),
    contact_order_over_sqrt_length = (col = "contact_order_over_sqrt_length", d = [0.7325347827051389, 1.0703505020915856, 1.2922797590340793, 1.5753821959445977]),
    radius_of_gyration_over_cube_root_length = (col = "radius_of_gyration_over_cube_root_length", d = [2.8571601663905595, 3.0451660334222233, 3.380207647237415, 3.996398936125884]),
    interchain_atom_prop = (col = "interchain_atom_prop", d = [0.07700198342597483, 0.13216181619044248, 0.19808429001635658, 0.31703558148141625]),
    helix_proportion = (col = "helix_proportion", d = [0.2172976323401412, 0.35179543193742807, 0.453972149186873, 0.5978494390234195]),
    sheet_proportion = (col = "sheet_proportion", d = [0.11064068180463674, 0.17411438399286255, 0.24727922780444928, 0.35310196114880565]),
    cys_prop = (col = "cys_prop", d = [0.00869702806824436, 0.01600465793304222, 0.024576411453660067, 0.038574819838914795]),
    pairedcys_ratio = (col = "pairedcys_ratio", d = [0.12580645161290324, 0.4358695652173913, 0.7471428571428571, 0.9227752639517346]),
    avg_temp = (col = "avg_temp", d = [48.8, 57.1, 69.8, 82.7]),
    scTM = (col = "scTM", d = [0.39088, 0.60834, 0.79817, 0.92681]),
    interchain_residue_prop = (col = "interchain_residue_prop", d = [0.09684589853010905, 0.18772309835961742, 0.2775269773269208, 0.4256217258921648]),
    length = (col = "length", d = [20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000]),
)


CHAIN_FEATS_64 = (
    genetax =   (col = "gene_superkingdom", d = Dict("Eukaryota" => 1, "Bacteria" => 2)),
    contact_order_over_sqrt_length = (col = "contact_order_over_sqrt_length", d = [0.7325347827051389, 1.0703505020915856, 1.2922797590340793, 1.5753821959445977]),
    radius_of_gyration_over_cube_root_length = (col = "radius_of_gyration_over_cube_root_length", d = [2.8571601663905595, 3.0451660334222233, 3.380207647237415, 3.996398936125884]),
    interchain_atom_prop = (col = "interchain_atom_prop", d = [0.07700198342597483, 0.13216181619044248, 0.19808429001635658, 0.31703558148141625]),
    helix_proportion = (col = "helix_proportion", d = [0.2172976323401412, 0.35179543193742807, 0.453972149186873, 0.5978494390234195]),
    sheet_proportion = (col = "sheet_proportion", d = [0.11064068180463674, 0.17411438399286255, 0.24727922780444928, 0.35310196114880565]),
    cys_prop = (col = "cys_prop", d = [0.00869702806824436, 0.01600465793304222, 0.024576411453660067, 0.038574819838914795]),
    pairedcys_ratio = (col = "pairedcys_ratio", d = [0.12580645161290324, 0.4358695652173913, 0.7471428571428571, 0.9227752639517346]),
    interchain_residue_prop = (col = "interchain_residue_prop", d = [0.09684589853010905, 0.18772309835961742, 0.2775269773269208, 0.4256217258921648]),
    symmetric = (col = "has_similar", d = Dict(false => 1, true => 2)),
    length = (col = "length", d = [30, 60, 100, 130, 190, 250, 350, 500]),
)

rowlookup(chainfeatures) = Dict(zip((zip(chainfeatures.pdb_id, chainfeatures.chain_id)), 1:length(chainfeatures.pdb_id)))
getrow(pdb_id, chain_id, rowmap) = get(rowmap, (pdb_id, chain_id), nothing)

swap_value(x,chain_id,featname,::Nothing) = x
swap_value(x,chain_id,featname,override) = get(get(override, chain_id, Dict()), featname, x)

rand_cats(v, p) = v .| (rand(length(v)) .< p)

"""
    featurizer(table::DataFrame, features::NamedTuple; all_mask_prob = 0, feat_mask_prob = 0)

Returns a function that will convert a PDB and chain to a feature vector.
`all_mask_prob` is the probability of masking all features for a chain.
`feat_mask_prob` is the probability of masking a feature for a chain, or a UnivariateDistribution that will be sampled for each chain, giving the probability of masking each feature for that chain.
`rand_cats_prob` is the probability of randomly setting an extra category for a feature to positve. This allows conditioning on a range of values.
`rand_cats_weight` is the number of chains that will have non-zero rand_cats_prob.

In the returned function, `override` is a dictionary of chain_id => dictionary of feature name => value.

Example:
chainfeatures = DLProteinFormats.load(PDBTable)
ff = featurizer(chainfeatures, features, all_mask_prob = 0.33, feat_mask_prob = Uniform(0,1))
ff("7A7B", "D", override = Dict(["D" => Dict(["gene_superkingdom" => "Eukaryota"])]))
"""
function featurizer(table::DataFrame, features::NamedTuple; all_mask_prob = 0, feat_mask_prob = 0, rand_cats_prob = 0, rand_cats_weight = 0)
    rl = rowlookup(table)
    cats = [(col = v.col, cat = DLProteinFormats.Categorizer(v.d)) for v in values(features)]
    feat_len = sum([c.cat.max + 1 for c in cats])
    function featurize(pdb_id, chain_id; override = nothing)
        r = getrow(pdb_id, chain_id, rl)
        if isnothing(r)
            return falses(feat_len)
        end
        feat_p = 0.0
        if feat_mask_prob isa Number
            feat_p = feat_mask_prob
        else
            feat_p = rand(feat_mask_prob)
        end
        rcw = rand() < rand_cats_weight ? rand_cats_prob : 0.0
        return (rand() > all_mask_prob) .* vcat([(rand() > feat_p) .* rand_cats(onehotcats(swap_value(table[r, cat.col], chain_id, cat.col, override), cat.cat), rcw) for cat in cats]...)
    end
    return featurize
end

pdbid_clean(s) = split(s, "_")[1]
chain_label_clean(s) = split(s, "-")[1]

"""
    broadcast_features(pdbids, chain_labels, batched_chainids, feature_func)

Broadcasts feature vectors for a set of PDBs and chain labels over a batched set of chain ids.
`pdbids` is a vector of PDB ids.
`chain_labels` is a vector of chain labels.
`batched_chainids` is a matrix of chain ids.
`feature_func` is a function that will convert a PDB and chain to a feature vector.

Example:
dat = DLProteinFormats.load(PDBSimpleFlatV2_500)
#Construct a little batch:
pdbids = [dat[10].name, dat[29].name]
chain_labels = [dat[10].chain_labels, dat[29].chain_labels]
chainids = zeros(Int, maximum([dat[29].len, dat[10].len]), 2)
chainids[1:dat[10].len,1] .= dat[10].chainids
chainids[1:dat[29].len,2] .= dat[29].chainids

bf = broadcast_features(pdbids, chain_labels, chainids, ff)
"""
function broadcast_features(pdbids, chain_labels, batched_chainids, feature_func)
    pdbids = pdbid_clean.(pdbids)
    chain_labels = [chain_label_clean.(cl) for cl in chain_labels]
    template_vec = feature_func("?", "?") .* 0
    feats = falses(length(template_vec), size(batched_chainids)...)
    for i in 1:length(pdbids)
        for j in 1:length(chain_labels[i])
            fv = feature_func(pdbids[i], chain_labels[i][j])
            feats[:,batched_chainids[:,i] .== j, i] .= fv
        end
    end
    return feats
end

export featurizer, broadcast_features
