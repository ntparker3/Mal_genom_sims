library(argparse)
library(dplyr)
library(recombuddy)
library(dcifer)

source("simulation_functions.R")   # for simulate_local_ancestors etc.

parser <- ArgumentParser()
parser$add_argument("--input",               required = TRUE)
parser$add_argument("--region_info",    required = TRUE)   # MAD4HATTER_Pf3D7_withPrimers.tsv
parser$add_argument("--primers",        required = TRUE) 
parser$add_argument("--output_rds",          required = TRUE)
parser$add_argument("--gen_pop_size",        type = "integer", default = 500)
parser$add_argument("--final_pop_size",      type = "integer", default = 500)
parser$add_argument("--final_pop_gen2_prop", type = "double",  default = 0.3)
parser$add_argument("--anc_ancestor_size", type = "double",  default = 100)
parser$add_argument("--anc_pop_alpha",  type = "double")
parser$add_argument("--anc_coi_r",      type = "double")
parser$add_argument("--anc_coi_p",      type = "double")
parser$add_argument("--anc_k_s",        type = "double")
parser$add_argument("--anc_max_k",      type = "integer")
parser$add_argument("--anc_max_coi",    type = "integer")
p <- parser$parse_args()

set.seed(12345)

MAD4HATTER_Pf3D7_region_info <- read_tsv(p$region_info)

MAD4HATTER_Pf3D7_region_info_diversity_targets <- MAD4HATTER_Pf3D7_region_info %>%
  filter(class == "Diversity")

MAD4HATTER_Pf3D7_primers <- read_tsv(p$primers)

raw_MAD4HATTER_Pf3D7_primers_diversity_targets <- MAD4HATTER_Pf3D7_primers %>%
  filter(target %in% MAD4HATTER_Pf3D7_region_info_diversity_targets$newName)

raw_MAD4HATTER_Pf3D7_region_info_diversity_targets_locs_for_sim = MAD4HATTER_Pf3D7_region_info_diversity_targets %>%
  select(`#chrom`, fullStart, fullStop, newName) %>%
  dplyr::rename(chrom = `#chrom`, start = fullStart, end = fullStop, target = newName)

ghana_data <- read_tsv(p$input)

ghana_diversity <- ghana_data %>%
  filter(target_name %in% MAD4HATTER_Pf3D7_region_info_diversity_targets$newName)

gen_params <- c(
    pop_alpha = p$anc_pop_alpha, coi_r = p$anc_coi_r, coi_p = p$anc_coi_p,
    k_s = p$anc_k_s, max_k = p$anc_max_k, max_coi = p$anc_max_coi
)

ancestors <- simulate_local_ancestors(
    og_data            = ghana_diversity,
    gen1_params        = gen_params, gen1_pop_size = p$gen_pop_size, gen1_ancestor_size = p$anc_ancestor_size,
    gen2_params        = gen_params, gen2_pop_size = p$gen_pop_size, gen2_ancestor_size = p$anc_ancestor_size,
    gen3_params        = gen_params, gen3_pop_size = p$gen_pop_size,
    final_pop_size     = p$final_pop_size,
    final_pop_gen_2_prop = p$final_pop_gen2_prop
)

# ancestors[[1]] is the final_pop dataframe your mass_sim_for_null expects
saveRDS(ancestors[[1]], p$output_rds)