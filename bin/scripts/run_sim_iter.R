library(argparse)
library(dplyr)
library(recombuddy)
library(dcifer)

source("simulation_functions.R")

parser <- ArgumentParser()
parser$add_argument("--ancestors_rds",  required = TRUE)
parser$add_argument("--region_info",    required = TRUE)  
parser$add_argument("--primers",        required = TRUE)   
parser$add_argument("--iter_id",        type = "integer", required = TRUE)
parser$add_argument("--seed",           type = "integer", required = TRUE)
parser$add_argument("--sim_pop_size",   type = "integer", default = 1000)
parser$add_argument("--num_samples",    type = "integer", default = 100)
parser$add_argument("--pop_alpha",      type = "double")
parser$add_argument("--coi_r",          type = "double")
parser$add_argument("--coi_p",          type = "double")
parser$add_argument("--k_s",            type = "double")
parser$add_argument("--max_k",          type = "integer")
parser$add_argument("--max_coi",        type = "integer")
parser$add_argument("--dcifer_output",  required = TRUE)
parser$add_argument("--summary_output", required = TRUE)
p <- parser$parse_args()

set.seed(p$seed)

# ── Reference data (same pattern as prep_ancestors.R) ────────────────────────
MAD4HATTER_Pf3D7_region_info <- read_tsv(p$region_info)

MAD4HATTER_Pf3D7_region_info_diversity_targets <- MAD4HATTER_Pf3D7_region_info %>%
  filter(class == "Diversity")

MAD4HATTER_Pf3D7_primers <- read_tsv(p$primers)

raw_MAD4HATTER_Pf3D7_primers_diversity_targets <- MAD4HATTER_Pf3D7_primers %>%
  filter(target %in% MAD4HATTER_Pf3D7_region_info_diversity_targets$newName)

raw_MAD4HATTER_Pf3D7_region_info_diversity_targets_locs_for_sim <- MAD4HATTER_Pf3D7_region_info_diversity_targets %>%
  select(`#chrom`, fullStart, fullStop, newName) %>%
  dplyr::rename(chrom = `#chrom`, start = fullStart, end = fullStop, target = newName)


ancestor_df <- readRDS(p$ancestors_rds)

sim_params <- c(
    pop_alpha = p$pop_alpha, coi_r = p$coi_r, coi_p = p$coi_p,
    k_s = p$k_s, max_k = p$max_k, max_coi = p$max_coi
)

# ── Run one simulation ────────────────────────────────────────────────────────
sim_pop <- do.call(sim_population,
    c(list(unique(ancestor_df$library_sample_name)), p$sim_pop_size, as.list(sim_params)))

intersected <- filter_and_intersect(sim_pop, ancestor_df, "Final_sim")

dsmp  <- formatDat(intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
coi   <- getCOI(dsmp, lrank = 2)
afreq <- calcAfreq(dsmp, coi, tol = 1e-5)

dcifer_result <- ibdDat(dsmp, coi, afreq,
    pval = TRUE, confint = TRUE, rnull = 0, alpha = 0.05, nr = 1e3)

# ── Full dcifer table (what you want to save per-sim) ────────────────────────
ibd_table <- as.data.frame.table(dcifer_result, responseName = "value") %>%
    rename(specimen_a = Var1, specimen_b = Var2, metric = Var3) %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(iter_id = p$iter_id)

gz_con <- gzfile(p$dcifer_output, "w")
write_tsv(ibd_table, gz_con)
close(gz_con)        # but it's saved per-iter as requested

# ── Per-sample mean/median relatedness (summary row) ─────────────────────────
mean_ibd <- ibd_table %>%
    group_by(specimen_a) %>%
    summarize(
        mean_relatedness   = mean(estimate,   na.rm = TRUE),
        median_relatedness = median(estimate, na.rm = TRUE),
        .groups = "drop"
    )

sampled_ibd <- mean_ibd %>%
    slice_sample(n = p$num_samples)

summary_row <- sampled_ibd %>%
    summarize(
        iter_id            = p$iter_id,
        seed               = p$seed,
        sim_pop_size       = p$sim_pop_size,
        pop_alpha          = p$pop_alpha,
        mean_rel_across_samples   = mean(mean_relatedness),
        median_rel_across_samples = median(mean_relatedness),
        sd_rel_across_samples     = sd(mean_relatedness)
    )

write_tsv(summary_row, p$summary_output)