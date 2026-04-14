# Function to calculate Jost's D for two populations
calculate_jost_d <- function(pop1_afreq, pop2_afreq) {

  # Get all unique alleles across both populations
  all_alleles <- unique(c(names(pop1_afreq), names(pop2_afreq)))

  # Initialize results
  d_values <- numeric(length(all_alleles))
  names(d_values) <- names(all_alleles)

  # Calculate D for each locus
  for (i in seq_along(all_alleles)) {
    locus_name <- names(all_alleles)[i]

    # Get allele frequencies for this locus
    p1 <- pop1_afreq[[i]]
    p2 <- pop2_afreq[[i]]

    # Get all alleles at this locus
    alleles_locus <- unique(c(names(p1), names(p2)))

    # Create aligned frequency vectors (0 for missing alleles)
    freq1 <- setNames(numeric(length(alleles_locus)), alleles_locus)
    freq2 <- setNames(numeric(length(alleles_locus)), alleles_locus)

    freq1[names(p1)] <- p1
    freq2[names(p2)] <- p2

    # Within-population expected heterozygosity
    Hs1 <- 1 - sum(freq1^2)
    Hs2 <- 1 - sum(freq2^2)
    Hs <- mean(c(Hs1, Hs2))

    # Total expected heterozygosity
    freq_total <- (freq1 + freq2) / 2
    Ht <- 1 - sum(freq_total^2)

    # Jost's D (with correction for 2 populations)
    if (Hs == 1) {
      d_values[i] <- NA  # Undefined when Hs = 1
    } else {
      d_values[i] <- ((Ht - Hs) / (1 - Hs)) * (2 / (2 - 1))
      # Simplifies to: d_values[i] <- 2 * (Ht - Hs) / (1 - Hs)
    }
  }

  return(d_values)
}

# Function to calculate Nei's FST for two populations
calculate_FST <- function(pop1_afreq, pop2_afreq) {

  # Get all unique alleles across both populations
  all_alleles <- unique(c(names(pop1_afreq), names(pop2_afreq)))

  # Initialize results
  fst_values <- numeric(length(all_alleles))
  names(fst_values) <- names(all_alleles)

  # Calculate D for each locus
  for (i in seq_along(all_alleles)) {
    locus_name <- names(all_alleles)[i]

    # Get allele frequencies for this locus
    p1 <- pop1_afreq[[i]]
    p2 <- pop2_afreq[[i]]

    # Get all alleles at this locus
    alleles_locus <- unique(c(names(p1), names(p2)))

    # Create aligned frequency vectors (0 for missing alleles)
    freq1 <- setNames(numeric(length(alleles_locus)), alleles_locus)
    freq2 <- setNames(numeric(length(alleles_locus)), alleles_locus)

    freq1[names(p1)] <- p1
    freq2[names(p2)] <- p2

    # Within-population expected heterozygosity
    Hs1 <- 1 - sum(freq1^2)
    Hs2 <- 1 - sum(freq2^2)
    Hs <- mean(c(Hs1, Hs2))

    # Total expected heterozygosity
    freq_total <- (freq1 + freq2) / 2
    Ht <- 1 - sum(freq_total^2)

    # Nei's FST
    fst_values[i] <- ((Ht - Hs) / Ht)

  }

  return(fst_values)
}


# Function to calculate Jost's D for two populations
calculate_pop_expected_het <- function(pop_afreq) {

  # Get all unique alleles across both populations
  all_alleles <- unique(names(pop_afreq))

  # Initialize results
  he_values <- numeric(length(all_alleles))
  names(he_values) <- names(all_alleles)

  # Calculate D for each locus
  for (i in seq_along(all_alleles)) {
    locus_name <- names(all_alleles)[i]

    # Get allele frequencies for this locus
    p1 <- pop_afreq[[i]]

    # Get all alleles at this locus
    alleles_locus <- unique(names(p1))

    # Create aligned frequency vectors (0 for missing alleles)
    freq1 <- setNames(numeric(length(alleles_locus)), alleles_locus)

    freq1[names(p1)] <- p1

    # Within-population expected heterozygosity
    Hs1 <- 1 - sum(freq1^2)

    he_values[i] <- Hs1
  }

  return(he_values)
}



###########
##### function to choose ancestor pop - from an allele table, choosing the samples with good coverage that have a dominant haplotype
###########

filter_good_coverage_dom_haplo <- function(data, real_data = F){

  if(real_data == F){
    initial_target_coverage_required = 0.90
    initial_sample_coverage_required = 0.85
  }else if(real_data == T){
    initial_target_coverage_required = 0.90
    initial_sample_coverage_required = 0.84
  }

  MAD4HATTER_diversity_per_sample_target_coverage = data %>%
    group_by(library_sample_name) %>%
    summarise(target_count = n_distinct(target_name)) %>%
    mutate(target_freq = target_count/n_distinct(data$target_name))%>%
    arrange(target_freq)

  MAD4HATTER_diversity_per_sample_target_coverage_filt = MAD4HATTER_diversity_per_sample_target_coverage %>%
    filter(target_freq >= initial_target_coverage_required)

  MAD4HATTER_diversity_sample_filt = data %>%
    filter(library_sample_name %in% MAD4HATTER_diversity_per_sample_target_coverage_filt$library_sample_name)

  MAD4HATTER_diversity_sample_filt_target_coverage = MAD4HATTER_diversity_sample_filt %>%
    group_by(target_name) %>%
    summarise(sample_count = n_distinct(library_sample_name)) %>%
    mutate(sample_freq = sample_count/n_distinct(MAD4HATTER_diversity_sample_filt$library_sample_name)) %>%
    arrange(sample_freq)


  MAD4HATTER_diversity_sample_filt_target_coverage_filt = MAD4HATTER_diversity_sample_filt_target_coverage %>%
    filter(sample_freq >= initial_sample_coverage_required)


  MAD4HATTER_diversity_sample_filt_target_filt = MAD4HATTER_diversity_sample_filt %>%
    filter(target_name %in% MAD4HATTER_diversity_sample_filt_target_coverage_filt$target_name)

  ##Now to filter to dominant haplotypes

  dominant_haplotype_freq_cut_off = 0.60

  MAD4HATTER_diversity_sample_filt_target_filt = MAD4HATTER_diversity_sample_filt_target_filt %>%
    group_by(library_sample_name, target_name) %>%
    mutate(total_reads_target_sample = sum(reads)) %>%
    mutate(freq = reads/total_reads_target_sample)

  # note, this only works for setting dom hap > 0.5001
  MAD4HATTER_diversity_sample_filt_target_filt_pass = MAD4HATTER_diversity_sample_filt_target_filt %>%
    mutate(domiant_hap = freq >= dominant_haplotype_freq_cut_off) %>%
    group_by(library_sample_name) %>%
    summarise(targets = n_distinct(target_name),
              dom_haps = sum(domiant_hap)) %>%
    filter(targets == dom_haps)

  MAD4HATTER_diversity_final_samples = MAD4HATTER_diversity_sample_filt_target_filt %>%
    filter(library_sample_name %in% MAD4HATTER_diversity_sample_filt_target_filt_pass$library_sample_name) %>%
    group_by(library_sample_name, target_name) %>%
    filter(reads == max(reads)) %>%
    ungroup()


  ##Now the simulation framework will allow for samples to have missing target info which allows to simulate
  ##deletions but likely given the target chosen here, the missing targets are just targets that didn't get
  ##successfully genotyped, so we will fill in the missing targets. Can either just fill in with the most freq
  ##haplotypes or try to be more complex and fill in with the most highly correlated haplotype with the other
  ##haplotypes. For now let's just fill in with just the most common genotype at that position

  MAD4HATTER_diversity_final_samples_hap_counts = MAD4HATTER_diversity_final_samples %>%
    group_by(target_name, seq) %>%
    summarise(n = n()) %>%
    group_by(target_name) %>%
    mutate(total = sum(n)) %>%
    mutate(freq = n/total)

  MAD4HATTER_diversity_final_samples_hap_counts_most_common = MAD4HATTER_diversity_final_samples_hap_counts %>%
    group_by(target_name) %>%
    slice_max(freq,with_ties = F)

  missing_targets_per_sample = as_tibble(expand.grid(unique(MAD4HATTER_diversity_final_samples$target_name),
                                                     unique(MAD4HATTER_diversity_final_samples$library_sample_name))) %>%
    dplyr::rename(target_name = Var1,
                  library_sample_name = Var2) %>%
    filter(!paste0(library_sample_name, "-", target_name) %in% paste0(MAD4HATTER_diversity_final_samples$library_sample_name, "-",
                                                                      MAD4HATTER_diversity_final_samples$target_name))

  missing_targets_per_sample = missing_targets_per_sample %>%
    left_join(MAD4HATTER_diversity_final_samples_hap_counts_most_common %>%
                ungroup() %>%
                select(target_name, seq)) %>%
    mutate(reads = 10)

  MAD4HATTER_diversity_final_samples_filed = MAD4HATTER_diversity_final_samples %>%
    bind_rows(MAD4HATTER_diversity_final_samples,
              missing_targets_per_sample)

  return(MAD4HATTER_diversity_final_samples_filed)
}


###########
##### function to filter only targets we kept during strain choosing and intersect with sim pop output
###########

filter_and_intersect <- function(sim_data, filter_good_coverage_dom_haplo_data, country){

  # filter MAD4HATTER info to just the targets we kept
  MAD4HATTER_Pf3D7_primers_diversity_targets = raw_MAD4HATTER_Pf3D7_primers_diversity_targets %>%
    filter(target %in% filter_good_coverage_dom_haplo_data$target_name)

  MAD4HATTER_Pf3D7_region_info_diversity_targets_locs_for_sim = raw_MAD4HATTER_Pf3D7_region_info_diversity_targets_locs_for_sim %>%
    filter(target %in% filter_good_coverage_dom_haplo_data$target_name)

  # add ancestral index to the allele table
  MAD4HATTER_final_samples_filed_for_sim = filter_good_coverage_dom_haplo_data %>%
    ungroup() %>%
    select(library_sample_name, target_name, seq) %>%
    dplyr::rename(ancestral_genotype = library_sample_name,
                  target = target_name)

  MAD4HATTER_final_samples_filed_for_sim2 <- MAD4HATTER_final_samples_filed_for_sim[!duplicated(MAD4HATTER_final_samples_filed_for_sim), ]

  MAD4HATTER_diversity_targets_population = recombuddy::intersect_panel_with_simulated_population(MAD4HATTER_Pf3D7_region_info_diversity_targets_locs_for_sim,
                                                                                                  sim_data)

  MAD4HATTER_diversity_targets_population_with_allele_data = MAD4HATTER_diversity_targets_population %>%
    left_join(MAD4HATTER_final_samples_filed_for_sim2, by = c("ancestral_genotype", "target"))

  MAD4HATTER_diversity_targets_population_with_allele_data <- MAD4HATTER_diversity_targets_population_with_allele_data %>%
    mutate(ancestral_genotype = paste0(country, "-", ancestral_genotype),
           simulated_sample = paste0(country, "-sim-", simulated_sample))

  return(MAD4HATTER_diversity_targets_population_with_allele_data)
}

############
############


###########
##### function to filter only targets we kept during strain choosing and intersect with sim pop output for real data
###########

filter_and_intersect_real_data <- function(sim_data, filter_good_coverage_dom_haplo_data, country){

  # filter MAD4HATTER info to just the targets we kept

  MAD4HATTER_Pf3D7_region_info_diversity_targets_locs_for_sim = filter_good_coverage_dom_haplo_data %>%
    mutate(chrom = str_extract(target_name, "^[^-]+"),
           start = as.numeric(str_extract(target_name, "(?<=-)[0-9]+(?=-)")),
           end = as.numeric(str_extract(target_name, "(?<=-)[0-9]+(?=[^0-9]*$)"))) %>%
    dplyr::select(chrom, start, end, target = target_name)

  # add ancestral index to the allele table
  MAD4HATTER_final_samples_filed_for_sim = filter_good_coverage_dom_haplo_data %>%
    ungroup() %>%
    select(library_sample_name, target_name, seq) %>%
    dplyr::rename(ancestral_genotype = library_sample_name,
                  target = target_name)

  MAD4HATTER_final_samples_filed_for_sim2 <- MAD4HATTER_final_samples_filed_for_sim[!duplicated(MAD4HATTER_final_samples_filed_for_sim), ]

  MAD4HATTER_diversity_targets_population = recombuddy::intersect_panel_with_simulated_population(MAD4HATTER_Pf3D7_region_info_diversity_targets_locs_for_sim,
                                                                                                  sim_data)

  MAD4HATTER_diversity_targets_population_with_allele_data = MAD4HATTER_diversity_targets_population %>%
    left_join(MAD4HATTER_final_samples_filed_for_sim2, by = c("ancestral_genotype", "target"))

  MAD4HATTER_diversity_targets_population_with_allele_data <- MAD4HATTER_diversity_targets_population_with_allele_data %>%
    mutate(ancestral_genotype = paste0(country, "-", ancestral_genotype),
           simulated_sample = paste0(country, "-sim-", simulated_sample))

  return(MAD4HATTER_diversity_targets_population_with_allele_data)
}

############
############




###################
###################
######## simple function to compare ibd between two pops, no mixing/imports
###################
###################


compare_pops_no_mixing_no_imports <- function(ancestors_pop1, params_pop1, size_pop1,
                                              ancestors_pop2, params_pop2, size_pop2,
                                              real_data = F){

  #### first pop

  first_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop1, real_data = T)

  sim_pop1 = do.call(sim_population, c(list(unique(first_pop_filtered$library_sample_name)),size_pop1,
                                       as.list(params_pop1)))

  if(real_data == F){
    pop1_intersected <- filter_and_intersect(sim_pop1, first_pop_filtered, "Population_1")
  }else if(real_data == T){
    pop1_intersected <- filter_and_intersect_real_data(sim_pop1, first_pop_filtered, "Population_1")
  }



  ### second pop

  second_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop2, real_data = T)

  sim_pop2 = do.call(sim_population, c(list(unique(second_pop_filtered$library_sample_name)),size_pop2,
                                       as.list(params_pop2)))

  if(real_data == F){
    pop2_intersected <- filter_and_intersect(sim_pop2, second_pop_filtered, "Population_2")
  }else if(real_data == T){
    pop2_intersected <- filter_and_intersect_real_data(sim_pop2, second_pop_filtered, "Population_2")
  }


  ### now combine pops

  combined_pops_interesected <- rbind(pop1_intersected, pop2_intersected) %>%
    mutate(population = substr(simulated_sample, 1, 12))


  ### seperate out pop1 and pop2/import so that can calculate ibd for same pop pairs
  pop_1_only_intersected <- filter(combined_pops_interesected, population == "Population_1")
  pop_2_only_intersected <- filter(combined_pops_interesected, population == "Population_2")


  combined_pops_dsmp <- formatDat(combined_pops_interesected, svar = "simulated_sample", lvar = "target", avar = "seq")
  combined_pops_coi   <- getCOI(combined_pops_dsmp, lrank = 2)
  combined_pops_afreq <- calcAfreq(combined_pops_dsmp, combined_pops_coi, tol = 1e-5)
  combined_pops_dcifer <- ibdDat(combined_pops_dsmp, combined_pops_coi, combined_pops_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                                 alpha = 0.05, nr = 1e3)

  combined_pop_ibd_table_diff_res_only <- as.data.frame.table(combined_pops_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0)) %>%
    filter(same_residence == 0)

  ### do the same for pop1 and pop2 alone
  pop1_dsmp <- formatDat(pop_1_only_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop1_coi   <- getCOI(pop1_dsmp, lrank = 2)
  pop1_afreq <- calcAfreq(pop1_dsmp, pop1_coi, tol = 1e-5)
  pop1_dcifer <- ibdDat(pop1_dsmp, pop1_coi, pop1_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                        alpha = 0.05, nr = 1e3)

  pop1_ibd_table <- as.data.frame.table(pop1_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))

  pop2_dsmp <- formatDat(pop_2_only_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop2_coi   <- getCOI(pop2_dsmp, lrank = 2)
  pop2_afreq <- calcAfreq(pop2_dsmp, pop2_coi, tol = 1e-5)
  pop2_dcifer <- ibdDat(pop2_dsmp, pop2_coi, pop2_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                        alpha = 0.05, nr = 1e3)

  pop2_ibd_table <- as.data.frame.table(pop2_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))


  all_ibd_tables <- rbind(combined_pop_ibd_table_diff_res_only, pop1_ibd_table, pop2_ibd_table)


  JD <- calculate_jost_d(pop1_afreq, pop2_afreq)
  FST <- calculate_FST(pop1_afreq, pop2_afreq)
  allele_freq_plot <- plot_af(pop1_afreq, pop2_afreq, 700)

  return(list(all_ibd_tables, JD, FST, allele_freq_plot))

}




###########
#### plot allele frequencies
###########

plot_af <- function(pop1_afreq, pop2_afreq, rows_to_plot){

  ### convert to dataframe for better manipulation
  pop1_afreq_df_values <- as.data.frame(unlist(pop1_afreq))[,1]
  pop1_afreq_df_names <- names(unlist(pop1_afreq))
  pop1_afreq_df_locus <- sub("\\..*", "", pop1_afreq_df_names)
  pop1_afreq_df_allele <- str_replace(pop1_afreq_df_names, "^[^.]*.", "")

  pop1_afreq_df <- data.frame(locus = pop1_afreq_df_locus,
                               allele = pop1_afreq_df_allele,
                               frequency = pop1_afreq_df_values,
                               country = "pop1")

  pop2_afreq_df_values <- as.data.frame(unlist(pop2_afreq))[,1]
  pop2_afreq_df_names <- names(unlist(pop2_afreq))
  pop2_afreq_df_locus <- sub("\\..*", "", pop2_afreq_df_names)
  pop2_afreq_df_allele <- str_replace(pop2_afreq_df_names, "^[^.]*.", "")


  pop2_afreq_df <- data.frame(locus = pop2_afreq_df_locus,
                                  allele = pop2_afreq_df_allele,
                                  frequency = pop2_afreq_df_values,
                                  country = "pop2")


  pop1_pop2_afreq_df <- rbind(pop2_afreq_df, pop1_afreq_df) %>%
    arrange(locus)


  ggplot((pop1_pop2_afreq_df[1:rows_to_plot,])) +
    geom_col(aes(x = allele, y = frequency, group = country, fill = country), position = "dodge") +
    facet_wrap(~locus, scales = "free_x", ncol = 10) +
    theme_minimal() +
    theme(axis.text.x = element_blank())
}




###############
###############
##### function to import from pop1 to pop2
###############
###############

imports_no_mixing_function <- function(ancestors_pop1, params_pop1, size_pop1,
                                       ancestors_pop2, params_pop2, size_pop2,
                                       num_imports, real_data = F){

  #### first pop

  first_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop1)

  sim_pop1 = do.call(sim_population, c(list(unique(first_pop_filtered$library_sample_name)),size_pop1,
                                       as.list(params_pop1)))

  if(real_data == F){
    pop1_intersected <- filter_and_intersect(sim_pop1, first_pop_filtered, "Population_1")
  }else if(real_data == T){
    pop1_intersected <- filter_and_intersect_real_data(sim_pop1, first_pop_filtered, "Population_1")
  }


  ### second pop

  second_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop2)

  sim_pop2 = do.call(sim_population, c(list(unique(second_pop_filtered$library_sample_name)),size_pop2,
                                       as.list(params_pop2)))

  if(real_data == F){
    pop2_intersected <- filter_and_intersect(sim_pop2, second_pop_filtered, "Population_2")
  }else if(real_data == T){
    pop2_intersected <- filter_and_intersect_real_data(sim_pop2, second_pop_filtered, "Population_2")
  }


  ### now import some number of cases from pop 1 into pop 2

  ##### now not combining all but all from cambodia and just 10 samples from ghana
  pop1_sampled_samples <- sample(unique(pop1_intersected$simulated_sample), num_imports)

  pop1_intersected_sampled <- pop1_intersected %>% filter(simulated_sample %in% pop1_sampled_samples) %>%
    mutate(simulated_sample = paste0("Import_", simulated_sample))

  pop1_intersected_not_sampled <- filter(pop1_intersected, !simulated_sample %in% pop1_sampled_samples)

  combined_pops_interesected <- rbind(pop1_intersected_sampled, pop1_intersected_not_sampled, pop2_intersected) %>%
    mutate(population = substr(simulated_sample, 1, 12))


  ### seperate out pop1 and pop2/import so that can calculate ibd for same pop pairs
  pop_1_only_intersected <- filter(combined_pops_interesected, population == "Population_1")
  pop_2_only_intersected <- filter(combined_pops_interesected, population %in% c("Population_2", "Import_Popul"))


  combined_pops_dsmp <- formatDat(combined_pops_interesected, svar = "simulated_sample", lvar = "target", avar = "seq")
  combined_pops_coi   <- getCOI(combined_pops_dsmp, lrank = 2)
  combined_pops_afreq <- calcAfreq(combined_pops_dsmp, combined_pops_coi, tol = 1e-5)
  combined_pops_dcifer <- ibdDat(combined_pops_dsmp, combined_pops_coi, combined_pops_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                                 alpha = 0.05, nr = 1e3)

  combined_pop_ibd_table_diff_res_only <- as.data.frame.table(combined_pops_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0)) %>%
    filter(same_residence == 0)

  ### do the same for pop1 and pop2 alone
  pop1_dsmp <- formatDat(pop_1_only_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop1_coi   <- getCOI(pop1_dsmp, lrank = 2)
  pop1_afreq <- calcAfreq(pop1_dsmp, pop1_coi, tol = 1e-5)
  pop1_dcifer <- ibdDat(pop1_dsmp, pop1_coi, pop1_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                        alpha = 0.05, nr = 1e3)

  pop1_ibd_table <- as.data.frame.table(pop1_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))

  pop2_dsmp <- formatDat(pop_2_only_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop2_coi   <- getCOI(pop2_dsmp, lrank = 2)
  pop2_afreq <- calcAfreq(pop2_dsmp, pop2_coi, tol = 1e-5)
  pop2_dcifer <- ibdDat(pop2_dsmp, pop2_coi, pop2_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                        alpha = 0.05, nr = 1e3)

  pop2_ibd_table <- as.data.frame.table(pop2_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))


  all_ibd_tables <- rbind(combined_pop_ibd_table_diff_res_only, pop1_ibd_table, pop2_ibd_table)

  all_ibd_tables <- all_ibd_tables %>%
    bind_rows(
      all_ibd_tables %>% rename(specimen_a = specimen_b, specimen_b = specimen_a,
                                population_a = population_b, population_b = population_a)  # flipped version
    )

  JD <- calculate_jost_d(pop1_afreq, pop2_afreq)
  FST <- calculate_FST(pop1_afreq, pop2_afreq)

  he1 <- calculate_pop_expected_het(pop1_afreq)
  he2 <- calculate_pop_expected_het(pop2_afreq)

  allele_freq_plot <- plot_af(pop1_afreq, pop2_afreq, 700)


  return(list(all_ibd_tables, JD, FST, he1, he2, allele_freq_plot))

}



#########
#### calculate KL
#########

calculate_KL <- function(table){
  table %>%
    group_by(specimen_a) %>%
    group_modify(~{

      # Create common breaks per id
      breaks <- seq(min(.x$estimate),
                    max(.x$estimate),
                    length.out = 100)

      # Compute histogram per category
      hists <- .x %>%
        group_by(same_residence) %>%
        summarise(
          probs = list(
            hist(estimate, breaks = breaks, plot = FALSE)$counts
          ),
          .groups = "drop"
        )

      # Convert counts to probabilities
      hists$probs <- map(hists$probs, ~ .x / sum(.x))

      # Assume two categories
      p <- hists$probs[[1]]
      q <- hists$probs[[2]]

      tibble(
        KL_pq = KL(rbind(p, q)),
        KL_qp = KL(rbind(q, p)),
        wasserstein1d(p, q)
      )
    }) %>%
    ungroup()

}





#############
####### simulation to weight different pop allele frequencies
#############
# Safe blending function that handles mismatched allele sets
blend_af_no_zeroes <- function(af1, af2, w) {
  # Get union of locus names
  all_loci <- union(names(af1), names(af2))

  blended <- lapply(all_loci, function(loc) {
    a <- af1[[loc]]
    b <- af2[[loc]]

    # Get union of allele names at this locus
    all_alleles <- union(names(a), names(b))

    # Expand both vectors to the full allele set, filling missing with 0
    a_full <- setNames(rep(0, length(all_alleles)), all_alleles)
    b_full <- setNames(rep(0, length(all_alleles)), all_alleles)

    a_full[names(a)] <- a
    b_full[names(b)] <- b

    ##weight by observed copies at this locus
    w1 <- n1[loc]
    w2 <- n1[loc]

    # Blend
    blended_locus <- (1 - w) * a_full + w * b_full

    # remvoe alleles with zero
    blended_locus <- blended_locus[blended_locus > 0]

    # Renormalize to sum to 1 (in case of floating point drift)
    blended_locus / sum(blended_locus)
  })

  names(blended) <- all_loci
  return(blended)
}

blend_af <- function(af1, af2, w) {
  # Get union of locus names
  all_loci <- union(names(af1), names(af2))

  blended <- lapply(all_loci, function(loc) {
    a <- af1[[loc]]
    b <- af2[[loc]]

    # Get union of allele names at this locus
    all_alleles <- union(names(a), names(b))

    # Expand both vectors to the full allele set, filling missing with 0
    a_full <- setNames(rep(0, length(all_alleles)), all_alleles)
    b_full <- setNames(rep(0, length(all_alleles)), all_alleles)

    a_full[names(a)] <- a
    b_full[names(b)] <- b

    # Blend
    blended_locus <- (1 - w) * a_full + w * b_full

    # Renormalize to sum to 1 (in case of floating point drift)
    blended_locus / sum(blended_locus)
  })

  names(blended) <- all_loci
  return(blended)
}





#### function to compare IBD across two pops under different allele frequencies
compare_pops_no_mixing_no_imports_different_af_weights <- function(ancestors_pop1, params_pop1, size_pop1,
                                                                   ancestors_pop2, params_pop2, size_pop2, weights){



  #### first pop

  first_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop1)

  sim_pop1 = do.call(sim_population, c(list(unique(first_pop_filtered$library_sample_name)),size_pop1,
                                       as.list(params_pop1)))

  pop1_intersected <- filter_and_intersect(sim_pop1, first_pop_filtered, "Population_1")



  ### second pop

  second_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop2)

  sim_pop2 = do.call(sim_population, c(list(unique(second_pop_filtered$library_sample_name)),size_pop2,
                                       as.list(params_pop2)))

  pop2_intersected <- filter_and_intersect(sim_pop2, second_pop_filtered, "Population_2")


  ### now combine pops

  combined_pops_interesected <- rbind(pop1_intersected, pop2_intersected) %>%
    mutate(population = substr(simulated_sample, 1, 12))


  combined_pops_dsmp <- formatDat(combined_pops_interesected, svar = "simulated_sample", lvar = "target", avar = "seq")
  combined_pops_coi   <- getCOI(combined_pops_dsmp, lrank = 2)

  pop1_dsmp <- formatDat(pop1_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop1_coi   <- getCOI(pop1_dsmp, lrank = 2)
  pop1_afreq <- calcAfreq(pop1_dsmp, pop1_coi, tol = 1e-5)

  pop2_dsmp <- formatDat(pop2_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop2_coi   <- getCOI(pop2_dsmp, lrank = 2)
  pop2_afreq <- calcAfreq(pop2_dsmp, pop2_coi, tol = 1e-5)


  full_results <- data.frame(
    weight = as.numeric(),
    specimen_a = as.character(),
    speciemn_b = as.character(),
    estimate = as.numeric(),
    p_value = as.numeric(),
    CI_lower = as.numeric(),
    CI_upper = as.numeric(),
    population_a = as.character(),
    population_b = as.character(),
    same_pop = as.numeric(),
    same_residence = as.numeric()
  )


  for(i in seq_along(weights)){
    w <- weights[i]

    af_blended <- blend_af(pop1_afreq, pop2_afreq, w = w)

    blended_dcifer <- ibdDat(combined_pops_dsmp, combined_pops_coi, af_blended, pval = TRUE, confint = TRUE, rnull = 0,
                             alpha = 0.05, nr = 1e3)

    combined_pop_ibd_table_diff_res_only <- as.data.frame.table(blended_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0)) %>%
      filter(same_residence == 0)

    af_blended_no_zeroes <- blend_af_no_zeroes(pop1_afreq, pop2_afreq, w = w)

    ### do the same for pop1 and pop2 alone
    if(w == 0){
      af_blended_reordered_for_pop1 <- mapply(function(blended, ref) {
        blended[names(ref)]  # reorder to match ref's allele order
      }, af_blended_no_zeroes, pop1_afreq[names(af_blended_no_zeroes)], SIMPLIFY = FALSE)
    }else if(w == 1){
      af_blended_reordered_for_pop1 <- af_blended
      }else{
      af_blended_reordered_for_pop1 = af_blended_no_zeroes
    }

    pop1_dcifer <- ibdDat(pop1_dsmp, pop1_coi, af_blended_reordered_for_pop1, pval = TRUE, confint = TRUE, rnull = 0,
                          alpha = 0.05, nr = 1e3)

    pop1_ibd_table <- as.data.frame.table(pop1_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))

    if(w == 1){
      af_blended_reordered_for_pop2 <- mapply(function(blended, ref) {
        blended[names(ref)]  # reorder to match ref's allele order
      }, af_blended_no_zeroes, pop2_afreq[names(af_blended_no_zeroes)], SIMPLIFY = FALSE)
    }else if(w == 0){
      af_blended_reordered_for_pop2 <- af_blended
    }else{
      af_blended_reordered_for_pop2 <- mapply(function(blended, ref) {
        blended[names(ref)]  # reorder to match ref's allele order
      }, af_blended_no_zeroes, pop2_afreq[names(af_blended_no_zeroes)], SIMPLIFY = FALSE)
    }

    pop2_dcifer <- ibdDat(pop2_dsmp, pop2_coi, af_blended_reordered_for_pop2, pval = TRUE, confint = TRUE, rnull = 0,
                          alpha = 0.05, nr = 1e3)

    pop2_ibd_table <- as.data.frame.table(pop2_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))


    all_ibd_tables <- rbind(combined_pop_ibd_table_diff_res_only, pop1_ibd_table, pop2_ibd_table) %>%
      mutate(weight = w) %>% relocate(weight)

    all_ibd_tables <- all_ibd_tables %>%
      bind_rows(
        all_ibd_tables %>% rename(specimen_a = specimen_b, specimen_b = specimen_a,
                                  population_a = population_b, population_b = population_a)  # flipped version
      )

    full_results <- rbind(full_results, all_ibd_tables)
  }


  JD <- calculate_jost_d(pop1_afreq, pop2_afreq)
  FST <- calculate_FST(pop1_afreq, pop2_afreq)

  return(list(full_results, JD, FST))

}



blend_af <- function(af1, af2, w, n1 = NULL, n2 = NULL) {
  all_loci <- union(names(af1), names(af2))

  blended <- lapply(all_loci, function(loc) {
    a <- af1[[loc]]
    b <- af2[[loc]]

    all_alleles <- union(names(a), names(b))

    a_full <- setNames(rep(0, length(all_alleles)), all_alleles)
    b_full <- setNames(rep(0, length(all_alleles)), all_alleles)

    a_full[names(a)] <- a
    b_full[names(b)] <- b

    # Use n1/n2 if provided, otherwise fall back to w as simple weight
    if (!is.null(n1) && !is.null(n2)) {
      w1 <- n1[loc] * (1 - w)
      w2 <- n2[loc] * w
    } else {
      w1 <- (1 - w)
      w2 <- w
    }

    blended_locus <- (w1 * a_full + w2 * b_full) / (w1 + w2)

    # Drop zeros and renormalize
    blended_locus <- blended_locus[blended_locus > 1e-10]
    blended_locus / sum(blended_locus)
  })

  names(blended) <- all_loci
  return(blended)
}




#### function to compare IBD across two pops under different allele frequencies
compare_pops_no_mixing_no_imports_different_af_weights <- function(ancestors_pop1, params_pop1, size_pop1,
                                                                   ancestors_pop2, params_pop2, size_pop2, weights){



  #### first pop

  first_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop1)

  sim_pop1 = do.call(sim_population, c(list(unique(first_pop_filtered$library_sample_name)),size_pop1,
                                       as.list(params_pop1)))

  pop1_intersected <- filter_and_intersect(sim_pop1, first_pop_filtered, "Population_1")



  ### second pop

  second_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop2)

  sim_pop2 = do.call(sim_population, c(list(unique(second_pop_filtered$library_sample_name)),size_pop2,
                                       as.list(params_pop2)))

  pop2_intersected <- filter_and_intersect(sim_pop2, second_pop_filtered, "Population_2")


  ### now combine pops

  combined_pops_interesected <- rbind(pop1_intersected, pop2_intersected) %>%
    mutate(population = substr(simulated_sample, 1, 12))


  combined_pops_dsmp <- formatDat(combined_pops_interesected, svar = "simulated_sample", lvar = "target", avar = "seq")
  combined_pops_coi   <- getCOI(combined_pops_dsmp, lrank = 2)

  pop1_dsmp <- formatDat(pop1_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop1_coi   <- getCOI(pop1_dsmp, lrank = 2)
  pop1_afreq <- calcAfreq(pop1_dsmp, pop1_coi, tol = 1e-5)

  pop2_dsmp <- formatDat(pop2_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop2_coi   <- getCOI(pop2_dsmp, lrank = 2)
  pop2_afreq <- calcAfreq(pop2_dsmp, pop2_coi, tol = 1e-5)

  pop1_n <- pop1_intersected %>%
    group_by(target, simulated_sample) %>%
    summarise(n_alleles = n_distinct(seq), .groups = "drop") %>%
    left_join(
      data.frame(simulated_sample = names(pop1_coi), coi = pop1_coi),
      by = "simulated_sample"
    ) %>%
    group_by(target) %>%
    summarise(total = sum(coi / n_alleles)) %>%
    { setNames(.$total, .$target) }

  pop2_n <- pop2_intersected %>%
    group_by(target, simulated_sample) %>%
    summarise(n_alleles = n_distinct(seq), .groups = "drop") %>%
    left_join(
      data.frame(simulated_sample = names(pop2_coi), coi = pop2_coi),
      by = "simulated_sample"
    ) %>%
    group_by(target) %>%
    summarise(total = sum(coi / n_alleles)) %>%
    { setNames(.$total, .$target) }


  full_results <- data.frame(
    weight = as.numeric(),
    specimen_a = as.character(),
    speciemn_b = as.character(),
    estimate = as.numeric(),
    p_value = as.numeric(),
    CI_lower = as.numeric(),
    CI_upper = as.numeric(),
    population_a = as.character(),
    population_b = as.character(),
    same_pop = as.numeric(),
    same_residence = as.numeric()
  )


  for(i in seq_along(weights)){
    w <- weights[i]

    af_blended <- blend_af(pop1_afreq, pop2_afreq, w = w)

    combined_rematch <- matchAfreq(combined_pops_dsmp, af_blended)
    combined_rematch_dsmp <- combined_rematch$dsmp
    combined_rematch_afreq <- combined_rematch$afreq

    blended_dcifer <- ibdDat(combined_rematch_dsmp, combined_pops_coi, combined_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                             alpha = 0.05, nr = 1e3)

    combined_pop_ibd_table_diff_res_only <- as.data.frame.table(blended_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0)) %>%
      filter(same_residence == 0)

    if(w == 0){
      pop1_rematch_dsmp <- pop1_dsmp
      pop1_rematch_afreq <- pop1_afreq
    }else{
      pop1_rematch <- matchAfreq(pop1_dsmp, af_blended)
      pop1_rematch_dsmp <- pop1_rematch$dsmp
      pop1_rematch_afreq <- pop1_rematch$afreq
    }

    pop1_dcifer <- ibdDat(pop1_rematch_dsmp, pop1_coi, pop1_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                          alpha = 0.05, nr = 1e3)

    pop1_ibd_table <- as.data.frame.table(pop1_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))

    if(w == 1){
      pop2_rematch_dsmp <- pop2_dsmp
      pop2_rematch_afreq <- pop2_afreq
    }else{
      pop2_rematch <- matchAfreq(pop2_dsmp, af_blended)
      pop2_rematch_dsmp <- pop2_rematch$dsmp
      pop2_rematch_afreq <- pop2_rematch$afreq
    }

    pop2_dcifer <- ibdDat(pop2_rematch_dsmp, pop2_coi, pop2_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                          alpha = 0.05, nr = 1e3)

    pop2_ibd_table <- as.data.frame.table(pop2_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))


    all_ibd_tables <- rbind(combined_pop_ibd_table_diff_res_only, pop1_ibd_table, pop2_ibd_table) %>%
      mutate(weight = w) %>% relocate(weight)

    all_ibd_tables <- all_ibd_tables %>%
      bind_rows(
        all_ibd_tables %>% rename(specimen_a = specimen_b, specimen_b = specimen_a,
                                  population_a = population_b, population_b = population_a)  # flipped version
      )

    full_results <- rbind(full_results, all_ibd_tables)
  }


  JD <- calculate_jost_d(pop1_afreq, pop2_afreq)
  FST <- calculate_FST(pop1_afreq, pop2_afreq)

  return(list(full_results, JD, FST))

}









#############
####### check against different weights for allele frequencies in an importation simulation
#############

###############
###############
##### function to import from pop1 to pop2
###############
###############

different_allele_weights_imports_no_mixing_function <- function(ancestors_pop1, params_pop1, size_pop1,
                                       ancestors_pop2, params_pop2, size_pop2,
                                       num_imports, weights){

  #### first pop

  first_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop1)

  sim_pop1 = do.call(sim_population, c(list(unique(first_pop_filtered$library_sample_name)),size_pop1,
                                       as.list(params_pop1)))

  pop1_intersected <- filter_and_intersect(sim_pop1, first_pop_filtered, "Population_1")


  ### second pop

  second_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop2)

  sim_pop2 = do.call(sim_population, c(list(unique(second_pop_filtered$library_sample_name)),size_pop2,
                                       as.list(params_pop2)))

  pop2_intersected <- filter_and_intersect(sim_pop2, second_pop_filtered, "Population_2")


  ### now import some number of cases from pop 1 into pop 2

  ##### now not combining all but all from cambodia and just 10 samples from ghana
  pop1_sampled_samples <- sample(unique(pop1_intersected$simulated_sample), num_imports)

  pop1_intersected_sampled <- pop1_intersected %>% filter(simulated_sample %in% pop1_sampled_samples) %>%
    mutate(simulated_sample = paste0("Import_", simulated_sample))

  pop1_intersected_not_sampled <- filter(pop1_intersected, !simulated_sample %in% pop1_sampled_samples)

  combined_pops_interesected <- rbind(pop1_intersected_sampled, pop1_intersected_not_sampled, pop2_intersected) %>%
    mutate(population = substr(simulated_sample, 1, 12))


  ### seperate out pop1 and pop2/import so that can calculate ibd for same pop pairs
  pop1_only_intersected <- filter(combined_pops_interesected, population == "Population_1")
  pop2_only_intersected <- filter(combined_pops_interesected, population == "Population_2")
  pop1_and_import_intersected <- filter(combined_pops_interesected, population %in% c("Population_1", "Import_Popul"))
  pop2_and_import_intersected <- filter(combined_pops_interesected, population %in% c("Population_2", "Import_Popul"))


  combined_pops_dsmp <- formatDat(combined_pops_interesected, svar = "simulated_sample", lvar = "target", avar = "seq")
  combined_pops_coi   <- getCOI(combined_pops_dsmp, lrank = 2)

  pop1_only_dsmp <- formatDat(pop1_only_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop1_only_coi   <- getCOI(pop1_only_dsmp, lrank = 2)
  pop1_only_afreq <- calcAfreq(pop1_only_dsmp, pop1_only_coi, tol = 1e-5)

  pop2_only_dsmp <- formatDat(pop2_only_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop2_only_coi   <- getCOI(pop2_only_dsmp, lrank = 2)
  pop2_only_afreq <- calcAfreq(pop2_only_dsmp, pop2_only_coi, tol = 1e-5)

  pop1_import_dsmp <- formatDat(pop1_and_import_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop1_import_coi   <- getCOI(pop1_import_dsmp, lrank = 2)
  pop1_import_afreq <- calcAfreq(pop1_import_dsmp, pop1_import_coi, tol = 1e-5)

  pop2_import_dsmp <- formatDat(pop2_and_import_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop2_import_coi   <- getCOI(pop2_import_dsmp, lrank = 2)
  pop2_import_afreq <- calcAfreq(pop2_import_dsmp, pop2_import_coi, tol = 1e-5)

  pop1_import_n <- pop1_and_import_intersected %>%
    group_by(target, simulated_sample) %>%
    summarise(n_alleles = n_distinct(seq), .groups = "drop") %>%
    left_join(
      data.frame(simulated_sample = names(pop1_coi), coi = pop1_coi),
      by = "simulated_sample"
    ) %>%
    group_by(target) %>%
    summarise(total = sum(coi / n_alleles)) %>%
    { setNames(.$total, .$target) }

  pop2_only_n <- pop2_only_intersected %>%
    group_by(target, simulated_sample) %>%
    summarise(n_alleles = n_distinct(seq), .groups = "drop") %>%
    left_join(
      data.frame(simulated_sample = names(pop2_coi), coi = pop2_coi),
      by = "simulated_sample"
    ) %>%
    group_by(target) %>%
    summarise(total = sum(coi / n_alleles)) %>%
    { setNames(.$total, .$target) }


  full_results <- data.frame(
    weight = as.numeric(),
    specimen_a = as.character(),
    speciemn_b = as.character(),
    estimate = as.numeric(),
    p_value = as.numeric(),
    CI_lower = as.numeric(),
    CI_upper = as.numeric(),
    population_a = as.character(),
    population_b = as.character(),
    same_pop = as.numeric(),
    same_residence = as.numeric()
  )


  for(i in seq_along(weights)){
    w <- weights[i]

    af_blended <- blend_af(pop1_only_afreq, pop2_import_afreq, w = w)

    combined_rematch <- matchAfreq(combined_pops_dsmp, af_blended)
    combined_rematch_dsmp <- combined_rematch$dsmp
    combined_rematch_afreq <- combined_rematch$afreq

    blended_dcifer <- ibdDat(combined_rematch_dsmp, combined_pops_coi, combined_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                             alpha = 0.05, nr = 1e3)

    combined_pop_ibd_table_diff_res_only <- as.data.frame.table(blended_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0)) %>%
      filter(same_residence == 0)

    if(w == 0){
      pop1_rematch_dsmp <- pop1_only_dsmp
      pop1_rematch_afreq <- pop1_only_afreq
      pop1_rematch_COI <- pop1_only_coi
    }else{
      pop1_rematch <- matchAfreq(pop1_only_dsmp, af_blended)
      pop1_rematch_dsmp <- pop1_rematch$dsmp
      pop1_rematch_afreq <- pop1_rematch$afreq
      pop1_rematch_COI <- pop1_only_coi
    }

    pop1_dcifer <- ibdDat(pop1_rematch_dsmp, pop1_rematch_COI, pop1_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                          alpha = 0.05, nr = 1e3)

    pop1_ibd_table <- as.data.frame.table(pop1_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))

    if(w == 1){
      pop2_rematch_dsmp <- pop2_import_dsmp
      pop2_rematch_afreq <- pop2_import_afreq
      pop2_rematch_coi <- pop2_import_coi
    }else{
      pop2_rematch <- matchAfreq(pop2_import_dsmp, af_blended)
      pop2_rematch_dsmp <- pop2_rematch$dsmp
      pop2_rematch_afreq <- pop2_rematch$afreq
      pop2_rematch_coi <- pop2_import_coi
    }

    pop2_dcifer <- ibdDat(pop2_rematch_dsmp, pop2_rematch_coi, pop2_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                          alpha = 0.05, nr = 1e3)

    pop2_ibd_table <- as.data.frame.table(pop2_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 12),
             population_b = substr(specimen_b, 1, 12),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))


    all_ibd_tables <- rbind(combined_pop_ibd_table_diff_res_only, pop1_ibd_table, pop2_ibd_table) %>%
      mutate(weight = w) %>% relocate(weight)

    all_ibd_tables <- all_ibd_tables %>%
      bind_rows(
        all_ibd_tables %>% rename(specimen_a = specimen_b, specimen_b = specimen_a,
                                  population_a = population_b, population_b = population_a)  # flipped version
      )

    full_results <- rbind(full_results, all_ibd_tables)
    print(w)
  }


  JD <- calculate_jost_d(pop1_only_afreq, pop2_import_afreq)
  FST <- calculate_FST(pop1_only_afreq, pop2_import_afreq)

  allele_freq_plot <- plot_af(pop1_only_afreq, pop2_afreq, 700)

  return(list(full_results, JD, FST, allele_freq_plot))

}








###############
###############
##### function to import from pop1 to pop2 with a 3rd population being compared to
###############
###############

three_pops_imports_no_mixing_function <- function(ancestors_pop1, params_pop1, size_pop1,
                                       ancestors_pop2, params_pop2, size_pop2,
                                       num_imports,
                                       ancestors_pop3, params_pop3, size_pop3){

  #### first pop

  first_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop1)

  sim_pop1 = do.call(sim_population, c(list(unique(first_pop_filtered$library_sample_name)),size_pop1,
                                       as.list(params_pop1)))

  pop1_intersected <- filter_and_intersect(sim_pop1, first_pop_filtered, "Population_1")


  ### second pop

  second_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop2)

  sim_pop2 = do.call(sim_population, c(list(unique(second_pop_filtered$library_sample_name)),size_pop2,
                                       as.list(params_pop2)))

  pop2_intersected <- filter_and_intersect(sim_pop2, second_pop_filtered, "Population_2")


  ### third pop

  third_pop_filtered <- filter_good_coverage_dom_haplo(ancestors_pop3)

  sim_pop3 = do.call(sim_population, c(list(unique(third_pop_filtered$library_sample_name)),size_pop3,
                                       as.list(params_pop3)))

  pop3_intersected <- filter_and_intersect(sim_pop3, third_pop_filtered, "Population_3")


  ### now import some number of cases from pop 1 into pop 2

  ##### now not combining all but all from cambodia and just 10 samples from ghana
  pop1_sampled_samples <- sample(unique(pop1_intersected$simulated_sample), num_imports)

  pop1_intersected_sampled <- pop1_intersected %>% filter(simulated_sample %in% pop1_sampled_samples) %>%
    mutate(simulated_sample = paste0("Import_", simulated_sample))

  pop1_intersected_not_sampled <- filter(pop1_intersected, !simulated_sample %in% pop1_sampled_samples)

  combined_1_2_pops_intersected <- rbind(pop1_intersected_sampled, pop1_intersected_not_sampled, pop2_intersected) %>%
    mutate(population = substr(simulated_sample, 1, 12))


  ### seperate out pop1 and pop2/import so that can calculate ibd for same pop pairs
  pop_1_only_intersected <- filter(combined_1_2_pops_intersected, population == "Population_1")
  pop_2_only_intersected <- filter(combined_1_2_pops_intersected, population %in% c("Population_2", "Import_Popul"))


  #### now make a combined pop2/imports and pop3 dataset. we don't care about IBD for pop1 vs pop3
  combined_2_3_pops_intersected <- rbind(pop1_intersected_sampled, pop2_intersected, pop3_intersected) %>%
    mutate(population = substr(simulated_sample, 1, 12))


  combined_1_2_pops_dsmp <- formatDat(combined_1_2_pops_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  combined_1_2_pops_coi   <- getCOI(combined_1_2_pops_dsmp, lrank = 2)
  combined_1_2_pops_afreq <- calcAfreq(combined_1_2_pops_dsmp, combined_1_2_pops_coi, tol = 1e-5)
  combined_1_2_pops_dcifer <- ibdDat(combined_1_2_pops_dsmp, combined_1_2_pops_coi, combined_1_2_pops_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                                 alpha = 0.05, nr = 1e3)

  combined_1_2_pop_ibd_table_diff_res_only <- as.data.frame.table(combined_1_2_pops_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0)) %>%
    filter(same_residence == 0)

  ### get AF for pop1 and pop 3
  pop1_dsmp <- formatDat(pop_1_only_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop1_coi   <- getCOI(pop1_dsmp, lrank = 2)
  pop1_afreq <- calcAfreq(pop1_dsmp, pop1_coi, tol = 1e-5)

  ### get AF for pop1 and pop 3
  pop3_dsmp <- formatDat(pop3_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop3_coi   <- getCOI(pop3_dsmp, lrank = 2)
  pop3_afreq <- calcAfreq(pop3_dsmp, pop3_coi, tol = 1e-5)

  ### do the same for pop2 alone
  pop2_dsmp <- formatDat(pop_2_only_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop2_coi   <- getCOI(pop2_dsmp, lrank = 2)
  pop2_afreq <- calcAfreq(pop2_dsmp, pop2_coi, tol = 1e-5)
  pop2_dcifer <- ibdDat(pop2_dsmp, pop2_coi, pop2_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                        alpha = 0.05, nr = 1e3)

  pop2_ibd_table <- as.data.frame.table(pop2_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))

  ### do the same for pop2/import and pop3
  combined_2_3_pops_dsmp <- formatDat(combined_2_3_pops_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  combined_2_3_pops_coi   <- getCOI(combined_2_3_pops_dsmp, lrank = 2)
  combined_2_3_pops_afreq <- calcAfreq(combined_2_3_pops_dsmp, combined_2_3_pops_coi, tol = 1e-5)
  combined_2_3_pops_dcifer <- ibdDat(combined_2_3_pops_dsmp, combined_2_3_pops_coi, combined_2_3_pops_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                                     alpha = 0.05, nr = 1e3)

  combined_2_3_pop_ibd_table_diff_res_only <- as.data.frame.table(combined_2_3_pops_dcifer, responseName = "value") %>%
    rename(
      specimen_a = Var1,
      specimen_b = Var2,
      metric = Var3) %>%
    pivot_wider(
      names_from = metric,
      values_from = value) %>%
    filter(specimen_a != specimen_b, !is.na(estimate)) %>%
    mutate(population_a = substr(specimen_a, 1, 12),
           population_b = substr(specimen_b, 1, 12),
           same_pop = ifelse(population_a == population_b, 1, 0),
           same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                     (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0)) %>%
    filter(same_residence == 0)


  all_ibd_tables <- rbind(combined_1_2_pop_ibd_table_diff_res_only, pop2_ibd_table, combined_2_3_pop_ibd_table_diff_res_only)

  all_ibd_tables <- all_ibd_tables %>%
    bind_rows(
      all_ibd_tables %>% rename(specimen_a = specimen_b, specimen_b = specimen_a,
                                population_a = population_b, population_b = population_a)  # flipped version
    )

  JD_1_2 <- calculate_jost_d(pop1_afreq, pop2_afreq)
  FST_1_2 <- calculate_FST(pop1_afreq, pop2_afreq)


  JD_2_3 <- calculate_jost_d(pop2_afreq, pop3_afreq)
  FST_2_3 <- calculate_FST(pop2_afreq, pop3_afreq)

  he1 <- calculate_pop_expected_het(pop1_afreq)
  he2 <- calculate_pop_expected_het(pop2_afreq)

  allele_freq_plot_1_2 <- plot_af(pop1_afreq, pop2_afreq, 700)
  allele_freq_plot_2_3 <- plot_af(pop2_afreq, pop3_afreq, 700)


  return(list(all_ibd_tables, JD_1_2, FST_1_2, JD_2_3, FST_2_3, he1, he2, allele_freq_plot_1_2, allele_freq_plot_2_3))

}



blend_af <- function(af1, af2, w, n1 = NULL, n2 = NULL) {
  all_loci <- union(names(af1), names(af2))

  blended <- lapply(all_loci, function(loc) {
    loc[1]
    a <- af1[[loc]]
    b <- af2[[loc]]

    all_alleles <- union(names(a), names(b))

    a_full <- setNames(rep(0, length(all_alleles)), all_alleles)
    b_full <- setNames(rep(0, length(all_alleles)), all_alleles)

    a_full[names(a)] <- a
    b_full[names(b)] <- b
    a_full

    # Use n1/n2 if provided, otherwise fall back to w as simple weight
    if (!is.null(n1) && !is.null(n2)) {
      w1 <- n1[loc] * (1 - w)
      w2 <- n2[loc] * w
    } else {
      w1 <- (1 - w)
      w2 <- w
    }

    blended_locus <- (w1 * a_full + w2 * b_full) / (w1 + w2)

  })

  names(blended) <- all_loci
  return(blended)
}



#### function to compare IBD across two pops under different allele frequencies
compare_pops_no_mixing_no_imports_different_af_weights_real_data <- function(pop1_data,
                                                                   pop2_data, weights){

  pop1_data = pop1_data %>% dplyr::select(simulated_sample = specimen_id, target = target_id, seq)
  pop2_data = pop2_data %>% dplyr::select(simulated_sample = specimen_id, target = target_id, seq)

  combined_data <- rbind(pop1_data, pop2_data)

  combined_pops_dsmp <- formatDat(combined_data, svar = "simulated_sample", lvar = "target", avar = "seq")
  combined_pops_coi   <- getCOI(combined_pops_dsmp, lrank = 2)

  pop1_dsmp <- formatDat(pop1_data, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop1_coi   <- getCOI(pop1_dsmp, lrank = 2)
  pop1_afreq <- calcAfreq(pop1_dsmp, pop1_coi, tol = 1e-5)

  pop2_dsmp <- formatDat(pop2_data, svar = "simulated_sample", lvar = "target", avar = "seq")
  pop2_coi   <- getCOI(pop2_dsmp, lrank = 2)
  pop2_afreq <- calcAfreq(pop2_dsmp, pop2_coi, tol = 1e-5)


  full_results <- data.frame(
    weight = as.numeric(),
    specimen_a = as.character(),
    speciemn_b = as.character(),
    estimate = as.numeric(),
    p_value = as.numeric(),
    CI_lower = as.numeric(),
    CI_upper = as.numeric(),
    population_a = as.character(),
    population_b = as.character(),
    same_pop = as.numeric(),
    same_residence = as.numeric()
  )


  for(i in seq_along(weights)){
    w <- weights[i]

    af_blended <- blend_af(pop1_afreq, pop2_afreq, w = w)

    combined_rematch <- matchAfreq(combined_pops_dsmp, af_blended)
    combined_rematch_dsmp <- combined_rematch$dsmp
    combined_rematch_afreq <- combined_rematch$afreq

    blended_dcifer <- ibdDat(combined_rematch_dsmp, combined_pops_coi, combined_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                             alpha = 0.05, nr = 1e3)

    combined_pop_ibd_table_diff_res_only <- as.data.frame.table(blended_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 2),
             population_b = substr(specimen_b, 1, 2),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0)) %>%
      filter(same_residence == 0)

    if(w == 0){
      pop1_rematch_dsmp <- pop1_dsmp
      pop1_rematch_afreq <- pop1_afreq
    }else{
      pop1_rematch <- matchAfreq(pop1_dsmp, af_blended)
      pop1_rematch_dsmp <- pop1_rematch$dsmp
      pop1_rematch_afreq <- pop1_rematch$afreq
    }

    pop1_dcifer <- ibdDat(pop1_rematch_dsmp, pop1_coi, pop1_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                          alpha = 0.05, nr = 1e3)

    pop1_ibd_table <- as.data.frame.table(pop1_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 2),
             population_b = substr(specimen_b, 1, 2),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))

    if(w == 1){
      pop2_rematch_dsmp <- pop2_dsmp
      pop2_rematch_afreq <- pop2_afreq
    }else{
      pop2_rematch <- matchAfreq(pop2_dsmp, af_blended)
      pop2_rematch_dsmp <- pop2_rematch$dsmp
      pop2_rematch_afreq <- pop2_rematch$afreq
    }

    pop2_dcifer <- ibdDat(pop2_rematch_dsmp, pop2_coi, pop2_rematch_afreq, pval = TRUE, confint = TRUE, rnull = 0,
                          alpha = 0.05, nr = 1e3)

    pop2_ibd_table <- as.data.frame.table(pop2_dcifer, responseName = "value") %>%
      rename(
        specimen_a = Var1,
        specimen_b = Var2,
        metric = Var3) %>%
      pivot_wider(
        names_from = metric,
        values_from = value) %>%
      filter(specimen_a != specimen_b, !is.na(estimate)) %>%
      mutate(population_a = substr(specimen_a, 1, 2),
             population_b = substr(specimen_b, 1, 2),
             same_pop = ifelse(population_a == population_b, 1, 0),
             same_residence = ifelse((population_a == "Import_Popul" & population_b == "Population_2") |
                                       (population_a == "Population_2" & population_b == "Import_Popul") | same_pop == 1, 1, 0))


    all_ibd_tables <- rbind(combined_pop_ibd_table_diff_res_only, pop1_ibd_table, pop2_ibd_table) %>%
      mutate(weight = w) %>% relocate(weight)

    all_ibd_tables <- all_ibd_tables %>%
      bind_rows(
        all_ibd_tables %>% rename(specimen_a = specimen_b, specimen_b = specimen_a,
                                  population_a = population_b, population_b = population_a)  # flipped version
      )

    full_results <- rbind(full_results, all_ibd_tables)

    print(w)
  }


  JD <- calculate_jost_d(pop1_afreq, pop2_afreq)
  FST <- calculate_FST(pop1_afreq, pop2_afreq)

  return(list(full_results, JD, FST))

}









##########
###### create function for getting local ancestors
##########

simulate_local_ancestors <- function(og_data, gen1_params, gen1_pop_size, gen1_ancestor_size,
                                     gen2_params, gen2_pop_size, gen2_ancestor_size,
                                     gen3_params, gen3_pop_size,
                                     final_pop_size, final_pop_gen_2_prop){


  og_dsmp <- formatDat(og_data, svar = "library_sample_name", lvar = "target_name", avar = "seq")
  og_coi   <- getCOI(og_dsmp, lrank = 1)
  og_afreq <- calcAfreq(og_dsmp,og_coi, tol = 1e-5)

  filtered_og <- filter_good_coverage_dom_haplo(og_data)

  ##### do first recombuddy sim population for gen1
  sim_pop_gen1 <- do.call(sim_population, c(list(unique(filtered_og$library_sample_name)),
                                            gen1_pop_size, as.list(gen1_params)))

  ## intersect with madhatter panel
  gen1_intersected <- filter_and_intersect(sim_pop_gen1, filtered_og, "Gen1")

  gen1_dsmp <- formatDat(gen1_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  gen1_coi   <- getCOI(gen1_dsmp, lrank = 2)
  gen1_afreq <- calcAfreq(gen1_dsmp, gen1_coi, tol = 1e-5)

  ## select 100 strains to make new population from
  gen1_new_ancestors = gen1_intersected %>%
    mutate(library_sample_name = paste0(simulated_sample, "-", within_sample_genotype),
           target_name = target) %>%
    select(library_sample_name, target_name, seq) %>%
    filter(library_sample_name %in% sample(unique(library_sample_name), gen1_ancestor_size))

  ## make second generation
  sim_pop_gen2 <- do.call(sim_population, c(list(unique(gen1_new_ancestors$library_sample_name)),
                                            gen2_pop_size, as.list(gen2_params)))



  ## intersect with madhatter panel
  gen2_intersected <- filter_and_intersect(sim_pop_gen2, gen1_new_ancestors, "Gen2")

  gen2_dsmp <- formatDat(gen2_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  gen2_coi   <- getCOI(gen2_dsmp, lrank = 2)
  gen2_afreq <- calcAfreq(gen2_dsmp, gen2_coi, tol = 1e-5)

  ## select 100 strains to make new population from
  gen2_new_ancestors = gen2_intersected %>%
    mutate(library_sample_name = paste0(simulated_sample, "-", within_sample_genotype),
           target_name = target) %>%
    select(library_sample_name, target_name, seq) %>%
    filter(library_sample_name %in% sample(unique(library_sample_name), gen2_ancestor_size))

  ## make third generation
  sim_pop_gen3 <- do.call(sim_population, c(list(unique(gen2_new_ancestors$library_sample_name)),
                                            gen3_pop_size, as.list(gen3_params)))

  ## intersect with madhatter panel
  gen3_intersected <- filter_and_intersect(sim_pop_gen3, gen2_new_ancestors, "Gen3")

  gen3_dsmp <- formatDat(gen3_intersected, svar = "simulated_sample", lvar = "target", avar = "seq")
  gen3_coi   <- getCOI(gen3_dsmp, lrank = 2)
  gen3_afreq <- calcAfreq(gen3_dsmp, gen3_coi, tol = 1e-5)


  ### sample some from gen 3 and some from gen 2 according to final pop size parameter and final pop prop gen 2 parameter
  sampled_gen3 <- gen3_intersected %>%
    mutate(library_sample_name = paste0(simulated_sample, "-", within_sample_genotype),
           target_name = target) %>%
    select(library_sample_name, target_name, seq) %>%
    filter(library_sample_name %in% sample(unique(library_sample_name), final_pop_size * (1-final_pop_gen_2_prop), replace = F))

  sampled_gen2 <- gen2_intersected %>%
    mutate(library_sample_name = paste0(simulated_sample, "-", within_sample_genotype),
           target_name = target) %>%
    select(library_sample_name, target_name, seq) %>%
    filter(library_sample_name %in% sample(unique(library_sample_name), final_pop_size * final_pop_gen_2_prop, replace = F))

  ### combine for final pop
  final_pop <- rbind(sampled_gen2, sampled_gen3)

  ## get afreq for final pop
  final_pop_dsmp <- formatDat(final_pop, svar = "library_sample_name", lvar = "target_name", avar = "seq")
  final_pop_coi   <- getCOI(final_pop_dsmp, lrank = 2)
  final_pop_afreq <- calcAfreq(final_pop_dsmp,final_pop_coi, tol = 1e-5)

  ### calculate distance measures between og pop and final pop
  final_og_JD <- mean(calculate_jost_d(og_afreq, final_pop_afreq), na.rm = T)
  final_og_FST <- mean(calculate_FST(og_afreq, final_pop_afreq), na.rm = T)

  ## plot afreq
  plot_final_og_af <- plot_af(og_afreq, final_pop_afreq, 700)

  return(list(final_pop, final_og_JD, final_og_FST, plot_final_og_af, og_afreq, gen1_afreq, gen2_afreq, gen3_afreq, final_pop_afreq))
}
