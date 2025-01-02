# File: mix_expt_examples.R
# Title: Example analysis of microbial mix experiments
# Authors: jeff smith, R. Fredrik Inglis
# Date: 2025-01
# Url: https://github.com/matryoshkev/mix-experiment-tools


# DESCRIPTION ==================================================================

# This script illustrates how one might calculate and visualize the fitness 
# effects of microbial interactions in mix experiments. 


# DEPENDENCIES =================================================================

# Made in R v4.2.1

# Install packages used by this script
install.packages("dplyr")

# Load packages
library(dplyr)  # General data-handling tools
source("mixexptr.R")


# DATASET 1: ABSOLUTE ABUNDANCE OF EACH STRAIN =================================

# Load data (ancestral and evolved Myxococcus)
data_smith2010 <-
	read.table("data_smith_2010.tsv", header = TRUE, sep = "\t") |>
	tibble() |>
	rename(
		initial_cells_evolved = "initial_cells_A",
		initial_cells_ancestral = "initial_cells_B",
		final_spores_evolved = "final_spores_A",
		final_spores_ancestral = "final_spores_B"
	) |>
	select(exptl_block, starts_with("initial"), starts_with("final"))

# Measured values are initial and final cell counts for each strain
data_smith2010

# Quick diagnostic report to visualize fitness results
report_mix_fitness(
	data = data_smith2010, 
	strain_names = c(A = "evolved", B = "ancestral"), 
	var_names = c(
		initial_count_A = "initial_cells_evolved",
		initial_count_B = "initial_cells_ancestral",
		final_count_A = "final_spores_evolved",
		final_count_B = "final_spores_ancestral"
	)
)
# Generated pdfs show:
# - Evolved strain decreases total-group fitness
# - Within-group fitness is frequency-dependent with functional form 
#   log(w_A/w_B) ~ log(n_A/n_B)

# Calculate fitness effects of mixing
fitness_smith2010 <- 
	calculate_mix_fitness(
		data = data_smith2010, 
		strain_names = c(A = "Evolved GVB206.3", B = "Ancestral GJV10"), 
		var_names = c(
			initial_count_A = "initial_cells_evolved",
			initial_count_B = "initial_cells_ancestral",
			final_count_A = "final_spores_evolved",
			final_count_B = "final_spores_ancestral"
		)
	)

# View fitness data
fitness_smith2010 |> select(
	exptl_block, initial_proportion_A, initial_ratio_A, starts_with("fitness")
)
# This is what we'd analyze statistically

# Plot fitness effects of mixing
dev.new(width = 6.25, height = 2.25, units = "in", noRStudioGD = TRUE)
fig_smith2010 <- 
	fitness_smith2010 |>
	plot_mix_fitness(
		# strain_names = c(A = "GVB206.3", B = "GJV10"), 
		strain_names = c(A = "evolved", B = "ancestral"), 
		mix_var = "initial_ratio_A"
	)
plot(fig_smith2010)
ggsave(
	"smith2010.pdf", 
	plot = fig_smith2010, 
	width = 6.25, 
	height = 2.25, 
	units = "in"
)


# DATASET 2: TOTAL ABUNDANCE + STRAIN FREQUENCY ================================

# Load data (antibiotic-resistant and sensitive E. coli)
data_Yurtsev2013 <-
	read.table("data_Yurtsev2013.tsv", header = TRUE, sep = "\t") |>
	tibble()

# Measured values are total cell density (OD_600) and strain frequency
data_Yurtsev2013 |> select(starts_with("initial"), starts_with("final"))

# Calculate fitness with mixexptr
fitness_Yurtsev2013 <- 
	calculate_mix_fitness(
		data = data_Yurtsev2013, 
		strain_names = c(A = "resistant", B = "sensitive"), 
		var_names = c(
			initial_count_total = "initial_number_total",
			initial_proportion_A = "initial_proportion_A",
			final_count_total = "final_number_total",
			final_proportion_A = "final_proportion_A"
		)
	) |>
	select(
		antibiotic, dilution, 
		initial_proportion_A, initial_ratio_A, 
		starts_with("fitness")
	)

# Plot fitness effects
dev.new(width = 6.25, height = 2.25, units = "in", noRStudioGD = TRUE)
fig_Yurtsev2013 <- 
	fitness_Yurtsev2013 |>
	plot_mix_fitness(
		strain_names = c(A = "resistant", B = "sensitive")
	)
plot(fig_Yurtsev2013)
	# - Total-group fitness not affected by mix frequency
	# - Within-group fitness is strongly frequency-dependent below some 
	#   threshold value

