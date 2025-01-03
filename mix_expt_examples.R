# File: mix_expt_examples.R
# Title: Example analysis of microbial mix experiments
# Authors: jeff smith, R. Fredrik Inglis
# Url: https://github.com/matryoshkev/mix-experiment-tools
# Date last modified: 2025-01

# Description: 
#   Illustrate how one might calculate and visualize the fitness effects of 
#   microbial interactions in mix experiments

# Dependencies: 
# install.packages("dplyr")  # Install packages used by this script
library(dplyr)  # General data-handling tools
source("mixexptr.R")  # Mix experiment tools


# DATASET 1: ABSOLUTE ABUNDANCE OF EACH STRAIN =================================

# Load data (ancestral and evolved Myxococcus)
data_smith2010 <-
	read.csv("data_smith2010.csv", comment.char = "#") |>
	tibble()

data_smith2010
# Data are raw colony counts and dilutions

# Calculate final cell number (surviving spores) from colony counts
data_smith2010 <- data_smith2010 |>
	mutate(
		final_spores_evolved = colonies_rif_1 * dilution_rif_1,
		final_spores_ancestral = colonies_kan_1 * dilution_kan_1,
	) |>
	select(exptl_block, starts_with("initial"), starts_with("final"))

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

# Calculate fitness effects of mixing using mixexptr::calculate_mix_fitness()
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

# Plot fitness effects of mixing (inside R)
dev.new(width = 6.25, height = 2.25, units = "in", noRStudioGD = TRUE)
fig_smith2010 <- 
	fitness_smith2010 |>
	plot_mix_fitness(
		strain_names = c(A = "evolved", B = "ancestral"), 
		mix_var = "initial_ratio_A"
	)
plot(fig_smith2010)


# DATASET 2: TOTAL ABUNDANCE + STRAIN FREQUENCY ================================

# Load data (antibiotic-resistant and sensitive E. coli)
data_Yurtsev2013 <-
	read.csv("data_Yurtsev2013.csv", comment.char = "#") |>
	tibble()

data_Yurtsev2013
# Measured values are total cell density (OD600) and strain frequency
# Experimental treatments are ampicillin and dilution

# Drop rows with negative strain frequencies (artifact from flow cytometry)
data_Yurtsev2013 <- data_Yurtsev2013 |>
	filter(fraction_resistant_initial > 0 & fraction_resistant_final > 0)

# mixexptr::calculate_mix_fitness() also works with total/fraction data
calculate_mix_fitness(
	data = data_Yurtsev2013, 
	strain_names = c(A = "resistant", B = "sensitive"), 
	var_names = c(
		initial_count_total = "OD_initial",
		initial_proportion_A = "fraction_resistant_initial",
		final_count_total = "OD_final",
		final_proportion_A = "fraction_resistant_final"
	)
) |>
select(
	ampicillin, dilution, 
	initial_proportion_A, initial_ratio_A, 
	starts_with("fitness")
)

# Or you could calculate fitness effects "by hand"
fitness_Yurtsev2013 <- 
	data_Yurtsev2013 |>
	mutate(
		fitness_total = OD_final / OD_initial, 
		fitness_AmpR = 
			(OD_final * fraction_resistant_final) / 
			(OD_initial * fraction_resistant_initial), 
		fitness_AmpS = 
			(OD_final * (1-fraction_resistant_final)) / 
			(OD_initial * (1-fraction_resistant_initial)), 
		fitness_ratio_AmpR_AmpS = 
			(fraction_resistant_final / (1-fraction_resistant_final)) /
			(fraction_resistant_initial / (1-fraction_resistant_initial))
	)

# Treatment levels
# sort(unique(data_Yurtsev2013$ampicillin))
# sort(unique(data_Yurtsev2013$dilution))

# Plot fitness effects for 100 ug/mL ampicillin, 100-fold dilution
fig_Yurtsev2013 <- 
	fitness_Yurtsev2013 |>
	filter(ampicillin == 100 & dilution == 100) |>
	# filter(ampicillin == 100 & dilution == 200) |>
	plot_mix_fitness(
		strain_names = c(A = "resistant", B = "sensitive")
	)
dev.new(width = 6.25, height = 2.25, units = "in", noRStudioGD = TRUE)
plot(fig_Yurtsev2013)
# Total-group fitness not affected by mix frequency
# Within-group fitness is strongly frequency-dependent below some threshold

# Compare to treatment with different ... 


