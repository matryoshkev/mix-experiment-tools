# File: mixexptr.R
# Title: Analyze microbial interactions in mix experiments
# Author: jeff smith 
# Date: 2025-01
# Url: https://github.com/matryoshkev/mix-experiment-tools


# DEPENDENCIES =================================================================

# Made in R v4.2.1
library(ggplot2)  # To make plots and figures
library(gtable)   # To combine subplots in one window
# library(grid)     # To plot gtable objects with grid.draw


# calculate_fitness.R ==========================================================

# Calculate recommended fitness measures
#
# Possible pop-state inputs:
#   (number_A, number_B)
#   (number_total, proportion_A)
#
# Calculates recommended fitness measures from data describing:
#   absolute abundance (cells or virions) of each strain, or
#   absolute abundance of total (sum of both strains) and 
#   proportion belonging to a focal strain
#
# Can be counts or densities (cfu/mL, OD_600). Have to be greater than zero. 
# Initial and final must be same units. 


#' Calculate fitness effects of mixing
#' 
#' [Description goes here]
#' 
#' @param data Data frame with ...
#' @param strain_names Named vector or list with names for ...
#' @param var_name Named vector or list with column names in data ...
#' 
#' @returns Object of same type as data with added columns: ...
#'   name_A, name_B,
#'   initial_proportion_A, initial_ratio_A,
#'   fitness_A, fitness_B, fitness_total, fitness_ratio_A
#' 
#' @export 
#' 
calculate_mix_fitness <- function(data, strain_names, var_names) {
	output <- data
	strain_names <- as.list(strain_names)
	var_names <- as.list(var_names)

	# Label which is A and which is B
	output$name_A <- strain_names$A
	output$name_B <- strain_names$B

	# Calculate fitness and initial frequency
	if (
		# Data is count_A, count_B
		!is.null(var_names$initial_count_A) & 
		!is.null(var_names$initial_count_B) &
		!is.null(var_names$final_count_A) & 
		!is.null(var_names$final_count_B)
	) {
		output <- fitness_from_counts(output, var_names)
	} else if (
		# Data is count_total, proportion_A
		!is.null(var_names$initial_count_total) & 
		!is.null(var_names$initial_proportion_A) &
		!is.null(var_names$final_count_total) & 
		!is.null(var_names$final_proportion_A)
	) {
		output <- fitness_from_freq(output, var_names)
	} else {
		stop("Cannot calculate fitness from data")
	}

	# Drop NaN from single-strain mixes etc
	output[sapply(output, is.nan)] <- NA
	# output[sapply(output, is.infinite)] <- NA

	# Warn about fitness zeroes
	if (any(c(output$fitness_A == 0, output$fitness_B == 0), na.rm = TRUE)) {
		warning(
			"Some fitness values equal to zero. Undefined on log scale.", 
			call. = FALSE
		)
	}

	# Return frame with calculated values added
	output
}

# Calculate fitness effects from count data
# 
fitness_from_counts <- function(data, var_names) {
	# Warn if invalid data
	for (var_name in var_names) {
		if (any(data[[var_name]] < 0, na.rm = TRUE)) {
			warning(
				var_name, " values < 0: not biologically meaningful", 
				call. = FALSE
			)
		}
	}

	# Calculate values
	initial_A <- data[[var_names$initial_count_A]]
	initial_B <- data[[var_names$initial_count_B]]
	final_A   <- data[[var_names$final_count_A]]
	final_B   <- data[[var_names$final_count_B]]
	output <- within(data, {
		initial_proportion_A <- initial_A / (initial_A + initial_B)
		initial_ratio_A      <- initial_A / initial_B
		fitness_A            <- final_A / initial_A
		fitness_B            <- final_B / initial_B
		fitness_total        <- (final_A + final_B) / (initial_A + initial_B)
		fitness_ratio_A      <- fitness_A / fitness_B
	})

	# Return frame with calculated values added
	output
}

# Calculate fitness effects from total count + proportion data
# 
fitness_from_freq <- function(data, var_names) {
	# Warn if invalid data
	# for (var_name in var_names)) {
	# for (var_name in c(name_initial_total, name_final_total)) {
		# if (any(data[[var_name]] < 0, na.rm = TRUE)) {
			# warning(
				# var_name, " values < 0: not biologically meaningful", 
				# call. = FALSE
			# )
		# }
	# }
	# for (var_name in c(name_initial_freq, name_final_freq)) {
		# if (any(c(data[[var_name]] < 0, data[[var_name]] > 1), na.rm = TRUE)) {
			# warning(var_name, " must be in range [0, 1]", call. = FALSE)
		# }
	# }

	# Calculate values
	initial_total <- data[[var_names$initial_count_total]]
	initial_freq  <- data[[var_names$initial_proportion_A]]
	final_total   <- data[[var_names$final_count_total]]
	final_freq    <- data[[var_names$final_proportion_A]]
	output <- within(data, {
		initial_proportion_A <- initial_freq
		initial_ratio_A <- initial_freq / (1-initial_freq)
		fitness_A <-
			(final_total * final_freq) / (initial_total * initial_freq)
		fitness_B <-
			(final_total * (1-final_freq)) / (initial_total * (1-initial_freq))
		fitness_total <- final_total / initial_total
		fitness_ratio_A <- 
			(final_freq / (1-final_freq)) / (initial_freq / (1-initial_freq))
	})

	# Return frame with calculated values
	output
}


# plot_fitness.R ===============================================================

#' Make pdf figures showing fitness effects of mixing
#'
#' [Description goes here] [strain, total group, within-group]
#' 
#' @param data [Input data]
#' @param strain_names [Strain names]
#' @param var_names [Measured vars]
#' 
#' @returns [Generates pdf plots (proportion and log-ratio)]
#' 
#' @export
#' 
report_mix_fitness <- function(data, strain_names, var_names) {
	data_fitness <- calculate_mix_fitness(data, strain_names, var_names)
	fig_proportion <- 
		plot_mix_fitness(
			data_fitness, 
			strain_names, 
			suppress_window = TRUE
		)
	fig_logratio <- 
		plot_mix_fitness(
			data_fitness, 
			strain_names, 
			mix_var = "initial_ratio_A", 
			suppress_window = TRUE
		)
	ggsave(
		"mix_fitness_1.pdf", 
		plot = fig_proportion, 
		width = 6.25, 
		height = 2.25, 
		units = "in"
	)
	ggsave(
		"mix_fitness_2.pdf", 
		plot = fig_logratio, 
		width = 6.25, 
		height = 2.25, 
		units = "in"
	)
	message("mixexptr: Results saved as mix_fitness_1.pdf, mix_fitness_2.pdf")
}

#' Plot fitness effects of mixing
#'
#' [Description goes here] [strain, total group, within-group]
#' 
#' @param data Data frame containing fitness measures
#' @param strain_names
#' @param var_names [Fitness vars]
#' @param x_scale [Optional string x-axis mixing scale]
#' @param suppress_plot [Optional don't make window with plot]
#'
#' @returns [ggplot object]
#'
#' @export
#' 
plot_mix_fitness <- function(
	data, 
	strain_names, 
	var_names = list(
		initial_proportion_A = "initial_proportion_A", 
		initial_ratio_A = "initial_ratio_A", 
		fitness_A = "fitness_A", 
		fitness_B = "fitness_B", 
		fitness_total = "fitness_total", 
		fitness_ratio_A = "fitness_ratio_A"
	), 
	mix_var = "initial_proportion_A", 
	suppress_window = TRUE
) {
	strain_names <- as.list(strain_names)
	var_names <- as.list(var_names)

	if (mix_var == "initial_ratio_A") {
		data <- subset(data, initial_proportion_A > 0 & initial_proportion_A < 1)
		# TODO: set x_limits to include 1
	}

	fig_output <- 
		combine_figures(
			plot_fitness(data, strain_names, var_names, mix_var), 
			plot_fitness_ratio(data, strain_names, var_names, mix_var), 
			units = 25, 
			middle = 17
		)

	if (!suppress_window) {
		# dev.new(width = 6.25, height = 2.25, units = "in", noRStudioGD = TRUE)
			# Inconsistently buggy...
		# grid.draw(fig_output)
		plot(fig_output)
	}
	fig_output
}

# Plot strain and total group fitness
# 
plot_fitness <- function(
	data, 
	strain_names, 
	var_names, 
	mix_var = "initial_proportion_A", 
	limits = NULL
	# drop_zero = FALSE
) {
	# TODO: Validate fitness data (> 0)
	# TODO: Limits for initial_ratio_A
	# TODO: Shared fitness scale

	figure_fitness_base() |>
	add_scale_fitness() |>
	add_scale_mix_var(strain_names, mix_var) %+%
	format_to_plot_fitness(data, strain_names, var_names)
}

# Plot within-group fitness ratio (fitness_A/fitness_B)
#   Will eventually be user-facing
#   TODO: drop_Inf = FALSE
#   TODO: drop_zero = FALSE
# 
plot_fitness_ratio <- function(
	data, 
	strain_names, 
	var_names, 
	mix_var = "initial_proportion_A", 
	limits = NULL
) {
	strain_names <- as.list(strain_names)
	var_names <- as.list(var_names)
	# TODO: Validate fitness data (> 0) ?
	# TODO: Limits for initial_ratio_A
	# TODO: Shared fitness scale (at least as much as strain & total group)

	figure_fitness_ratio_base() |>
	add_scale_fitness_ratio(strain_names) |>
	add_scale_mix_var(strain_names, mix_var) %+%
	subset(data, !is.na(fitness_ratio_A))	
}

# Combine subfigures: strain/total + within-group
#   Eventually user-facing?
# 
combine_figures <- function(fig1, fig2, units, middle) {
	# widths <- c(1, 1)
	# width1 <- widths[[1]]
	# width2 <- widths[[2]]
	pdf(file = NULL)  
		# So ggplotGrob() doesn't create blank plot window
		# see https://github.com/tidyverse/ggplot2/issues/809
	fig_output <- gtable::gtable_add_grob(
		gtable::gtable(
			widths = unit(rep(1, units), "null"), 
			heights = unit(rep(1, 1), "null")
		), 
		grobs = list(ggplotGrob(fig1), ggplotGrob(fig2)), 
		l = c(1, middle), # Left extents
		r = c(middle - 1, units),  # Right extents
		t = c(1, 1), # Top extents
		b = c(1, 1)  # Bottom extents
	)
	dev.off()  # end of workaround
	fig_output
}

# Reshape data to plot fitness_A, fitness_B, fitness_total
# 
format_to_plot_fitness <- function(data, strain_names, var_names) {	
	output <- as.data.frame(data)  # So reshape() doesn't choke on tibbles
	output$initial_proportion_A <- output[[var_names$initial_proportion_A]]
	output$initial_ratio_A      <- output[[var_names$initial_ratio_A]]
	output$fitness_A            <- output[[var_names$fitness_A]]
	output$fitness_B            <- output[[var_names$fitness_B]]
	output$fitness_total        <- output[[var_names$fitness_total]]
	output <- output |>
		subset(select = c(
			"initial_proportion_A", "initial_ratio_A", 
			"fitness_A", "fitness_B", "fitness_total"
		)) |>
		reshape(
			direction = "long", 
			varying = c("fitness_A", "fitness_B", "fitness_total"), 
			v.names = c("fitness"), 
			times = c(strain_names$A, strain_names$B, "Total group"), 
			timevar = "strain"
		) |>
		subset(!is.na(fitness))
	output$strain <- factor(output$strain, 
		levels = c(strain_names$A, strain_names$B, "Total group")
	)
	output$my_facet <- output$strain == "Total group"
	output
}


# figure_elements.R ============================================================

figure_base <- function() {
	ggplot() +
	theme(
		text                 = element_text(size = 9), 
		legend.title         = element_blank(), 
		legend.background    = element_blank(), 
		legend.direction     = "horizontal", 
		legend.justification = c(0.5, 0.15), 
		legend.position      = c(0.5, 1),
		strip.text           = element_blank(),
		strip.background     = element_blank()
	) + 
	ggtitle("")  # Space for legend, keeps plot height consistent across figs
}

# Base elements for fitness figure (strain A, strain B, total group)
# 
figure_fitness_base <- function() {
	figure_base() +
	aes(color = strain, fill = strain) + 
	scale_fill_manual(values  = c("tan",  "lightsteelblue",  gray(0.65))) + 
	scale_color_manual(values = c("tan4", "lightsteelblue4", gray(0.1))) + 
		# Strain A, strain B, total group
	geom_point(shape = 21) + 
	geom_point(shape = 1) +
	facet_wrap(~ my_facet, nrow = 1)
}

# Base elements for fitness ratio figure
# 
figure_fitness_ratio_base <- function() {
	figure_base() +
	geom_point(shape = 21, color = gray(0.1), fill = gray(0.65)) + 
	geom_point(shape = 1, color = gray(0.1))
}

# Add x-axis scale: initial proportion A or initial ratio A/B
# 
add_scale_mix_var <- function(input_fig, strain_names, mix_var) {
	if (mix_var == "initial_proportion_A") {
		output_fig <- input_fig |>
			add_scale_initial_proportion(strain_names)
	} else if (mix_var == "initial_ratio_A") {
		output_fig <- input_fig |>
			add_scale_initial_ratio(strain_names)
	} else {
		stop(
			"mix_var must be 'initial_proportion_A' or 'initial_ratio_A'", 
			"given: ", mix_var
		)
	}
	output_fig
}

# Add x-axis scale: initial proportion A
# 
add_scale_initial_proportion <- function(input_fig, strain_names) {
	breaks <- seq(0, 1, by = 0.2)
	input_fig + 
		aes(x = initial_proportion_A) +
		scale_x_continuous(
			name = paste("Initial proportion", strain_names$A), 
			limits = c(0, 1), 
			breaks = breaks, 
			minor_breaks = breaks
		)
}

# Add x-axis scale: initial ratio A/B (log10)
#   TODO: limits always include 1
#   TODO: minor_breaks
# 
add_scale_initial_ratio <- function(input_fig, strain_names, limits = NULL) {
	breaks <- 10^c(-10:10)
	# if (!is.null(limits)) {
		# # Set minor breaks based on range (1e1 or 1e2)
	# } 
	input_fig + 
		aes(x = initial_ratio_A) +
		scale_x_log10(
			name = bquote(
				"Initial ratio" ~ .(strain_names$A)/.(strain_names$B)
			),
			breaks = breaks, 
			labels = scales::label_log()
			# labels = log_labels(breaks)
			# limits = limits.initial.ratio, 
			# minor_breaks = breaks.initial.ratio.minor,
		)
	# geom_vline(xintercept = 1, color = "white", linewidth = 0.8)
}

# Add y-axis scale: fitness (strain A, strain B, total group)
#   TODO: limits (min 10-fold range, shared scale across subplots, always include 1)
#   TODO: breaks
#   TODO: minor_breaks
#   TODO: log labels with "1"
#   TODO: linear_labels when appropriate
# 
add_scale_fitness <- function(input_fig) {
	input_fig + 
	aes(y = fitness) +
	scale_y_log10(
		name = "Wrightian fitness\n (final no. / initial no.)", 
		labels = scales::label_log(), 
		# labels = log_labels()
		# limits       = limits.fitness, 
		# breaks       = breaks.fitness,
		# minor_breaks = breaks.fitness.minor, 
	) 
	# geom_hline(yintercept = 1, color = "white", linewidth = 0.8)
}

# Add y-axis scale: Within-group fitness ratio A/B
#   TODO: limits (min 10-fold range, shared scale across subplots)
#   TODO: breaks
#   TODO: labels
# 
add_scale_fitness_ratio <- function(input_fig, strain_names) {
	input_fig + 
	aes(y = fitness_ratio_A) +
	scale_y_log10(
		# name = bquote("Fitness ratio" ~ .(strain_names$A)/.(strain_names$B)),
		name = paste("Fitness ratio\n", strain_names$A, "/", strain_names$B),
		# limits       = limits.fitness, 
		# breaks       = breaks.fitness,
		# minor_breaks = breaks.fitness.minor, 
		# labels       = labels.fitness
	) 
	# geom_hline(yintercept = 1, color = "white", linewidth = 0.8)
}

# get_fitness_limits() ? 
# get_fitness_breaks() ? 
# get_fitness_labels() ? 

# Linear-scale axis labels: 0, 1, or x
# 
linear_labels <- function(x) {
	if (x == 1) { 
		as.expression(1) 
	} else if (x == 0) { 
		as.expression(0)
	} else { 
		as.expression(x)
	}
}
linear_labels <- Vectorize(linear_labels)

# Log-scale axis labels: 10^x or 1
# 
log_labels <- function(x) {
	if (x == 1) { 
		as.expression(1) 
	} else {
		as.expression(bquote(10^.(log10(x))))
	}
}
log_labels <- Vectorize(log_labels)



# ==============================================================================

