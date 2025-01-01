# TODO: mix-experiment-tools

## Example analysis (mix_expt_examples.R)

- Myxo example: Calculate cfu from colony counts & dilution
- Show how to calculate fitness measures "by hand"
- For more info see: smith j & Inglis RF (2025) Best practices for analyzing the fitness effects of microbial interactions. bioRxiv ... DOI: ...

Approach  
- Model best practices described in paper
- Use this file to figure out clearest/easiest syntax to design for
- Only brief comments: let code be foreground
- Comments about what the data show (e.g. x-proportion vs x-logratio)

Use cases  
- Quick diagnostic tool to visualize data (preliminary analysis)
- Starting materials for new lab member
- Supplemental figure in paper

Be sure to include  
- Dataset with number A + number B
- Dataset with total + proportion
- `report_mix_fitness()` = quick diagnostic plots
- `calculate_mix_fitness()` 
- Calculate fitness "by hand", then into `plot_mix_fitness()`
- Both x-scales: proportion, log-ratio

Open questions  
- Flesh out description? 
- Describe math in script? 


## Data files

- Both tsv and csv examples
- Should be formatted well, but not same `var_names` as script

smith 2010  
- Just colony counts and dilution

Yurtsev 2013
- Reformat: simplify, follow original more closely
- csv
- Re-use author language in description
- ampicillin, culture_id, OD_initial, OD_final, fraction_gfp_initial, fraction_gfp_final, replicate, dilution


## Other files
- Flesh out README.md
- spreadsheet template
- tsv, csv templates (with comment fields)


## Package: mixexptr.R

- Skeleton package
- GitHub repository

### Later scope
- plot fitness measures separately
  . plot_strain_fitness()
  . plot_multilevel_fitness()
    . plot_total_group_fitness()
    . plot_within_group_fitness()
- custom color, fill
- custom point shape
- custom theme
- custom limits
- custom breaks
- combine figures (after adding lines for fitted models, for example)
- log-transformed values not 10^x

