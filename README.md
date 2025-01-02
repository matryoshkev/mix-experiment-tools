# Digital tools to analyze microbial interactions in mix experiments

A common experimental design for measuring microbial interactions is to mix together two different microbes (different genotypes, for example) and measure how their behavior and reproductive success depends on mix frequency. How do they perform together compared to performing separately? 

Here we offer a set of digital tools to help analyze mix experiments using the R environment for statistical computing. Aimed at researchers and students with basic R experience, we hope these will be useful as: 
- starting materials for new lab member
- a quick diagnostic tool to visualize data in preliminary analyses
- example code to modify and reuse

## R code

`mix_expt_examples.R`
: Example analysis of data from mix experiments

`mixexptr.R`
: Library of R functions to calculate and plot fitness effects in mix experiments

## Data handling

`data_smith2010.csv`
: Example dataset in archive-quality csv format. Metadata embedded as comment lines beginning with `#`. 

`data_smith2010.xslx`
: Example dataset in spreadsheet format. Follows best practices described by [Broman \& Woo (2017)](https://doi.org/10.1080/00031305.2017.1375989). 

## Further reading

- Broman KW and Woo KH (2017) Data organization in spreadsheets. The American Statistician 72:2--10. [https://doi.org/10.1080/00031305.2017.1375989](https://doi.org/10.1080/00031305.2017.1375989)
- [Dryad | Good data practices](https://datadryad.org/stash/best_practices)

