# Tools to analyze microbial mix experiments

A common experimental design for measuring microbial interactions is to mix together two different microbes (different genotypes, for example) and measure how their behavior and reproductive success depends on mix frequency. How do they perform together compared to performing separately? 

Here we offer a set of digital tools to help analyze mix experiments using the R environment for statistical computing. Aimed at researchers and students with basic R experience, we hope these will be useful as: 
- a quick diagnostic tool to visualize data in preliminary analyses
- example code to modify and reuse
- starting materials for new lab member

These materials aim for best practice in data storage and analysis to enable open, transparent, repeatable science. 

## Files

`mix_expt_examples.R`
: Example analyses of data from mix experiments

`mixexptr.R`
: Library of R functions to calculate and plot fitness effects in mix experiments

`data_smith2010.csv`
: Example dataset where measured values are cell counts for each strain. Archive-quality csv format. Metadata embedded as comment lines beginning with `#`. 

`data_smith2010.xslx`
: Example dataset in spreadsheet format. Follows best practices described by [Broman \& Woo (2017)](https://doi.org/10.1080/00031305.2017.1375989). 

`data_Yurtsev2013.csv`
: Example dataset where measured values are total density and strain frequency


## Further reading

- smith j and Inglis RF (2021) Evaluating kin and group selection as tools for quantitative analysis of microbial data. Proceedings B 288:20201657. [https://doi.org/10.1098/rspb.2020.1657](https://doi.org/10.1098/rspb.2020.1657)
- Broman KW and Woo KH (2017) Data organization in spreadsheets. The American Statistician 72:2-10. [https://doi.org/10.1080/00031305.2017.1375989](https://doi.org/10.1080/00031305.2017.1375989)
- [Dryad | Good data practices](https://datadryad.org/stash/best_practices)

