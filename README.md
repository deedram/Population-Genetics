# Population-Genetics
------

### fst.R

#### Overview
This code can be run from either the command line or within R. It generates Fst (fixation index) values single sites of a sample when provided a base comparison sample. It also provides total fst for each sample as well as cds values for Zika Virus, West Nile Virus, Chikungunya Virus, and Dengue Virus Type 2. Other viruses can be added easily as long as the coding region for the virus is known. 

Fst is a measure of population differentiation due to genetic structure. It is an extremely common statistic for population genetics. There are a few different ways to measure Fst. In this case, a method-of-moments estimator that does not rely on assumptions about the shape of the sampling distribution is used. It is an extension of a Fst estimator proposed by Reynolds et al. (1983).

#### Reference
The specific calculations were used from Fumagalli, M., Vieira, F. G., Korneliussen, T. S., Linderoth, T., Huerta-Sánchez, E., Albrechtsen, A., & Nielsen, R. (2013). Quantifying population genetic differentiation from next-generation sequencing data. Genetics, 195(3), 979–992. https://doi.org/10.1534/genetics.113.154740

#### Code Usage
##### Within R
This code can be used in R if the first 4 lines are not used. The following need to be defined prior to running the code:
```
infile <- "filepath to inputfile containing samples"
outfile <- "filepath to outputfile with new fst values"
input <- "name of base sample"
```
The infile must provide the following columns:
  
-"sample" contains the sample names for the base and all desired samples   
-"pos" contains the genome positions   
-"af" contains the allele frequency   
  
##### From Command Line
After loading the R module in the command line, the code can easily be run with the following:
```
Rscript fst.R "filepath to inputfile containing samples" "filepath to outputfile with new fst values" "name of base sample"
```
