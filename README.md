
# Custom analysis with PLINK R plugin
_*If like me, you thought this would be great but hadn't actually got around to figuring out how to use it, here is a script to play with and some setup instructions._

<br>

It is possible to call R from PLINK. This facility allows you to keep genotype and phenotype data in PLINK binary format and perform a custom analysis. Below is an example of how this facility can be used to retrieve model fit statistics.

More information for PLINK's **R Plugin functions** is available in the [1.07](http://zzz.bwh.harvard.edu/plink/rfunc.shtml) and [1.9](https://www.cog-genomics.org/plink/1.9/rserve) documentation, including details for changing port, host, socket. 

<br>


## Getting started

First, you will need to install the development version of [PLINK](https://www.cog-genomics.org/plink/1.9/), and the latest version of [R](https://cran.r-project.org/). Open R and install relevant packages. `Rserve` is required; `broom` and a couple of `tidyverse` packages are needed for the specific example below. Make a note of the `Rserve` installation location printed by `install.packages`. You will need to point to it later.

To copy the R script, clone this repository.

```{bash}
git clone https://github.com/kenhanscombe/plink-custom-r.git
```

<br>


## Retrieve model fit statistics

In an R script (e.g. `plink_custom_analysis.R`), define a custom function. This script defines a pseudo R-squared for a logistic regression analysis, and uses the `broom` functions `glance` and `tidy` to collect fit statistics. (Note: Before changing anything to suit your needs, see the **Details** section at the end.)

```{r}
Rplink <- function(PHENO, GENO, CLUSTER, COVAR) {
  
    library(tidyverse)
    library(broom)

    pseudo_rsq <- function(model){
        dev <- model$deviance
        null_dev <- model$null.deviance
        model_n <- length(model$fitted.values)
        r2_cox_snell <- 1 - exp(-(null_dev - dev) / model_n)
        r2_nagelkerke <- r2_cox_snell / (1 - (exp(-(null_dev / model_n))))
        r2_nagelkerke
    }

    func <- function(snp) { 
        m <- glm(PHENO == 2 ~ COVAR + snp, family = "binomial")
        rsq <- pseudo_rsq(m)
        glance_m  <- glance(m) %>% unlist(.[1, ])
        tidy_m <- tidy(m) %>% select(-term) %>% tail(n = 1) %>% unlist()
        summary_m <- c(tidy_m, glance_m, rsq)
        c(length(summary_m), summary_m)
    }

    apply(GENO, 2, func)
}
```

<br>


To run the custom analysis, first start `Rserve` (supply the full path to `R CMD`). All data input and filtering flags to PLINK remain the same. Simply add `--R [R script filename]` to the PLINK call. The results of the custom analysis are written to `plink.auto.R` by default (As usual, you can change the file stem `plink` with `--out`).

```{bash}
R CMD /full/path/to/Rserve

plink \
--bfile {prefix} \
--pheno [filename] \
--covar [filename] \
(other optional filters ...)
--logistic \
--R custom_plink_analysis.R
```

**NB.** In the above example we're collecting model fit statistics from a logistic regresion (using the excellent package `broom`). `--logistic` is an optional sanity check. Compare `plink.assoc.logistic` to `plink.auto.R` for effect size, signed statistic, and p-value. (Adding a header to the `plink --R` output helps. See **Output** section below)

<br>


### Output

For each SNP in your analysis (i.e., each row in the output `plink.auto.R`), PLINK combines the vector of outputs `summary_m`, with the 4 values for CHR, SNP, BP, and A1. The R read command below adds a header to the custom output. You could of course do this in a bash one-liner, but if you're going to use R to visualize your association results and model fit statistics, you can add column names on reading in the data.

```{bash}
library(tidyverse)

# These col_names correspond to the custom analysis above.
custom_plink_result <- read_table2(
  "plink.auto.R",
  col_names = c("chr", "snp", "bp", "a1", "estimate", "std_error", "statistic",
    "p_value", "null_deviance", "df_null", "logLik", "aic", "bic", "deviance",
    "df_residual", "pseudo_rsq"),
  col_type = cols(snp = col_character())
)

# If your snp names are in CHR:BP format (or contain a colon), readr::read_table2 interprets them as time data. If you use a readr function to read in your results, supply the correct type as above.

```

<br>


### Multi-SNP model

If you want to inspect the overall model fit of a multi-SNP model, or compare the relative fit of multiple genetic variants (e.g. your 3 favourite SNPs), against a null model (e.g. 10 PCs), you cannot include the SNPs with the `--condition` flag. PLINK's `--R` always runs the analysis defined in `Rplink`. There are a couple of workarounds. 
One solution is to add the SNPs to the covariate file. First, convert the 3 SNPs to a 0/1/2 count of the reference allele with `--recode A`. The recoded SNPs appear in the last 3 columns of `plink.raw`. Add these 3 columns to the covariate file. Next, edit the function call in Rplink to not include snps (i.e., delete `+ snp` ) then run your custom analysis once with `--covar-number 1-10` (null), and a second time with `--covar-number 1-13`. Compare the 2 models.

<br>
<br>


***
### Details (summarised from PLINK 1.07 and 1.9 documentation)

For a sample of size n, genotyped at l genetic variants, including c covariates, all genotypes, phenotypes, covariates and cluster membership are accessible within the custom R script as:

<br>

**PHENO** A vector of phenotypes of length n.

**GENO** An n x l matrix of genotypes.

**CLUSTER** A vector of cluster membership codes of length n.

**COVAR** An n x c matrix of covariates.

<br>

The R script defines a function `Rplink`, with obligatory header, and return value, as follows, 

```{r}
Rplink <- function(PHENO, GENO, CLUSTER, COVAR) {
    
    # A function f is applied to the columns of GENO (i.e. to each genetic variant) and
    # must return a numeric vector v, combined with its length.
    f <- function(s) {
        
        # Function body
      
        c(length(v), v)
    }
    
    apply(GENO, 2, f)
}
```
