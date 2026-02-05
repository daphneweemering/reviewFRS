# Statistical analysis of the ALSFRS-R in randomized controlled clinical trials for amyotrophic lateral sclerosis: A systematic review

This research repository contains all material related to our study *"Statistical analysis of the ALSFRS-R in randomized controlled clinical trials for amyotrophic lateral sclerosis: A systematic review"*

## What's the project about? # 
Disability rating scales, like the revised amyotrophic lateral sclerosis functional rating scale (ALSFRS-R), play a pivotal role in clinical trials by assessing how experimental treatments affect the daily lives of people with neurodegenerative diseases. Limited guidance on analyzing these scales may contribute to the high failure rate of clinical trials. Using amyotrophic lateral sclerosis (ALS) as a case study, we aim to systematically review how disability rating scales have been analyzed in clinical trials, and how these approaches influence
the validity and the precision of trial results.

## What is in the repository? #

| Folder/file | Content |
|:-----------|:-----------------------------------------------------------------------|
| figures    | Holds the figures that are available in the manuscript and supplementary material. |
| scripts    | Contains analysis scripts:<br> • [`clean.CEF.R`](scripts/0.%20clean.CEF.R) – Cleans data from the Ceftriaxone ALS trial.<br> • [`clean.RT001.R`](scripts/0.%20clean.RT001.R) – Cleans data from the RT001 ALS trial.<br> • [`TE.R`](scripts/1.%20TE.R) – Simulates permuted treatment assignment and adds treatment effects.<br> • [`studyadj.R`](scripts/2.%20studyadj.R) – Adjusts data for study-based modifications.<br> • [`perm.CEF.R`](scripts/3.%20perm.CEF.R) – Fits models to the simulated Ceftriaxone ALS trial data.<br> • [`perm.RT001.R`](scripts/3.%20perm.RT001.R) – Fits models to the simulated RT001 ALS trial data.<br> • [`sim.JM.R`](scripts/4.%20sim.JM.R) – Simulates data using the joint modeling framework.<br> • [`models.JM.R`](scripts/5.%20models.JM.R) – Fits models repeatedly to the simulated data created in [`sim.JM.R`](scripts/4.%20sim.JM.R) to obtain false-positives rate and power.<br> • [`example.sim.R`](scripts/6.%20example.sim.R) – Creates a single simulated dataset using [`sim.JM.R`](scripts/4.%20sim.JM.R) and fits different models as an example of the analyses that we ran on the Ceftriaxone ALS trial data.<br> |

## Example of the analyses that were done in this study #
In our study, we used the Ceftriaxone ALS clinical trial dataset to guide the simulation of realistic trial scenarios with covariates. Because the original dataset cannot be shared, we provide simulated clinical trial data generated using a joint modeling framework. For a detailed description of the simulation method, see the eMethods section of the manuscript, or refer to the simulation code in [`sim.JM.R`](scripts/4.%20sim.JM.R). These simulated datasets allow users to reproduce the analyses from our study. The models can be run using these simulated data, but only include the baseline covariate (baseline score) as in our original analysis. Additional covariates are not included in the provided simulated data because simulating them while preserving the realistic dependencies observed in the original trial is not feasible.

For these analyses, R(studio) is needed.

#### Required files
- [`sim.JM.R`](scripts/4.%20sim.JM.R): Contains the function `JMD()` to simulate clinical trial data.  
- [`example.sim.R`](scripts/6.%20example.sim.R): Runs the example analyses on the simulated data.  

Make sure both files are in the same folder.

#### Step 1: Set working directory

```r
# Check current working directory
getwd()

# Change working directory if needed
setwd("path/to/where/files/are")
``` 

#### Step 2: Source the simulation code
```r
# Load the simulation code
source("4. sim.JM.R")
```
#### Step 3: Simulate a trial dataset
```r
D <- JMD(
  N = 500,        # number of patients
  Nm = 7,         # trial measurements
  Nm.pre = 4,     # pre-baseline measurements
  B2 = 0, B4 = 0, # no treatment effect on function or survival
  g1 = 0,         
  FUT = 12,       # months follow-up
  DO.TRT = 0.1,   # dropout for treatment
  DO.PLB = 0.1    # dropout for placebo
)
```
- D is a list containing 15 datasets, each in a different format or with a different missing data strategy (e.g., long, wide, LOCF, complete cases).
- A short description of each dataset is included in the script.

#### Step 4: Fit the models
- `est()` extracts effect estimates, p-values and z-scores.
- Lines 137–269 fit multiple models (m1, m2, …, m30).
- For each model, a short description is given of: the outcome, the analysis approach, the missing data strategy for death, the missing data strategy for non-death missingness, and adjustment for baseline score (respectively).

#### Step 5. Explore results
- Results can be explored using the `est()` function for the different models. For example: `est(m1)` gives the effect estimate, p-value and z-score for m1. 
- Use these outputs to understand how different modeling choices and missing data strategies affect the results in simulated trial scenarios.

#### Step 6. Visual display
- The last lines of code create a figure similar to Figure 4 in the article, showing z-scores for each model. All models should be available in the global environment to create the figure. 

## Machine and package information
```
─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.0 (2024-04-24)
 os       macOS 15.6.1
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Amsterdam
 date     2025-09-08
 rstudio  2025.05.1+513 Mariposa Orchid (desktop)
 pandoc   NA

─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package      * version   date (UTC) lib source
 abind          1.4-5     2016-07-21 [1] CRAN (R 4.4.0)
 backports      1.5.0     2024-05-23 [1] CRAN (R 4.4.0)
 boot           1.3-30    2024-02-26 [2] CRAN (R 4.4.0)
 broom          1.0.7     2024-09-26 [1] CRAN (R 4.4.1)
 cachem         1.1.0     2024-05-16 [1] CRAN (R 4.4.0)
 car            3.1-2     2023-03-30 [1] CRAN (R 4.4.0)
 carData        3.0-5     2022-01-06 [1] CRAN (R 4.4.0)
 cellranger     1.1.0     2016-07-27 [1] CRAN (R 4.4.0)
 checkmate      2.3.1     2023-12-04 [1] CRAN (R 4.4.0)
 cli            3.6.4     2025-02-13 [1] CRAN (R 4.4.1)
 coda           0.19-4.1  2024-01-31 [1] CRAN (R 4.4.0)
 codetools      0.2-20    2024-03-31 [2] CRAN (R 4.4.0)
 coin         * 1.4-3     2023-09-27 [1] CRAN (R 4.4.0)
 colorspace     2.1-1     2024-07-26 [1] CRAN (R 4.4.1)
 cowplot      * 1.1.3     2024-01-22 [1] CRAN (R 4.4.0)
 data.table   * 1.17.0    2025-02-22 [1] CRAN (R 4.4.1)
 DBI            1.2.2     2024-02-16 [1] CRAN (R 4.4.0)
 devtools       2.4.5     2022-10-11 [1] CRAN (R 4.4.0)
 digest         0.6.36    2024-06-23 [1] CRAN (R 4.4.0)
 dplyr        * 1.1.4     2023-11-17 [1] CRAN (R 4.4.0)
 ellipsis       0.3.2     2021-04-29 [1] CRAN (R 4.4.0)
 emmeans      * 1.10.3    2024-07-01 [1] CRAN (R 4.4.0)
 estimability   1.5.1     2024-05-12 [1] CRAN (R 4.4.0)
 fastmap        1.2.0     2024-05-15 [1] CRAN (R 4.4.0)
 forcats      * 1.0.0     2023-01-29 [1] CRAN (R 4.4.0)
 foreach        1.5.2     2022-02-02 [1] CRAN (R 4.4.0)
 fs             1.6.4     2024-04-25 [1] CRAN (R 4.4.0)
 generics       0.1.3     2022-07-05 [1] CRAN (R 4.4.0)
 ggplot2      * 3.5.1     2024-04-23 [1] CRAN (R 4.4.0)
 ggpubr       * 0.6.0     2023-02-10 [1] CRAN (R 4.4.0)
 ggsankey     * 0.0.99999 2025-03-05 [1] Github (davidsjoberg/ggsankey@b675d0d)
 ggsignif       0.6.4     2022-10-13 [1] CRAN (R 4.4.0)
 glmnet         4.1-8     2023-08-22 [1] CRAN (R 4.4.0)
 glue           1.8.0     2024-09-30 [1] CRAN (R 4.4.1)
 gtable         0.3.6     2024-10-25 [1] CRAN (R 4.4.1)
 haven        * 2.5.4     2023-11-30 [1] CRAN (R 4.4.0)
 hms            1.1.3     2023-03-21 [1] CRAN (R 4.4.0)
 htmltools      0.5.8.1   2024-04-04 [1] CRAN (R 4.4.0)
 htmlwidgets    1.6.4     2023-12-06 [1] CRAN (R 4.4.0)
 httpuv         1.6.15    2024-03-26 [1] CRAN (R 4.4.0)
 iterators      1.0.14    2022-02-05 [1] CRAN (R 4.4.0)
 JM           * 1.5-2     2022-08-08 [1] CRAN (R 4.4.0)
 jomo           2.7-6     2023-04-15 [1] CRAN (R 4.4.0)
 later          1.3.2     2023-12-06 [1] CRAN (R 4.4.0)
 lattice        0.22-6    2024-03-20 [2] CRAN (R 4.4.0)
 libcoin        1.0-10    2023-09-27 [1] CRAN (R 4.4.0)
 lifecycle      1.0.4     2023-11-07 [1] CRAN (R 4.4.0)
 lme4         * 1.1-35.3  2024-04-16 [1] CRAN (R 4.4.0)
 lubridate    * 1.9.3     2023-09-27 [1] CRAN (R 4.4.0)
 magrittr       2.0.3     2022-03-30 [1] CRAN (R 4.4.0)
 MASS         * 7.3-60.2  2024-04-05 [1] local
 Matrix       * 1.7-0     2024-03-22 [2] CRAN (R 4.4.0)
 matrixStats    1.3.0     2024-04-11 [1] CRAN (R 4.4.0)
 memoise        2.0.1     2021-11-26 [1] CRAN (R 4.4.0)
 mice         * 3.16.0    2023-06-05 [1] CRAN (R 4.4.0)
 miceadds     * 3.17-44   2024-01-09 [1] CRAN (R 4.4.0)
 mime           0.12      2021-09-28 [1] CRAN (R 4.4.0)
 miniUI         0.1.1.1   2018-05-18 [1] CRAN (R 4.4.0)
 minqa          1.2.7     2024-05-20 [1] CRAN (R 4.4.0)
 mitml          0.4-5     2023-03-08 [1] CRAN (R 4.4.0)
 mitools        2.4       2019-04-26 [1] CRAN (R 4.4.1)
 mmrm         * 0.3.11    2024-03-05 [1] CRAN (R 4.4.0)
 modeltools     0.2-23    2020-03-05 [1] CRAN (R 4.4.0)
 multcomp       1.4-26    2024-07-18 [1] CRAN (R 4.4.0)
 munsell        0.5.1     2024-04-01 [1] CRAN (R 4.4.0)
 mvtnorm        1.2-5     2024-05-21 [1] CRAN (R 4.4.0)
 nlme         * 3.1-164   2023-11-27 [2] CRAN (R 4.4.0)
 nloptr         2.0.3     2022-05-26 [1] CRAN (R 4.4.0)
 nnet           7.3-19    2023-05-03 [2] CRAN (R 4.4.0)
 pan            1.9       2023-12-07 [1] CRAN (R 4.4.0)
 pbmcapply    * 1.5.1     2022-04-28 [1] CRAN (R 4.4.0)
 pillar         1.10.1    2025-01-07 [1] CRAN (R 4.4.1)
 pkgbuild       1.4.4     2024-03-17 [1] CRAN (R 4.4.0)
 pkgconfig      2.0.3     2019-09-22 [1] CRAN (R 4.4.0)
 pkgload        1.3.4     2024-01-16 [1] CRAN (R 4.4.0)
 plotrix      * 3.8-4     2023-11-10 [1] CRAN (R 4.4.0)
 profvis        0.3.8     2023-05-02 [1] CRAN (R 4.4.0)
 promises       1.3.0     2024-04-05 [1] CRAN (R 4.4.0)
 purrr        * 1.0.4     2025-02-05 [1] CRAN (R 4.4.1)
 R6             2.6.1     2025-02-15 [1] CRAN (R 4.4.1)
 rbibutils      2.2.16    2023-10-25 [1] CRAN (R 4.4.0)
 Rcpp           1.0.12    2024-01-09 [1] CRAN (R 4.4.0)
 Rdpack         2.6       2023-11-08 [1] CRAN (R 4.4.0)
 readr        * 2.1.5     2024-01-10 [1] CRAN (R 4.4.0)
 readxl       * 1.4.3     2023-07-06 [1] CRAN (R 4.4.0)
 remotes        2.5.0     2024-03-17 [1] CRAN (R 4.4.1)
 rlang          1.1.5     2025-01-17 [1] CRAN (R 4.4.1)
 rpart          4.1.23    2023-12-05 [2] CRAN (R 4.4.0)
 rstatix      * 0.7.2     2023-02-01 [1] CRAN (R 4.4.0)
 rstudioapi     0.16.0    2024-03-24 [1] CRAN (R 4.4.0)
 sandwich       3.1-1     2024-09-15 [1] CRAN (R 4.4.1)
 scales         1.3.0     2023-11-28 [1] CRAN (R 4.4.0)
 sessioninfo  * 1.2.2     2021-12-06 [1] CRAN (R 4.4.0)
 shape          1.4.6.1   2024-02-23 [1] CRAN (R 4.4.0)
 shiny          1.8.1.1   2024-04-02 [1] CRAN (R 4.4.0)
 stringi        1.8.4     2024-05-06 [1] CRAN (R 4.4.0)
 stringr      * 1.5.1     2023-11-14 [1] CRAN (R 4.4.0)
 survival     * 3.5-8     2024-02-14 [2] CRAN (R 4.4.0)
 TH.data        1.1-2     2023-04-17 [1] CRAN (R 4.4.0)
 tibble       * 3.2.1     2023-03-20 [1] CRAN (R 4.4.0)
 tidyr        * 1.3.1     2024-01-24 [1] CRAN (R 4.4.0)
 tidyselect     1.2.1     2024-03-11 [1] CRAN (R 4.4.0)
 tidyverse    * 2.0.0     2023-02-22 [1] CRAN (R 4.4.0)
 timechange     0.3.0     2024-01-18 [1] CRAN (R 4.4.0)
 tzdb           0.4.0     2023-05-12 [1] CRAN (R 4.4.0)
 urlchecker     1.0.1     2021-11-30 [1] CRAN (R 4.4.0)
 usethis        2.2.3     2024-02-19 [1] CRAN (R 4.4.0)
 vctrs          0.6.5     2023-12-01 [1] CRAN (R 4.4.0)
 withr          3.0.2     2024-10-28 [1] CRAN (R 4.4.1)
 xtable         1.8-4     2019-04-21 [1] CRAN (R 4.4.0)
 zoo          * 1.8-12    2023-04-13 [1] CRAN (R 4.4.0)
```

## Access and permissions
The creation and maintenance of this research repository are the responsibilities of the author (Daphne Weemering). This archive is completely open access, and can be accessed for an indefinite period.

## Contact 
You can reach me via email at <d.n.weemering@umcutrecht.nl> or <dnweemering@gmail.com> for questions regarding this repository. For correspondence regarding the article, contact Dr. Ruben van Eijk (corresponding author) at <r.p.a.vaneijk-2@umcutrecht.nl>.
