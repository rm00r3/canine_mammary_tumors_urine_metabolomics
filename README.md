# Mammary Urine Analysis (Refactored)

This project refactors the four model scripts into a single modular pipeline with a top-level driver.

## Project Layout

- `analysis_driver.R`
- `R/` (shared helper modules)
- `outputs/model1/`
- `outputs/model2/`
- `outputs/model3/`
- `outputs/model4/`

## Required Input Files

By default, the driver expects:

- `~/Desktop/Mammary_urine/Mammary_urine/data_raw/newjunemsms_allmodes_msms_data_new.csv`
- `~/Desktop/Mammary_urine/Mammary_urine/data_raw/Meta_all_fixed_test.csv`

You can override these with options before sourcing the driver:

- `options(mammary.data_file = "/absolute/path/to/data.csv")`
- `options(mammary.meta_file = "/absolute/path/to/meta.csv")`

## How To Run

From a fresh R session:

```r
setwd("/Users/Moorebid/Desktop/Mammary_urine/Mammary_urine_cleaned")
source("analysis_driver.R")
```

Choose model(s) using `mammary.models`:

- Model 1 only:

```r
options(mammary.models = "model1")
source("analysis_driver.R")
```

- Model 2 only:

```r
options(mammary.models = "model2")
source("analysis_driver.R")
```

- Model 3 only:

```r
options(mammary.models = "model3")
source("analysis_driver.R")
```

- Model 4 only:

```r
options(mammary.models = "model4")
source("analysis_driver.R")
```

- All models:

```r
options(mammary.models = "all")
source("analysis_driver.R")
```

Optional output root override:

```r
options(mammary.output_root = "/absolute/path/to/outputs")
```

## Output Behavior

- Original output basenames are preserved within each model folder.
- Each model writes to its own directory so files are not overwritten across models.
- Freedman-Lane permutation testing is retained.
- PLS-DA is removed from this refactored pipeline.
