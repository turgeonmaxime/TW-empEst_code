# TW-empEst code

Running the code in this repository will create datasets and figures that will be stored in `cache` and `figures`, respectively.

To generate the datasets, run the following command from the root directory:

```r
Rscript distribution_approx/generate_data.R
```

Then, to generate the figures, run the following command (still from the root directory):

```r
Rscript distribution_approx/analyse_results.R
```
