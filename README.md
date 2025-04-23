# PTDBoot

PTDBoot is an R package containing methods conducting leveraging large scale machine learning predictions in downstream statistical inference when a small number of labelled calibration samples are available. The package has functions that implement the Predict-Then-Debias bootstrap algorithms and its extensions from the paper referenced below. 

# Installation and tutorial

You can install the development version of PTDBoot directly from GitHub: 

```r
# install.packages("devtools") # if you don't have it already
devtools::install_github("DanKluger/PTDBoot")
```

👉 [**View the tutorial here**](https://dankluger.github.io/PTDBootTutorial/Tutorial.html)

# Citations

If you use PTDBoot in your work, please cite the following paper:

> Dan M. Kluger, Kerri Lu, Tijana Zrnic, Sherrie Wang, Stephen Bates (2025). Prediction-Powered Inference with Imputed Covariates and Nonuniform Sampling.  
> 2501.18577 [stat.ME]  
> [https://arxiv.org/abs/2504.12345](https://arxiv.org/abs/2501.18577)](https://arxiv.org/abs/2501.18577)

# Questions and Reporting Bugs

Please report any bugs or direct any questions to dkluger@mit.edu.
