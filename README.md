# Introduction

The code is for the paper "Differentially Private Bootstrap: New Privacy Analysis and Inference Strategies" by Zhanyu Wang, Guang Cheng, and Jordan Awan. https://arxiv.org/abs/2210.06140

# Run codes and draw figures
We use `R` to perform all computation and visualization.

To reproduce the table 2 in our paper, please run `Table2.r` and the results are shown in the `txt` files in `./results/`.

To reproduce the corresponding figures in our paper, please run the files `Figure*.r`, and the results will be `pdf` files in the current folder `./`. The figure files should be the same as the files in `./figures/`. For those experiments taking longer time to run, e.g., Figure 5, we split their procedure into two steps and provide intermediate results saved in the folders `./results*`. To save time, one can skip the first step and use the files `Figure*.r` (the second step) to draw the figures.
For how to reproduce the intermediate results in `./results*`, i.e., the `csv` files, see the details in the next section.

# Obtain intermediate results

## Figure 5: Coverage check of private CI with μ-GDP where μ = 1, 0.5, 0.3, 0.1.
Please run `coverage_check.r`, and the results will be in `./results/`.

## Figure 6: Results of 90% CIs for the slope parameters in logistic regression and quantile regression.
Please run `logistic_regression.R` and `quantile_regression`, and the results will be in `./results_logistic/` and `./results_quantile/` correspondingly. 

## Figure 7: Results for various choices of B and μ-GDP for the population mean inference.
Please run `clamp_normal_mean.r`, and the results will be in `./results_clamp_normal_mean/`.

## Figure 8 and 9: Results for various choices of B and μ-GDP for logistic regression with c = 1 and 0.01.
Please run `logistic_regression.r`, and the results will be in `./results_logistic/`.

## Figure 10 and 11: Results for various choices of B and μ-GDP for quantile regression with c = 1 and 0.01.
Please run `quantile_regression.r`, and the results will be in `./results_quantile/`.

## License
Copyright 2023 Zhanyu Wang, Guang Cheng, and Jordan Awan

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

