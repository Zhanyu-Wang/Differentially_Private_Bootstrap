# Introduction

The code is for the paper "Differentially Private Bootstrap: New Privacy Analysis and Inference Strategies" by Zhanyu Wang, Guang Cheng, and Jordan Awan. https://arxiv.org/abs/2210.06140

# Run codes and draw figures
We use `R` and `Python` to perform all computation and visualization.

## Figure 4: Comparison between the privacy profile δ(ε) of composition computed from asymptotics (Theorem 16) and numerical evaluation (Proposition 14). 
Please run `Figure4.sh`, and the figures will be in `./figures/Figure4/`.

## Other figures 
Please run the files `Figure*.r`, and the results will be `pdf` files in the current folder `./`. The figure files should be the same as the files in `./figures/`. For those experiments taking longer time to run, e.g., Figure 6, we split their procedure into two steps and provide intermediate results saved in the folders `./results*`. To save time, one can skip the first step and use the files `Figure*.r` (the second step) to draw the figures.
For how to reproduce the intermediate results in `./results*`, i.e., the `csv` files, see the details in the next section.

# Obtain intermediate results

## Figure 6: Coverage check of private CI with μ-GDP where μ = 1, 0.5, 0.3, 0.1.
Please run `coverage_check.r`, and the results will be in `./results/`.

## Figure 7: Coverage and width of CIs for the population mean with different privacy guarantees. 
Please run `Figure7-compute.r`, and the results will be in `./results/`.

## Figure 8: Results of 90% CIs for the slope parameters in logistic regression and quantile regression.
Please run `logistic_regression.r` and `quantile_regression.r`, and the results will be in `./results_logistic/` and `./results_quantile/` correspondingly. 

## Supplementary Figure 1: Results for various choices of B and μ-GDP for the population mean inference.
Please run `clamp_normal_mean.r`, and the results will be in `./results_clamp_normal_mean/`.

## Supplementary Figure 2 and 3: Results for various choices of B and μ-GDP for logistic regression with c = 1 and 0.01.
Please run `logistic_regression.r`, and the results will be in `./results_logistic/`.

## Supplementary Figure 4 and 5: Results for various choices of B and μ-GDP for quantile regression with c = 1 and 0.01.
Please run `quantile_regression.r`, and the results will be in `./results_quantile/`.

## Supplementary Figure 6: Comparison between the privacy profile δ(ε) of composition computed from asymptotics (Theorem 16) and numerical evaluation (Proposition 14).
Please run `Figure4.sh` (since this Figure shows more of the settings used in Figure 4), and the figures will be in `./figures/Figure4/`.

## License
Copyright 2024 Zhanyu Wang, Guang Cheng, and Jordan Awan

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

