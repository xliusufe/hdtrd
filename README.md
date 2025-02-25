# hdtest
R package "hdtrd" for calculating p-value of the test statistic for relevant difference in high-dimensional generalized linear regression models and its application to transfer learning. In the paper Liu (2024), we propose novel statistics to test relevant difference of two high dimensional coefficients. The proposed method can serve as the high dimensional transfer learning.

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/hdtrd")

    # or
    #install.packages("remotes")
    library(remotes)
    remotes::install_github("xliusufe/hdtrd")    

# Usage

   - [x] [hdtrd-manual.pdf](https://github.com/xliusufe/hdtrd/blob/master/inst/hdtrd-manual.pdf) ---------- Details of the usage of the package.
# Example
    library(hdtrd)

    ## testing for relevant difference 
    data(simulData_test_gauss)
    pvals <- pvalrd(data = datahb)
    pvals


    ## testing for transferibability
    data(simulData_trans_gauss)
    pvals <- pvaltrans(target = dataset[[1]], source = dataset[-1])
    pvals

    ## estimation for transfer learning
    data(simulData_trans_gauss)
    fit <- utrans(target = dataset[[1]], source = dataset[-1], idtrans = seq(5))
    fit$beta[1:10]


# References
Cui, H., Guo, W. and Zhong, W. (2018). Test for high-dimensional regression coefficients using refitted cross-validation variance estimation. The Annals of Statistics, 46, 958-988.

Chen, Z., Cheng, V. X. and Liu, X. (2024). Hypothesis testing on high dimensional quantile regression. Journal of Econometrics.

Chen, J., Li, Q., and Chen, H. Y. (2022). Testing generalized linear models with highdimensional nuisance parameters. Biometrika, 110. 83-99.

Guo, B.and Chen, S. X. (2016). Tests for high dimensional generalized linear models. Journal of the Royal Statistical Society, Series B, 78, 1079-1102.

Karoui, N, E. (2008) Spectrum estimation for large dimensional covariance matrices using random matrix theory. The Annals of Statistics, 36(6), 2757-2790.

Kong, W. and Valiant, G. (2017). Spectrum estimation from samples. Annals of Statistics. 45, 2218-2247.

Li, S., Cai, T. T., and Li, H. (2022a). Transfer Learning in Large-Scale Gaussian Graphical Models with False Discovery Rate Control. Journal of the American Statistical Association, 118, 2171-2183.

Li, S., Cai, T. T., and Li, H. (2022b). Transfer Learning for High-Dimensional Linear Regression: Prediction, Estimation and Minimax Optimality. Journal of the Royal Statistical Society Series B, 84, 149-173.

Liu, S. (2024). Unified Transfer Learning Models for High-Dimensional Linear Regression. Proceedings of The 27th International Conference on Artificial Intelligence and Statistics, PMLR. 238, 1036-1044.

Liu, X. (2024). Testing relevant difference in high-dimensional linear regression with applications to detect transferability. Manuscript.

Liu, X., Zheng, S. and Feng, X. (2020). Estimation of error variance via ridge regression. Biometrika. 107, 481-488.

Tian, Y. and Feng, Y. (2023) Transfer Learning Under High-Dimensional Generalized Linear Models. Journal of the American Statistical Association, 118, 2684-2697.

Yang, W., Guo, X. and Zhu, L. (2023). Score function-based tests for ultrahigh-dimensional linear models. arXiv:2212.08446.

Zhang, X. and Cheng, G. (2017). Simultaneous inference for high-dimensional linear models. Journal of the American Statistical Association, 112, 757-768.

Van den Meersche, K., Soetaert, K., and Van Oevelen, D. (2009). xsample(): An R Function for Sampling Linear Inverse Problems. Journal of Statistical Software, Code Snippets, 30, 1-15.

# Development
This R package is developed by Xu Liu (liu.xu@sufe.edu.cn).
