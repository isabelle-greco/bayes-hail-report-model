# Figures

Each figure in the paper and each supplementary figure in the supplementary materials has here its own independent notebook, labeled accordingly.
Each notebook is designed to be run using the same R installation as previously (`NCI-data-analysis/2023.02` via project `dk92`).
As a reminder, the figures in the paper are:

* [Figure 1](fig1_map.ipynb): Map of the study area in south-east Queensland and north-east New South Wales. The study area (blue) encompasses the Queensland state capital (Brisbane) as well as several smaller urban areas indicated in the map by points sized according to their population (Australian Bureau of Statistics 2021). The location of the Mount Stapylton radar is also highlighted (red triangle) along with its 135 km range ring highlighted (pink).

* [Figure 2](fig2_sim_func.ipynb): Hail and reporting functions estimated by each of the three experiments on one synthetic dataset. The functions used to generate the data are given by the black dashed line and the 50% and 95% credible intervals are shaded. The corresponding values of the LOOCV ELPD estimate are given in each figure.

* [Figure 3](fig3_sim_params.ipynb): The prior (dashed) and posterior (solid) distributions of the four parameters of one LHHR simulation experiment. The top row gives the two parameters of the hail function and the bottom those of the reporting function in this simulation.

* [Figure 4](fig4_hail_function.ipynb): Above, the hail function for the best performing model as determined by the ELPD. Below, the distribution of the value of MESH at which severe hail becomes more likely than not with the different prior perturbations. The dashed line indicates a 30 mm MESH threshold. References to strength in the prior names refer to informativeness of the prior (e.g. weakly informative, strongly informative, etc.).

* [Figure 5](fig5_reporting_rate.ipynb): Left, the posterior mean reporting rate across the region. Coastlines and the New South Wales-Queensland state border are shown in red. Right, the reporting function from the best fitting model to illustrate the uncertainty and possibility of a plateau.

* [Figure 6](fig6_climatology.ipynb): The posterior mean estimate of the average number of annual hail days over the region with a 95\% credible interval (CI), in comparison with three threshold-based estimates of the same quantity. Again the coastlines and state border are shown in red with three regional centres referenced in other regional climatologies marked in white.

Further, the supplementary figures are:

* [Supplementary Figure 1](supp_fig01_pop_reports.ipynb): The population density of Australia from the 2016 Census (top) and the number of reports in the SSA from Jan 2010 - Apr 2016 (bottom). For clarity, only 1 km$^2$ grid cells with density over 1 person km$^{-2}$ are shown. Note the similarities between the two.

* [Supplementary Figure 2](supp_fig02_emp_evidence.ipynb): An example of the applicability of the inverse logit function when considering the approximate empirical probability of hail as a function of population density. A 30 mm MESH threshold was used in this figure to define hail events. Results for other MESH thresholds are similar qualitatively.

* [Supplementary Figure 3](supp_fig03_bulk_sims.ipynb): The distribution of posterior mean hail and reporting functions estimated from 100 simulations of synthetic data. The true probability is indicated by the dotted lines in each figure and 50% and 95% credible intervals are also shown as more and less opaque shaded regions respectively.

* [Supplementary Figure 4](supp_fig04_sim_params.ipynb): The prior (dashed) and posterior (solid) distributions of the four parameters of one HHLR experiment. The top row gives the two parameters of the hail function and the bottom of the reporting function.

* [Supplementary Figure 5](supp_fig05_elpd_diff.ipynb): The difference in pointwise contributions to the ELPD between the best performing model and the model with the same priors on the Yeo-Johnson parameter but a filter on the expected probability of hail. The best performing model is superior to the filtered model when the difference is positive.

* [Supplementary Figure 6](supp_fig06_thresholds.ipynb): The optimal threshold when all cells above a given population density are considered, as determined by the Heidke skill score. Every threshold possible given the data was attempted though only the largest nine are shown here.

* [Supplementary Figure 7](supp_fig07_logistic.ipynb): Fitting a logistic regression model to all cells above a given population density. The points used are jittered and shown along with the fitted curve. The dashed lines indicate the optimal threshold as per Supplementary Figure~\ref{fig-supp-hss} and its corresponding hail probability under the logistic model. Note the differences in the fitted curve between the last three plots and the remaining plots.

* [Supplementary Figure 8](supp_fig08_expected_error.ipynb): The expected model residuals (difference in total expected less total observed reports), stratified by weekdays and weekends. Note the out of sample cross-validation predictions are used. The location of the radar is marked in yellow to highlight the lack of range dependence in errors. Note the lighter cells west of Brisbane and north-west of the radar site perhaps showing under-prediction over the highway to Toowoomba, particularly on weekends.


* [Supplementary Figure 9](supp_fig09_sim_priors.ipynb): The priors utilised in the LHHR and HHLR experiments. The opacity indicates 50%, 95%, and 99% credible intervals in order of decreasing opacity.

* [Supplementary Figure 10](supp_fig10_yjnarrow_report.ipynb): The ten different prior report functions when the strongly informative prior was applied to the Yeo-Johnson transform parameters.

* [Supplementary Figure 11](supp_fig11_yjnarrow_hail.ipynb): The ten different prior hail functions when the strongly informative prior was applied to the Yeo-Johnson transform parameters.

* [Supplementary Figure 12](supp_fig12_yjnormal_report.ipynb): The ten different prior report functions when the weakly informative prior was applied to the Yeo-Johnson transform parameters.

* [Supplementary Figure 13](supp_fig13_yjnormal_hail.ipynb): The ten different prior hail functions when the weakly informative prior was applied to the Yeo-Johnson transform parameters.

* [Supplementary Figure 14](supp_fig14_yjwide_report.ipynb): The ten different prior report functions when the very weakly informative prior was applied to the Yeo-Johnson transform parameters.

* [Supplementary Figure 15](supp_fig15_yjwide_hail.ipynb): The ten different prior hail functions when the very weakly informative prior was applied to the Yeo-Johnson transform parameters.