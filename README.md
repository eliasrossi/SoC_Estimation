# SoC_Estimation

Algorithms for SOC Estimation of Ion-Lithium Batteries (https://doi.org/10.1016/j.est.2023.107677).

Function func_EKF performs the estimation using the Extended Kalman Filter

Function func_MHSE performs the estimation using the Moving-Horizon State Estimation. The function costfcn2 is used inside the MHSE for evaluating the cost function of the constrained optimization problem.
