using BayesianMixtures
using CSV
using DataFrames
using DelimitedFiles
using Plots
using Statistics

cd()
cd("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/")

data_mixture = CSV.read("data_MFM.csv", DataFrame)

length_data = size(data_mixture)[1]

x = data_mixture[1:length_data,1]

n_total = 50000

n_burn = 1000

# n_total = 1000000  # total number of MCMC sweeps to run
# log_pk = "k -> log(k in (1:30) ? 1/30 : 0)"   # log prior on k is Uniform{1,...,30}
# options = BayesianMixtures.options("Normal","MFM",x,n_total; n_keep=5000,log_pk=log_pk,t_max=30)

options = BayesianMixtures.options("Normal","MFM",x,n_total,n_burn = n_burn, n_keep=5000, t_max = 20)

result = BayesianMixtures.run_sampler(options)

cd()
cd("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/")

# cd()
# cd("Dropbox/DSS_MIX/DSS_mixture_02_13_2023/dssmix_functions_11_2023/MFM_julia_test_09_14_2023/")

temp = result.keepers
result_weight = result.N[:,temp]

#writedlm("result_weight_MFM.txt", result_weight)
#writedlm("result_t_MFM.txt", result.t)
#writedlm("result_theta_MFM.txt", result.theta)

writedlm("output_MFM_DPM/result_theta_MFM.txt", result.theta)
writedlm("output_MFM_DPM/result_weight_MFM.txt", result_weight)
writedlm("output_MFM_DPM/result_t_MFM.txt", result.t)

#--------------------------------------#

min_x = minimum(x) - 0.5*sqrt(var(x))
max_x = maximum(x) + 0.5*sqrt(var(x))

# min.y<-min(y.data)-0.5*sqrt(var(y.data))
# max.y<-max(y.data)+0.5*sqrt(var(y.data))

# min_x = 10
# max_x = 40

xs = range(min_x,stop=max_x,length=500)
fs = Float64[BayesianMixtures.density_estimate(xi,result) for xi in xs]
writedlm("output_MFM_DPM/result_dens_MFM.txt", fs)

# n_total = 1000000  # total number of MCMC sweeps to run
# log_pk = "k -> log(k in (1:30) ? 1/30 : 0)"   # log prior on k is Uniform{1,...,30}
# options = BayesianMixtures.options("Normal","MFM",x,n_total; n_keep=5000,log_pk=log_pk,t_max=30)
