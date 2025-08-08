using BayesianMixtures
using CSV
using DataFrames
using DelimitedFiles
using LinearAlgebra
using Statistics
B = BayesianMixtures

#cd()
#cd("Dropbox/DSS_MIX/DSS_mixture_02_13_2023/functions_paper_8_10_2023/")

cd()
cd("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/")

data_mixture = CSV.read("data_MFM_mult.csv", DataFrame)
length_data = size(data_mixture)[1]
dim_data = size(data_mixture)[2]

Mat_temp = Matrix(data_mixture)
x = [zeros(dim_data)::Array{Float64,1} for i = 1:length_data]

for i in 1:length_data
  x[i] = Mat_temp[i,:]
end

t_max = 30
n_burn = 5000 # 5000
#n_total = 20*n_burn
n_total = 60000
n_keep = 1000

options = B.options("MVN","MFM",x,n_total, n_keep=n_keep, n_burn=n_burn, t_max=t_max)
result = B.run_sampler(options)

weights_matrix = result.N[:,result.keepers]
number_fact = size(weights_matrix)[1]

Z = zeros(number_fact,dim_data*n_keep)
for j in 1:n_keep
  for c in findall(weights_matrix[:,j].>0)
    Z[c,(dim_data*j-(dim_data-1)):(dim_data*j)] = result.theta[c,j].m
  end
end

M = zeros(number_fact,(dim_data^2)*n_keep)
for j in 1:n_keep
  for c in findall(weights_matrix[:,j].>0)
    M[c,(dim_data^2*j-((dim_data^2)-1)):((dim_data^2)*j)] = reshape(inv(result.theta[c,j]._R),((dim_data^2),1))
  end
end

#cd()
#cd("Dropbox/DSS_MIX/DSS_mixture_02_13_2023/functions_paper_8_10_2023/")

cd()
cd("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/")

n = length_data
writedlm("output_MFM_DPM_mult/result_weight_mult.txt", weights_matrix/n)

writedlm("output_MFM_DPM_mult/result_mult_mu.txt", Z)

writedlm("output_MFM_DPM_mult/result_mult_sigma2.txt", M)

writedlm("output_MFM_DPM_mult/result_t_mult.txt", result.t)

writedlm("output_MFM_DPM_mult/result_mult_class_MFM.txt", result.z)


#x = result.options.x
#z = result.z[:,1000][:]
#B.plot_clusters(x,z; markersize=5)
