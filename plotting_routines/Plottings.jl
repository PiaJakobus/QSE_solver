a = Network_qse.extract_partition_function()

qse = load("../qse_tables/vary_R/QSE_table.jld")["data"]
nse = load("../NSE_higher_QSE/logScale/NSE_table.jld")["data"]
params_qse = load("../NSE_higher_QSE/logScale/QSE_params.jld")
params_nse = load("../NSE_higher_QSE/logScale/NSE_params.jld")
