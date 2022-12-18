module Network_qse
using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end



export testing
#"""
# https://docs.julialang.org/en/v1/manual/unicode-input/
#"""
__precompile__(false)

using Distributed

include("dependencies.jl")
include("IO.jl")
include("Constants.jl")
include("DataTypes.jl")
include("Boltzmann.jl")
include("Tools.jl")

function solve_NSE(yrange::Array{Float64,1}, trange::Array{Float64,1}, rrange::Array{Float64,1})
    srange = Array{Float64, 1}(undef, size(trange, 1))
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 4}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1))
    tmp = Network_qse.initial_guess(a)
    c12 = find_el("C12", a)
    print(c12)
    sp = Network_qse.StepParameter(-50,50,30)

    for (i,y) in enumerate(yrange)
        for (j, t) in enumerate(trange)
            for (k, r) in enumerate(rrange)
                any(isnan, tmp) ? tmp = initial_guess(a) : nothing
                th = Network_qse.ThermoProperties(t, r, y, -5.0)
                ff = Network_qse.Func(2, x -> Network_qse.nse_condition(x, th, a), x -> Network_qse.df_nse_condition(x,th,a), false)
                any(isnan.(tmp)) ? tmp = Network_qse.initial_guess(a) : nothing
                tmp = Network_qse.MultiNewtonRaphson(tmp, ff, th, a, sp)
                res[:,i,j,k] = Network_qse.x_i(tmp, th, a)
                srange[j] = sum(res[:,i,j,k][c12:end])
                println(">>>> ", j, " ", " sum ",sum(res[:,i,j,k]))
            end
        end
    end
    save("output/NSE_table.jld", "data", res)
    save("output/NSE_params.jld", "trange", trange, "yrange", yrange, "rrange", rrange, "srange", srange)
    open("output/NSE_README.txt"; write=true) do f
        write(f, "# Netwon Raphson parametrization: $(sp)\n")
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# load parameters with load(\"data.jld\")[\"yrange\"]\n")
        write(f, "# y-range $(size(yrange)), T-range $(size(trange)), rho-range $(size(rrange)),  s-range $(size(srange))\n")
        writedlm(f, [yrange, "\n",  trange, "\n", rrange, "\n", srange])
    end
    return res, srange
end


function parallel_NSE(i::Int64, yrange::Array{Float64,1}, trange::Array{Float64,1}, rrange::Array{Float64,1})
    srange = Array{Float64, 1}(undef, size(trange, 1))
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 4}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1))
    #res = Array{Float64, 4}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1))
    tmp = Network_qse.initial_guess(a)
    c12 = find_el("C12", a)
    sp = Network_qse.StepParameter(-50,50,30)
    y = yrange[i]
    #println(nprocs())
    #@distributed
    #for (i,y) in enumerate(yrange)
        #addprocs(2)
        for (j, t) in enumerate(trange)
            #println("-----")
            for (k, r) in enumerate(rrange)
                any(isnan, tmp) ? tmp = initial_guess(a) : nothing
                th = Network_qse.ThermoProperties(t, r, y, -5.0)
                ff = Network_qse.Func(2, x -> Network_qse.nse_condition(x, th, a), x -> Network_qse.df_nse_condition(x,th,a), false)
                any(isnan.(tmp)) ? tmp = Network_qse.initial_guess(a) : nothing
                tmp = Network_qse.MultiNewtonRaphson(tmp, ff, th, a, sp)
                res[:,i,j,k] = Network_qse.x_i(tmp, th, a)
                srange[j] = sum(res[:,i,j,k][c12:end])
                println(">>>> ", j, " ", " sum ",sum(res[:,i,j,k]))
                #save("./testing_out", "data", i)
            end
        end
    #end
    save("./NSE_table.jld", "data", res)
    save("./NSE_params.jld", "trange", trange, "yrange", yrange, "rrange", rrange, "srange", srange)
    open("./NSE_README.txt"; write=true) do f
        write(f, "# Netwon Raphson parametrization: $(sp)\n")
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# load parameters with load(\"data.jld\")[\"yrange\"]\n")
        write(f, "# y-range $(size(yrange)), T-range $(size(trange)), rho-range $(size(rrange)),  s-range $(size(srange))\n")
        writedlm(f, [yrange, "\n",  trange, "\n", rrange, "\n", srange])
    end
    return res, srange
end


function parallel_QSE(ARGS, yrange::Array{Float64,1}, trange::Array{Float64,1}, rrange::Array{Float64,1}, srange::Array{Float64,1})
    #i = parse(Int64, ARGS[1]) + 1 
    k  = parse(Int64, ARGS[1]) + 1 
    #l = parse(Int64, ARGS[1]) + 1 
    #j = parse(Int64, ARGS[1]) + 1 
    a  = Network_qse.extract_partition_function()
    res = Array{Float64, 5}(undef, size(a,1), size(yrange, 1), 1, 1, size(srange,1))
    #res = Array{Float64, 5}(undef, size(a,1), size(yrange, 1), size(trange,1), 1, size(srange,1))
    #res_mu = Array{Float64, 5}(undef, 3, size(yrange, 1), size(trange, 1), size(rrange, 1), size(srange, 1))
    sp = Network_qse.StepParameter(-10,10,50)
    #y = yrange[i]
    #t = trange[j]
    r = rrange[k]
    #q = srange[l]
    l_t = size(trange,1)
    l_y = size(yrange,1)
    l_r = size(rrange,1)
    l_s = size(srange,1)
    #for (k, r) in enumerate(rrange)
    	tmp = Network_qse.qse_initial_guess(a, trange[1], rho = rrange[1])
        tmp1 = tmp 
        for (j, t) in enumerate(trange)
#            mkdir("output/r$(k-1)")
#            mkdir("output/r$(k-1)/t$j")
            for  (l, q) in enumerate(srange)
                tmp = tmp1
    		for (i,y) in enumerate(yrange)
		    print("ye",y,i,"\n")
                    th = Network_qse.ThermoProperties(t, r, y, q)
                    ff = Network_qse.Func(3, x -> Network_qse.qse_condition(x, th, a), x -> Network_qse.df_qse_condition(x,th,a), false)
                    any(isnan.(tmp)) ? tmp = Network_qse.qse_initial_guess(a, t) : nothing
                    tmp = Network_qse.MultiNewtonRaphson(tmp, ff, th, a, sp)
                    x_nse, x_qse = Network_qse.x_i_QSE(tmp, th, a)
                    res[:,i,1,1,l] = vcat(x_nse,x_qse)
                    #println("PID", rpad(getpid(),10), "T: ", rpad(string(t),20), "rho: ",rpad(string(r),20), "sum X_cl: ", rpad(string(srange[l]),20), "y_e: ", rpad(string(y),20), "mu: ", rpad(string.(tmp),20))
                   if (i == 1) && (!any(isnan.(tmp))) 
                        tmp1 = tmp
                    end 
                end
            end
     		    print("\n ------ new temp directory -------- ", k-1, "\n")
    		    save("./output/r$(k-1)/t$j/QSE_table.jld", "data", res)
     	            save("./output/r$(k-1)/t$j/params.jld", "trange", trange, "t-value", t, "yrange", yrange, "srange", srange, "rrange", rrange[k])
    	            open("./output/r$(k-1)/t$j/QSE_README.txt"; write=true) do f
        	    	write(f, "# Netwon Raphson parametrization: $(sp)\n")
                    	#write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
                    	write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
                    	write(f, "# load parameters with load(\"data.jld\")[\"yrange\"]\n")
                    	write(f, "# y-range $(size(yrange)), T-range $(size(trange)), rho-range $(size(rrange)),  s-range $(size(srange))\n")
                    	writedlm(f, [yrange, "\n",  t, "\n", r, "\n", srange])
                    end
 
        end
    #end
    println("PID", rpad(getpid(),10), "rho: ",rpad(string(r),20), " mu ", tmp)
    println("\n\n --- finished ---")
end

#parallel_QSE(ARGS, collect(LinRange(0.5,0.497,5)), collect(LinRange(4.6e9, 4.4e9, 5)), collect(LinRange(1e8, 7.12e8, 16)), collect(LinRange(-1.5, -1.75, 5)));
#parallel_QSE(ARGS, collect(LinRange(0.5,0.47,75)), collect(LinRange(5e9, 3e9, 75)), collect(LinRange(1e8, 1e9, 64)), collect(LinRange(-1.5, -4, 75)));


function solve_QSE(yrange::Array{Float64,1}, trange::Array{Float64,1}, rrange::Array{Float64,1}, srange::Array{Float64,1})
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 5}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1), size(srange, 1))
    res_mu = Array{Float64, 5}(undef, 3, size(yrange, 1), size(trange, 1), size(rrange, 1), size(srange, 1))
    tmp = Network_qse.qse_initial_guess(a,trange[1])
    sp = Network_qse.StepParameter(-10,10,50)
    l_t = size(trange,1)
    l_y = size(yrange,1)
    l_r = size(rrange,1)
    l_s = size(srange,1)
    for (l, q) in enumerate(srange)
        tmp = Network_qse.qse_initial_guess(a,trange[1]; rho = rrange[1])
        for (j, t) in enumerate(trange)
            for (k, r) in enumerate(rrange), (i,y) in enumerate(yrange)
                    th = Network_qse.ThermoProperties(t, r, y, q)
                    ff = Network_qse.Func(3, x -> Network_qse.qse_condition(x, th, a), x -> Network_qse.df_qse_condition(x,th,a), false)
                    any(isnan.(tmp)) ? tmp = Network_qse.qse_initial_guess(a,t) : nothing
                    tmp = Network_qse.MultiNewtonRaphson(tmp, ff, th, a, sp)
                    x_nse, x_qse = Network_qse.x_i_QSE(tmp, th, a)
                    res_mu[:,i,j,k,l] = tmp
                    res[:,i,j,k,l] = vcat(x_nse,x_qse)
                    println("sum X: ",sum(res[:,i,j,k,l]), "   sum X_cl: ", sum(x_qse))
                    println("y:    ", y, "   T:    ", t, " ", " r:     ",r, "   sum X_cl: ", 1.0 - 10.0^(srange[l]))
            end
        end
    end
    save("output/QSE_table.jld", "data", res)
    save("output/QSE_params.jld", "trange", trange, "yrange", yrange, "rrange", rrange, "srange", 1.0 .- 10.0.^(srange))
    open("output/QSE_README.txt"; write=true) do f
        write(f, "# Netwon Raphson parametrization: $(sp)\n")
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# load parameters with load(\"data.jld\")[\"yrange\"]\n")
        write(f, "# y-range $(size(yrange)), T-range $(size(trange)), rho-range $(size(rrange)),  s-range $(size(srange))\n")
        writedlm(f, [yrange, "\n",  trange, "\n", rrange, "\n", 1.0 .- 10.0.^(srange)])
    end
    return res, 1.0 .- 10.0.^(srange)
end

















end
