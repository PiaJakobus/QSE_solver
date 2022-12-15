
a = Network_qse.extract_partition_function()

qse = load("/home/../../media/user1/My Passport/Network_qse/src/output/r1/t1/QSE_table.jld")["data"]
nse = load("..//Network_qse_run1/NSE_higher_QSE/logScale/NSE_table.jld")["data"]
params_qse = load("../Network_qse_run1/NSE_higher_QSE/logScale/QSE_params.jld")
params_nse = load("../Network_qse_run1/NSE_higher_QSE/logScale/NSE_params.jld")

qse = load("/home/../../media/user1/My Passport/Network_qse/src/output/r1/t1/QSE_table.jld")["data"]
params_qse = load("/home/../../media/user1/My Passport/Network_qse/src/output/r1/t1/params.jld")

params_qse = load("../Network_qse/y0/QSE_params.jld")
qse = load("../Network_qse/QSE_table.jld")["data"]



qse = load("../Network_qse_run1/qse_tables/vary_R/QSE_table.jld")["data"]
qse = load("../Network_qse_run1/NSE_higher_QSE/error/QSE_table.jld")["data"]
nse_alt = load("../Network_qse_run1/qse_tables/vary_R/NSE_table.jld")["data"]
params_qse = load("../Network_qse_run1/qse_tables/vary_R/QSE_params.jld")
params_nse_alt = load("../Network_qse_run1/qse_tables/vary_R/NSE_params2.jld")
params_qse = load("../Network_qse_run1/NSE_higher_QSE/error/QSE_params.jld")



print(qse[Network_qse.find_el("Fe56", a),1,:,1,2])

Plots.plot(cl_qse, qse[Network_qse.find_el("H1", a),10,1,1,:],yscale =:log10)
Plots.plot!(cl_qse, qse[Network_qse.find_el("1n", a),10,1,1,:],yscale =:log10)
Plots.plot!(cl_qse, qse[3,10,1,1,:],yscale =:log10)
Plots.plot!(cl_qse, qse[4,10,1,1,:],yscale =:log10)

cl_qse = params_qse["srange"]
#cl_nse = params_nse["srange"]
trange = params_qse["trange"]
yrange = params_qse["yrange"]
rrange = params_qse["rrange"]


# Julia: fl__ -> x_i, ye, t, rho, s
AexNse(temp::Int64) = sum([nse[i,1,temp,1] * a[i].A for i in Network_qse.find_el("C12", a):size(a,1)])
ZexNse(temp::Int64) = sum([nse[i,1,temp,1,1] * a[i].Z for i in Network_qse.find_el("C12", a):size(a,1)])


AexQse(ye) = sum([qse[i,ye,1,1,5] * a[i].A for i in Network_qse.find_el("C12", a):size(a,1)])
ZexQse(s) = sum([qse[i,10,1,1,s] * a[i].Z for i in Network_qse.find_el("C12", a):size(a,1)])

Plots.plot(yrange,ZexQse.(1:size(yrange, 1)),color=:red)
Plots.plot!(cl_qse,ZexQse.(1:size(cl_qse,1), size(cl_qse,1)),color=:red)
Plots.plot!(yrange,AexQse.(1:size(trange, 1), 2),color=:blue)
Plots.plot!(yrange,AexQse.(1:size(trange, 1), size(cl_qse,1)),color=:blue)


p = Plots.plot(trange, cl_nse, label = "NSE", axis = "T₉", yaxis = "Clustersize", xaxis=:flip,legend = :bottomright, yscale =:log10)
hline!([cl_qse[1], cl_qse[end]], fillrange=[cl_qse[1] cl_qse[end]], alpha = 0.2,label = "QSE", axis = "T₉", yaxis = "Clustersize", xaxis=:flip)
Plots.savefig(p,"../Network_qse_run1/NSE_higher_QSE/logScale/clustersize_new.pdf")


p = Plots.plot(trange, AexNse.(1:size(trange, 1)), label = "<A> NSE", linestyle = :dash, c = :blue,axis = "T₉", yaxis = "<Z>, <A>", xaxis=:flip, legend=:bottomleft)
plot!(trange, AexQse.(1:size(trange, 1), size(cl_qse,1)), label = "<A> QSE", c = :blue, axis = "T₉", yaxis = "<Z>, <A>", xaxis=:flip, legend=:bottomleft)
plot!(trange, ZexNse.(1:size(trange, 1)), label = "<Z> NSE", linestyle = :dash, c = :red, axis = "T₉", yaxis = "<Z>, <A>", xaxis=:flip, legend=:bottomleft)
plot!(trange, ZexQse.(1:size(trange, 1), size(cl_qse,1)), label = "<Z> QSE", c = :red, axis = "T₉", yaxis = "<Z>, <A>", xaxis=:flip, legend=:bottomleft)
plot!([trange trange], [ZexQse.(1:size(trange, 1), 3) ZexQse.(1:size(trange, 1), size(cl_qse,1))],  label = :false, fillrange=[ZexQse.(1:size(trange, 1), size(cl_qse,1)) ZexQse.(1:size(trange, 1), 2)], fillalpha=:0.1, color =:red)
plot!([trange trange], [AexQse.(1:size(trange, 1), 2) AexQse.(1:size(trange, 1), size(cl_qse,1))],  label = :false, fillrange=[AexQse.(1:size(trange, 1), size(cl_qse,1)) AexQse.(1:size(trange, 1), 2)], fillalpha=:0.1, color =:blue)
Plots.savefig("../Network_qse_run1/NSE_higher_QSE/logScale12/compare_below.png")
Plots.savefig("../confirmation/average.pdf")

cl_qse[end]




Plots.plot(trange, qse[Network_qse.find_el("Ni56", a), 1,:,1,end], xlabel = "T₉", ylabel = "Xᵢ",
    c = "blue", dt = "--", label = "QSE Ni62", xaxis =:flip,  title = L"X_{\mathrm{QSE}}=%$(cl_qse[1])",
        yaxis = :log, ylims = (1e-10, 1))
Plots.plot(trange, nse[Network_qse.find_el("Ni56", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", label = "NSE Ni62", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
Plots.savefig("Ni56.png")



ind = 6
for ind in collect(1:100)
    println(ind)
    Plots.plot(trange, qse[Network_qse.find_el("H1", a), 1,:,1,ind], title = L"X_{\mathrm{NSE}}\in[%$(round(minimum(cl_nse),digits=4)),%$(round(maximum(cl_nse),digits=4))]", titleloc = :middle, xlabel = "T₉", ylabel = "Xᵢ", linestyle = :dash, color = "red", label = L"\mathrm{p\,\,QSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
    annotate!(trange[40],1, Plots.text(L"X_{\mathrm{QSE}}=%$(round(cl_qse[ind],digits=6))"))
    #annotate!(trange[40], 0.9, Plots.)
    plot!(trange, nse[Network_qse.find_el("H1", a), 1,:,1], titleloc=:right, xlabel = "T₉", ylabel = "Xᵢ", color = "red", label = L"\mathrm{p\,\,NSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
    #Plots.savefig("H1.png")
    plot!(trange, qse[Network_qse.find_el("1n", a), 1,:,1,ind], xlabel = "T₉", ylabel = "Xᵢ", linestyle = :dash, color = "blue",label = L"\mathrm{n\,\,QSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
    Plots.plot!(trange, nse[Network_qse.find_el("1n", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", color = "blue", label = L"\mathrm{n\,\,NSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
    Plots.savefig("../confirmation/png/alphas$ind.pdf")
end


minimum(cl_nse)


plot!(trange, qse[Network_qse.find_el("1n", a), 1,:,1,end], xlabel = "T₉", ylabel = "Xᵢ", linestyle = :dash, color = "blue",label = L"\mathrm{n\,\,QSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
Plots.plot!(trange, nse[Network_qse.find_el("1n", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", color = "blue", label = L"\mathrm{n\,\,NSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))

plot!(trange, qse[Network_qse.find_el("He4", a), 1,:,1,end], xlabel = "T₉", ylabel = "Xᵢ", linestyle = :dash, color = "green",label = L"\mathrm{He-4\,\,QSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
Plots.plot!(trange, nse[Network_qse.find_el("He4", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", color = "green", label = L"\mathrm{He-4\,\,NSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))

Plots.plot!(trange, qse[Network_qse.find_el("C12", a), 1,:,1,end], xlabel = "T₉", ylabel = "Xᵢ", linestyle = :dash, color = "green",label = L"\mathrm{He-4\,\,QSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
Plots.plot!(trange, nse[Network_qse.find_el("C12", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", color = "green", label = L"\mathrm{He-4\,\,NSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))



plot!(trange, qse[Network_qse.find_el("Si28", a), 1,:,1,end], xlabel = "T₉", ylabel = "Xᵢ", linestyle = :dash, color = "green", label = L"{}^{28}\mathrm{Si}\,\,\mathrm{QSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
plot!(trange, nse[Network_qse.find_el("Si28", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", color = "green", label = L"{}^{28}\mathrm{Si}\,\,\mathrm{NSE}", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
Plots.savefig("../Network_qse_run1/NSE_higher_QSE/logScale12/compare_end.png")




function species_str(a_i)
    if a_i.name == "1n"
        x = "n+p"
        y = " "
    elseif a_i.name == "H1"
        x = "n+p"
        y = " "
    elseif a_i.name == "H2"
        x = "H"
        y = "2"
    elseif a_i.name == "H3"
        x = "H"
        y = "3"
    elseif a_i.name == "He3"
        x = "He"
        y = "3"
    elseif a_i.name == "He4"
        x = "He"
        y = "4"
    elseif length(a_i.name) == 4.0
        y = a_i.name[3:4]
        x = a_i.name[1:2]
    else
        y = a_i.name[2:3]
        x = a_i.name[1]
    end
    return x, y
end


mutable struct NSE end

@recipe function f(::NSE, n::Integer, el::Integer, range::Vector, xaxis::String,
                    yVal::Array, customcolor = :green; add_marker = true)
    linecolor   --> customcolor
    seriestype  :=  :path
    markershape --> (add_marker ? :star5 : :none)
    #ytickfont -> font(20, "Courier")
    #markercolor --> "green"
    label --> false
    #xguide --> "Yₑ"
    xguide --> xaxis
    yguide --> "Xᵢ"
    #xtickfont --> font(50, "Courier")
    #ytickfont --> font(50, "Courier")
    yaxis --> :log
    linewidth = 4
    ylims --> (1e-4, 1)
    size --> (1920,1080)
    #linewidth = 4
    delete!(plotattributes, :true)
    #trange, all[el,n,:,1] # particle - Y - : - rho
    #yrange, all[el,:,n,1] # particle - : - tem - rho
    range, yVal
end



#col = shuffle(shuffle(Base.range(colorant"green", stop=colorant"red", length=size(a,1))))
#col = cgrad(:roma, rev = true, size(a,1), categorical = false, scale = :exp)
#col = cgrad(:tab20c, size(a,1), categorical = true, scale = :linear)
##col = shuffle(shuffle(Base.range(HSV(0,1,1), stop=HSV(-360,1,1),length=size(a,1))))
#col = shuffle(shuffle(Base.range(HSV(0,1,1), stop=HSV(-360,1,1),length=size(a,1))))

col = distinguishable_colors(size(a,1), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
shuffle!(col)



animY = @animate for i ∈ reverse(1:size(cl_qse,1)-1) # timeframes T
    println(i)
    Plots.plot()

    for (k, el) in zip(Iterators.countfrom(3), a[3:end])#for (k, el) in enumerate(a[3:end])
        #println(k, el)
        ind = argmax(filter!(!isnan, nse[k,1,:,1]))
        if nse[k,1,:,1][ind] > 0.01
            plot!(trange, nse[k, 1, :, 1], label = nothing, color = "grey", yaxis = :log)
            x,y = species_str(a[k])

            annotate!(trange[ind], 1.2*nse[k,1,:,1][ind], Plots.text(L"{}^{%$y}\!\textrm{%$x}", 25, "grey"))
        end
    end
    ind = argmax(filter!(!isnan, qse[1,1,:,1,i]))
    ind = argmax(filter!(!isnan, qse[1,1,:,1,i]))
    Plots.plot!(trange, qse[1,1,:,1, i]+qse[2,1,:,1, i], c = "red", label = "QSE p + n",ls=:dashdot,legendfontsize=12)
    Plots.plot!(trange, nse[1,1,:,1]+nse[2,1,:,1], c = "grey", label = "NSE p + n", ls=:dashdot,legendfontsize=12)
    ##annotate!(trange[50], 1.2, Plots.text(L"\log\,(1 - \mathrm{X}_{\mathrm{Cl}}) \;\; = ", 19, "black"))
    #annotate!(trange[50], 0.9, Plots.text(L"\overline{\mathrm{X}}_{\mathrm{NSE}} \;\;\;\;\; = ", 19, "black"))
    #annotate!(trange[50], 0.7, Plots.text(L"\mathrm{X}_\mathrm{Cl}/\overline{\mathrm{X}}_\mathrm{NSE} = ", 19, "black"))
    annotate!(trange[37], 1.2, Plots.text("X = "*lpad(rpad(string(round(cl_qse[i]; sigdigits = 7, base = 10)),4), 1), 19, "black"))
    #annotate!(trange[40], 0.9, Plots.text(lpad(rpad(string(round(mean(cl_nse[i]); sigdigits = 2, base = 10)),4), 4), 19, "black"))
    #annotate!(trange[40], 0.7, Plots.text(lpad(rpad(string(round(cl_qse[i]/mean(cl_nse[i]); sigdigits = 2, base = 10)),4), 4), 19, "black"))

    annotate!(trange[90], 1.2, Plots.text(string(round(rrange[1]; sigdigits = 2, base = 10))*" g/cm³", 19, "black"))
    #annotate!(trange[end - 5], 1.2*(qse[1,1,:,1, i]+qse[2,1,:,1, i]), Plots.text(L"n + p", 15, "red"))
    #annotate!(trange[end - 5], 1.2*(nse[1,1,:,1]+nse[2,1,:,1]),Plots.text(L"n + p", 15, "grey"))
    for (k, el) in zip(Iterators.countfrom(3), a[3:end]) #enumerate(a[3:end])
        ind = argmax(filter!(!isnan, qse[k,1,:,1,i]))
        if qse[k,1,:,1,i][ind] > 0.01
            #plot!(nse,i, k, yrange, "Yₑ", all[k,:,i,1], colT[rand(1:30)], marker = (:circle,2),
            Plots.plot!(trange, qse[k,1,:,1, i],
            yaxis = :log,
            lw = 4,
            label = false,
            ylims = (1e-3, 1.5),
            yticks = ([1e-3, 1e-2, 1e-1, 1], ["10⁻³", "10⁻²", "10⁻¹", "1"]),
            size = (1920,1380),
            xlabel = "T₉ [K]",
            ylabel = "Xᵢ",
            c = col[k],
            #marker = (:circle,2),
            #markercolor = "white",
            linewidth = 2,
            guidefont=font(23),
            xtickfont = font(16),
            ytickfont = font(16),
            thickness_scaling = 1.3,
            margin=5Plots.mm)
            x,y = species_str(a[k])
            if qse[k,1,:,1, i][ind] > 0.001 && qse[k,1,:,1, i][ind] < 0.01
                annotate!(trange[ind], 1.2*qse[k,1,:,1, i][ind],
                Plots.text(L"{}^{%$y}\!\textrm{%$x}", 15, col[k]))
            else
                annotate!(trange[ind], 1.2*qse[k,1,:,1, i][ind],
                Plots.text(L"{}^{%$y}\!\textrm{%$x}", 20, col[k]))
            end
        end
    end
end

gif(animY, "../Network_qse_run1/NSE_higher_QSE/logScale12/x_cl_evol_2.mp4", fps = 7)
