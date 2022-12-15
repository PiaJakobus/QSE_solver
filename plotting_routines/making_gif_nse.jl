using Colors
using Random 
using RecipesBase
using LaTeXStrings

a = Network_qse.extract_partition_function()
yrange = collect(LinRange(0.5,0.47, 64))
trange = collect(LinRange(6e9, 3e9, 64))
rrange = collect(LinRange(1e7, 1e7, 1))


#all = Network_qse.testing(yrange, trange, rrange)

#save("../Network_qse_run1/data.jld", "data", all)
all = load("../Network_qse_run1/nse_tables/150150.jld")["data"]




function species_str(a_i)
    if a_i.name == "1n"
        x = "n"
        y = "1"
    elseif a_i.name == "H1"
        x = "H"
        y = "1"
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

species_str(a[55])

length("hi")

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
    ylims --> (1e-2, 1)
    size --> (1920,1080)
    #linewidth = 4
    delete!(plotattributes, :true)
    #trange, all[el,n,:,1] # particle - Y - : - rho
    #yrange, all[el,:,n,1] # particle - : - tem - rho
    range, yVal
end




nse = NSE()
elem = [1,2,3, 6,762,764,772,774,660,661,664,710,563,565,470,997,884]

col = shuffle(shuffle(Base.range(colorant"green", stop=colorant"red", length=size(a,1))))
col = cgrad(:roma, rev = true, size(a,1), categorical = false, scale = :exp)
col = cgrad(:tab20c, size(a,1), categorical = true, scale = :linear)
col = shuffle(shuffle(Base.range(HSV(0,1,1), stop=HSV(-360,1,1),length=size(a,1))))
#col = shuffle(shuffle(Base.range(HSV(0,1,1), stop=HSV(-360,1,1),length=size(a,1))))

col = distinguishable_colors(size(a,1), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

shuffle!(col)

animY = @animate for i ∈ reverse(2:size(trange,1)) # timeframes T
    Plots.plot()
    annotate!(yrange[30], 1.2, Plots.text(lpad(string(round(trange[i]; sigdigits = 2, base = 10)),5)*" K", 22, "black"))
    annotate!(yrange[10], 1.2, Plots.text(string(round(rrange[1]; sigdigits = 2, base = 10))*" g/cm³", 22, "black"))
    for (k, el) in enumerate(a)
        ind = argmax(nse_table[k,:,i,1])
        if nse_table[k,:,i,1][ind] > 0.001
            #plot!(nse_table,i, k, yrange, "Yₑ", nse_table[k,:,i,1], colT[rand(1:30)], marker = (:circle,2),
            Plots.plot!(yrange, nse_table[k,:,i,1],
            yaxis = :log,
            lw = 4,
            label = false,
            ylims = (1e-3, 1.5),
            yticks = ([1e-3, 1e-2, 1e-1, 1], ["10⁻³", "10⁻²", "10⁻¹", "1"]),
            size = (1920,1380),
            xlabel = "Yₑ",
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
            if nse_table[k,:,i,1][ind] > 0.001 && nse_table[k,:,i,1][ind] < 0.01
                annotate!(yrange[ind], 1.2*nse_table[k,:,i,1][ind],
                Plots.text(L"{}^{%$y}\!\textrm{%$x}", 15, col[k]))
            else
                annotate!(yrange[ind], 1.2*nse_table[k,:,i,1][ind],
                Plots.text(L"{}^{%$y}\!\textrm{%$x}", 20, col[k]))
            end
        end
    end
end

gif(animY, "output/y_evol_2.gif", fps = 12)




animT = @animate for i ∈ reverse(1:size(yrange,1)) # timeframes Y
    Plots.plot()
    annotate!(trange[30], 1.2,Plots.text(rpad(string(round(yrange[i]; sigdigits = 3, base = 10)), 5, '0')*" Yₑ", 22, "black"))
    annotate!(trange[10], 1.2, Plots.text(string(round(rrange[1]; sigdigits = 2, base = 10))*" g/cm³", 22, "black"))
    for (k,el) in enumerate(a)
        ind = argmax(nse_table[k,i,:,1])
        if nse_table[k,i,:,1][ind] > 0.001
            #plot!(nse,i, k, yrange, "Yₑ", nse_table[k,:,i,1], colT[rand(1:30)], marker = (:circle,2),
            plot!(trange, nse_table[k,i,:,1],
            yaxis = :log,
            alpha = 1,
            fontfamily = :Courier,
            lw = 4,
            label = false,
            ylims = (1e-6, 1.5),
            yticks = ([1e-3, 1e-2, 1e-1, 1], ["10⁻³", "10⁻²", "10⁻¹", "1"]),
            size = (1920,1380),
            xlabel = "T [K]",
            ylabel = "Xᵢ",
            c = col[k],
            #c = colT[rand(1:30)],
            #marker = (:circle,2),
            #markercolor = "white",
            linewidth = 2,
            guidefont=font(23),
            xtickfont = font(16),
            ytickfont = font(16),
            thickness_scaling = 1.3,
            margin=5Plots.mm)
            x,y = species_str(a[k])
            if nse_table[k,i,:,1][ind] > 0.001 && nse_table[k,i,:,1][ind] < 0.01
                annotate!(trange[ind], 1.2*nse_table[k,i,:,1][ind],
                Plots.text(L"{}^{%$y}\!\textrm{%$x}", 15, col[k]))
            else
                annotate!(trange[ind], 1.2*nse_table[k,i,:,1][ind],
                Plots.text(L"{}^{%$y}\!\textrm{%$x}", 20, col[k]))
            end
        end
    end
end



gif(animT, "output/temperature_evol2.gif", fps = 9)
