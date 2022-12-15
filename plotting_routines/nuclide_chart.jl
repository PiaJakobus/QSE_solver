
using LaTeXStrings
a = Network_qse.extract_partition_function()


isops = map(i -> a[i].A, 1:size(a,1))
nuclei = map(i -> a[i].Z, 1:size(a,1))




p = Plots.scatter(isops, map(i -> a[i].Eb, 1:size(a,1)) ./ isops, yaxis=:flip, legend=:false,xlabel="A", ylabel="Binding energy per nucleon [MeV]",alpha=:0.3)
scatter!([isops[Network_qse.find_el("Fe56", a)]], [a[Network_qse.find_el("Fe56", a)].Eb / isops[Network_qse.find_el("Fe56", a)]],c=:red,alpha=:0.6,subplot=1)
scatter!([isops[Network_qse.find_el("Ni62", a)]], [a[Network_qse.find_el("Ni62", a)].Eb / isops[Network_qse.find_el("Ni62", a)]],c=:red,alpha=:0.6,subplot=1)
scatter!([isops[Network_qse.find_el("Fe58", a)]], [a[Network_qse.find_el("Fe58", a)].Eb / isops[Network_qse.find_el("Fe58", a)]],c=:red,alpha=:0.6,subplot=1)
scatter!([isops[Network_qse.find_el("Ni56", a)]], [a[Network_qse.find_el("Ni56", a)].Eb / isops[Network_qse.find_el("Ni56", a)]],c=:red,alpha=:0.6,subplot=1)
scatter!([isops[Network_qse.find_el("Si28", a)]], [a[Network_qse.find_el("Si28", a)].Eb / isops[Network_qse.find_el("Si28", a)]],c=:blue,alpha=:1.0,subplot=1)


l = Plots.scatter!(
    isops[1:1:end], (map(i -> a[i].Eb, 1:size(a,1)) ./ isops)[1:1:end], yaxis=:flip,
    inset = (1, bbox(0.05, 0.05, 0.5, 0.25, :bottom, :right)),
    #ticks = nothing,
    #yaxis=:flip,
    alpha = 0.7,
    ylims = (-8.9, -8.6), xlims = (50,70),
    #yaxis=:flip,
    subplot = 2,
    legend=:false,
    bg_inside = nothing
)


#Plots.quiver!([50],[-2],quiver=([50],[-1]))
Plots.plot!([isops[Network_qse.find_el("Fe56", a)] - 1,140],[a[Network_qse.find_el("Fe56", a)].Eb / isops[Network_qse.find_el("Fe56", a)],-0.8],c=:grey,alpha=:0.8)
Plots.plot!([isops[Network_qse.find_el("Ni62", a)] + 1,250],[a[Network_qse.find_el("Ni62", a)].Eb / isops[Network_qse.find_el("Ni62", a)],-2],c=:grey,alpja=:0.2)
#Plots.plot!([isops[Network_qse.find_el("Ni56", a)] + 10,250],[a[Network_qse.find_el("Ni56", a)].Eb / isops[Network_qse.find_el("Ni56", a)],-2],c=:grey,alpja=:0.2)
#Plots.plot!([isops[Network_qse.find_el("Fe58", a)] + 10,250],[a[Network_qse.find_el("Fe58", a)].Eb / isops[Network_qse.find_el("Fe58", a)],-2],c=:grey,alpja=:0.2)



scatter!([isops[Network_qse.find_el("Fe56", a)]], [a[Network_qse.find_el("Fe56", a)].Eb / isops[Network_qse.find_el("Fe56", a)]],c=:red,subplot=2)
scatter!([isops[Network_qse.find_el("Ni62", a)]], [a[Network_qse.find_el("Ni62", a)].Eb / isops[Network_qse.find_el("Ni62", a)]],c=:red,subplot=2)
scatter!([isops[Network_qse.find_el("Fe58", a)]], [a[Network_qse.find_el("Fe58", a)].Eb / isops[Network_qse.find_el("Fe58", a)]],c=:red,subplot=2)
scatter!([isops[Network_qse.find_el("Ni56", a)]], [a[Network_qse.find_el("Ni56", a)].Eb / isops[Network_qse.find_el("Ni56", a)]],c=:red,subplot=2)
scatter!([isops[Network_qse.find_el("Si28", a)]], [a[Network_qse.find_el("Si28", a)].Eb / isops[Network_qse.find_el("Si28", a)]],c=:red,subplot=2)


annotate!((isops[Network_qse.find_el("Fe56", a)], a[Network_qse.find_el("Fe56", a)].Eb / isops[Network_qse.find_el("Fe56", a)]-0.05, Plots.text("Fe56", 9, "red")),subplot=2)
annotate!((isops[Network_qse.find_el("Ni62", a)]+0.4, a[Network_qse.find_el("Ni62", a)].Eb / isops[Network_qse.find_el("Ni62", a)]-0.06, Plots.text("Ni62", 9, "red")),subplot=2)
annotate!((isops[Network_qse.find_el("Fe58", a)]+1, a[Network_qse.find_el("Fe58", a)].Eb / isops[Network_qse.find_el("Fe58", a)]-0.05, Plots.text("Fe58", 9, "red")),subplot=2)
#annotate!((isops[Network_qse.find_el("Ni56", a)]+0.5, a[Network_qse.find_el("Ni56", a)].Eb / isops[Network_qse.find_el("Ni56", a)]-0.03, Plots.text("Ni56", 9, "red")),subplot=2)
annotate!((52.5, -8.82, Plots.text("Ni56", 9, "red")),subplot=2)
Plots.plot!([52.5,56],[-8.79,a[Network_qse.find_el("Ni56", a)].Eb / isops[Network_qse.find_el("Ni56", a)]],c=:red,alpha=:0.5,subplot=2)

Plots.plot!([a[Network_qse.find_el("Si28", a)].A,65],[a[Network_qse.find_el("Si28", a)].Eb / isops[Network_qse.find_el("Si28", a)], -4],c=:darkblue,alpha=:1.0)
annotate!((70, -3.8, Plots.text("Si28", 9, "darkblue")))


Plots.savefig("plots/nuclides.pdf")
