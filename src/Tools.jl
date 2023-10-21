
"""
    logsumexp(arr)

max + log(sum(exp{arr - max}))
instead of
log(sum(exp(arr)))
"""
function logsumexp(arr::Array{Float64,2})
    max = maximum(arr)
    dx = arr .- max
    sumexp = sum(exp.(dx))
    return max + log.(sumexp)
end

function find_el(el::String, ap::Any)
    filter(i -> (ap[i].name == el), 1:size(ap,1))[1]
end

function inv_3x3(m)
    det = m[1,1] * (m[2,2] * m[3,3] - m[2,3] * m[3,2]) -
          m[1,2] * (m[2,1] * m[3,3] - m[2,3] * m[3,1]) +
          m[1,3] * (m[2,1] * m[3,2] - m[2,2] * m[3,1])
    m⁻¹ = (1.0/det) .*
        [m[2,2]*m[3,3]-m[2,3]*m[3,2] m[1,3]*m[3,2]-m[1,2]*m[3,3] m[1,2]*m[2,3]-m[1,3]*m[2,2];
         m[2,3]*m[3,1]-m[2,1]*m[3,3] m[1,1]*m[3,3]-m[1,3]*m[3,1] m[1,3]*m[2,1]-m[1,1]*m[2,3];
         m[2,1]*m[3,2]-m[2,2]*m[3,1] m[1,2]*m[3,1]-m[1,1]*m[3,2] m[1,1]*m[2,2]-m[1,2]*m[2,1]]
    #any(isnan.(m⁻¹)) ? m⁻¹ = pinv(m) : nothing
    return m⁻¹
end


function inv_2x2(df::Array{Float64,2})
    detInv = 1.0 / (df[1,1] * df[2,2] - df[1,2] * df[2,1])
    inv = ([df[2,2] -df[1,2]; -df[2,1] df[1,1]]) .* detInv
end

"""
    MultiNewtonRaphson(x::Vector, func::Any, th::Any, ap::Any, props::Any)
props contains parameters for NR. For NSE a min/max value of [-50,50] and
alpha = 30 is recommended. For QSE min/max = [-10,10] and alpha = 50
seems to run stable
"""
function MultiNewtonRaphson(x::Array{Float64,1}, func::Func, th::ThermoProperties, ap::Array{AtomicProperties, 1}, props::StepParameter;linsys=false)
    zaehler = 0
    ϵ = 1.0
    while (abs(ϵ) > 1e-4) && (zaehler < 10000)
        f = func.f(x)
        α = props.alpha
        maxi = props.max
        mini = props.min
        step = 1.0
        # can be more stable at low T to solve linear system instead of inverting matrix ...
        if linsys == true
            step = solve(LinearProblem(func.df(x),func.f(x)))
        else
            inv = func.inv(x)
            step =  inv * f
        end
        x = x .- min.(1, zaehler/α) * max.(min.(step, maxi), mini)
        ϵ = LinearAlgebra.norm2(f)
        zaehler += 1
    end
    println("res: ",x)
    return x

"""
    restmass(xi::Array{Float64,1}, ap::Array{Network_qse.AtomicProperties, 1},rho)
computes restmass for fixed rho,ye,temp
if array of different thermoprops desired use sth like:
map(k->restmass(qse[:,1,k,1,1],a,1e8),collect(1:length(trange)))
"""
function restmass(xi::Array{Float64,1}, ap::Array{Network_qse.AtomicProperties, 1},rho)
    return sum((rho .* xi ./ m_by) .* ((map(k->k.M,ap) .* meverg ./ map(k->k.A,ap)) .+ map(k->k.Z,a) .* m_e .* c^2))
end

end
