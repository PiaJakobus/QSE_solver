using Distributed


@everywhere include("./src/Network_qse.jl")


function run_parallel_nse(yrange, trange, rrange)
    all_results = pmap(1:size(yrange, 1)) do index
        result = Network_qse.testing(index, yrange, trange, rrange)
    end
    return all_results
end

function run_parallel_qse(yrange, trange, rrange, srange)
    all_results = pmap(1:size(yrange, 1)) do index
        result = Network_qse.testing_QSE(index, yrange, trange, rrange, srange)
    end
    return all_results
end
