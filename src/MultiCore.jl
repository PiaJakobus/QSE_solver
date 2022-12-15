# start julia specifying how many cores to make available
##$ julia -p 4

# julia multicore start command in a shell script
##$ julia -p $(nproc)
# http://stackoverflow.com/questions/2314750/how-to-assign-the-output-of-a-shell-command-to-a-variable

# alternatively add local cores dynamically (e.g. after starting julia single-core)
##$ julia
addprocs_local(3)

# check the number of cores made available
nprocs()

# total number of cores (with hyperthreading)
#  http://stackoverflow.com/questions/27931026/obtain-the-number-of-cpu-cores-in-julia
CPU_CORES

# make functions from a source file available to all processes
require("Network_qse.jl")
# altermnatively: (e.g. in interactive work)
# @everywhere macro executes a statement on all processes
@everywhere id = myid()
@everywhere function exec()
  # do something
end
