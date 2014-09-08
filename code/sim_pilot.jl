using StatsBase
using JSON

include("utils/graph_generators.jl")
include("utils/web_generators.jl")
include("utils/simulation_functions.jl")

OUT_PATH = "../outputs/pilot/"

voids = true
while voids
   global A = type_one(25, 30, 0.2)
   nx = sum(A, 1) .> 0
   ny = sum(A. 2) .> 0
   voids = !(all(nx) & all(ny))
end

for treshold in linspace(0.1, 0.4, 3)
   for replicate in 1:10
      filename = OUT_PATH*"geom_"*string(treshold)*"_rep"*string(replicate)*".json"
      G = geometric(80, treshold)
      simulate(G,A, output_file=filename)
      filename = OUT_PATH*"er_"*string(treshold)*"_rep"*string(replicate)*".json"
      G = ER(80, treshold)
      simulate(G,A, output_file=filename)
   end
end
