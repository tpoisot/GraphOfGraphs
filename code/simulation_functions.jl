using StatsBase
using JSON

include("web_generators.jl")
include("graph_generators.jl")

function migration(sp, p)
   return sp[1:length(sp)][rand(length(sp)).<p]
end

function extinction(sp, p)
   return migration(sp, 1-p)
end

function keep(n, k)
   nn = typeof(n[1])[]
   for i in 1:length(n)
      if k[i]
         append!(nn, Int64[n[i]])
      end
   end
   return nn
end

function simulate(G,A; init_occupancy=0.2, migration_rate=0.3, extinction_rate=0.1, timesteps=100, output_file="sim.json")
   n_patches = size(G)[1]
   n_red = size(A)[2]
   n_blue = size(A)[1]


   output = {
      "inputs"=>{
         "G"=>G, "A"=>A,
         "c"=>migration_rate, "e"=>extinction_rate, "io"=>init_occupancy
      },"timepoints"=>Array(Array{Dict{ASCIIString, Array{Int64, 1}}, 1}, timesteps+1)
      }

   ### init data
   landscape = Array(Dict{ASCIIString, Array{Int64, 1}}, n_patches)
   for i=1:n_patches
      if rand() < init_occupancy
         landscape[i] = {"red"=>[1:n_red], "blue"=>[1:n_blue]}
      else
         landscape[i] = {"red"=>[], "blue"=>[]}
      end
   end

   output["timepoints"][1] = landscape

   ### simulation
   for current_time in 1:timesteps
      new_landscape = landscape
      ## MIGRATION
      for patch in 1:n_patches
         current_patch = landscape[patch]
         neighbours = [1:n_patches]'[G[patch,:]]
         for node_color in ["blue", "red"]
            if length(current_patch[node_color]) > 0
               migrants = migration(current_patch[node_color], migration_rate)
               for m in migrants
                  append!(new_landscape[sample(neighbours)][node_color], [m])
               end
            end
         end
      end
      ## MERGE NON-UNIQUE SPECIES AND DO EXTINCTIONS
      for patch in 1:n_patches
         new_landscape[patch]["red"] = extinction(unique(new_landscape[patch]["red"]), extinction_rate)
         new_landscape[patch]["blue"] = extinction(unique(new_landscape[patch]["blue"]), extinction_rate)
      end
      ## DO NETWORK-BASED EXTINCTIONS
      for patch in 1:n_patches
         global current_patch = new_landscape[patch]
         if length(current_patch["red"]) == 0
            new_landscape[patch]["blue"] = Int64[]
         elseif length(current_patch["blue"]) == 0
            new_landscape[patch]["red"] = new_landscape[patch]["red"]
         else
            global current_A  = A[current_patch["blue"],current_patch["red"]]
            p_blue = sum(current_A, 2) .> 0
            p_red = sum(current_A, 1) .> 0
            current_patch["red"] = keep(current_patch["red"], p_red)
            current_patch["blue"] = keep(current_patch["blue"], p_blue')
         end
      end
      ## UPDATE THE LANDSCAPE
      landscape = new_landscape
      output["timepoints"][current_time + 1] = landscape
   end
   ### WRITE
   JSON.print(open(output_file, "w"), output)
   return output
end
