function ER(N,p)
   G = rand((N, N)) .< p
   return G & !eye(Bool, N)
end

function geometric(N,d)
   x = rand(N)
   y = rand(N)
   dmat = sqrt((x .- x').^2 .+ (y .- y').^2)
   G = dmat .< d
   return G & !eye(Bool, N)
end
