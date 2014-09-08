library(jsonlite)
library(igraph)
library(stringr)
library(plyr)
library(stringr)
library(RColorBrewer)

make_unipartite <- function(x)
{
   nsp <- NROW(x) + NCOL(x)
   A <- matrix(0, ncol=nsp, nrow=nsp)
   colnames(A) <- rownames(A) <- c(1:nsp)
   A[1:NROW(x),NROW(x)+c(1:NCOL(x))] = x
   return(A)
}

output_path <- "../outputs/pilot/"
sim_output <- str_c(output_path, list.files(path=output_path))

SIM <- fromJSON(sim_output[42])
G <- graph.adjacency(SIM$inputs$G)
A <- graph.adjacency(make_unipartite(SIM$inputs$A))
endpoint <- SIM$timepoints[length(SIM$timepoints)]

lgraphs <- list()
for(i in c(1:vcount(G)))
{
   reds <- endpoint[[1]]$red[[i]]
   blues <- endpoint[[1]]$blue[[i]]
   if (length(reds) > 0 & length(blues) > 0)
      lgraphs[[as.character(i)]] <- graph.adjacency(make_unipartite(SIM$inputs$A[reds, blues]))
}

### bipartite motifs
s36 <- list("1"=c(2,3))
s6 <- list("1"=3, "2"=3)
s76 <- list("1"=c(3,4),"2"=4)
s2184 <- list("1"=c(2,3,4))
s14 <- list("1"=4, "2"=4, "3"=4)
s204 <- list("1"=c(3,4),"2"=c(3,4))
motifs <- llply(list("36"=s36, "6"=s6, "76"=s76, "2184"=s2184, "14"=s14, "204"=s204), graph.adjlist)

count_motifs = function(G) laply(motifs, function(x) graph.count.subisomorphisms.vf2(G, x))
n_motifs <- ldply(lgraphs, count_motifs)
colnames(n_motifs) <- c("V", str_c("M",names(motifs)))

### Correct the distribution of motifs
for(i in c(1:nrow(n_motifs)))
{
   tA <- lgraphs[[as.character(n_motifs$V[i])]]
   n_blue <- sum(degree(tA, mode="out")>0)
   n_red <- vcount(tA) - n_blue
   Correct <- c(choose(n_blue,1)*choose(n_red,2),
                choose(n_blue,2)*choose(n_red,1),
                choose(n_blue,2)*choose(n_red,2),
                choose(n_blue,1)*choose(n_red,3),
                choose(n_blue,3)*choose(n_red,1),
                choose(n_blue,2)*choose(n_red,2)
                )
   c_count = unlist(n_motifs[i,c(2:7)] / Correct)
   c_count[is.nan(c_count)] = 0.0
   c_count[c_count>1] = 1
   n_motifs[i,c(2:7)] = c_count
}

#### PLOT

lay = layout.fruchterman.reingold(G)

par(mfcol=c(3,2))
for(m in names(motifs))
{
   mcount = rep(0, vcount(G))
   mcount[as.numeric(n_motifs$V)] = n_motifs[[str_c("M",m)]]
   plot(G, layout=lay, vertex.size=mcount*15, vertex.label=NA, vertex.color="palegreen", edge.arrow.size=0, vertex.frame.color=NA, edge.color="#ccccccbb")
   title(m)
}

