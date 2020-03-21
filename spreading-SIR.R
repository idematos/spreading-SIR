rm(list=ls()) # clear all variables
set.seed(101) # start the random seed

library(igraph) # Load the igraph package

#### MODEL PARAMETERS ####
N = 200; # number of nodes
av.dg = 8; # average degree
m = av.dg/2; # parameter of the BA model
q = 0.5 # rewiring probability in the WS model
p = av.dg/N # probability in the ER model

## MODELS ###
# BA network
BA <- barabasi.game(N, m = av.dg/2, directed = FALSE)
# # ER network
ER <- erdos.renyi.game(N,p, type =c("gnp"))
# # WS network
WS <- sample_smallworld(dim=1,size=N, nei = av.dg/2, p = q)

######## READ FROM FILE ####
#net <- read.table("test-star.txt")
#G.Real <- graph.data.frame(net, directed=FALSE)
#plot(G)
# Zachary network
#G <- make_graph("Zachary")

######### Dynamics ############
# E = gsize(G) # number of edges
## Parameters
# beta = 0.2 # probability of transmission
mu = 0.05  # Probability of recovery


A = as_adjacency_matrix(BA)
x = eigen(A)
lambda.max = max(x$values)

lambdaBA = list()
rhoBA = list()
for(C in 0:10){
	beta = C*mu/lambda.max
	# print(lambda.max)
	# states: S:0 I:1 R:2

	tmp = 0
	for(i in 0:N){
		vstates = matrix(0,nrow = N, ncol = 1) # node states
		vstates[i] = 1 # i is infected
		Ninf = list() # number of infected nodes
		vt = list() # vt stores the time step
		t = 1 # time step
		Ninfected = 1 # number of infected nodes
		while(Ninfected > 0){ # while there are infected nodes
			vinfected = which(vstates %in% 1) #f ind the positions in vstates with value 1
			# try to infect all the nodes in step t
			for(j in vinfected){
				ng = neighbors(BA, j) # find the neighbours of j
				for(k in ng){ # for each neighbor
					if(vstates[k] != 2){ # if k is not recovered
						if(runif(1,0,1) < beta){ #try to infect k
							vstates[k] = 1
						}
					}
					# try to recover j
					if(runif(1,0,1) < mu){
						vstates[j] = 2
					}
				}
			}
			Ninfected = enlgth(which(vstates %in% 1)) # number of infected nodes
			Ninf[t] = Ninfected/N # number of infected nodes at time t
			vt[t] = t # time step
			# cat(t, Ninfected/N, "\n")
			t = t + 1
		}
		# plot(vt, Ninf, xlab="t", ylab="% Infected", col = 'black',pch = 21,  bg = "red", type="b")
		# fraction of recovery nodes when the disease starts in i
		rho.i = length(which(vstates %in% 2))/N
		tmp = tmp + rho.i
	}
	# cat("lambda:", beta/mu, "rho:", tmp/N, "\n")
	lambdaBA[C] = mu/beta
	rhoBA[C] = tmp/N
}

A = as_adjacency_matrix(ER)
x = eigen(A)
lambda.max = max(x$values)

lambdaER = list()
rhoER = list()
for(C in 0:10){
	beta = C*mu/lambda.max
	# print(lambda.max)
	# states: S:0 I:1 R:2

	tmp = 0
	for(i in 0:N){
		vstates = matrix(0,nrow = N, ncol = 1) # node states
		vstates[i] = 1 # i is infected
		Ninf = list() # number of infected nodes
		vt = list() # vt stores the time step
		t = 1 # time step
		Ninfected = 1 # number of infected nodes
		while(Ninfected > 0){ # while there are infected nodes
			vinfected = which(vstates %in% 1) # find the positions in vstates with value 1
			# try to infect all the nodes in step t
			for(j in vinfected){
				ng = neighbors(ER, j) # find the neighbours of j
				for(k in ng){ # for each neighbor
					if(vstates[k] != 2){# if k is not recovered
						if(runif(1,0,1) < beta){ # try to infect k
							vstates[k] = 1
						}
					}
					# try to recover j
					if(runif(1,0,1) < mu){
						vstates[j] = 2
					}
				}
			}
			Ninfected = length(which(vstates %in% 1)) # number of infected nodes
			Ninf[t] = Ninfected/N # number of infected nodes at time t
			vt[t] = t # time step
			# cat(t, Ninfected/N, "\n")
			t = t + 1
		}
		# plot(vt, Ninf, xlab="t", ylab="% Infected", col = 'black',pch = 21,  bg = "red", type="b")
		# fraction of recovery nodes when the disease starts in i
		rho.i = length(which(vstates %in% 2))/N
		tmp = tmp + rho.i
	}
	# cat("lambda:", beta/mu, "rho:", tmp/N, "\n")
	lambdaER[C] = mu/beta
	rhoER[C] = tmp/N
}

A = as_adjacency_matrix(WS)
x = eigen(A)
lambda.max = max(x$values)

lambdaWS = list()
rhoWS = list()
for(C in 0:10){
	beta = C*mu/lambda.max
	# print(lambda.max)
	# states: S:0 I:1 R:2

	tmp = 0
	for(i in 0:N){
		vstates = matrix(0,nrow = N, ncol = 1) # node states
		vstates[i] = 1 #i is infected
		Ninf = list() # number of infected nodes
		vt = list() # vt stores the time step
		t = 1 # time step
		Ninfected = 1 # number of infected nodes
		while(Ninfected > 0){ # while there are infected nodes
			vinfected = which(vstates %in% 1) # find the positions in vstates with value 1
			# try to infect all the nodes in step t
			for(j in vinfected){
				ng = neighbors(WS, j) # find the neighbours of j
				for(k in ng){ # for each neighbor
					if(vstates[k] != 2){ # if k is not recovered
						if(runif(1,0,1) < beta){ # try to infect k
							vstates[k] = 1
						}
					}
					# try to recover j
					if(runif(1,0,1) < mu){
						vstates[j] = 2
					}
				}
			}
			Ninfected = length(which(vstates %in% 1)) # number of infected nodes
			Ninf[t] = Ninfected/N # number of infected nodes at time t
			vt[t] = t # time step
			# cat(t, Ninfected/N, "\n")
			t = t + 1
		}
		# plot(vt, Ninf, xlab="t", ylab="% Infected", col = 'black',pch = 21,  bg = "red", type="b")
		# fraction of recovery nodes when the disease starts in i
		rho.i = length(which(vstates %in% 2))/N
		tmp = tmp + rho.i
	}
	# cat("lambda:", beta/mu, "rho:", tmp/N, "\n")
	lambdaWS[C] = mu/beta
	rhoWS[C] = tmp/N
}

plot(lambdaBA, rhoBA, xlab="", ylab="", type="l", col="blue", lwd=3)
lines(lambdaER, rhoER, col="red", lwd=3)
lines(lambdaWS, rhoWS, col="green", lwd=3)

legend("topleft",
      c("BA","ER", "WS"),
      fill=c("blue","red", "green")
)
