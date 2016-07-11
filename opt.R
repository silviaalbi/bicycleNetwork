####### Combine solutions of the gennetic

combine <-function (matinitsols,elite,tampob,tamelite)
{
  listasolsmat<-list()
  
  for (i in 1:tamelite) # we do not touch the elite
  {
    listasolsmat[[i]]<-matinitsols[[elite[i]]]
  }
  
  for (i in 1:(tampob-tamelite)) # The rest are combinations
  {
    sel1<-matinitsols[[sample(1:tamelite,1)]]
    sel2<-matinitsols[[sample(1:tamelite,1)]]
    
    newcomb<-sel1
    
    for (j in 1:nrow(sel1))
    {
      if (runif(1)>0.75) newcomb[j,]<-sel2[j,] #we randomly mix the rows
      
      if (runif(1)>0.8) 
        {
        seli<-sample(1:nrow(sel1),1)
        newcomb[j,seli]<- !newcomb[j,seli] # mutate also some values
      }
      
    }
    
    listasolsmat[[i+tamelite]]<-newcomb
  }
  combine<-listasolsmat
}


############# Genetic

genetic<- function(STEPS=1000,dreal,LIMRED,mpesos)

{
  tampob<-100
  tamelite<-10
  initsols<-list()
  matinitsols<-list()
  valcost<-rep(0,tampob)
  
  for (i in 1:tampob)
  {
 #   print(i)
    k<-Inf
   while (k==Inf)
   {
   initsols[[i]]<-montecarlo(100,dreal,LIMRED,mpesos)
    matinitsols[[i]]<-initsols[[i]]$msol
    valcost[i]<- initsols[[i]]$costsol
    k<-valcost[i]
   }
  }
  elite<-order(valcost)[1:tamelite] ### elements if the elite
  
  
  ###MAIN LOOP
  for (i in 1:STEPS)
  {
  #  print(paste("step",i,"best",valcost[elite[1]]))
  
    newborns<-combine(matinitsols,elite,tampob,tamelite)
    costinitsols<-list()
    for (i in 1:tampob)
    {
      costinitsols[[i]]<-costsol(newborns[[i]],LIMRED,dreal,mpesos)    
      valcost[i]<- costinitsols[[i]]$cost
      matinitsols[[i]]<-newborns[[i]]
    #  print(paste(i,valcost[i]))
    }
    elite<-order(valcost)[1:tamelite] ### elements if the elite
    
  }
  sel<-elite[1]
  genetic<-list("msol"=matinitsols[[sel]],"costsol"=costinitsols[[sel]]$cost,"igb"=costinitsols[[sel]]$gr,"distred"=costinitsols[[sel]]$distred)
}
  
  
  
#### annealing

annealing<- function(dreal,LIMRED,mpesos,pct=0.015,Temp=10000,cool=0.95)
{
  
  psubst<-0.1
  k<-0
  #### find a feasible solution using a montecarlo method
  
  initsol<-montecarlo(1000,dreal,LIMRED,mpesos)
  #print("initsol")
  #print(initsol$costsol)
  
  
  costini<-costsol(initsol$msol,LIMRED,dreal,mpesos)
  
  numcols<-ncol(initsol$msol)
  
  while (Temp>0.1)
  {
    #choose direction
    ##dir<-sample(1:2,1)
    dir <-1
    if(runif(1)<psubst) dir<-0 ### most of cases we add 
    if (dir==1) #add  a connection
    {
      newsol<-initsol$msol
      i <-sample(1:numcols,1)## row randomly selected
      if (sum(newsol[,i]==FALSE)>0) # solo si hay elementos a False
      {
        sel<-which(newsol[,i]==FALSE)[sample(1:sum(newsol[,i]==FALSE),1)] ### We select the position to changue
        newsol[sel,i]<-TRUE # we convert randomly one FALSE in TRUE
        
      }
    }
    else #remove a connection
    {
      newsol<-initsol$msol
      i <-sample(1:numcols,1)## row randomly selected
      if (sum(newsol[,i]==TRUE)>0) # solo si hay elementos a False
      {
        sel<-which(newsol[,i]==TRUE)[sample(1:sum(newsol[,i]==TRUE),1)] ### We select the position to changue
        newsol[sel,i]<-FALSE # we convert randomly one FALSE in TRUE
      }
      
    }
  #  print(sum(newsol==TRUE))
  #  print(sum(initsol$msol==TRUE))
    
    costnew<-costsol(newsol,LIMRED,dreal,mpesos)   
    ea<-costini$cost
    eb<-costnew$cost
    
   
    
    k<-k+1
    p<-exp((-(eb-ea)/Temp))## prob of not selecting the best sol
   # p<-0.05
    #print(paste("bestsol",k,Temp,p))
    #print(ea)
    
    if ((eb<ea) ||  runif(1)<p)
    {  
      initsol<-list("msol"=newsol,"costsol"=costnew$cost,"igb"=costnew$gr,"distred"=costnew$distred)
      costini<-costnew
    }
    # Decrease the temperature
    Temp=Temp*cool
    
  }
  annealing<-list("msol"=newsol,"costsol"=costnew$cost,"igb"=costnew$gr,"distred"=costnew$distred)
}
  
  ### Get a list of neighbours from a solution
  
  getlistneigh<-function(msol)
  {
    ### We simply add a TRUE for row   
    lsol<-list()
    iters<-ncol(msol)
    for (i in 1:iters)
    {
      newsol<-msol
      if(sum(newsol[,i]==FALSE)>0) msol[which(newsol[,i]==FALSE)[sample(1:sum(newsol[,i]==FALSE),1)],i]<-TRUE # we convert randomly one FALSE in TRUE
      #newsol<-lower.tri(newsol)*newsol
      lsol[[i]]   <-newsol
    }
    getlistneigh<-lsol
  }
  
  #### computes the best solution from a list of solutions
  
  getbestsol<-function(listneigh,dreal,LIMRED,mpesos)
  {
    mincost<-Inf
    iters<-length(listneigh)
    for (i in 1:iters)
    {
      dob<-costsol(listneigh[[i]],LIMRED,dreal,mpesos)
      if (mincost>dob$cost) 
      {
        mincost<-dob$cost
        igb<-dob$gr
        distredmin<-dob$distred
        msol<-listneigh[[i]]
        
      }
    }
    if (mincost==Inf) getbestsol<-list("msol"=NULL,"costsol"=Inf,"igb"=NULL,"distred"=NULL)
    else getbestsol<-list("msol"=msol,"costsol"=mincost,"igb"=igb,"distred"=distredmin)
  }
  
  
  
  
  ### hillclimb method
  
  hillclimb<- function(dreal,LIMRED,mpesos,pct=0.015)
  {
    
    minsol<-Inf
    #### find a feasible solution using a montecarlo method
    
    initsol<-montecarlo(1000,dreal,LIMRED,mpesos)
   # print("initsol")
  #  print(initsol$costsol)
    
    k<-1
    current<-initsol
    while(1)
    {
      listneigh<-getlistneigh(current$msol)  
      bestsol<-getbestsol(listneigh,dreal,LIMRED,mpesos)
      if(bestsol$costsol>=current$costsol) break
      current<-bestsol
     # print(paste("bestsol",k))
    #  print(current$costsol)
      k<-k+1
    }
    
    
    hillclimb<-current
  }
  
  
  
  ### Montecarlo method
  
  montecarlo<- function(STEPS,dreal,LIMRED,mpesos,pct=0.015)
    
  {
    costminsol<-Inf
    
    for (i in 1:STEPS)
    {
      sol<-matrix(runif(N*N),N,N) 
      sol<-sol<pct  #Solucion al azar , True indica persencia arista entre nodos , a menor valor menos aristas
      sol<-lower.tri(sol)*sol
      
      dob<-costsol(sol,LIMRED,dreal,mpesos)
      #  print (dob$cost)
      if (costminsol>dob$cost) 
      {
        costminsol<-dob$cost
        igb<-dob$gr
        distredmin<-dob$distred
        msol<-sol
        
      }
    }
    if(costminsol==Inf) montecarlo<-list("msol"=sol,"costsol"=costminsol,"igb"=NULL,"distred"=Inf)
    else montecarlo<-list("msol"=msol,"costsol"=costminsol,"igb"=igb,"distred"=distredmin)
  }
  
  
  #### computes the cost of a solution
  costsol <- function(sol,LIMRED,dreal,mpesos)
  {  
    distred<-sum(dreal*sol)## distancia de la red, debemos comprobar que sea menor que LIMRED para que sea válida
    
    if (distred>LIMRED)
    {
      #print("solucion larga")
      costsol<-list("cost"= Inf,"gr"=NULL)
    }
    else
    {
      ig <- graph.adjacency(dreal*sol, mode="undirected", weighted=TRUE)
      
      # plot(ig, edge.label=round(E(ig)$weight, 3)) #para mostrar la red
      
      cl<-clusters(ig)
      
      if (cl$no>1) 
      {
        #print("solucion no valida")
        #numero de clusters. Debe ser 1, si es mayor indica que las estaciones no están comunicadas
        costsol<-list("cost"= Inf,"gr"=NULL)
      }
      else
      {
        distancia_sol <-shortest.paths(ig, mode="out")##calcula distancias optimas entre todas las estaciones de la red
        
        objetivo<-distancia_sol*mpesos## objetivo a minimizar
        print(objetivo)
        dob<-sum(objetivo)
        costsol<-list("cost"= dob,"gr"=ig,"distred"=distred)
      }
    }
  }