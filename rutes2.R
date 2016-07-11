
options( java.parameters = "-Xmx6g" )
.lib<- c("igraph","geosphere", "dplyr")



.inst <- .lib %in% installed.packages()
if (length(.lib[!.inst])>0) install.packages(.lib[!.inst])
lapply(.lib, require, character.only=TRUE)
#setwd("C:/Users/Lidia/Dropbox/HackForGood2016/bicis")
setwd("/home/cesar/Dropbox/research/ECML2016")

source("opt.R")
MSEL<-2 # 1 Montecarlo, 2 Hill, 3 Annealing,4 Genetic
NITER<-1
NETFACT<-70 #Factor tha divide the total length of the network and it div

estaciones<-read.csv("datos/estaciones_valenbisi.csv")
index <- with(estaciones, order(estaciones$NUM_STATION))
estaciones<-estaciones[index,]

N<-nrow(estaciones) #estaciones


#### CODIGO NANDO

#weights by use
stationUse <- tbl_df(read.csv("datos/summaryBikes.csv"))
stationUse.station <- group_by(stationUse,station)
useByStation <- summarise(stationUse.station, use=mean(mean.day)) 
max<-max(useByStation$use)
min<-min(useByStation$use)
useByStation$weight <- 1-(useByStation$use - min)/( max - min )




for(i in 1:N){
  estaciones[i,"peso"]<-(estaciones[i,"STANDS"]*10)/40
  #estaciones[i,"peso"] <- filter(useByStation, station == i)[[1,"weight"]] 
  
}
max2 <- max(estaciones$peso)
min2 <- min(estaciones$peso)
estaciones$peso <- (estaciones$peso - min)/( max - min )

w1 = 0.5
w2 = 0.5
if (sum(useByStation$station == estaciones$NUM_STATION)==nrow(useByStation)){
  estaciones$pondWeight <- w1 * useByStation$weight + w2 * estaciones$peso
}else{
  estaciones$pondWeight <- useByStation$weight
}

pesos<-estaciones$pondWeight #pesos entre 1 y 10 (aqui deberiamos poner los pesos reales) #Nando, pesos (ponderados) normalizados por uso de las estaciones en 2015 + 
#pesos<-runif(N,0,10)
mpesos<-matrix(0,N,N) # indica la importancia entre par de estaciones
for (i in 1:N)#### seguro que se puede hacer sin bucles
  for (j in 1:N)
  {
    mpesos[i,j]<-pesos[i]+pesos[j]
  }
mpesos<-lower.tri(mpesos)*mpesos #triangulamos por abajo


dreal<-matrix(0,N,N) # indica la distancia entre par de estaciones
for (i in 1:N)#### seguro que se puede hacer sin bucles
  for (j in 1:N)
  {
    distancia<-distHaversine(c(as.numeric(estaciones[i,"LON"]),as.numeric(estaciones[i,"LAT"])),c(as.numeric(estaciones[j,"LON"]),as.numeric(estaciones[j,"LAT"])))
    dreal[i,j]<-distancia
  }
dreal<-lower.tri(dreal)*dreal #triangulamos por abajo


LIMRED<-sum(dreal)/NETFACT # tamaño de la red en distancia , ponemos este limite para evitar tener una red infinita
MAXD<-N*(N-1)/2 #tamaño de la matriz triangular

rescost<-c(1:NITER)
reskm<-c(1:NITER)
restmp<-c(1:NITER)
for (i in 1:NITER)
{
  
  set.seed(i)
  if (MSEL==1)
  {
    k=Inf
    STEPS<-1000
    tiempo<-system.time(while (k==Inf)
    {
      mtsol<-montecarlo(STEPS,dreal,LIMRED,mpesos)
      k<-mtsol$costsol
    #  print("kk")
    #  print(k)
    })
 #   mtsol<-montecarlo(STEPS,dreal,LIMRED,mpesos))
    print (paste("sol minima Montecarlo",i))
    print (mtsol$costsol)
    print("km")
    print(mtsol$distred)
    print(tiempo)
  }
  
  if (MSEL==2)
  {
    tiempo<-system.time(mtsol<-hillclimb(dreal,LIMRED,mpesos))
    costminsol<-mtsol$costsol
    print ("sol minima Hillclimb")
    print (mtsol$costsol)
    print("km")
    print(mtsol$distred)
    print(tiempo)
  }
  
  if (MSEL==3)
  {
    tiempo<-system.time(mtsol<-annealing(dreal,LIMRED,mpesos,Temp=100000,cool=0.999))
    print ("sol minima annealing")
    print (mtsol$costsol)
    print(tiempo)
  }
  
  if (MSEL==4)
  {
    tiempo<-system.time(mtsol<-genetic(300,dreal,LIMRED,mpesos))
    print ("sol minima genetic")
    print (mtsol$costsol)
    print(tiempo)
  }
  rescost[i]<-mtsol$costsol
  restmp[i]<-tiempo[1]
  reskm[i]<-mtsol$distred
}

totres<-list("cost"=rescost,"tmp"=restmp,"km"=reskm)
titfile<-paste("results","met",MSEL,"iter",NITER,"NETFACT",NETFACT,".RData",sep="_")
  
  save(totres,file=titfile)

#costminsol<-mtsol$costsol
#igb<-mtsol$igb
#distredmin<-mtsol$distred


#print("distancia red total")
#print(sum(dreal))
#print ("sol minima")
#print (costminsol)
#print("km")
#print(distredmin)
#plot(igb, edge.label=round(E(igb)$weight, 3)) #para mostrar la red


### Nando 2
