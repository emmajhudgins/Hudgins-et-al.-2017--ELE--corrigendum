#### Dispersal model from Hudgins et al. 2017 - UPDATED JULY 4, 2018
#### "Predicting the spread of all invasive forest pests in the United States". Ecol. Lett.
#### Written by Emma Hudgins

rm(list=ls()) 
library(pdist)
setwd('~/Downloads/RData')# replace with path to the data folder

#Read in Data
data<-read.csv('countydatanorm_march.csv', stringsAsFactors = FALSE) #spatial predictor data, scaled to mean 0, standard deviation of 1 with scale() function
data2<-read.csv('datanorm.csv', stringsAsFactors = FALSE) #species-level data (including life history traits)
gen<-matrix(0,3372, 78)
for (sppp in 1:75)
{gen[,sppp]<-c(which(data[,paste(data2$COLUMN_NAM[sppp])]>0), rep(0,(3372-length(which(data[,paste(data2$COLUMN_NAM[sppp])]>0)))))} #determine where species are present
rr<-which(data2$YEAR>=10) #determine which species have been present for 10 years
n_spp<-length(rr)
data2<-data2[rr,] 
gen<-gen[,rr]
FIA<-read.csv('FIAcodes_notypos.csv', stringsAsFactors = FALSE) #host identities for each pest
FIA<-FIA[rr,] # only examine pests present for 10 years
FIA2<-read.csv('FIA_march.csv', stringsAsFactors = FALSE) #host presences
fia<-list()
FIA$FIA<-as.character(FIA$FIA)
fia<-strsplit(FIA$FIA, split=", ")
currpopden<-as.matrix(read.csv("currpopden_march.csv", stringsAsFactors = FALSE)) #population density in each decade, scaled to mean 0, sd 1
sources<-as.list(read.csv('Psources_notypos.csv')[,1]) #identities of centroids of host ranges, used as proxies for initial introduction site of each pest
body<-read.csv('body_pred.csv', stringsAsFactors = F)
body$size[is.na(body$Source)==FALSE]<-scale(body$size[is.na(body$Source)==FALSE], center=TRUE) #body size of each pest


#Include only species with known host distributions, get host distribution info
L<-rep(0,n_spp)
prez<-matrix(0,3372,n_spp)

for (spp in 1:n_spp) #Find out where hosts are present by linking FIA data with list of known hosts for each pest
{ 
  cols<-paste("FIA", fia[[spp]], sep="_")
  pres<-rep(0, 3372)
  for (q in 1:length(cols))
  {
    if (cols[q] %in% colnames(FIA2))
    {
      nn<-which(FIA2[,cols[q]]!=0)
      pres[nn]=pres[nn]+1
    }
  }
  prez[,spp]<-c(which(pres!=0),rep(0, length(which(pres==0))))
  L[spp]<-length(which(pres!=0))
}
good<-which(L!=0)
prez<-prez[,good]
L<-L[good]
data2<-data2[good,]
gen<-gen[,good]
host.density2<-read.csv("hostvol_notypos.csv", stringsAsFactors = FALSE) #host density at each site
fia<-fia[good]
n_spp=length(data2[,1]) #number of species in cleaned dataset
prez2<-matrix(0,3372,64)
for (sppp in 1:64) #limit species presences recorded to areas where hosts are known to be present
{prez2[,sppp]<-c(intersect(prez[which(prez[,sppp]!=0),sppp], gen[which(gen[,sppp]!=0),sppp]), rep(0,(3372-length(intersect(prez[which(prez[,sppp]!=0),sppp], gen[which(gen[,sppp]!=0),sppp])))))}

### Finished data cleaning and loading

# Create distance matrix between grid cells
Tr1<-function(x)
{
  sqrt((data$X_coord-data$X_coord[x])^2+(data$Y_coord-data$Y_coord[x])^2)
}
dists<-sapply(1:3372, Tr1)
T1<-exp(-dists/50000) #take exponential of distance matrix in terms of # of grid cells
YEARS<-data2$YEAR #time since discovery for each species

LLfit=function(par) #Fitting function that is to be optimized - includes dispersal simulations
  {
    glm_pars<-c(1:23)
    pars<-rep(0,23)
    pars[c(1,21,22,4,18,20,8)]<-par #input variables you want to fit (1 is intercept, 22 is threshold for presence/to send out propagules, 21 is growth rate)
    par<-pars
    

    # Pest life history trait predictors
    constpD=par[2]*data2$NumHosts
    constpD[which(data2$Guild=="Borers")]=constpD[which(data2$Guild=="Borers")]-par[6]
    constpD[which(data2$Guild=="Defoliators")]=constpD[which(data2$Guild=="Defoliators")]-par[6]
    constpD[which(data2$Guild=="Pathogens")]=constpD[which(data2$Guild=="Pathogens")]+par[6]
    constpD[which(data2$Guild=="Suckers")]=constpD[which(data2$Guild=="Suckers")]-par[6]
    constpD[which(data2$Continent=="Asia")]=constpD[which(data2$Continent=="Asia")]+par[7]
    constpD[which(data2$Continent=="Australia")]=constpD[which(data2$Continent=="Australia")]-par[7]
    constpD[which(data2$Continent=="Eurasia")]=constpD[which(data2$Continent=="Eurasia")]+par[7]
    constpD[which(data2$Continent=="Europe")]=constpD[which(data2$Continent=="Europe")]+par[7]
    constpD[which(data2$Continent=="North America")]=constpD[which(data2$Continent=="North America")]-par[7]
    constpD[which(data2$Continent=="UNK")]=constpD[which(data2$Continent=="UNK")]-par[7]
    constpD[which(body$strong==1)]=constpD[which(body$strong==1)]+par[11]
    constpD[which(body$strong==0)]=constpD[which(body$strong==0)]-par[11]
    constpD[which(body$size!=0)]=constpD[which(body$size!=0)]+par[14]*body$size[which(body$size!=0)]
    constpD[which(body$size==0)]=constpD[which(body$size==0)]+par[23]
    constpD=matrix(rep(constpD),3372,64, byrow=TRUE)
    

    
    #Predictors of dispersal out of a site
    constpD2<-matrix(rep(par[9]*data[,18]+par[10]*data[,16]+par[16]*data[,19]+par[17]*data[,20]+par[18]*data[,21]),3372,64)+par[19]*host.density2
    constpD<-as.numeric(constpD)+constpD2
    
    #Predictors of dispersal into a site
    constpD3<-matrix(rep(par[4]*data[,19]+par[12]*data[,20]+par[3]*data[,21]+par[13]*data[,18]+par[15]*data[,16]),3372,64)+par[5]*host.density2
    
    Pfull<<-matrix(0, 3372, n_spp) #initialize matrix of final pest locations

    all_spp<-which(1:64!=3) # remove ALB (see corrigendum)

    #Dispersal Simulations
    for (spp in all_spp)
    {
    
      YEAR=YEARS[spp]

      #Pest Parameters
      Psource=sources[[spp]]
      Discovery<-2009-YEAR
      rem<-Discovery%/%10
      Discovery<-rem*10 #convert years since discovery to a multiple of 10

      T2<-T1[prez[1:L[spp],spp],prez[1:L[spp],spp]] #limit distance matrix between sites to sites within host range
      vecP<-rep(0,L[spp]) #pest presences at current timestep
      for (rrr in 1:length(Psource))
      {vecP[which(prez[,spp]==Psource[rrr])]=1} #add maximum propagule pressure at source location
      r0<-par[22] #growth rate
      if (par[20]==0) #if human population not a predictor, this doesn't need to be calculated within the time loop
      {zzz<-matrix(rep(constpD3[prez[1:L[spp],spp],spp], L[spp]), nrow=L[spp], ncol=L[spp], byrow=TRUE)}  #move host density?

      for (time in 1:(YEAR/10))
      {
        
        Pnext<-rep(0,L[spp]) #next timestep's pest distribution
        column<-(((Discovery+10*(time-1))-1790)/10)+1 #figure out column of human population matrix
        
        qq<-matrix(rep(constpD[prez[which(vecP>=par[21]),spp],spp]+par[8]*currpopden[prez[which(vecP>=par[21]),spp],column], L[spp]), nrow=length(which(vecP>=par[21])), ncol=L[spp]) #update predictors of dispersal out to include current human population density at this timestep
        if (par[20]!=0) 
        {zzz<-matrix(rep(constpD3[prez[1:L[spp],spp],spp]+par[20]*currpopden[prez[1:L[spp],spp],column], L[spp]), nrow=L[spp], ncol=L[spp], byrow=TRUE)} #update predictors fo dispersal in to include current human population density at this timestep
        qq<-(2*par[1]*exp((zzz[which(vecP>=par[21]),]+qq)))/(1+exp(zzz[which(vecP>=par[21]),]+qq)) #scale probability of dispersal logistically, combine predictors of dispersal into and out of each site
        qq<-T2[which(vecP>=par[21]),]^qq #multiply by distance matrix
        
        if (length(which(vecP>=par[21]))>1){qq<-qq/rowSums(qq)} #scale to maximum propagule pressure that can leave each site
        if (length(which(vecP>=par[21]))==1){qq<-qq/sum(qq)}
        qq[which(qq<0.001)]=0 #round to 3 decimal places
        Pnext=(vecP[which(vecP>=par[21])])%*%(qq) #matrix multiplication to derive propagule pressure moving between all sites at this timestep
        Pnext[which(Pnext>=par[21])]=Pnext[which(Pnext>=par[21])]*r0 #growth of propagule pressure
        Pnext[which(Pnext>=1)]<-1 #limit to max propagule pressure of 1
        vecP=Pnext
        vecP[which(prez[,spp]==Psource)]=1 #maintain maximum propagule pressure at source location
      }
      Pfull[,spp]<<-c(prez[which(vecP>=par[21]),spp], rep(0, 3372-length(which(vecP>=par[21])))) #record locations of final pest presence
    }
    
    deviation=function(spp) #MET calculation (see supplementary material for theoretical background)
    {
        dii<-2*sum(dist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), upper=F))
        dij<-sum(pdist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), cbind(data$X_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]], data$Y_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]]))@dist)
        djj<-2*sum(dist(cbind(data$X_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]], data$Y_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]])))
        return(((1/(length(which(Pfull[,spp]!=0))*length(which(prez2[,spp]!=0))))*dij)-((1/(2*length(which(Pfull[,spp]!=0))^2))*dii)-((1/(2*length(which(prez2[,spp]!=0))^2))*djj))
    }
    return (sum(tapply(all_spp, all_spp, deviation))) #MET score returned is summed across all species and is in meters
}
#Fit best fitting paramter values for dispersal model (minimize MET)
m3<-optim(par=c(1.47513873697857,0.000125427515111294,1.31547456769312, -0.571318199516109, 14.1309612456393, -0.165979641909258, 0.239137626880927), fn=LLfit, control=list(trace=100, maxit=1000,parscale=c(0.1,0.0001,1,1,1,1,1))) #parscale tunes optimizer based on relative size of paramters, maxit sets max iterations, trace prints out fitting process
#Given that the optimizer sometimes gets stuck in local minima, I use this loop to restart the optimization until the MET value does not change
yy<-m3$value
xx<-1000000
while (yy!=xx)
  {
    yy<-m3$value
    m3<-optim(par=m3$par, fn=LLfit, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1,1,1,1,1)))
    xx<-m3$value
}
#save outputted locations of spread and parameters
write.table(Pfull, file=paste("presences", i,sep='.'))
write.table(m3$par, file=paste("par",i,sep='.'))
