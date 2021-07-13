rm(list = ls())


#### Set current working directory using setwd("")

library(deSolve)
library(DEoptim)
library(pracma)
library(Hmisc)
library(scales)

################################
# Model Function for within host viral dynamics
################################

SARS_COV_model <- function(t,x,params){
  with(as.list(x),{   
    
    ddt_S = -beta*V*S
    ddt_I = beta*V*S - delta*I^k*I - m*E^r*I/(E^r+E50^r) 
    ddt_V = p*I-c*V
    
    ddt_M1 = w*I*M1 - q*M1
    ddt_M2 = q*(M1-M2)
    ddt_E = q*M2- de*E
    
    der <- c(ddt_S,ddt_I,ddt_V,ddt_E,ddt_M1,ddt_M2)
    
    list(der)
  })       
}


################################   
# AUC calculation function
################################  
AUC <- function(t,y) {
  y1=na.omit(y)
  mn=length(y1)
  val_AUC=(trapz(t[seq(1,mn,by=1)],na.omit(y)))
  return(val_AUC)
}

################################
# Main to calculate number of superspreader events, time of superspreader events and other in the new variant emergence phenomenon 
################################

library(deSolve)

################################
# 1. Read Population parameters for heterogeneity in viral dynamics
################################

Parameters = read.csv("SimulatedParameters.txt",header=TRUE)

################# 
####### determine which Re your want to simulate data for
################

ReInput = 1.2

if(ReInput==1){
  no_contact_per_day_Re = 2.3
  ind=1
}else if(ReInput==1.2){
  no_contact_per_day_Re = 3.1
  ind=2
}else if(ReInput==1.4){
  no_contact_per_day_Re = 3.5
  ind=3
} else if(ReInput==1.6){
  no_contact_per_day_Re = 3.75
  ind=4
}else if(ReInput==1.8){
  no_contact_per_day_Re = 4.0
  ind=5
}else if(ReInput==2.0){
  no_contact_per_day_Re = 5.0
  ind=6
}else if(ReInput==2.2){
  no_contact_per_day_Re = 5.5
  ind=7
}

total_pop=10^6 # Assume a million people in the city

## simulate COVID-19 epidemics for 150 days 
no_days=150

### seed to the pandemic (could be 1, 10, 100 infections)
seed_number=1

infection_threshold=1000 ## stop simulation when total number of infections since the start of the pandemic reaches 1000 infections
SSE_threshold=5 ## Superspreader event threshold of 5
SSE_threshold2=10## Superspreader event threshold of 10
SSE_threshold3=20## Superspreader event threshold of 20
SSE_threshold4=40## Superspreader event threshold of 40
SSE_threshold5=60## Superspreader event threshold of 60
SSE_threshold6=80## Superspreader event threshold of 80


no_runs=1000 ## number of replicates as the simulation is stochastic
max_replicates=no_runs

#### Pre-define matrix for storage
II_metrics = matrix(0,nrow=no_days*no_runs,ncol=13)   # for 45 days time, incidence, recovered and instance
colnames(II_metrics)=c("time","Incidence","Recovered","Instance","SuperSpreader","Cumulative_infections",
                       "SSE10","SSE20","SSE40","SSE60","SSE80","VL","contacts")

instance_index=1 ###counter for each replicate

for (instances in seq(1,no_runs,by=1)){     ### running simulation for each replicate 
  
  ## tt_total_COVID is total COVID infections at all times
  tt_total_COVID = matrix(0,nrow=no_days,ncol=1)  # we run the epidemics for a maximum of no_days 
  tt_total_COVID = seq(1,length(tt_total_COVID),by=1)
  
  II_total_COVID = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  II_total_COVID[1]=seed_number ## seed the epidemic with one infection
  
  SSE_total_COVID = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID2 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID3 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID4 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID5 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID6 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 

  ### first_SSE_event_size is the first SSE event size
  first_SSE_event_size =  matrix(0,nrow=no_days,ncol=1)  # we run the epidemics for a maximum of no_days 
  ### II_recover_COVID is the naturally recovered (such that they have no to low viral load for transmission) 
  II_recover_COVID = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  
  Cumulative_infections = matrix(0,nrow=no_days,ncol=1)
  # print(II_total_COVID)
  
  tcurrent=0 ## starting time point
  
  t_options=seq(1,no_days,by=1) # we run the epidemics for a maximum of no_days 
  
  dispersion_options<-c(40) ## dispersion rho parameter value

  no_contact_per_day_final_ref <- rep(0,no_days) ## predefine array representing number of contact an individual will have on a daily basis and initiate it with 0
  V_ref <- rep(0,no_days) # predefine array representing viral loads will have on a daily basis and initiate it with 0
  
  for (tt1 in t_options) {   # simulate for 150 days
   
    Number_of_Active_infections <- c() ## number of active infections on that day
    Days_in_which_infection_will_recover<-c() ## how many days in which that particular individual will recover
    Day_at_which_infection_infection_will_recover<-c() ## On what day a particular individual will recover
    
    Number_of_Active_infections = II_total_COVID[tt1] ## how many active infections do we have at this time point 
    
    #### we have to simulate each active infection viral loads so we use a for loop
    
    if (Number_of_Active_infections>0){
      first_loop=seq(1,Number_of_Active_infections,by=1)
      for (yu in first_loop){ ## simulating epidemics originating because of each active infection at time t

        tmin=min(Parameters$tzero[c(seq(1,max_replicates,by=1))])
        ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
        tzero=0
        beta=Parameters$beta[ijxxx]
        delta=Parameters$delta[ijxxx]
        k=Parameters$k[ijxxx]
        p=Parameters$p[ijxxx]
        m=Parameters$m[ijxxx]
        w=Parameters$w[ijxxx]
        E50=Parameters$E50[ijxxx]
        r=Parameters$r[ijxxx]
        q=Parameters$q[ijxxx]
        de=Parameters$de[ijxxx]
        c=Parameters$c[ijxxx]
        
        # Initial conditions
        S_0 = 1e7
        I_0 = 1
        V_0 = p*I_0/c
        E_0 = 0
        M1_0 = 1
        M2_0 = 0
        
        dT=1 # 
        init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
        t.out <- seq(0,30,by=dT)
        params=c()
        out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
        ## out$V is the simulated viral load
        
        ###### this particular infections will also recover based on viral suppression duration
        Days_in_which_infection_will_recover= length(which(out$V>(10^2))) # out$V >100 means viral load detectable 
        Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1 ## +1 as the 0 cannot be input to array element
        II_recover_COVID[Day_at_which_infection_infection_will_recover]=II_recover_COVID[Day_at_which_infection_infection_will_recover]+1

        ### Now to generate new infections from that particular and to know when
        AUC_factor=7 # tzero	alpha	betat	no_contact	rho
        alpha=0.8
        no_contact_per_day=no_contact_per_day_Re ##
        dispersion_options= c(40) ## value for rho

        ## simulate number of exposed contact rate for that individual at all time points
        no_contact_per_day_final<-c()
        no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                        scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
        ### determine prob of transmission according to viral loads
        Prob_V=(pmax(0,(out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(out$V)))^alpha)
 
        time_V_final=t.out[is.finite(Prob_V)] # only for finite probabilities Prob_V, we record next generation of infections
        Prob_V_final=Prob_V[is.finite(Prob_V)] # only for finite probabilities Prob_V, we record next generation of infections

        V_final=(out$V[is.finite(Prob_V)])
        
        
        if(tt1==1){ #### recording number of exposed contacts of that person and we do that on the first of pandemic (seed day)
          
          no_contact_per_day_final_ref[c(seq(1,length(no_contact_per_day_final),by=1))]= round(no_contact_per_day_final)
          
          log_Vfinal= round(log10(V_final),digits=1)
          log_Vfinal[is.na(log_Vfinal)]=0
          log_Vfinal[log_Vfinal<0]=0
          
          V_ref[c(seq(1,length(log_Vfinal),by=1))] = log_Vfinal
          
          
        }
        
        
        ## check for successful transmissions (since transmission risk is not successful transmission)
        random_numbers= runif(n=length(Prob_V_final),min=0,max=1)  # uniform distribution on the interval from min to max.
        ## account for teh depletion in susceptibles
        susc_prob=(total_pop-sum((II_total_COVID[c(seq(1,tt1,by=1))]))-sum(II_recover_COVID[c(seq(1,tt1,by=1))]))/total_pop
        no_contact_per_day_final_person=c()
        no_contact_per_day_final_person=no_contact_per_day_final[is.finite(Prob_V)]
        
        ## determine next generation of infections
        IND_Success=Prob_V_final>random_numbers ## days on whicj successful transmission occured
        number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob # successful transmissions (incidence) occured to these many new individuals
        time_new_infections_actual=time_V_final[Prob_V_final>random_numbers] ## days on which successful transmission occured
     
        if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){

          total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
          ks_options=seq(1,length(number_of_new_infections_actual),by=1)
          for (ks in ks_options) { ## re-write number_of_new_infections_actual as total_number_of_new_infections over time
            if(ks==1){ 
              total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
            } else {
              total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
            }
            
          }

          total_number_of_new_infections=round(total_number_of_new_infections) ## rounding off as we cannot simulate for floating numbers
          ### remove any duplicate entries 
          time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
          total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]

          
          if(any(total_number_of_new_infections == 0)){   #### if no successful transmission occurred during 150 days
            total_number_of_new_infections==0
            ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
          } else{ #### if successful transmission occurred on any day
            ind=seq(1,length(total_number_of_new_infections),by=1)
          }
          
          index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]
          time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]

          # calculating number of new infections at those time points when successful transmission occured
          cumulative_incidence_temp=total_number_of_new_infections[index_where_new_infections_happen]
   
          ### finding incidence every day
          incidence<-vector()
          y_length=seq(1,length(cumulative_incidence_temp),by=1)
          for (yid in y_length){
            if(yid==1){
              incidence[yid]=cumulative_incidence_temp[yid]
            } else{
              incidence[yid]=cumulative_incidence_temp[yid]-cumulative_incidence_temp[yid-1]  ################ this is imp
            }
          }
          
          
          ####### Now we want to update II_total_COVID  with new infections art future time points
          if (length(time_when_new_infection_is_added)>0){ 
            fl_explore=seq(1,length(time_when_new_infection_is_added),by=1)
          }
          else {
            fl_explore=c()
          }
          
          # First add first generation new infections
          for (ijktt in fl_explore){
            tr=time_when_new_infection_is_added[ijktt]
            tr_actual= 1 + tr + tt1  # actual time at which new infections spread is relative to the current time (tt1)
            ## 1 here reflects time from infection to start of the viral replication
            Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time
            ## update II_total_COVID
             II_total_COVID[tr_actual]=II_total_COVID[tr_actual]+Ir
            
             
             ## detemine if these were superspreader events of varying sizes
            Sr = Ir>SSE_threshold
            SSE_total_COVID[tr_actual] =  SSE_total_COVID[tr_actual] +  Sr
            
            Sr2 = Ir>SSE_threshold2
            SSE_total_COVID2[tr_actual] =  SSE_total_COVID2[tr_actual] +  Sr2
            
            Sr3 = Ir>SSE_threshold3
            SSE_total_COVID3[tr_actual] =  SSE_total_COVID3[tr_actual] +  Sr3
            
            Sr4 = Ir>SSE_threshold4
            SSE_total_COVID4[tr_actual] =  SSE_total_COVID4[tr_actual] +  Sr4
            
            Sr5 = Ir>SSE_threshold5
            SSE_total_COVID5[tr_actual] =  SSE_total_COVID5[tr_actual] +  Sr5
            
            Sr6 = Ir>SSE_threshold6
            SSE_total_COVID6[tr_actual] =  SSE_total_COVID6[tr_actual] +  Sr6

          }
    
        } 
        
      }  ### simulating epidemics for each active infection
      
    }  ### running for all active infections at a particular time
    
    
   #### Stop simulaton is cumulative number of infections at any time reaches 1000 infections
    if(sum(II_total_COVID)>=infection_threshold){  
      break
    }
    
  }  ##### for all time points
  
  
  
  Cumulative_infections[tt1]=sum(II_total_COVID)
  
  ### Store final output
  
  II_metrics[seq(1+no_days*(instance_index-1),no_days*instance_index,by=1),]=cbind(tt_total_COVID[seq(1,no_days,by=1)],
                                                                                   II_total_COVID[seq(1,no_days,by=1)],
                                                                                   II_recover_COVID[seq(1,no_days,by=1)],
                                                                                   rep(instance_index,no_days),
                                                                                   SSE_total_COVID,
                                                                                   Cumulative_infections,
                                                                                   SSE_total_COVID2,
                                                                                   SSE_total_COVID3,
                                                                                   SSE_total_COVID4,
                                                                                   SSE_total_COVID5 ,
                                                                                   SSE_total_COVID6,
                                                                                   V_ref , 
                                                                                   no_contact_per_day_final_ref)
  
  
  instance_index=instance_index+1
}       


## write the output
write.csv(II_metrics,"COVID_5_inst_to_100_infections_Re_general.csv",row.names=FALSE)





##################################################################
################################# Plotting the results
##################################################################

infection_threshold=1000
no_runs=1000

nrows = 3 
ncolumns =2
par(mfrow = c(nrows,ncolumns))
par(mfrow = c(nrows,ncolumns),     
    oma = c(0, 0, 0, 0), # 
    mar = c(3.5, 3.5, 0.0, 0.0), 
    mgp = c(2.5, 1, 0),    
    xpd = FALSE)

library(ggplot2)
library(egg)



A_inst = read.csv("COVID_5_inst_to_100_infections_Re_general.csv",header=TRUE)
 
  A_inst=as.data.frame(A_inst)
  A_inst$Instance=as.factor(A_inst$Instance)
  A_inst$time=as.numeric(A_inst$time)
  A_inst$Incidence=as.numeric(A_inst$Incidence)
  A_inst$Recovered=as.numeric(A_inst$Recovered)
  

  #######################################################
  ####################### PREPARING DATA 
  #######################################################
  
  ## times at which we reach infection_threshold
  Dist_time = A_inst$time[A_inst$Cumulative_infections>infection_threshold]
  
  Burn_out = 100*(no_runs-length(Dist_time))/no_runs
  print("Burn out rate")
  print(Burn_out)

  
  ### instances when infection takes off and reaches 1000 infections
  Inst_infection_yes  = unique(A_inst$Instance[A_inst$Cumulative_infections>=infection_threshold])
  Inst_infection_no  = setdiff(seq(1,no_runs,by=1),Inst_infection_yes)
  
  inds_yes <- which(A_inst$Instance %in% Inst_infection_yes)
  inds_no <- which(A_inst$Instance %in% Inst_infection_no)
  
  Data_infection_yes=A_inst[inds_yes,]
  Data_infection_no=A_inst[inds_no,]
  
  #max(Data_infection_no$Cumulative_infections)
  
  SSE_events_infection_yes= Data_infection_yes$time[Data_infection_yes$SuperSpreader>1]
  SSE_events_infection_no= Data_infection_no$time[Data_infection_no$SuperSpreader>1]
  
  #######################################################
  ####################### PLotting Time to reach 1000 infections
  #######################################################
  
  breaks_dist = c(seq(0,500,by=1))
  
  hist(Dist_time, #
       freq= TRUE, # breaks=breaks) #,
       xlab = "Days to reach 1000 infections",
       ylab = "Frequency",
       main="",
       #main=paste("theta",  theta),
       xlim=c(0,100),
       ylim=c(0,15),
       #freq=TRUE,
       #ylim = c(0,0.5),
       breaks=breaks_dist,
       col=scales::alpha('red',.5), # Filling Color
       border = "red",
       right = FALSE,
       xaxt="n",yaxt="n") # "hotpink"
  
  axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
  axis(2,at=c(seq(0,35,by=5)),labels=c(seq(0,35,by=5)),
       cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)
  
  text(73,12,as.expression(bquote(paste("Median")==.(round(median(Dist_time),digits=1)))),
       cex=1.0,col="red")
  text(80,14,as.expression(bquote(paste("Burn out (%)")==.(round(Burn_out,digits=1)))),
       cex=1.0,col="red")
  
  
  rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
  rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
  rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
  rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
  rug(x = c(seq(0,15,by=5))+1, ticksize = -0.04, side = 2)
  rug(x = c(seq(0,15,by=5))+2, ticksize = -0.04, side = 2)
  rug(x = c(seq(0,15,by=5))+3, ticksize = -0.04, side = 2)
  rug(x = c(seq(0,15,by=5))+4, ticksize = -0.04, side = 2)
  
  #######################################################
  ####################### PLotting across 1000 repetitions during the time to reach 1000 infections
  #######################################################
  
    breaks_dist = c(seq(0,500,by=1))
    
    hist(SSE_events_infection_yes, #
         freq= TRUE, # breaks=breaks) #,
         xlab = "Timing of super-spreader events during the time to reach 1000 infection in those trajectories where infections reaches 1000 infections",
         ylab = "Frequency",
         main="",
         #main=paste("theta",  theta),
         xlim=c(0,100),
         ylim=c(0,500),
         #freq=TRUE,
         #ylim = c(0,0.5),
         breaks=breaks_dist,
         col=scales::alpha('red',.5), # Filling Color
         border = "red",
         right = FALSE,
         xaxt="n",yaxt="n") # "hotpink"
    
    
    
    text(80,1510,as.expression(bquote(paste("Old strain, burn out (%)")==.(round(Burn_out,digits=1)))),
         cex=1.0,col="red")

    
    axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
    axis(2,at=c(seq(0,1500,by=100)),labels=c(seq(0,1500,by=100)),
         cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)
    
    rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
    rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
    rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
    rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
    rug(x = c(seq(0,1500,by=100))+50, ticksize = -0.04, side = 2)
    

  
  
  
  
