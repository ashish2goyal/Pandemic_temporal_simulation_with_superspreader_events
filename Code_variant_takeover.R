rm(list = ls())

## please go through COVID_to_1000_infections_Re_general.R for a more commented code to understand basic model and code
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

total_pop=10^6 # Assume a million people in the city

## simulate COVID-19 epidemics for 150 days 
no_days=150

### seed to the pandemic (could be 1, 10, 100, 1000 infections)
seed_number=1000

infection_threshold=100000 ## stop simulation when total number of infections since the start of the pandemic reaches 100000 infections
SSE_threshold=5 ## Superspreader event threshold of 5
SSE_threshold2=10## Superspreader event threshold of 10
SSE_threshold3=20## Superspreader event threshold of 20
SSE_threshold4=40## Superspreader event threshold of 40
SSE_threshold5=60## Superspreader event threshold of 60
SSE_threshold6=80## Superspreader event threshold of 80


no_runs=100 ## number of replicates as the simulation is stochastic

#### Pre-define matrix for storage
II_metrics = matrix(0,nrow=no_days*no_runs,ncol=11)   # for no_days days time, record incidence, recovered and instance
colnames(II_metrics)=c("time","Total_Incidence","Instance", "Incidence_clade1","Incidence_clade2",
                       "Incidence_clade3","Incidence_clade4","Incidence_clade5","Incidence_clade6",
                       "Incidence_clade7","Incidence_clade8")

instance_index=1###counter for each replicate

II_metrics_old <- c() ## initialize a temporary array

for (instances in seq(1,no_runs,by=1)){    ### running simulation for each replicate 
  
  ## tt_total_COVID is total COVID infections at all times
  tt_total_COVID = matrix(0,nrow=no_days,ncol=1)  # we run the epidemics for a maximum of no_days
  tt_total_COVID = seq(1,length(tt_total_COVID),by=1)
  
  ## tt_total_COVID_clade1 is total COVID infections for clade 1 (Re=0.8) at all times
  tt_total_COVID_clade1 = matrix(0,nrow=no_days,ncol=1)  # for 45 days
  tt_total_COVID_clade1 = seq(1,length(tt_total_COVID_clade1),by=1)
  ## tt_total_COVID_clade2 is total COVID infections for clade 2 (Re=1.0) at all times
  tt_total_COVID_clade2 = tt_total_COVID_clade1
  ## tt_total_COVID_clade3 is total COVID infections for clade 3 (Re=1.2) at all times
  tt_total_COVID_clade3 = tt_total_COVID_clade1
  ## tt_total_COVID_clade4 is total COVID infections for clade 4 (Re=1.4) at all times
  tt_total_COVID_clade4 = tt_total_COVID_clade1
  ## tt_total_COVID_clade5 is total COVID infections for clade 5 (Re=1.6) at all times
  tt_total_COVID_clade5 = tt_total_COVID_clade1
  ## tt_total_COVID_clade6 is total COVID infections for clade 6 (Re=1.8) at all times
  tt_total_COVID_clade6 = tt_total_COVID_clade1
  ## tt_total_COVID_clade7 is total COVID infections for clade 7 (Re=2.0) at all times
  tt_total_COVID_clade7 = tt_total_COVID_clade1
  ## tt_total_COVID_clade8 is total COVID infections for clade 8 (Re=2.2) at all times
  tt_total_COVID_clade8 = tt_total_COVID_clade1
  
  clades_options=c(seq(1,8,by=1)) ## for 8 pre-defined clades
  
  II_total_COVID_clade1 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days
  II_total_COVID_clade1[1]=seed_number ### starting with which clade -- we start with clade 1 of Re=0.8
  
  II_total_COVID_clade2 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days
  II_total_COVID_clade2[1]=0 ### at t=0, clade 2 is absent
  
  II_total_COVID_clade3 = matrix(0,nrow=no_days,ncol=1)  # we run the epidemics for a maximum of no_days
  II_total_COVID_clade3[1]=0 ### at t=0, clade 3 is absent
  
  II_total_COVID_clade4 = matrix(0,nrow=no_days,ncol=1)  # we run the epidemics for a maximum of no_days
  II_total_COVID_clade4[1]=0 ### at t=0, clade 4 is absent
  
  II_total_COVID_clade5 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  II_total_COVID_clade5[1]=0 ### at t=0, clade 5 is absent
  
  II_total_COVID_clade6 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days
  II_total_COVID_clade6[1]=0 ### at t=0, clade 6 is absent
  
  II_total_COVID_clade7 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days
  II_total_COVID_clade7[1]=0 ### at t=0, clade 7 is absent
  
  II_total_COVID_clade8 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days
  II_total_COVID_clade8[1]=0 ### at t=0, clade 8 is absent
  
  ### total infections will be cumulative of all clades
  II_total_COVID = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days
  II_total_COVID = II_total_COVID_clade1 + II_total_COVID_clade2+II_total_COVID_clade3+
    II_total_COVID_clade4+II_total_COVID_clade5+II_total_COVID_clade6+II_total_COVID_clade7+II_total_COVID_clade8

  
  
  SSE_total_COVID = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID2 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID3 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID4 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID5 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  SSE_total_COVID6 = matrix(0,nrow=no_days,ncol=1)   # we run the epidemics for a maximum of no_days 
  
  ### first_SSE_event_size is the first SSE event size
  first_SSE_event_size =  matrix(0,nrow=no_days,ncol=1)  # we run the epidemics for a maximum of no_days 
  ### II_recover_COVID_clade1 is the naturally recovered for clade 1 (such that they have no to low viral load for transmission): at t=0. they are zero
  II_recover_COVID_clade1 = matrix(0,nrow=no_days,ncol=1)   # for 45 days
  ### II_recover_COVID_clade2 is the naturally recovered for clade 2
  II_recover_COVID_clade2 = II_recover_COVID_clade1
  ### II_recover_COVID_clade3 is the naturally recovered for clade 3
  II_recover_COVID_clade3 = II_recover_COVID_clade1
  ### II_recover_COVID_clade4 is the naturally recovered for clade 4
  II_recover_COVID_clade4 = II_recover_COVID_clade1
  ### II_recover_COVID_clade5 is the naturally recovered for clade 5
  II_recover_COVID_clade5 = II_recover_COVID_clade1
  ### II_recover_COVID_clade6 is the naturally recovered for clade 6
  II_recover_COVID_clade6 = II_recover_COVID_clade1
  ### II_recover_COVID_clade7 is the naturally recovered for clade 7
  II_recover_COVID_clade7 = II_recover_COVID_clade1
  ### II_recover_COVID_clade8 is the naturally recovered for clade 8
  II_recover_COVID_clade8 = II_recover_COVID_clade1
  
  ### total recover will be cumulative of all clades
  II_recover_COVID = matrix(0,nrow=no_days,ncol=1)   # for 45 days
  II_recover_COVID = II_recover_COVID_clade1 + II_recover_COVID_clade2+II_recover_COVID_clade3+
    II_recover_COVID_clade4+II_recover_COVID_clade5+II_recover_COVID_clade6+II_recover_COVID_clade7++II_recover_COVID_clade8

  
  Cumulative_infections = matrix(0,nrow=no_days,ncol=1)
  # print(II_total_COVID)
  
  tcurrent=0 ## starting time point
  
  t_options=seq(1,no_days,by=1) # we run the epidemics for a maximum of no_days 
  
  dispersion_options<-c(40) ## dispersion rho parameter value

  
  no_contact_per_day_final_ref <- rep(0,no_days) ## predefine array representing number of contact an individual will have on a daily basis and initiate it with 0
  V_ref <- rep(0,no_days) # predefine array representing viral loads will have on a daily basis and initiate it with 0
  
  person_ID=2 ### Assigning starting infected individual number the ID of 2
  
  
  for (tt1 in t_options) {   # simulate for 150 days
    
    Number_of_Active_infections_clade1 <- c() ## number of active infections of clade 1 on that day
    Number_of_Active_infections_clade2 <- c() ## number of active infections of clade 2 on that day
    Number_of_Active_infections_clade3 <- c() ## number of active infections of clade 3 on that day
    Number_of_Active_infections_clade4 <- c() ## number of active infections of clade 4 on that day
    Number_of_Active_infections_clade5 <- c() ## number of active infections of clade 5 on that day
    Number_of_Active_infections_clade6 <- c() ## number of active infections of clade 6 on that day
    Number_of_Active_infections_clade7 <- c() ## number of active infections of clade 7 on that day
    Number_of_Active_infections_clade8 <- c() ## number of active infections of clade 8 on that day
    
    Days_in_which_infection_will_recover<-c() ## how many days in which that particular individual will recover
    Day_at_which_infection_infection_will_recover<-c() ## On what day a particular individual will recover
    
    Number_of_Active_infections_clade1 = II_total_COVID_clade1[tt1] ## how many active infections of clade 1 do we have at this time point 
    Number_of_Active_infections_clade2 = II_total_COVID_clade2[tt1] ## how many active infections of clade 2 do we have at this time point 
    Number_of_Active_infections_clade3 = II_total_COVID_clade3[tt1] ## how many active infections of clade 3 do we have at this time point 
    Number_of_Active_infections_clade4 = II_total_COVID_clade4[tt1] ## how many active infections of clade 4 do we have at this time point 
    Number_of_Active_infections_clade5 = II_total_COVID_clade5[tt1] ## how many active infections of clade 5 do we have at this time point 
    Number_of_Active_infections_clade6 = II_total_COVID_clade6[tt1] ## how many active infections of clade 6 do we have at this time point 
    Number_of_Active_infections_clade7 = II_total_COVID_clade7[tt1] ## how many active infections of clade 7 do we have at this time point 
    Number_of_Active_infections_clade8 = II_total_COVID_clade8[tt1] ## how many active infections of clade 8 do we have at this time point 

    
    #### we have to simulate each active infection viral loads in all 8 clades so we use a for loop
    
    for(clades in clades_options) {  
      
      if(clades==1){ ## if infection belong to clade 1
        Number_of_Active_infections=II_total_COVID_clade1[tt1]

        if (Number_of_Active_infections>0){
          first_loop=seq(1,Number_of_Active_infections,by=1)

          for (yu in first_loop){ 
      
            tzero=0 
            ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
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
            
            dT=1 # this has to be less than 1
            init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
            t.out <- seq(0,30,by=dT)
            params=c()
            out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
            
            
            ###### thsi particular infections will also recover based on viral suppression duration
            Days_in_which_infection_will_recover= length(which(increase_VL*out$V>(10^2))) # out$V >100 means viral load detectable 
            Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1
            II_recover_COVID_clade1[Day_at_which_infection_infection_will_recover]=II_recover_COVID_clade1[Day_at_which_infection_infection_will_recover]+1
            
            
            ### Now to get new infections and when
            AUC_factor=7 # tzero	alpha	betat	no_contact	rho
            alpha=0.8
            
            if(clades==1){
              no_contact_per_day_Re = 1.1
            }else if(clades ==2){
              no_contact_per_day_Re = 2.3
            }else if(clades ==3){
              no_contact_per_day_Re = 3.1
            } else if(clades ==4){
              no_contact_per_day_Re = 3.5
            }else if(clades ==5){
              no_contact_per_day_Re = 3.75
            }else if(clades ==6){
              no_contact_per_day_Re = 4.0
            }else if(ReInput==7){
              no_contact_per_day_Re = 5.0
            }else if(ReInput==8){
              no_contact_per_day_Re = 5.5
            }
            
            no_contact_per_day=no_contact_per_day_Re ##
            dispersion_options= c(40) ## value for rho
            
            no_contact_per_day=1.1 ## clade 1 default
            original_strain_re=0.8  ## clade 1 default
            
            ## simulate number of exposed contact rate for that individual at all time points
            no_contact_per_day_final<-c()
            no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                            scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 

            Prob_V=(pmax(0,(increase_VL*out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(increase_VL*out$V)))^alpha)
            
            time_V_final=t.out[is.finite(Prob_V)]
            Prob_V_final=Prob_V[is.finite(Prob_V)]
            V_final=(increase_VL*out$V[is.finite(Prob_V)])
            
            
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
            
            number_of_new_infections_ideal_case= dT*no_contact_per_day_final_person*Prob_V_final*susc_prob
            IND_Success=Prob_V_final>random_numbers
            number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob
            time_new_infections_actual=time_V_final[Prob_V_final>random_numbers]
 
            
            if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){

              total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
              ks_options=seq(1,length(number_of_new_infections_actual),by=1)
              for (ks in ks_options) {
                if(ks==1){ 
                  total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                } else {
                  total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                }
                
              }

              total_number_of_new_infections=round(total_number_of_new_infections)
              time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
              
              total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]

              if(any(total_number_of_new_infections == 0)){                    #### new
                total_number_of_new_infections==0
                ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
              } else{
                ind=seq(1,length(total_number_of_new_infections),by=1)
              }
              
              index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]
              
              time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]
              
              time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added-tzero  ################ this is imp  (if there is incubation period for presymptomatic and symptomatic incidence)
              # calculating number of new infections at these time points
              
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
              
              ### this time new infections could be of new clade and thus, for that we define next_strain_mat for those time points at which successful transmission occured
              next_strain_mat <- c()
              
              if (length(time_when_new_infection_is_added_relative_to_tzero)>0){ 

                fl_explore=seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1)
                
                for (op in c(seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1))){
                  no_contact_per_day_options = c(1.1, 2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) ## theta options for new variants
                  
                  ### we assume that 1/100 chance that an infected person will generate a new variant at the of successful transmission and it will also be dominant for transmission 
                  if(runif(1)>0.99){
                    x <- 1:length(no_contact_per_day_options)
                    x1 = sample(x)[1]
                    no_contact_per_day=no_contact_per_day_options[x1] ## new variant gets a new Re
                    new_strain =1
                    if(no_contact_per_day==1.1){
                      new_strain_re=0.8
                    }else if(no_contact_per_day==2.3){
                      new_strain_re=1.0
                    } else if(no_contact_per_day==3.1){
                      new_strain_re=1.2
                    }else if(no_contact_per_day==3.5){
                      new_strain_re=1.4
                    }else if(no_contact_per_day==3.75){
                      new_strain_re=1.6
                    }else if(no_contact_per_day==4.0){
                      new_strain_re=1.8
                    }else if(no_contact_per_day==5.0){
                      new_strain_re=2.0
                    } else {
                      new_strain_re=2.2
                    }

                    next_strain_mat[op] = new_strain_re
                  } else {
                    new_strain =0
                    if(no_contact_per_day==1.1){                      
                      new_strain_re=0.8                     
                      }  else if(no_contact_per_day==2.3){
                      new_strain_re=1.0
                    } else if(no_contact_per_day==3.1){
                      new_strain_re=1.2
                    }else if(no_contact_per_day==3.5){
                      new_strain_re=1.4
                    }else if(no_contact_per_day==3.75){
                      new_strain_re=1.6
                    }else if(no_contact_per_day==4.0){
                      new_strain_re=1.8
                    }else if(no_contact_per_day==5.0){
                      new_strain_re=2.0
                    } else {
                      new_strain_re=2.2
                    }
                    #print("new_strain_re")
                    #print(new_strain_re)
                    next_strain_mat[op] = new_strain_re
                  }
                }
                
                
              }   else {
                fl_explore=c()
              }
              
              
              # First add first generation new infections
              for (ijktt in fl_explore){
                
                
                incubation_period_sample=rnorm(1, mean = 4, sd = 1)
                incubation_period_sample=0 ## we do not account for incubation period and tzero
                tr=time_when_new_infection_is_added_relative_to_tzero[ijktt]
                tr_actual= 1 + tr + tt1 + incubation_period_sample # actual time at which new infections spread is relative to the current time (tt1)
                ## 1 reflects time from infection to viral replication
                Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time

                if(next_strain_mat[ijktt]==0.8){
                  II_total_COVID_clade1[tr_actual]=II_total_COVID_clade1[tr_actual]+Ir
                } else if(next_strain_mat[ijktt]==1){
                  II_total_COVID_clade2[tr_actual]=II_total_COVID_clade2[tr_actual]+Ir
                } else if(next_strain_mat[ijktt]==1.2 ) {
                  II_total_COVID_clade3[tr_actual]=II_total_COVID_clade3[tr_actual]+Ir
                } else if(next_strain_mat[ijktt]==1.4 ) {
                  II_total_COVID_clade4[tr_actual]=II_total_COVID_clade4[tr_actual]+Ir
                } else if(next_strain_mat[ijktt]==1.6 ) {
                  II_total_COVID_clade5[tr_actual]=II_total_COVID_clade5[tr_actual]+Ir
                } else if(next_strain_mat[ijktt]==1.8 ) {
                  II_total_COVID_clade6[tr_actual]=II_total_COVID_clade6[tr_actual]+Ir
                } else if(next_strain_mat[ijktt]==2.0 ) {
                  II_total_COVID_clade7[tr_actual]=II_total_COVID_clade7[tr_actual]+Ir
                } else if(next_strain_mat[ijktt]==2.2 ) {
                  II_total_COVID_clade8[tr_actual]=II_total_COVID_clade8[tr_actual]+Ir
                }
                
                II_total_COVID_clade1 <- II_total_COVID_clade1[!is.na(II_total_COVID_clade1)]
                II_total_COVID_clade2 <- II_total_COVID_clade2[!is.na(II_total_COVID_clade2)]
                II_total_COVID_clade3 <- II_total_COVID_clade3[!is.na(II_total_COVID_clade3)]
                II_total_COVID_clade4 <- II_total_COVID_clade4[!is.na(II_total_COVID_clade4)]
                II_total_COVID_clade5 <- II_total_COVID_clade5[!is.na(II_total_COVID_clade5)]
                II_total_COVID_clade6 <- II_total_COVID_clade6[!is.na(II_total_COVID_clade6)]
                II_total_COVID_clade7 <- II_total_COVID_clade7[!is.na(II_total_COVID_clade7)]
                II_total_COVID_clade8 <- II_total_COVID_clade8[!is.na(II_total_COVID_clade8)]
                
                II_total_COVID= II_total_COVID_clade1 + II_total_COVID_clade2+
                  II_total_COVID_clade3 + II_total_COVID_clade4+
                  II_total_COVID_clade5 + II_total_COVID_clade6+
                  II_total_COVID_clade7+ II_total_COVID_clade8
                
                
              }
              
              
              
            } 
            
            
            person_ID = person_ID+1  ## update ID for the next infected individual
            
          }  ### simulating epidemics for each active infetion
          
 
        }  ### running fro all active infections at a particular time
        
        
        #### We repeat what we did for clade 1 for all clades below
        
      } else if(clades==2){
          Number_of_Active_infections=II_total_COVID_clade2[tt1]
          
          if (Number_of_Active_infections>0){
            first_loop=seq(1,Number_of_Active_infections,by=1)

            for (yu in first_loop){ ## simulating epidemics originating because of each active infection at time t in that clade
              
              ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
              #ijxxx=1
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
              
              dT=1 # this has to be less than 1
              init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
              t.out <- seq(0,30,by=dT)
              params=c()
              out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
              
              
              ###### thsi particular infections will also recover based on viral suppression duration
              Days_in_which_infection_will_recover= length(which(increase_VL*out$V>(10^2))) # out$V >100 means viral load detectable -
              Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1
              II_recover_COVID_clade2[Day_at_which_infection_infection_will_recover]=II_recover_COVID_clade2[Day_at_which_infection_infection_will_recover]+1

              ### Now to get new infections and when
              
              AUC_factor=7 # tzero	alpha	betat	no_contact	rho
              alpha=0.8
              
              if(clades==1){
                no_contact_per_day_Re = 1.1
              }else if(clades ==2){
                no_contact_per_day_Re = 2.3
              }else if(clades ==3){
                no_contact_per_day_Re = 3.1
              } else if(clades ==4){
                no_contact_per_day_Re = 3.5
              }else if(clades ==5){
                no_contact_per_day_Re = 3.75
              }else if(clades ==6){
                no_contact_per_day_Re = 4.0
              }else if(ReInput==7){
                no_contact_per_day_Re = 5.0
              }else if(ReInput==8){
                no_contact_per_day_Re = 5.5
              }
              
              no_contact_per_day=no_contact_per_day_Re ##
              dispersion_options= c(40) ## value for rho
              
              
              no_contact_per_day=2.3 ## clade 2 default
              original_strain_re=1.0  ## clade 2 default
              
              ## different contact rate at each time point
              no_contact_per_day_final<-c()
              no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                              scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
             
              Prob_V=(pmax(0,(increase_VL*out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(increase_VL*out$V)))^alpha)
              
             
              time_V_final=t.out[is.finite(Prob_V)]
              Prob_V_final=Prob_V[is.finite(Prob_V)]

              V_final=(increase_VL*out$V[is.finite(Prob_V)])
              
              
              if(tt1==1){
                no_contact_per_day_final_ref[c(seq(1,length(no_contact_per_day_final),by=1))]= round(no_contact_per_day_final)
                
                log_Vfinal= round(log10(V_final),digits=1)
                log_Vfinal[is.na(log_Vfinal)]=0
                log_Vfinal[log_Vfinal<0]=0
                
                V_ref[c(seq(1,length(log_Vfinal),by=1))] = log_Vfinal
                
                
              }
              
           
              ## check for successful transmissions (since transmission risk is not successful transmission)
              random_numbers= runif(n=length(Prob_V_final),min=0,max=1)  # uniform distribution on the interval from min to max.
              susc_prob=(total_pop-sum((II_total_COVID[c(seq(1,tt1,by=1))]))-sum(II_recover_COVID[c(seq(1,tt1,by=1))]))/total_pop
              no_contact_per_day_final_person=c()
              no_contact_per_day_final_person=no_contact_per_day_final[is.finite(Prob_V)]
              
              number_of_new_infections_ideal_case= dT*no_contact_per_day_final_person*Prob_V_final*susc_prob

              IND_Success=Prob_V_final>random_numbers
              number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob
              time_new_infections_actual=time_V_final[Prob_V_final>random_numbers]

              
              if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){

                total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                for (ks in ks_options) {
                  if(ks==1){ 
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                  } else {
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                  }
                  
                }

                total_number_of_new_infections=round(total_number_of_new_infections)

                time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
                
                total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]

                if(any(total_number_of_new_infections == 0)){                    #### new
                  total_number_of_new_infections==0
                  ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
                } else{
                  ind=seq(1,length(total_number_of_new_infections),by=1)
                }

                index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]

                time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]

                time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added-tzero  ################ this is imp  (add incubation period)

                # caluclating number of new infections at these time points
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

                next_strain_mat <- c()
                
                ####### Add this data to II_total_COVID 
                if (length(time_when_new_infection_is_added_relative_to_tzero)>0){ 
                
                  fl_explore=seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1)
                  
                  for (op in c(seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1))){
                    ## c(2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) for 1.0, 1.2 , 1.4, 1.6, 1.8, 2.0 and 2.2 Re values
                    no_contact_per_day_options = c(1.1, 2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) ## options for new varitnas
                    
                    ### 1/100 chance that an infected person will generate a new variant and it will be dominant for transmission 
                    if(runif(1)>0.99){
                      x <- 1:length(no_contact_per_day_options)
                      x1 = sample(x)[1]
                      no_contact_per_day=no_contact_per_day_options[x1] ## new variant gets a new Re
                      new_strain =1
                      if(no_contact_per_day==1.1){                      
                        new_strain_re=0.8                     
                        } else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }
                      next_strain_mat[op] = new_strain_re
                    } else {
                      new_strain =0
                      if(no_contact_per_day==1.1){                      
                        new_strain_re=0.8                     
                        }   else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }
                      next_strain_mat[op] = new_strain_re
                    }
                  }
                  
                  
                }   else {
                  fl_explore=c()
                }
   
                
                # First add first generation new infections
                for (ijktt in fl_explore){
                  
                  
                  incubation_period_sample=rnorm(1, mean = 4, sd = 1)
                  incubation_period_sample=0 ## we do not account for incubation period and tzero
                  tr=time_when_new_infection_is_added_relative_to_tzero[ijktt]
                  tr_actual= 1 + tr + tt1 + incubation_period_sample # actual time at which new infections spread is relative to the current time (tt1)
                  ## 1 reflects time from infection to viral replication
                  Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time

                  if(next_strain_mat[ijktt]==0.8){
                    II_total_COVID_clade1[tr_actual]=II_total_COVID_clade1[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1){
                    II_total_COVID_clade2[tr_actual]=II_total_COVID_clade2[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.2 ) {
                    II_total_COVID_clade3[tr_actual]=II_total_COVID_clade3[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.4 ) {
                    II_total_COVID_clade4[tr_actual]=II_total_COVID_clade4[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.6 ) {
                    II_total_COVID_clade5[tr_actual]=II_total_COVID_clade5[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.8 ) {
                    II_total_COVID_clade6[tr_actual]=II_total_COVID_clade6[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.0 ) {
                    II_total_COVID_clade7[tr_actual]=II_total_COVID_clade7[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.2 ) {
                    II_total_COVID_clade8[tr_actual]=II_total_COVID_clade8[tr_actual]+Ir
                  }
                  
                  II_total_COVID_clade1 <- II_total_COVID_clade1[!is.na(II_total_COVID_clade1)]
                  II_total_COVID_clade2 <- II_total_COVID_clade2[!is.na(II_total_COVID_clade2)]
                  II_total_COVID_clade3 <- II_total_COVID_clade3[!is.na(II_total_COVID_clade3)]
                  II_total_COVID_clade4 <- II_total_COVID_clade4[!is.na(II_total_COVID_clade4)]
                  II_total_COVID_clade5 <- II_total_COVID_clade5[!is.na(II_total_COVID_clade5)]
                  II_total_COVID_clade6 <- II_total_COVID_clade6[!is.na(II_total_COVID_clade6)]
                  II_total_COVID_clade7 <- II_total_COVID_clade7[!is.na(II_total_COVID_clade7)]
                  
                  II_total_COVID_clade8 <- II_total_COVID_clade8[!is.na(II_total_COVID_clade8)]
                  
                  II_total_COVID= II_total_COVID_clade1 + II_total_COVID_clade2+
                    II_total_COVID_clade3 + II_total_COVID_clade4+
                    II_total_COVID_clade5 + II_total_COVID_clade6+
                    II_total_COVID_clade7+ II_total_COVID_clade8
                  
                }

                
              } 
              
             
              person_ID = person_ID+1 ## update infected individual ID
              
            }  ### simulating epidemics for each active infetion
   
          }  ### running fro all active infections at a particular time
          
        } else if(clades==3){ ## ding the same for clade 3
          Number_of_Active_infections=II_total_COVID_clade3[tt1]

          
          if (Number_of_Active_infections>0){
            first_loop=seq(1,Number_of_Active_infections,by=1)

            for (yu in first_loop){ ## simulating epidemics originating because of each active infection at time t in that clade
              
              ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
              #ijxxx=1
              tzero=4+Parameters$tzero[ijxxx]   # +4 becuase mean is 4 days to convert tzero to days since onset of symptoms (when is this person infectious)
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
              
              dT=1 # this has to be less than 1
              init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
              t.out <- seq(0,30,by=dT)
              params=c()
              out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
              
              
              ###### thsi particular infections will also recover based on viral suppression duration
              Days_in_which_infection_will_recover= length(which(increase_VL*out$V>(10^2))) # out$V >100 means viral load detectable 

              Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1

              II_recover_COVID_clade3[Day_at_which_infection_infection_will_recover]=II_recover_COVID_clade3[Day_at_which_infection_infection_will_recover]+1

              ### Now to get new infections and when
              AUC_factor=7 # tzero	alpha	betat	no_contact	rho
              alpha=0.8
              
              if(clades==1){
                no_contact_per_day_Re = 1.1
              }else if(clades ==2){
                no_contact_per_day_Re = 2.3
              }else if(clades ==3){
                no_contact_per_day_Re = 3.1
              } else if(clades ==4){
                no_contact_per_day_Re = 3.5
              }else if(clades ==5){
                no_contact_per_day_Re = 3.75
              }else if(clades ==6){
                no_contact_per_day_Re = 4.0
              }else if(ReInput==7){
                no_contact_per_day_Re = 5.0
              }else if(ReInput==8){
                no_contact_per_day_Re = 5.5
              }
              
              no_contact_per_day=no_contact_per_day_Re ##
              dispersion_options= c(40) ## value for rho

              no_contact_per_day=3.1 ## clade 3 default
              original_strain_re=1.2  ## clade 3 default
       
              ## different contact rate at each time point
              no_contact_per_day_final<-c()
              no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                              scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
           
              Prob_V=(pmax(0,(increase_VL*out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(increase_VL*out$V)))^alpha)
              
              time_V_final=t.out[is.finite(Prob_V)]
              Prob_V_final=Prob_V[is.finite(Prob_V)]

              V_final=(increase_VL*out$V[is.finite(Prob_V)])
              
              
              if(tt1==1){
                no_contact_per_day_final_ref[c(seq(1,length(no_contact_per_day_final),by=1))]= round(no_contact_per_day_final)
                
                log_Vfinal= round(log10(V_final),digits=1)
                log_Vfinal[is.na(log_Vfinal)]=0
                log_Vfinal[log_Vfinal<0]=0
                
                V_ref[c(seq(1,length(log_Vfinal),by=1))] = log_Vfinal
                
                
              }
              
           
              #######################################
              ################ Another way of estimating actual transmission risk to succeffsul transmision
              #######################################
              ## check for successful transmissions (since transmission risk is not successful transmission)
              random_numbers= runif(n=length(Prob_V_final),min=0,max=1)  # uniform distribution on the interval from min to max.
              susc_prob=(total_pop-sum((II_total_COVID[c(seq(1,tt1,by=1))]))-sum(II_recover_COVID[c(seq(1,tt1,by=1))]))/total_pop
              no_contact_per_day_final_person=c()
              no_contact_per_day_final_person=no_contact_per_day_final[is.finite(Prob_V)]
              
              number_of_new_infections_ideal_case= dT*no_contact_per_day_final_person*Prob_V_final*susc_prob
              IND_Success=Prob_V_final>random_numbers
              number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob
              time_new_infections_actual=time_V_final[Prob_V_final>random_numbers]
              
              
              if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){
                
                total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                for (ks in ks_options) {
                  if(ks==1){ 
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                  } else {
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                  }
                  
                }
 
                total_number_of_new_infections=round(total_number_of_new_infections)

                time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
                
                total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]

                if(any(total_number_of_new_infections == 0)){                    #### new
                  total_number_of_new_infections==0
                  ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
                } else{
                  ind=seq(1,length(total_number_of_new_infections),by=1)
                }
                
                index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]

                time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]

                time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added-tzero  ################ this is imp  (add incubation period)

                # caluclating number of new infections at these time points
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
                
                next_strain_mat <- c()
                
                ####### Add this data to II_total_COVID 
                if (length(time_when_new_infection_is_added_relative_to_tzero)>0){ 
      
                  fl_explore=seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1)
                  
                  for (op in c(seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1))){
                    ## c(2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) for 1.0, 1.2 , 1.4, 1.6, 1.8, 2.0 and 2.2 Re values
                    no_contact_per_day_options = c(1,1, 2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) ## options for new varitnas
                    
                    ### 1/100 chance that an infected person will generate a new variant and it will be dominant for transmission 
                    if(runif(1)>0.99){
                      x <- 1:length(no_contact_per_day_options)
                      x1 = sample(x)[1]
      
                      no_contact_per_day=no_contact_per_day_options[x1] ## new variant gets a new Re
                      new_strain =1
                      if(no_contact_per_day==1.1){                      
                        new_strain_re=0.8                    
                      } else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }

                      next_strain_mat[op] = new_strain_re
                    } else {
                      new_strain =0
                      if(no_contact_per_day==1.1){                      
                        new_strain_re=0.8                    
                      } else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }

                      next_strain_mat[op] = new_strain_re
                    }
                  }
                  
                  
                }   else {
                  fl_explore=c()
                }
                

                # First add first generation new infections
                for (ijktt in fl_explore){
                  
                  
                  incubation_period_sample=rnorm(1, mean = 4, sd = 1)
                  incubation_period_sample=0 ## we do not account for incubation period and tzero
                  tr=time_when_new_infection_is_added_relative_to_tzero[ijktt]
                  tr_actual= 1 + tr + tt1 + incubation_period_sample # actual time at which new infections spread is relative to the current time (tt1)
                  ## 1 reflects time from infection to viral replication
                  Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time
     

                  if(next_strain_mat[ijktt]==0.8){
                    II_total_COVID_clade1[tr_actual]=II_total_COVID_clade1[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1){
                    II_total_COVID_clade2[tr_actual]=II_total_COVID_clade2[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.2 ) {
                    II_total_COVID_clade3[tr_actual]=II_total_COVID_clade3[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.4 ) {
                    II_total_COVID_clade4[tr_actual]=II_total_COVID_clade4[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.6 ) {
                    II_total_COVID_clade5[tr_actual]=II_total_COVID_clade5[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.8 ) {
                    II_total_COVID_clade6[tr_actual]=II_total_COVID_clade6[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.0 ) {
                    II_total_COVID_clade7[tr_actual]=II_total_COVID_clade7[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.2 ) {
                    II_total_COVID_clade8[tr_actual]=II_total_COVID_clade8[tr_actual]+Ir
                  }
                  
                  II_total_COVID_clade1 <- II_total_COVID_clade1[!is.na(II_total_COVID_clade1)]
                  II_total_COVID_clade2 <- II_total_COVID_clade2[!is.na(II_total_COVID_clade2)]
                  II_total_COVID_clade3 <- II_total_COVID_clade3[!is.na(II_total_COVID_clade3)]
                  II_total_COVID_clade4 <- II_total_COVID_clade4[!is.na(II_total_COVID_clade4)]
                  II_total_COVID_clade5 <- II_total_COVID_clade5[!is.na(II_total_COVID_clade5)]
                  II_total_COVID_clade6 <- II_total_COVID_clade6[!is.na(II_total_COVID_clade6)]
                  II_total_COVID_clade7 <- II_total_COVID_clade7[!is.na(II_total_COVID_clade7)]
                  II_total_COVID_clade8 <- II_total_COVID_clade8[!is.na(II_total_COVID_clade8)]
                  
                  II_total_COVID= II_total_COVID_clade1 + II_total_COVID_clade2+
                    II_total_COVID_clade3 + II_total_COVID_clade4+
                    II_total_COVID_clade5 + II_total_COVID_clade6+
                    II_total_COVID_clade7+ II_total_COVID_clade8
                  
                  
                }
                
    
              } 
              
              
              person_ID = person_ID+1 
              
            }  ### simulating epidemics for each active infetion
            
  
          }  ### running fro all active infections at a particular time
          
        }else if(clades==4){ ## Doing the above for clade 4
          Number_of_Active_infections=II_total_COVID_clade4[tt1]

          if (Number_of_Active_infections>0){
            first_loop=seq(1,Number_of_Active_infections,by=1)

            for (yu in first_loop){ ## simulating epidemics originating because of each active infection at time t in that clade

              
              
              ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
              #ijxxx=1
              tzero=4+Parameters$tzero[ijxxx]   # +4 becuase mean is 4 days to convert tzero to days since onset of symptoms (when is this person infectious)
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
              
              dT=1 # this has to be less than 1
              init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
              t.out <- seq(0,30,by=dT)
              params=c()
              out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
              
              
              ###### thsi particular infections will also recover based on viral suppression duration
              # Plotting_VL()
              Days_in_which_infection_will_recover= length(which(increase_VL*out$V>(10^2))) # out$V >100 means viral load detectable 

              Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1

              II_recover_COVID_clade4[Day_at_which_infection_infection_will_recover]=II_recover_COVID_clade4[Day_at_which_infection_infection_will_recover]+1

              
              ### Now to get new infections and when
              
              AUC_factor=7 # tzero	alpha	betat	no_contact	rho
              alpha=0.8
              
              if(clades==1){
                no_contact_per_day_Re = 1.1
              }else if(clades ==2){
                no_contact_per_day_Re = 2.3
              }else if(clades ==3){
                no_contact_per_day_Re = 3.1
              } else if(clades ==4){
                no_contact_per_day_Re = 3.5
              }else if(clades ==5){
                no_contact_per_day_Re = 3.75
              }else if(clades ==6){
                no_contact_per_day_Re = 4.0
              }else if(ReInput==7){
                no_contact_per_day_Re = 5.0
              }else if(ReInput==8){
                no_contact_per_day_Re = 5.5
              }
              
              no_contact_per_day=no_contact_per_day_Re ##
              dispersion_options= c(40) ## value for rho
              
              
              no_contact_per_day=3.5 ## clade 4 default
              original_strain_re=1.4  ## clade 4 default
              
             
              ## different contact rate at each time point
              no_contact_per_day_final<-c()
              no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                              scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
            
              Prob_V=(pmax(0,(increase_VL*out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(increase_VL*out$V)))^alpha)
             
              time_V_final=t.out[is.finite(Prob_V)]
              Prob_V_final=Prob_V[is.finite(Prob_V)]
   
              V_final=(increase_VL*out$V[is.finite(Prob_V)])
              
              
              if(tt1==1){
                no_contact_per_day_final_ref[c(seq(1,length(no_contact_per_day_final),by=1))]= round(no_contact_per_day_final)
                
                log_Vfinal= round(log10(V_final),digits=1)
                log_Vfinal[is.na(log_Vfinal)]=0
                log_Vfinal[log_Vfinal<0]=0
                
                V_ref[c(seq(1,length(log_Vfinal),by=1))] = log_Vfinal
                
                
              }
              
              
          
              ## check for successful transmissions (since transmission risk is not successful transmission)
              random_numbers= runif(n=length(Prob_V_final),min=0,max=1)  # uniform distribution on the interval from min to max.
              susc_prob=(total_pop-sum((II_total_COVID[c(seq(1,tt1,by=1))]))-sum(II_recover_COVID[c(seq(1,tt1,by=1))]))/total_pop

              no_contact_per_day_final_person=c()
              no_contact_per_day_final_person=no_contact_per_day_final[is.finite(Prob_V)]
              
              number_of_new_infections_ideal_case= dT*no_contact_per_day_final_person*Prob_V_final*susc_prob
              
              IND_Success=Prob_V_final>random_numbers
             
              number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob
              time_new_infections_actual=time_V_final[Prob_V_final>random_numbers]

              
              if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){

                total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                for (ks in ks_options) {
                  if(ks==1){ 
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                  } else {
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                  }
                  
                }

                total_number_of_new_infections=round(total_number_of_new_infections)

                time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
                
                total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]

                if(any(total_number_of_new_infections == 0)){                    #### new
                  total_number_of_new_infections==0
                  ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
                } else{
                  ind=seq(1,length(total_number_of_new_infections),by=1)
                }

                index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]

                time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]


                time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added-tzero  ################ this is imp  (add incubation period)

                # caluclating number of new infections at these time points
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

                next_strain_mat <- c()
                
                ####### Add this data to II_total_COVID 
                if (length(time_when_new_infection_is_added_relative_to_tzero)>0){ 
                 
                  fl_explore=seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1)
                  
                  for (op in c(seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1))){
                    ## c(2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) for 1.0, 1.2 , 1.4, 1.6, 1.8, 2.0 and 2.2 Re values
                    no_contact_per_day_options = c(1,1, 2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) ## options for new varitnas
                    
                    ### 1/100 chance that an infected person will generate a new variant and it will be dominant for transmission 
                    if(runif(1)>0.99){
                      x <- 1:length(no_contact_per_day_options)
                      x1 = sample(x)[1]
                      no_contact_per_day=no_contact_per_day_options[x1] ## new variant gets a new Re
                      new_strain =1
                      if(no_contact_per_day==1.1){                       new_strain_re=0.8                    
                      }                       else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }

                      next_strain_mat[op] = new_strain_re
                    } else {
                      new_strain =0
                      if(no_contact_per_day==1.1){                       new_strain_re=0.8                     
                      }                      else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }

                      next_strain_mat[op] = new_strain_re
                    }
                  }
                  
                  
                }   else {
                  fl_explore=c()
                }
                
                
                # First add first generation new infections
                for (ijktt in fl_explore){
                  
                  
                  incubation_period_sample=rnorm(1, mean = 4, sd = 1)
                  incubation_period_sample=0 ## we do not account for incubation period and tzero
                  #print("Inside loop")
                  tr=time_when_new_infection_is_added_relative_to_tzero[ijktt]
                  tr_actual= 1 + tr + tt1 + incubation_period_sample # actual time at which new infections spread is relative to the current time (tt1)
                  ## 1 reflects time from infection to viral replication
                  Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time
                  #print("actual time at which new infections spread is relative to the current time (tt1)")

                  if(next_strain_mat[ijktt]==0.8){
                    II_total_COVID_clade1[tr_actual]=II_total_COVID_clade1[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1){
                    II_total_COVID_clade2[tr_actual]=II_total_COVID_clade2[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.2 ) {
                    II_total_COVID_clade3[tr_actual]=II_total_COVID_clade3[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.4 ) {
                    II_total_COVID_clade4[tr_actual]=II_total_COVID_clade4[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.6 ) {
                    II_total_COVID_clade5[tr_actual]=II_total_COVID_clade5[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.8 ) {
                    II_total_COVID_clade6[tr_actual]=II_total_COVID_clade6[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.0 ) {
                    II_total_COVID_clade7[tr_actual]=II_total_COVID_clade7[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.2 ) {
                    II_total_COVID_clade8[tr_actual]=II_total_COVID_clade8[tr_actual]+Ir
                  }
                  
                  II_total_COVID_clade1 <- II_total_COVID_clade1[!is.na(II_total_COVID_clade1)]
                  II_total_COVID_clade2 <- II_total_COVID_clade2[!is.na(II_total_COVID_clade2)]
                  II_total_COVID_clade3 <- II_total_COVID_clade3[!is.na(II_total_COVID_clade3)]
                  II_total_COVID_clade4 <- II_total_COVID_clade4[!is.na(II_total_COVID_clade4)]
                  II_total_COVID_clade5 <- II_total_COVID_clade5[!is.na(II_total_COVID_clade5)]
                  II_total_COVID_clade6 <- II_total_COVID_clade6[!is.na(II_total_COVID_clade6)]
                  II_total_COVID_clade7 <- II_total_COVID_clade7[!is.na(II_total_COVID_clade7)]
                  
                  II_total_COVID_clade8 <- II_total_COVID_clade8[!is.na(II_total_COVID_clade8)]
                  
                  II_total_COVID= II_total_COVID_clade1 + II_total_COVID_clade2+
                    II_total_COVID_clade3 + II_total_COVID_clade4+
                    II_total_COVID_clade5 + II_total_COVID_clade6+
                    II_total_COVID_clade7+ II_total_COVID_clade8
                  
                  
                }
                
        
              } 
              
              
              person_ID = person_ID+1 
              
            }  ### simulating epidemics for each active infetion
  
            
          }  ### running fro all active infections at a particular time
          
        }else if(clades==5){
          Number_of_Active_infections=II_total_COVID_clade5[tt1]

          if (Number_of_Active_infections>0){
            first_loop=seq(1,Number_of_Active_infections,by=1)

            for (yu in first_loop){ ## simulating epidemics originating because of each active infection at time t in that clade
  
              
              ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
              tzero=4+Parameters$tzero[ijxxx]   # +4 becuase mean is 4 days to convert tzero to days since onset of symptoms (when is this person infectious)
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
              
              dT=1 # this has to be less than 1
              init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
              t.out <- seq(0,30,by=dT)
              params=c()
              out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
              
              
              ###### thsi particular infections will also recover based on viral suppression duration
              Days_in_which_infection_will_recover= length(which(increase_VL*out$V>(10^2))) # out$V >100 means viral load detectable 
              Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1
              II_recover_COVID_clade5[Day_at_which_infection_infection_will_recover]=II_recover_COVID_clade5[Day_at_which_infection_infection_will_recover]+1
              
              
              ### Now to get new infections and when
              AUC_factor=7 # tzero	alpha	betat	no_contact	rho
              alpha=0.8
              
              if(clades==1){
                no_contact_per_day_Re = 1.1
              }else if(clades ==2){
                no_contact_per_day_Re = 2.3
              }else if(clades ==3){
                no_contact_per_day_Re = 3.1
              } else if(clades ==4){
                no_contact_per_day_Re = 3.5
              }else if(clades ==5){
                no_contact_per_day_Re = 3.75
              }else if(clades ==6){
                no_contact_per_day_Re = 4.0
              }else if(ReInput==7){
                no_contact_per_day_Re = 5.0
              }else if(ReInput==8){
                no_contact_per_day_Re = 5.5
              }
              
              no_contact_per_day=no_contact_per_day_Re ##
              dispersion_options= c(40) ## value for rho
              
              no_contact_per_day=3.75 ## clade 5 default
              original_strain_re=1.6  ## clade 5 default

              
              ## different contact rate at each time point
              no_contact_per_day_final<-c()
              no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                              scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 

              Prob_V=(pmax(0,(increase_VL*out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(increase_VL*out$V)))^alpha)
              
              # Prob_V=1-exp(-1*out$V/(AUC_max))  -- another way of calculating infectivity from viral laods (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005481)
              
              time_V_final=t.out[is.finite(Prob_V)]
              Prob_V_final=Prob_V[is.finite(Prob_V)]

              V_final=(increase_VL*out$V[is.finite(Prob_V)])
              
              
              if(tt1==1){
                no_contact_per_day_final_ref[c(seq(1,length(no_contact_per_day_final),by=1))]= round(no_contact_per_day_final)
                
                log_Vfinal= round(log10(V_final),digits=1)
                log_Vfinal[is.na(log_Vfinal)]=0
                log_Vfinal[log_Vfinal<0]=0
                
                V_ref[c(seq(1,length(log_Vfinal),by=1))] = log_Vfinal
                
                
              }

    
              ## check for successful transmissions (since transmission risk is not successful transmission)
              random_numbers= runif(n=length(Prob_V_final),min=0,max=1)  # uniform distribution on the interval from min to max.
              susc_prob=(total_pop-sum((II_total_COVID[c(seq(1,tt1,by=1))]))-sum(II_recover_COVID[c(seq(1,tt1,by=1))]))/total_pop
              no_contact_per_day_final_person=c()
              no_contact_per_day_final_person=no_contact_per_day_final[is.finite(Prob_V)]
              
              number_of_new_infections_ideal_case= dT*no_contact_per_day_final_person*Prob_V_final*susc_prob

              IND_Success=Prob_V_final>random_numbers

              number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob
              time_new_infections_actual=time_V_final[Prob_V_final>random_numbers]

              if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){

                total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                for (ks in ks_options) {
                  if(ks==1){ 
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                  } else {
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                  }
                  
                }
                
                total_number_of_new_infections=round(total_number_of_new_infections)
                
                time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
                
                total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]
                
                if(any(total_number_of_new_infections == 0)){                    #### new
                  total_number_of_new_infections==0
                  ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
                } else{
                  ind=seq(1,length(total_number_of_new_infections),by=1)
                }
                
                index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]

                time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]

                time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added-tzero  ################ this is imp  (add incubation period)

                # caluclating number of new infections at these time points
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
                
                next_strain_mat <- c()
                
                ####### Add this data to II_total_COVID 
                if (length(time_when_new_infection_is_added_relative_to_tzero)>0){ 
                  
                  fl_explore=seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1)
                  
                  for (op in c(seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1))){
                    ## c(2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) for 1.0, 1.2 , 1.4, 1.6, 1.8, 2.0 and 2.2 Re values
                    no_contact_per_day_options = c(1.1, 2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) ## options for new varitnas
                    
                    ### 1/100 chance that an infected person will generate a new variant and it will be dominant for transmission 
                    if(runif(1)>0.99){
                      x <- 1:length(no_contact_per_day_options)
                      x1 = sample(x)[1]
                      no_contact_per_day=no_contact_per_day_options[x1] ## new variant gets a new Re
                      new_strain =1
                      if(no_contact_per_day==1.1){                       new_strain_re=0.8                    
                      }    else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }
                      next_strain_mat[op] = new_strain_re
                    } else {
                      new_strain =0
                      if(no_contact_per_day==1.1){                       new_strain_re=0.8                     
                      }                else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }
                      next_strain_mat[op] = new_strain_re
                    }
                  }
                  
                  
                }   else {
                  fl_explore=c()
                }
                
                
                # First add first generation new infections
                for (ijktt in fl_explore){
                  
                  
                  incubation_period_sample=rnorm(1, mean = 4, sd = 1)
                  incubation_period_sample=0 ## we do not account for incubation period and tzero
                  tr=time_when_new_infection_is_added_relative_to_tzero[ijktt]
                  tr_actual= 1 + tr + tt1 + incubation_period_sample # actual time at which new infections spread is relative to the current time (tt1)
                  ## 1 reflects time from infection to viral replication
                  Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time

                  if(next_strain_mat[ijktt]==0.8){
                    II_total_COVID_clade1[tr_actual]=II_total_COVID_clade1[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1){
                    II_total_COVID_clade2[tr_actual]=II_total_COVID_clade2[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.2 ) {
                    II_total_COVID_clade3[tr_actual]=II_total_COVID_clade3[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.4 ) {
                    II_total_COVID_clade4[tr_actual]=II_total_COVID_clade4[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.6 ) {
                    II_total_COVID_clade5[tr_actual]=II_total_COVID_clade5[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.8 ) {
                    II_total_COVID_clade6[tr_actual]=II_total_COVID_clade6[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.0 ) {
                    II_total_COVID_clade7[tr_actual]=II_total_COVID_clade7[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.2 ) {
                    II_total_COVID_clade8[tr_actual]=II_total_COVID_clade8[tr_actual]+Ir
                  }
                  
                  II_total_COVID_clade1 <- II_total_COVID_clade1[!is.na(II_total_COVID_clade1)]
                  II_total_COVID_clade2 <- II_total_COVID_clade2[!is.na(II_total_COVID_clade2)]
                  II_total_COVID_clade3 <- II_total_COVID_clade3[!is.na(II_total_COVID_clade3)]
                  II_total_COVID_clade4 <- II_total_COVID_clade4[!is.na(II_total_COVID_clade4)]
                  II_total_COVID_clade5 <- II_total_COVID_clade5[!is.na(II_total_COVID_clade5)]
                  II_total_COVID_clade6 <- II_total_COVID_clade6[!is.na(II_total_COVID_clade6)]
                  II_total_COVID_clade7 <- II_total_COVID_clade7[!is.na(II_total_COVID_clade7)]
                  
                  II_total_COVID_clade8 <- II_total_COVID_clade8[!is.na(II_total_COVID_clade8)]
                  
                  II_total_COVID= II_total_COVID_clade1 + II_total_COVID_clade2+
                    II_total_COVID_clade3 + II_total_COVID_clade4+
                    II_total_COVID_clade5 + II_total_COVID_clade6+
                    II_total_COVID_clade7+ II_total_COVID_clade8
                  
                }
                
       
              } 
              
              
              person_ID = person_ID+1 
              
            }  ### simulating epidemics for each active infetion
            
     
          }  ### running fro all active infections at a particular time
          
        }else if(clades==6){
          Number_of_Active_infections=II_total_COVID_clade6[tt1]

          if (Number_of_Active_infections>0){
            first_loop=seq(1,Number_of_Active_infections,by=1)

            for (yu in first_loop){ ## simulating epidemics originating because of each active infection at time t in that clade

              ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
              #ijxxx=1
              tzero=4+Parameters$tzero[ijxxx]   # +4 becuase mean is 4 days to convert tzero to days since onset of symptoms (when is this person infectious)
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
              
              dT=1 # this has to be less than 1
              init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
              t.out <- seq(0,30,by=dT)
              params=c()
              out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
              
              
              ###### thsi particular infections will also recover based on viral suppression duration

              Days_in_which_infection_will_recover= length(which(increase_VL*out$V>(10^2))) # out$V >100 means viral load detectable 

              Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1

              II_recover_COVID_clade6[Day_at_which_infection_infection_will_recover]=II_recover_COVID_clade6[Day_at_which_infection_infection_will_recover]+1

              ### Now to get new infections and when
              
              AUC_factor=7 # tzero	alpha	betat	no_contact	rho
              alpha=0.8
              
              if(clades==1){
                no_contact_per_day_Re = 1.1
              }else if(clades ==2){
                no_contact_per_day_Re = 2.3
              }else if(clades ==3){
                no_contact_per_day_Re = 3.1
              } else if(clades ==4){
                no_contact_per_day_Re = 3.5
              }else if(clades ==5){
                no_contact_per_day_Re = 3.75
              }else if(clades ==6){
                no_contact_per_day_Re = 4.0
              }else if(ReInput==7){
                no_contact_per_day_Re = 5.0
              }else if(ReInput==8){
                no_contact_per_day_Re = 5.5
              }
              
              no_contact_per_day=no_contact_per_day_Re ##
              dispersion_options= c(40) ## value for rho
              
              no_contact_per_day=4.0 ## clade 6 default
              original_strain_re=1.8  ## clade 6 default
              

              ## different contact rate at each time point
              no_contact_per_day_final<-c()
              no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                              scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 

              Prob_V=(pmax(0,(increase_VL*out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(increase_VL*out$V)))^alpha)

              time_V_final=t.out[is.finite(Prob_V)]
              Prob_V_final=Prob_V[is.finite(Prob_V)]

              V_final=(increase_VL*out$V[is.finite(Prob_V)])
              
              
              if(tt1==1){
                no_contact_per_day_final_ref[c(seq(1,length(no_contact_per_day_final),by=1))]= round(no_contact_per_day_final)
                
                log_Vfinal= round(log10(V_final),digits=1)
                log_Vfinal[is.na(log_Vfinal)]=0
                log_Vfinal[log_Vfinal<0]=0
                
                V_ref[c(seq(1,length(log_Vfinal),by=1))] = log_Vfinal
                
                
              }
    
              ## check for successful transmissions (since transmission risk is not successful transmission)
              random_numbers= runif(n=length(Prob_V_final),min=0,max=1)  # uniform distribution on the interval from min to max.
              susc_prob=(total_pop-sum((II_total_COVID[c(seq(1,tt1,by=1))]))-sum(II_recover_COVID[c(seq(1,tt1,by=1))]))/total_pop

              no_contact_per_day_final_person=c()
              no_contact_per_day_final_person=no_contact_per_day_final[is.finite(Prob_V)]
              
              number_of_new_infections_ideal_case= dT*no_contact_per_day_final_person*Prob_V_final*susc_prob

              IND_Success=Prob_V_final>random_numbers
              number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob
              time_new_infections_actual=time_V_final[Prob_V_final>random_numbers]

              
              if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){

                total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                for (ks in ks_options) {
                  if(ks==1){ 
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                  } else {
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                  }
                  
                }

                total_number_of_new_infections=round(total_number_of_new_infections)

                time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
                
                total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]

                if(any(total_number_of_new_infections == 0)){                    #### new
                  total_number_of_new_infections==0
                  ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
                } else{
                  ind=seq(1,length(total_number_of_new_infections),by=1)
                }
     
                index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]
    
                time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]
                
                time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added-tzero  ################ this is imp  (add incubation period)

                # caluclating number of new infections at these time points
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

                next_strain_mat <- c()
                
                ####### Add this data to II_total_COVID 
                if (length(time_when_new_infection_is_added_relative_to_tzero)>0){ 
                  
                  fl_explore=seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1)
                  
                  for (op in c(seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1))){
                    ## c(2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) for 1.0, 1.2 , 1.4, 1.6, 1.8, 2.0 and 2.2 Re values
                    no_contact_per_day_options = c(1.1, 2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) ## options for new varitnas
                    
                    ### 1/100 chance that an infected person will generate a new variant and it will be dominant for transmission 
                    if(runif(1)>0.99){
                      x <- 1:length(no_contact_per_day_options)
                      x1 = sample(x)[1]
                      no_contact_per_day=no_contact_per_day_options[x1] ## new variant gets a new Re
                      new_strain =1
                      if(no_contact_per_day==1.1){                       new_strain_re=0.8                     
                      }else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }
                      next_strain_mat[op] = new_strain_re
                    } else {
                      new_strain =0
                      if(no_contact_per_day==1.1){                       new_strain_re=0.8                     
                      }else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }
                      next_strain_mat[op] = new_strain_re
                    }
                  }
                  
                  
                }   else {
                  fl_explore=c()
                }
                
                # First add first generation new infections
                for (ijktt in fl_explore){
                  
                  
                  incubation_period_sample=rnorm(1, mean = 4, sd = 1)
                  incubation_period_sample=0 ## we do not account for incubation period and tzero
                  tr=time_when_new_infection_is_added_relative_to_tzero[ijktt]
                  tr_actual= 1 + tr + tt1 + incubation_period_sample # actual time at which new infections spread is relative to the current time (tt1)
                  ## 1 reflects time from infection to viral replication
                  Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time

                  if(next_strain_mat[ijktt]==0.8){
                    II_total_COVID_clade1[tr_actual]=II_total_COVID_clade1[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1){
                    II_total_COVID_clade2[tr_actual]=II_total_COVID_clade2[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.2 ) {
                    II_total_COVID_clade3[tr_actual]=II_total_COVID_clade3[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.4 ) {
                    II_total_COVID_clade4[tr_actual]=II_total_COVID_clade4[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.6 ) {
                    II_total_COVID_clade5[tr_actual]=II_total_COVID_clade5[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.8 ) {
                    II_total_COVID_clade6[tr_actual]=II_total_COVID_clade6[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.0 ) {
                    II_total_COVID_clade7[tr_actual]=II_total_COVID_clade7[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.2 ) {
                    II_total_COVID_clade8[tr_actual]=II_total_COVID_clade8[tr_actual]+Ir
                  }
                  
                  II_total_COVID_clade1 <- II_total_COVID_clade1[!is.na(II_total_COVID_clade1)]
                  II_total_COVID_clade2 <- II_total_COVID_clade2[!is.na(II_total_COVID_clade2)]
                  II_total_COVID_clade3 <- II_total_COVID_clade3[!is.na(II_total_COVID_clade3)]
                  II_total_COVID_clade4 <- II_total_COVID_clade4[!is.na(II_total_COVID_clade4)]
                  II_total_COVID_clade5 <- II_total_COVID_clade5[!is.na(II_total_COVID_clade5)]
                  II_total_COVID_clade6 <- II_total_COVID_clade6[!is.na(II_total_COVID_clade6)]
                  II_total_COVID_clade7 <- II_total_COVID_clade7[!is.na(II_total_COVID_clade7)]
                  
                  II_total_COVID_clade8 <- II_total_COVID_clade8[!is.na(II_total_COVID_clade8)]
                  
                  II_total_COVID= II_total_COVID_clade1 + II_total_COVID_clade2+
                    II_total_COVID_clade3 + II_total_COVID_clade4+
                    II_total_COVID_clade5 + II_total_COVID_clade6+
                    II_total_COVID_clade7+ II_total_COVID_clade8
                  
                }
           
                
              } 

              person_ID = person_ID+1 
              
            }  ### simulating epidemics for each active infetion
            

          }  ### running fro all active infections at a particular time
          
        }else if(clades==7){
          Number_of_Active_infections=II_total_COVID_clade7[tt1]

          if (Number_of_Active_infections>0){
            first_loop=seq(1,Number_of_Active_infections,by=1)

            for (yu in first_loop){ ## simulating epidemics originating because of each active infection at time t in that clade

              
              ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
              #ijxxx=1
              tzero=4+Parameters$tzero[ijxxx]   # +4 becuase mean is 4 days to convert tzero to days since onset of symptoms (when is this person infectious)
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
              
              dT=1 # this has to be less than 1
              init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
              t.out <- seq(0,30,by=dT)
              params=c()
              out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
              
              
              ###### thsi particular infections will also recover based on viral suppression duration
              # Plotting_VL()
              Days_in_which_infection_will_recover= length(which(increase_VL*out$V>(10^2))) # out$V >100 means viral load detectable 
              
              Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1

              II_recover_COVID_clade7[Day_at_which_infection_infection_will_recover]=II_recover_COVID_clade7[Day_at_which_infection_infection_will_recover]+1

              ### Now to get new infections and when
              
              AUC_factor=7 # tzero	alpha	betat	no_contact	rho
              alpha=0.8
              
              if(clades==1){
                no_contact_per_day_Re = 1.1
              }else if(clades ==2){
                no_contact_per_day_Re = 2.3
              }else if(clades ==3){
                no_contact_per_day_Re = 3.1
              } else if(clades ==4){
                no_contact_per_day_Re = 3.5
              }else if(clades ==5){
                no_contact_per_day_Re = 3.75
              }else if(clades ==6){
                no_contact_per_day_Re = 4.0
              }else if(ReInput==7){
                no_contact_per_day_Re = 5.0
              }else if(ReInput==8){
                no_contact_per_day_Re = 5.5
              }
              
              no_contact_per_day=no_contact_per_day_Re ##
              dispersion_options= c(40) ## value for rho
              
              no_contact_per_day=5.0 ## clade 7 default
              original_strain_re=2.0  ## clade 7 default

              
              ## different contact rate at each time point
              no_contact_per_day_final<-c()
              no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                              scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
  
              Prob_V=(pmax(0,(increase_VL*out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(increase_VL*out$V)))^alpha)
              

              time_V_final=t.out[is.finite(Prob_V)]
              Prob_V_final=Prob_V[is.finite(Prob_V)]

              V_final=(increase_VL*out$V[is.finite(Prob_V)])
              
              
              if(tt1==1){
                no_contact_per_day_final_ref[c(seq(1,length(no_contact_per_day_final),by=1))]= round(no_contact_per_day_final)
                
                log_Vfinal= round(log10(V_final),digits=1)
                log_Vfinal[is.na(log_Vfinal)]=0
                log_Vfinal[log_Vfinal<0]=0
                
                V_ref[c(seq(1,length(log_Vfinal),by=1))] = log_Vfinal
                
                
              }
              

              ## check for successful transmissions (since transmission risk is not successful transmission)
              random_numbers= runif(n=length(Prob_V_final),min=0,max=1)  # uniform distribution on the interval from min to max.
              susc_prob=(total_pop-sum((II_total_COVID[c(seq(1,tt1,by=1))]))-sum(II_recover_COVID[c(seq(1,tt1,by=1))]))/total_pop
              no_contact_per_day_final_person=c()
              no_contact_per_day_final_person=no_contact_per_day_final[is.finite(Prob_V)]
              
              number_of_new_infections_ideal_case= dT*no_contact_per_day_final_person*Prob_V_final*susc_prob
              IND_Success=Prob_V_final>random_numbers
              number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob
              time_new_infections_actual=time_V_final[Prob_V_final>random_numbers]

              if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){

                total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                for (ks in ks_options) {
                  if(ks==1){ 
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                  } else {
                    total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                  }
                  
                }

                total_number_of_new_infections=round(total_number_of_new_infections)
 
                time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
                
                total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]
     
                if(any(total_number_of_new_infections == 0)){                    #### new
                  total_number_of_new_infections==0
                  ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
                } else{
                  ind=seq(1,length(total_number_of_new_infections),by=1)
                }

                index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]

                time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]

                time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added-tzero  ################ this is imp  (add incubation period)

                # caluclating number of new infections at these time points
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

                next_strain_mat <- c()
                
                ####### Add this data to II_total_COVID 
                if (length(time_when_new_infection_is_added_relative_to_tzero)>0){ 

                  fl_explore=seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1)
                  
                  for (op in c(seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1))){
                    ## c(2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) for 1.0, 1.2 , 1.4, 1.6, 1.8, 2.0 and 2.2 Re values
                    no_contact_per_day_options = c(1.1,2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) ## options for new varitnas
                    
                    ### 1/100 chance that an infected person will generate a new variant and it will be dominant for transmission 
                    if(runif(1)>0.99){
                      x <- 1:length(no_contact_per_day_options)
                      x1 = sample(x)[1]
                      no_contact_per_day=no_contact_per_day_options[x1] ## new variant gets a new Re
                      new_strain =1
                      if(no_contact_per_day==1.1){                       new_strain_re=0.8                     
                      }else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }
                      next_strain_mat[op] = new_strain_re
                    } else {
                      new_strain =0
                      if(no_contact_per_day==1.1){                       new_strain_re=0.8                     
                      }else if(no_contact_per_day==2.3){
                        new_strain_re=1.0
                      } else if(no_contact_per_day==3.1){
                        new_strain_re=1.2
                      }else if(no_contact_per_day==3.5){
                        new_strain_re=1.4
                      }else if(no_contact_per_day==3.75){
                        new_strain_re=1.6
                      }else if(no_contact_per_day==4.0){
                        new_strain_re=1.8
                      }else if(no_contact_per_day==5.0){
                        new_strain_re=2.0
                      } else {
                        new_strain_re=2.2
                      }

                      next_strain_mat[op] = new_strain_re
                    }
                  }
                  
                  
                }   else {
                  fl_explore=c()
                }
                
                # First add first generation new infections
                for (ijktt in fl_explore){
                  
                  
                  incubation_period_sample=rnorm(1, mean = 4, sd = 1)
                  incubation_period_sample=0 ## we do not account for incubation period and tzero
                  tr=time_when_new_infection_is_added_relative_to_tzero[ijktt]
                  tr_actual= 1 + tr + tt1 + incubation_period_sample # actual time at which new infections spread is relative to the current time (tt1)
                  ## 1 reflects time from infection to viral replication
                  Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time

                  if(next_strain_mat[ijktt]==0.8){
                    II_total_COVID_clade1[tr_actual]=II_total_COVID_clade1[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1){
                    II_total_COVID_clade2[tr_actual]=II_total_COVID_clade2[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.2 ) {
                    II_total_COVID_clade3[tr_actual]=II_total_COVID_clade3[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.4 ) {
                    II_total_COVID_clade4[tr_actual]=II_total_COVID_clade4[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.6 ) {
                    II_total_COVID_clade5[tr_actual]=II_total_COVID_clade5[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==1.8 ) {
                    II_total_COVID_clade6[tr_actual]=II_total_COVID_clade6[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.0 ) {
                    II_total_COVID_clade7[tr_actual]=II_total_COVID_clade7[tr_actual]+Ir
                  } else if(next_strain_mat[ijktt]==2.2 ) {
                    II_total_COVID_clade8[tr_actual]=II_total_COVID_clade8[tr_actual]+Ir
                  }
                  
                  II_total_COVID_clade1 <- II_total_COVID_clade1[!is.na(II_total_COVID_clade1)]
                  II_total_COVID_clade2 <- II_total_COVID_clade2[!is.na(II_total_COVID_clade2)]
                  II_total_COVID_clade3 <- II_total_COVID_clade3[!is.na(II_total_COVID_clade3)]
                  II_total_COVID_clade4 <- II_total_COVID_clade4[!is.na(II_total_COVID_clade4)]
                  II_total_COVID_clade5 <- II_total_COVID_clade5[!is.na(II_total_COVID_clade5)]
                  II_total_COVID_clade6 <- II_total_COVID_clade6[!is.na(II_total_COVID_clade6)]
                  II_total_COVID_clade7 <- II_total_COVID_clade7[!is.na(II_total_COVID_clade7)]
                  II_total_COVID_clade8 <- II_total_COVID_clade8[!is.na(II_total_COVID_clade8)]
                  
                  II_total_COVID= II_total_COVID_clade1 + II_total_COVID_clade2+
                    II_total_COVID_clade3 + II_total_COVID_clade4+
                    II_total_COVID_clade5 + II_total_COVID_clade6+
                    II_total_COVID_clade7+ II_total_COVID_clade8
                  
                  
                }
                

              } 

              person_ID = person_ID+1 
              
            }  ### simulating epidemics for each active infetion

            
          }  ### running fro all active infections at a particular time
 
          
        } else if(clades==8){
          Number_of_Active_infections=II_total_COVID_clade8[tt1]
          
                        if (Number_of_Active_infections>0){
                          first_loop=seq(1,Number_of_Active_infections,by=1)

                          for (yu in first_loop){ ## simulating epidemics originating because of each active infection at time t in that clade
    
                            ijxxx=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
                            #ijxxx=1
                            tzero=4+Parameters$tzero[ijxxx]   # +4 becuase mean is 4 days to convert tzero to days since onset of symptoms (when is this person infectious)
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
                            
                            dT=1 # this has to be less than 1
                            init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
                            t.out <- seq(0,30,by=dT)
                            params=c()
                            out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
                            
                            
                            ###### thsi particular infections will also recover based on viral suppression duration
                            # Plotting_VL()
                            Days_in_which_infection_will_recover= length(which(increase_VL*out$V>(10^2))) # out$V >100 means viral load detectable 

                            Day_at_which_infection_infection_will_recover= Days_in_which_infection_will_recover + tt1 +1

                            II_recover_COVID_clade8[Day_at_which_infection_infection_will_recover]=II_recover_COVID_clade8[Day_at_which_infection_infection_will_recover]+1
                            
                            ### Now to get new infections and when
                            
                            AUC_factor=7 # tzero	alpha	betat	no_contact	rho
                            alpha=0.8
                            
                            if(clades==1){
                              no_contact_per_day_Re = 1.1
                            }else if(clades ==2){
                              no_contact_per_day_Re = 2.3
                            }else if(clades ==3){
                              no_contact_per_day_Re = 3.1
                            } else if(clades ==4){
                              no_contact_per_day_Re = 3.5
                            }else if(clades ==5){
                              no_contact_per_day_Re = 3.75
                            }else if(clades ==6){
                              no_contact_per_day_Re = 4.0
                            }else if(ReInput==7){
                              no_contact_per_day_Re = 5.0
                            }else if(ReInput==8){
                              no_contact_per_day_Re = 5.5
                            }
                            
                            no_contact_per_day=no_contact_per_day_Re ##
                            dispersion_options= c(40) ## value for rho
                            
                            no_contact_per_day=5.5 ## clade 8 default
                            original_strain_re=2.2  ## clade 8 default
                       
                            ## different contact rate at each time point
                            no_contact_per_day_final<-c()
                            no_contact_per_day_final=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[1]),
                                                            scale=dispersion_options[1]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
    
                            Prob_V=(pmax(0,(increase_VL*out$V)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(increase_VL*out$V)))^alpha)
                            
                 
                            time_V_final=t.out[is.finite(Prob_V)]
                            Prob_V_final=Prob_V[is.finite(Prob_V)]
                    
                            V_final=(increase_VL*out$V[is.finite(Prob_V)])
                            
                            
                            if(tt1==1){
                              no_contact_per_day_final_ref[c(seq(1,length(no_contact_per_day_final),by=1))]= round(no_contact_per_day_final)
                              
                              log_Vfinal= round(log10(V_final),digits=1)
                              log_Vfinal[is.na(log_Vfinal)]=0
                              log_Vfinal[log_Vfinal<0]=0
                              
                              V_ref[c(seq(1,length(log_Vfinal),by=1))] = log_Vfinal
                              
                              
                            }

                            ## check for successful transmissions (since transmission risk is not successful transmission)
                            random_numbers= runif(n=length(Prob_V_final),min=0,max=1)  # uniform distribution on the interval from min to max.
                            susc_prob=(total_pop-sum((II_total_COVID[c(seq(1,tt1,by=1))]))-sum(II_recover_COVID[c(seq(1,tt1,by=1))]))/total_pop
                            no_contact_per_day_final_person=c()
                            no_contact_per_day_final_person=no_contact_per_day_final[is.finite(Prob_V)]
                            
                            number_of_new_infections_ideal_case= dT*no_contact_per_day_final_person*Prob_V_final*susc_prob
                           
                            IND_Success=Prob_V_final>random_numbers
                            number_of_new_infections_actual= dT*no_contact_per_day_final_person[IND_Success]*Prob_V_final[IND_Success]*susc_prob
                            time_new_infections_actual=time_V_final[Prob_V_final>random_numbers]

                            if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){

                              total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                              ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                              for (ks in ks_options) {
                                if(ks==1){ 
                                  total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                                } else {
                                  total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                                }
                                
                              }
                              
                              total_number_of_new_infections=round(total_number_of_new_infections)

                              time_new_infections_actual=time_new_infections_actual[!duplicated(total_number_of_new_infections)]
                              
                              total_number_of_new_infections=total_number_of_new_infections[!duplicated(total_number_of_new_infections)]
  
                              if(any(total_number_of_new_infections == 0)){                    #### new
                                total_number_of_new_infections==0
                                ind=seq(1+length(total_number_of_new_infections[total_number_of_new_infections==0]),length(total_number_of_new_infections),by=1)
                              } else{
                                ind=seq(1,length(total_number_of_new_infections),by=1)
                              }
        
                              index_where_new_infections_happen=ind[!duplicated(total_number_of_new_infections[ind])]
      
                              time_when_new_infection_is_added=time_new_infections_actual[index_where_new_infections_happen]
                        
                              time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added-tzero  ################ this is imp  (add incubation period)
                  
                              # caluclating number of new infections at these time points
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
                              
                     
                              next_strain_mat <- c()
                              
                              ####### Add this data to II_total_COVID 
                              if (length(time_when_new_infection_is_added_relative_to_tzero)>0){ 
                                
                                fl_explore=seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1)
                                
                                    for (op in c(seq(1,length(time_when_new_infection_is_added_relative_to_tzero),by=1))){
                                          ## c(2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) for 1.0, 1.2 , 1.4, 1.6, 1.8, 2.0 and 2.2 Re values
                                          no_contact_per_day_options = c(1.1, 2.3, 3.1, 3.5, 3.75 , 4.0 , 5.0, 5.5) ## options for new varitnas
                                          
                                          ### 1/100 chance that an infected person will generate a new variant and it will be dominant for transmission 
                                          if(runif(1)>0.99){
                                                      x <- 1:length(no_contact_per_day_options)
                                                      x1 = sample(x)[1]
                                                      no_contact_per_day=no_contact_per_day_options[x1] ## new variant gets a new Re
                                                      new_strain =1
                                                      if(no_contact_per_day==1.1){                      
                                                        new_strain_re=0.8                     
                                                      }else if(no_contact_per_day==2.3){
                                                        new_strain_re=1.0
                                                      } else if(no_contact_per_day==3.1){
                                                        new_strain_re=1.2
                                                      }else if(no_contact_per_day==3.5){
                                                        new_strain_re=1.4
                                                      }else if(no_contact_per_day==3.75){
                                                        new_strain_re=1.6
                                                      }else if(no_contact_per_day==4.0){
                                                        new_strain_re=1.8
                                                      }else if(no_contact_per_day==5.0){
                                                        new_strain_re=2.0
                                                      } else {
                                                        new_strain_re=2.2
                                                      }
                                                      
                                            next_strain_mat[op] = new_strain_re
                                          } else {
                                                      new_strain =0
                                                    if(no_contact_per_day==1.1){                   
                                                      new_strain_re=0.8                     
                                                    }else if(no_contact_per_day==2.3){
                                                      new_strain_re=1.0
                                                    } else if(no_contact_per_day==3.1){
                                                      new_strain_re=1.2
                                                    }else if(no_contact_per_day==3.5){
                                                      new_strain_re=1.4
                                                    }else if(no_contact_per_day==3.75){
                                                      new_strain_re=1.6
                                                    }else if(no_contact_per_day==4.0){
                                                      new_strain_re=1.8
                                                    }else if(no_contact_per_day==5.0){
                                                      new_strain_re=2.0
                                                    } else {
                                                      new_strain_re=2.2
                                                    }
                                            next_strain_mat[op] = new_strain_re
                                          }
                                  }
                                
                                
                              }   else {
                                fl_explore=c()
                              }
                              

                              # First add first generation new infections
                              for (ijktt in fl_explore){
                                
                                
                                incubation_period_sample=rnorm(1, mean = 4, sd = 1)
                                incubation_period_sample=0 ## we do not account for incubation period and tzero
        
                                tr=time_when_new_infection_is_added_relative_to_tzero[ijktt]
                                tr_actual= 1 + tr + tt1 + incubation_period_sample # actual time at which new infections spread is relative to the current time (tt1)
                                ## 1 reflects time from infection to viral replication
                                Ir=incidence[ijktt] # these are the number of new infections caused by each infection moving forward in time
                              
                                      if(next_strain_mat[ijktt]==1){
                                        II_total_COVID_clade1[tr_actual]=II_total_COVID_clade1[tr_actual]+Ir
                                      } else if(next_strain_mat[ijktt]==1.2 ) {
                                        II_total_COVID_clade2[tr_actual]=II_total_COVID_clade2[tr_actual]+Ir
                                      } else if(next_strain_mat[ijktt]==1.4 ) {
                                        II_total_COVID_clade3[tr_actual]=II_total_COVID_clade3[tr_actual]+Ir
                                      } else if(next_strain_mat[ijktt]==1.6 ) {
                                        II_total_COVID_clade4[tr_actual]=II_total_COVID_clade4[tr_actual]+Ir
                                      } else if(next_strain_mat[ijktt]==1.8 ) {
                                        II_total_COVID_clade5[tr_actual]=II_total_COVID_clade5[tr_actual]+Ir
                                      } else if(next_strain_mat[ijktt]==2.0 ) {
                                        II_total_COVID_clade6[tr_actual]=II_total_COVID_clade6[tr_actual]+Ir
                                      } else if(next_strain_mat[ijktt]==2.2 ) {
                                       
                                        II_total_COVID_clade7[tr_actual]=II_total_COVID_clade7[tr_actual]+Ir
                                      }
                            
                                
                                II_total_COVID_clade1 <- II_total_COVID_clade1[!is.na(II_total_COVID_clade1)]
                                II_total_COVID_clade2 <- II_total_COVID_clade2[!is.na(II_total_COVID_clade2)]
                                II_total_COVID_clade3 <- II_total_COVID_clade3[!is.na(II_total_COVID_clade3)]
                                II_total_COVID_clade4 <- II_total_COVID_clade4[!is.na(II_total_COVID_clade4)]
                                II_total_COVID_clade5 <- II_total_COVID_clade5[!is.na(II_total_COVID_clade5)]
                                II_total_COVID_clade6 <- II_total_COVID_clade6[!is.na(II_total_COVID_clade6)]
                                II_total_COVID_clade7 <- II_total_COVID_clade7[!is.na(II_total_COVID_clade7)]
                                II_total_COVID_clade8 <- II_total_COVID_clade8[!is.na(II_total_COVID_clade8)]
                                
                                II_total_COVID= II_total_COVID_clade1 + II_total_COVID_clade2+
                                  II_total_COVID_clade3 + II_total_COVID_clade4+
                                  II_total_COVID_clade5 + II_total_COVID_clade6+
                                  II_total_COVID_clade7+ II_total_COVID_clade8
  
                                
                              }
                            
                            } 

                            
                            person_ID = person_ID+1 
                            
                          }  ### simulating epidemics for each active infetion
    
                        }  ### running fro all active infections at a particular time
                        
        }  ## for clade 8
        
     
   }  #### for all clades
    
    if(sum(II_total_COVID)>=infection_threshold){  ### if 10^5 number of cumulative infections occured then no need to run more
      break
    }

  }  ##### for all time points
  

  Cumulative_infections[tt1]=sum(II_total_COVID)
  
  
  II_metrics[seq(1+no_days*(instance_index-1),no_days*instance_index,by=1),]=cbind(tt_total_COVID[seq(1,no_days,by=1)],
                                                                                   II_total_COVID[seq(1,no_days,by=1)],
                                                                                   rep(instance_index,no_days),
                                                                                   II_total_COVID_clade1,
                                                                                   II_total_COVID_clade2,
                                                                                   II_total_COVID_clade3,
                                                                                   II_total_COVID_clade4,
                                                                                   II_total_COVID_clade5,
                                                                                   II_total_COVID_clade6 ,
                                                                                   II_total_COVID_clade7,
                                                                                   II_total_COVID_clade8)
  
  write.csv(II_metrics,"Variants_overtake_Re_general.csv",row.names=FALSE) #
  

  instance_index=instance_index+1

}       













#######################################################
############## Plotting figure
#######################################################





library(deSolve)
library(DEoptim)
library(pracma)
library(Hmisc)
library(scales)
library(ggplot2)
library(egg)
library(pracma)
library(RColorBrewer)



infection_threshold=100000
no_runs=1000
single_plot=0
breaks_dist = c(seq(0,200,by=1))
breaks_dist2 = c(seq(0,9,by=1))
##################################################################
##################### Re=0.8
##################################################################

nrows = 3
ncolumns =3
par(mfrow = c(nrows,ncolumns))
par(mfrow = c(nrows,ncolumns),     
    oma = c(0, 0, 0, 0), # 
    mar = c(3.5, 3.5, 0.0, 0.0), 
    mgp = c(2.5, 1, 0),    
    xpd = FALSE)



jBrewColors <- brewer.pal(n = 8, name = "Spectral")
#display.brewer.pal(n = 8, name = 'Spectral')

A_inst = read.csv("Variants_overtake_Re_general.csv",header=TRUE)  #


S2_taking_S1<- c()
S3_taking_S1<- c()
S4_taking_S1<- c()
S5_taking_S1<- c()
S6_taking_S1<- c()
S7_taking_S1<- c()
S8_taking_S1<- c()

Dominant_strain<- c()


ijk=1
ijxx=1
counter1=1

counter_burout=0

while(ijxx<101){  # ijk<100
  
  A_traj<- c()
  ind_traj<- c()
  
  ind_traj = which(A_inst$Instance==ijxx)
  
  A_traj = A_inst[ind_traj,]
  
  indexA=length(A_traj$time)
  cumulative_incidence_clade1 <- vector()
  cumulative_incidence_clade2 <- vector()
  cumulative_incidence_clade3 <- vector()
  cumulative_incidence_clade4 <- vector()
  cumulative_incidence_clade5 <- vector()
  cumulative_incidence_clade6 <- vector()
  cumulative_incidence_clade7 <- vector()
  cumulative_incidence_clade8 <- vector()
  
  for (iA in c(seq(1,indexA,by=1))) { 
    #print(iA)
    if(iA==1){ 
      cumulative_incidence_clade1[iA]=A_traj$Incidence_clade1[iA]
      cumulative_incidence_clade2[iA]=A_traj$Incidence_clade2[iA]
      cumulative_incidence_clade3[iA]=A_traj$Incidence_clade3[iA]
      cumulative_incidence_clade4[iA]=A_traj$Incidence_clade4[iA]
      cumulative_incidence_clade5[iA]=A_traj$Incidence_clade5[iA]
      cumulative_incidence_clade6[iA]=A_traj$Incidence_clade6[iA]
      cumulative_incidence_clade7[iA]=A_traj$Incidence_clade7[iA]
      cumulative_incidence_clade8[iA]=A_traj$Incidence_clade8[iA]
    } else {
      cumulative_incidence_clade1[iA]=A_traj$Incidence_clade1[iA]+cumulative_incidence_clade1[iA-1] 
      cumulative_incidence_clade2[iA]=A_traj$Incidence_clade2[iA]+cumulative_incidence_clade2[iA-1]
      cumulative_incidence_clade3[iA]=A_traj$Incidence_clade3[iA]+cumulative_incidence_clade3[iA-1]
      cumulative_incidence_clade4[iA]=A_traj$Incidence_clade4[iA]+cumulative_incidence_clade4[iA-1]
      cumulative_incidence_clade5[iA]=A_traj$Incidence_clade5[iA]+cumulative_incidence_clade5[iA-1]
      cumulative_incidence_clade6[iA]=A_traj$Incidence_clade6[iA]+cumulative_incidence_clade6[iA-1]
      cumulative_incidence_clade7[iA]=A_traj$Incidence_clade7[iA]+cumulative_incidence_clade7[iA-1]
      cumulative_incidence_clade8[iA]=A_traj$Incidence_clade8[iA]+cumulative_incidence_clade8[iA-1]
    }
  }
  
  Incidence_max = which.max(c(max(A_traj$Incidence_clade1),max(A_traj$Incidence_clade2),max(A_traj$Incidence_clade3),
                              max(A_traj$Incidence_clade4),max(A_traj$Incidence_clade5),max(A_traj$Incidence_clade6),
                              max(A_traj$Incidence_clade7),max(A_traj$Incidence_clade8)))
  
  Incidence_max_number=max(c(max(A_traj$Incidence_clade1),max(A_traj$Incidence_clade2),max(A_traj$Incidence_clade3),
                             max(A_traj$Incidence_clade4),max(A_traj$Incidence_clade5),max(A_traj$Incidence_clade6),
                             max(A_traj$Incidence_clade7),max(A_traj$Incidence_clade8)))
  
  
  if(Incidence_max==1){
    strain_last_time=  A_traj$time[A_traj$Incidence_clade1 == max(A_traj$Incidence_clade1)]
  } else if (Incidence_max==2){
    strain_last_time=  A_traj$time[A_traj$Incidence_clade2 == max(A_traj$Incidence_clade2)]
  } else if (Incidence_max==3){
    strain_last_time=  A_traj$time[A_traj$Incidence_clade3 == max(A_traj$Incidence_clade3)]
  } else if (Incidence_max==4){
    strain_last_time=  A_traj$time[A_traj$Incidence_clade4 == max(A_traj$Incidence_clade4)]
  } else if (Incidence_max==5){
    strain_last_time=  A_traj$time[A_traj$Incidence_clade5 == max(A_traj$Incidence_clade5)]
  } else if (Incidence_max==6){
    strain_last_time=  A_traj$time[A_traj$Incidence_clade6 == max(A_traj$Incidence_clade6)]
  } else if (Incidence_max==7){
    strain_last_time=  A_traj$time[A_traj$Incidence_clade7 == max(A_traj$Incidence_clade7)]
  } else if (Incidence_max==8){
    strain_last_time=  A_traj$time[A_traj$Incidence_clade8 == max(A_traj$Incidence_clade8)]
  }
  
  print("strain_last_time")
  print(strain_last_time)
  
  if(Incidence_max_number==1000 & Incidence_max==1){ # means t=0 was maximum and the infection did not take off
    Dominant_strain[ijxx] = 0
  } else{
    Dominant_strain[ijxx] = Incidence_max
  }
  
  
  
  if(Incidence_max_number==1000){
    counter_burout=counter_burout+1
  }
  
  
  cumulative_incidence = cumulative_incidence_clade1 + cumulative_incidence_clade2 + cumulative_incidence_clade3 +
    cumulative_incidence_clade4 + cumulative_incidence_clade5+
    cumulative_incidence_clade6+ cumulative_incidence_clade7 + cumulative_incidence_clade8
  
  time_to_10000_infections_temp = A_traj$time[cumulative_incidence>infection_threshold]
  if(isempty(time_to_10000_infections_temp)){
    time_to_10000_infections_temp1=0
  } else {
    time_to_10000_infections_temp1=time_to_10000_infections_temp[1] 
  }
  
  time_to_10000_infections_temp1[1]=strain_last_time[1]
  
  print("time_to_10000_infections_temp1")
  print(time_to_10000_infections_temp1)
  
  S2_taking_S1_temp1 = A_traj$time[cumulative_incidence_clade2>cumulative_incidence_clade1]
  
  if(isempty(S2_taking_S1_temp1)){
    S2_taking_S1[ijk]=0
    #  print(S2_taking_S1)
  } else {
    S2_taking_S1[ijk]=S2_taking_S1_temp1[1]
  }
  
  S3_taking_S1_temp1 = A_traj$time[cumulative_incidence_clade3>cumulative_incidence_clade1]
  
  if(isempty(S3_taking_S1_temp1)){
    S3_taking_S1[ijk]=0
  } else {
    S3_taking_S1[ijk]=S3_taking_S1_temp1[1]
  }
  
  S4_taking_S1_temp1 = A_traj$time[cumulative_incidence_clade4>cumulative_incidence_clade1]
  
  if(isempty(S4_taking_S1_temp1)){
    S4_taking_S1[ijk]=0
  } else {
    S4_taking_S1[ijk]=S4_taking_S1_temp1[1]
  }
  
  S5_taking_S1_temp1 = A_traj$time[cumulative_incidence_clade5>cumulative_incidence_clade1]
  
  if(isempty(S5_taking_S1_temp1)){
    S5_taking_S1[ijk]=0
  } else {
    S5_taking_S1[ijk]=S5_taking_S1_temp1[1]
  }
  
  
  S6_taking_S1_temp1 = A_traj$time[cumulative_incidence_clade6>cumulative_incidence_clade1]
  
  if(isempty(S6_taking_S1_temp1)){
    S6_taking_S1[ijk]=0
  } else {
    S6_taking_S1[ijk]=S6_taking_S1_temp1[1]
  }
  
  
  S7_taking_S1_temp1 = A_traj$time[cumulative_incidence_clade7>cumulative_incidence_clade1]
  
  if(isempty(S7_taking_S1_temp1)){
    S7_taking_S1[ijk]=0
  } else {
    S7_taking_S1[ijk]=S7_taking_S1_temp1[1]
  }
  
  S8_taking_S1_temp1 = A_traj$time[cumulative_incidence_clade8>cumulative_incidence_clade1]
  
  if(isempty(S8_taking_S1_temp1)){
    S8_taking_S1[ijk]=0
  } else {
    S8_taking_S1[ijk]=S8_taking_S1_temp1[1]
  }
  
  
  colr1 = "darkgoldenrod4"
  #colr1= jBrewColors[1]
  colr2 = "black"
  #colr2= jBrewColors[2]
  colr3 = "blue"
  #colr3= jBrewColors[3]
  colr4 = "Dark green"
  #colr4= jBrewColors[4]
  colr5 = "red"
  #colr5= jBrewColors[5]
  colr6 = "orange"
  #colr6= jBrewColors[6]
  colr7 = "sky blue"
  #colr7= jBrewColors[7]
  colr8 = "Violet"
  
  
  time_to_10000_infections_temp1 = 100
  
  plot(A_traj$time[c(seq(1,time_to_10000_infections_temp1,by=1))],
       A_traj$Incidence_clade1[c(seq(1,time_to_10000_infections_temp1,by=1))] , 
       col = colr1 , lwd=1.5 , lty=1 ,
       xlim=c(0,50),ylim=c(1,10000),  bty="n", type="l", log="y", yaxt="n" , #xaxt="n" , 
       cex.main=1,cex.axis=1,cex.lab=1,
       xlab="Time (days)",
       ylab="Incidence")
  
  lines(A_traj$time[c(seq(1,time_to_10000_infections_temp1,by=1))],
        A_traj$Incidence_clade2[c(seq(1,time_to_10000_infections_temp1,by=1))] , 
        col = colr2 , lwd=1.5 , lty=1)
  lines(A_traj$time[c(seq(1,time_to_10000_infections_temp1,by=1))],
        A_traj$Incidence_clade3[c(seq(1,time_to_10000_infections_temp1,by=1))] ,
        col = colr3 , lwd=1.5 , lty=1)
  lines(A_traj$time[c(seq(1,time_to_10000_infections_temp1,by=1))],
        A_traj$Incidence_clade4[c(seq(1,time_to_10000_infections_temp1,by=1))] ,
        col = colr4 , lwd=1.5 , lty=1)
  lines(A_traj$time[c(seq(1,time_to_10000_infections_temp1,by=1))],
        A_traj$Incidence_clade5[c(seq(1,time_to_10000_infections_temp1,by=1))] , 
        col = colr5 , lwd=1.5 , lty=1)
  lines(A_traj$time[c(seq(1,time_to_10000_infections_temp1,by=1))],
        A_traj$Incidence_clade6[c(seq(1,time_to_10000_infections_temp1,by=1))] , 
        col = colr6 , lwd=1.5 , lty=1)
  lines(A_traj$time[c(seq(1,time_to_10000_infections_temp1,by=1))],
        A_traj$Incidence_clade7[c(seq(1,time_to_10000_infections_temp1,by=1))] , 
        col = colr7 , lwd=1.5 , lty=1)
  lines(A_traj$time[c(seq(1,time_to_10000_infections_temp1,by=1))],
        A_traj$Incidence_clade8[c(seq(1,time_to_10000_infections_temp1,by=1))] , 
        col = colr8 , lwd=1.5 , lty=1)
  
  legend(x=30,y=10000,c("Re=0.8","Re=1.0","Re=1.2","Re=1.4","Re=1.6","Re=1.8","Re=2.0","Re=2.2"),
         col=c("darkgoldenrod4","black","blue","Dark green","red","orange","sky blue","violet"),
         lty=c(1,1,1,1,1,1,1,1),bty = "n",cex = 0.8)
  
  axis(2,at=10^c(seq(0,4,by=1)),labels=expression("10"^"0","10"^"1","10"^"2","10"^"3","10"^"4"),
       cex.axis=0.8,font=1,lwd=1.5,las=1)#,las=1)
  
  ijxx=ijxx+1    
  
  
}

nrows = 1
ncolumns =1
par(mfrow = c(nrows,ncolumns))
par(mfrow = c(nrows,ncolumns),     
    oma = c(0, 0, 0, 0), # 
    mar = c(3.5, 3.5, 0.0, 0.0), 
    mgp = c(2.5, 1, 0),    
    xpd = FALSE)

plot.new()

nrows = 3
ncolumns =3
par(mfrow = c(nrows,ncolumns))
par(mfrow = c(nrows,ncolumns),     
    oma = c(0, 0, 0, 0), # 
    mar = c(3.5, 3.5, 0.0, 0.0), 
    mgp = c(2.5, 1, 0),    
    xpd = FALSE)

hist(S2_taking_S1, #
     freq= FALSE, # breaks=breaks) #,
     xlab = "Time of strain 2 cumulative infections  takeover",
     ylab = "Frequency",
     main="",
     #main=paste("theta",  theta),
     xlim=c(0,150),
     ylim=c(0,1),
     #freq=TRUE,
     #ylim = c(0,0.5),
     breaks=breaks_dist,
     col=scales::alpha('black',.5), # Filling Color
     border = "black",
     right = FALSE,
     xaxt="n",yaxt="n")

axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
axis(2,at=c(seq(0,1,by=0.2)),labels=c(seq(0,1,by=0.2)),
     cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)

text(73,0.8,as.expression(bquote(paste("Takeover (%)")==.(round(0.01*length(S2_taking_S1[S2_taking_S1>0]),digits=2)))),
     cex=1.0,col="red")


rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
rug(x = c(seq(0,1,by=0.2))+0.1, ticksize = -0.04, side = 2)



hist(S3_taking_S1, #
     freq= FALSE, # breaks=breaks) #,
     xlab = "Time of strain 3 cumulative infections  takeover",
     ylab = "Frequency",
     main="",
     #main=paste("theta",  theta),
     xlim=c(0,150),
     ylim=c(0,1),
     #freq=TRUE,
     #ylim = c(0,0.5),
     breaks=breaks_dist,
     col=scales::alpha('black',.5), # Filling Color
     border = "black",
     right = FALSE,
     xaxt="n",yaxt="n")

axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
axis(2,at=c(seq(0,1,by=0.2)),labels=c(seq(0,1,by=0.2)),
     cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)

text(73,0.8,as.expression(bquote(paste("Takeover (%)")==.(round(0.01*length(S3_taking_S1[S3_taking_S1>0]),digits=2)))),
     cex=1.0,col="red")


rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
rug(x = c(seq(0,1,by=0.2))+0.1, ticksize = -0.04, side = 2)


hist(S4_taking_S1, #
     freq= FALSE, # breaks=breaks) #,
     xlab = "Time of strain 4 cumulative infections takeover",
     ylab = "Frequency",
     main="",
     #main=paste("theta",  theta),
     xlim=c(0,150),
     ylim=c(0,1),
     #freq=TRUE,
     #ylim = c(0,0.5),
     breaks=breaks_dist,
     col=scales::alpha('black',.5), # Filling Color
     border = "black",
     right = FALSE,
     xaxt="n",yaxt="n")

axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
axis(2,at=c(seq(0,1,by=0.2)),labels=c(seq(0,1,by=0.2)),
     cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)

text(73,0.8,as.expression(bquote(paste("Takeover (%)")==.(round(0.01*length(S4_taking_S1[S4_taking_S1>0]),digits=2)))),
     cex=1.0,col="red")


rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
rug(x = c(seq(0,1,by=0.2))+0.1, ticksize = -0.04, side = 2)


hist(S5_taking_S1, #
     freq= FALSE, # breaks=breaks) #,
     xlab = "Time of strain 5 cumulative infections takeover",
     ylab = "Frequency",
     main="",
     #main=paste("theta",  theta),
     xlim=c(0,150),
     ylim=c(0,1),
     #freq=TRUE,
     #ylim = c(0,0.5),
     breaks=breaks_dist,
     col=scales::alpha('black',.5), # Filling Color
     border = "black",
     right = FALSE,
     xaxt="n",yaxt="n")

axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
axis(2,at=c(seq(0,1,by=0.2)),labels=c(seq(0,1,by=0.2)),
     cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)

text(73,0.8,as.expression(bquote(paste("Takeover (%)")==.(round(0.01*length(S5_taking_S1[S5_taking_S1>0]),digits=2)))),
     cex=1.0,col="red")


rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
rug(x = c(seq(0,1,by=0.2))+0.1, ticksize = -0.04, side = 2)


hist(S6_taking_S1, #
     freq= FALSE, # breaks=breaks) #,
     xlab = "Time of strain 6 cumulative infections takeover",
     ylab = "Frequency",
     main="",
     #main=paste("theta",  theta),
     xlim=c(0,150),
     ylim=c(0,1),
     #freq=TRUE,
     #ylim = c(0,0.5),
     breaks=breaks_dist,
     col=scales::alpha('black',.5), # Filling Color
     border = "black",
     right = FALSE,
     xaxt="n",yaxt="n")

axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
axis(2,at=c(seq(0,1,by=0.2)),labels=c(seq(0,1,by=0.2)),
     cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)

text(73,0.8,as.expression(bquote(paste("Takeover (%)")==.(round(0.01*length(S6_taking_S1[S6_taking_S1>0]),digits=2)))),
     cex=1.0,col="red")


rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
rug(x = c(seq(0,1,by=0.2))+0.1, ticksize = -0.04, side = 2)


hist(S7_taking_S1, #
     freq= FALSE, # breaks=breaks) #,
     xlab = "Time of strain 7 cumulative infections takeover",
     ylab = "Frequency",
     main="",
     #main=paste("theta",  theta),
     xlim=c(0,150),
     ylim=c(0,1),
     #freq=TRUE,
     #ylim = c(0,0.5),
     breaks=breaks_dist,
     col=scales::alpha('black',.5), # Filling Color
     border = "black",
     right = FALSE,
     xaxt="n",yaxt="n")

axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
axis(2,at=c(seq(0,1,by=0.2)),labels=c(seq(0,1,by=0.2)),
     cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)

text(73,0.8,as.expression(bquote(paste("Takeover (%)")==.(round(0.01*length(S7_taking_S1[S7_taking_S1>0]),digits=2)))),
     cex=1.0,col="red")


rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
rug(x = c(seq(0,1,by=0.2))+0.1, ticksize = -0.04, side = 2)


hist(S8_taking_S1, #
     freq= FALSE, # breaks=breaks) #,
     xlab = "Time of strain 8 cumulative infections takeover",
     ylab = "Frequency",
     main="",
     #main=paste("theta",  theta),
     xlim=c(0,150),
     ylim=c(0,1),
     #freq=TRUE,
     #ylim = c(0,0.5),
     breaks=breaks_dist,
     col=scales::alpha('black',.5), # Filling Color
     border = "black",
     right = FALSE,
     xaxt="n",yaxt="n")

axis(1,at=c(seq(0,150,by=25)),labels=c(seq(0,150,by=25)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
axis(2,at=c(seq(0,1,by=0.2)),labels=c(seq(0,1,by=0.2)),
     cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)

text(73,0.8,as.expression(bquote(paste("Takeover (%)")==.(round(0.01*length(S8_taking_S1[S8_taking_S1>0]),digits=2)))),
     cex=1.0,col="red")


rug(x = c(seq(0,150,by=25))+5, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+10, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+15, ticksize = -0.02, side = 1)
rug(x = c(seq(0,150,by=25))+20, ticksize = -0.02, side = 1)
rug(x = c(seq(0,1,by=0.2))+0.1, ticksize = -0.04, side = 2)



hist(Dominant_strain, #
     freq= FALSE, # breaks=breaks) #,
     xlab = "Dominant strain",
     ylab = "Frequency",
     main="",
     #main=paste("theta",  theta),
     xlim=c(0,9),
     ylim=c(0,1),
     #freq=TRUE,
     #ylim = c(0,0.5),
     breaks=breaks_dist2,
     col=scales::alpha('black',.5), # Filling Color
     border = "black",
     right = FALSE,
     xaxt="n",
     yaxt="n"
)

axis(1,at=c(seq(0,9,by=1))+0.5,labels=c(seq(0,9,by=1)),cex.axis=1.0,font=1,lwd=1,las=2)#,las=1)
axis(2,at=c(seq(0,1,by=0.2)),labels=c(seq(0,1,by=0.2)),
     cex.axis=1.0,font=1,lwd=1,las=1)#,las=1)

text(4, 0.7, "Dominant Strain corresponds to burnouts", cex=0.8)

rug(x = c(seq(0,1,by=0.2))+0.1, ticksize = -0.04, side = 2)



print(length(Dominant_strain[Dominant_strain>1]))



print(100*length(Dominant_strain[Dominant_strain>1])/(100 - counter_burout))





