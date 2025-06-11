# Load deSolve library
library(deSolve)
library(ggpubr)


library(deSolve)
library(readxl)

# Define the parameters
parameters <- c(a = 0.025,   # birth  rate  (1.34/year)
                b = 0.016,  # death  rate (0.836yyear)
                gamma = 0.00076,   # density-dependence in mortality (0.0397km2/yr)
                sigma = 0.25, #incubation  period  (1/sigma = 0.13  year)
                alpha = 0.5, #rabies-induced mortality (66.36/year)
                beta = 0.021, #transmission rate (33.25 year)
                rho = 1- 0.20)  # proportion  of  raccoons  that  develop  natural  immunity(50.20)


parameters <- c(a = 0.013,   # birth  rate  (1.34/year)
                b = 0.0076,  # death  rate (0.836yyear)
                gamma = 0.00007,   # density-dependence in mortality (0.0397km2/yr)
                sigma = 0.24, #incubation  period  (1/sigma = 0.13  year)
                alpha = 0.5, # rabies-inducedmortality (66.36/year)
                beta = 0.031, #transmission rate (33.25 year)
                rho = 1 - 0.)  # proportion  of  raccoons  that  develop  natural  immunity(50.20)

# Define the function for the differential equations
lotka_volterra <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    
    dX <- (a * (X + I)) - ((b + (beta * Y) + (gamma * N)) * X)
      
    dH1 <- (rho * beta * X * Y) - (b + sigma + (gamma * N)) * H1
    
    dH2  <- ((1- rho) * beta * X * Y) - (b + sigma + (gamma * N)) * H2
    
    dY <- sigma*H1 - (b + alpha + (gamma * N)) * Y
    
    dI <- sigma*H2 - (b + (gamma * N)) * I
    
    N = dX + dH1 + dH2 + dY + dI
                                 
    
    return(list(c(dX, dH1, dH2, dY, dI, N)))
  })
}






# Set initial conditions
initial_state <- c(X = 64, H1 = 0, H2 = 0, Y = 0.008, I = 0, N = 64+0.008)  # initial number of rabbits and foxes

# Set time points for which to solve the ODEs
times <- seq(0, 2000, by = 1)  # from 0 to 200, in increments of 0.1

# Solve the ODEs
output <- ode(y = initial_state, times = times, func = lotka_volterra, parms = parameters)






plot(times,output[,4]*162,type="l",col="blue",ylab="Proportion")

# lines(times,output[,3],col="orange")
# lines(times,output[,4],col="red")  
# lines(times,output[,5],col="green")
# legend(300,0.7,legend=c("S","E","I","R"),col=c("blue","orange","red","green"), lty=1, cex=0.8)


sum(output[,4][1:100]*162)



#####################
############################


library(deSolve)




parameters <- c(a = 0.013,   # birth  rate  (1.34/year)
                b = 0.0076,  # death  rate (0.836yyear)
                gamma = 0.00007,   # density-dependence in mortality (0.0397km2/yr)
                sigma = 0.24, #incubation  period  (1/sigma = 0.13  year)
                alpha = 0.5, # rabies-induced mortality (66.36/year)
                beta = 0.031, #transmission rate (33.25 year)
                rho = 1 - 0.20, # proportion  of  raccoons  that  develop  natural  immunity(50.20)
                b_h	<-  0.000326923,  #(Human_birth/52)/1000, #Human birth rate	0.000326923
                lambda_h <- 0, 	#Human loss of vaccination immunity rate	0
                m_h	<- (1/72)/52, #Human mortality rate	0.000274725
                v_h	<- 0.969, # (REQUIRES INPUT)Human vaccine efficacy	0.969
                alpha_h	<- 0.0, #Human prophylactic rate	0
                beta_dh <- 0.0000156, # (REQUIRES INPUT) Dog human transmission rate	0.0000156
                P10 <- 0.95, 	#(REQUIRES INPUT)#PEP vaccination rate	25%
                
                mu_h <- 	(1/10)*7, # (REQUIRES INPUT)Inverse of average infective period, rabid human mortality rate	1.0000
                
                p_ExptoInf <-	0.015,   #Rates of developing rabies from exposure (bites), per week
                p_ExptoNoInf <- 0.101) #Rates of not developing rabies from exposure (bites), per week)  

# parameters <- c(a = 0.026,   # birth  rate  (1.34/year) Childs et al 2000
#                 b = 0.016,  # death  rate (0.836yyear) Childs et al 2000
#                 gamma = 0.00076,   # density-dependence in mortality (0.0397km2/yr) Childs et al 2000
#                 sigma = 0.16, #inverse incubation  period  (1/sigma = 0.16) Rabies econ,  0.13  Childs et al 2000
#                 alpha = 1.2, # rabies-induced mortality (66.36/year) Childs et al 2000
#                 beta = 0.64, #0.64, #transmission rate (33.25 year) 
#                 rho = (0.98), # 1- proportion  of  raccoons  that  develop  natural  immunity 
#                 
#                 ### Raccoon to human transmission
#                 
#                 b_h	<-  0.26, #0.000326923,  #(Human_birth/52)/1000, #Human birth rate	0.000326923
#                 lambda_h <- 0, 	#Human loss of vaccination immunity rate	0
#                 m_h	<- (1/72)/52, #Human mortality rate	0.000274725
#                 v_h	<- 0.969, # (REQUIRES INPUT)Human vaccine efficacy	0.969
#                 alpha_h	<- 0.0, #Human prophylactic rate	0
#                 beta_dh <- 0.00034, # (REQUIRES INPUT) Dog human transmission rate	0.0000156
#                 P10 <- 0.80, 	#(REQUIRES INPUT)#PEP vaccination rate	25%
#                 
#                 mu_h <- 	(1/10)*7, # (REQUIRES INPUT)Inverse of average infective period, rabid human mortality rate	1.0000
#                 
#                 p_ExptoInf <-	0.025, #0.025,   #Rates of developing rabies from exposure (bites), per week
#                 p_ExptoNoInf <- 0.097) #0.097) #Rates of not developing rabies from exposure (bites), per week)  
# 

#R0*((sigma+b)*(b+0.75))/(sigma*0.45*dX)

# 1.5*((0.13+0.016)*(0.016+0.75))/(0.16*0.45*64)
# 
# 
# # Define the function for the differential equations
# lotka_volterra <- function(time, state, parameters) {
#   with(as.list(c(state, parameters)), {
#     
#     
#     dX <- (a * (X + I)) - ((b + (beta * Y) + (gamma * N)) * X)
#     
#     dH1 <- (rho * beta * X * Y) - (b + sigma + (gamma * N)) * H1
#     
#     dH2  <- ((1- rho) * beta * X * Y) - (b + sigma + (gamma * N)) * H2
#     
#     dY <- sigma*H1 - (b + alpha + (gamma * N)) * Y
#     
#     dI <- sigma*H2 - (b + (gamma * N)) * I
#     
#     N = dX + dH1 + dH2 + dY + dI
#     
#     dSh <- (b_h*(Sh+Eh+Rh))+(lambda_h*Rh) -(m_h*Sh)-(v_h*alpha_h*Sh)-(beta_dh*Sh*dH2) +(Eh*p_ExptoNoInf)# change Rh 
#     
#     dEh <- (beta_dh*Sh*dH2) -(m_h*Eh) -(Eh*p_ExptoInf*P10*v_h)-(Eh*p_ExptoInf*(1-P10*v_h)) -(Eh*p_ExptoNoInf)
#     
#     dExp <- (beta_dh*Sh*dH2) 
#     
#     dIh <- (Eh*p_ExptoInf*(1-P10*v_h))-(m_h*Ih)-(mu_h*Ih)
#     
#     dDh <- (Eh*p_ExptoInf*(1-P10*v_h)) 
#     
#     dRh <- (Eh*p_ExptoInf*P10*v_h)+(v_h*alpha_h*Sh)-(m_h*Rh)-(lambda_h*Rh)
#     
#     dNh = dSh+dEh+dIh+dDh+dRh
#     
#     
#     return(list(c(dX, dH1, dH2, dY, dI, N, dSh, dEh, dExp, dIh, dDh, dRh, dNh)))
#   })
# }
# 
# #=R12+(βdh_DogHumanTransmission*Q12*J12)-(mh_HumanDeath*R12)-(R12*p_ExptoInf*P10_PEPVaccination*vh_HumanVaccineEfficacy)-(R12*p_ExptoInf*(1-P10_PEPVaccination*vh_HumanVaccineEfficacy))-(R12*p_ExptoNoInf)
# 
# # Set initial conditions
# initial_state <- c(X = 10, H1 = 0, H2 = 0, Y = 0.008, I = 0, N = 30+0.008,
#                    Sh = 750, Eh = 0, Exp = 0,   Ih = 0, Dh = 0, Rh = 0, Nh = 750)  # initial number of rabbits and foxes
# 
# # Set time points for which to solve the ODEs
# times <- seq(0, 4000, by = 1)  # from 0 to 200, in increments of 0.1
# 
# # Solve the ODEs
# output <- ode(y = initial_state, times = times, func = lotka_volterra, parms = parameters)
# 
# output
#   
# 
# 
# #sum(output[,"Dh"][1:200]*162)
# 
# sum(output[,"Y"][1:100]*162)
# sum(output[,"Exp"][1:100]*162)
# #sum(output[,"Eh"][1:250]*162)
# 
# 
# rac_deaths <- as.data.frame(output) %>%
#               ggplot(aes(x=times, y=Y * 162)) +
#               geom_line() + xlim(1, 4000) + ylim(0, 25)
# 
# 
# hum_exp   <- as.data.frame(output) %>%
#   ggplot(aes(x=times, y=cumsum(Exp* 162))) +
#   geom_line() + xlim(1, 4000)
# 
# 
# ggarrange(rac_deaths, hum_exp, ncol = 2, labels = c("Raccoon deaths", "Cumulative human exposures"))



#########################

Ss <- 60

Sh <- 100

ddp <- (0.026 - 0.0076)/Ss*(1+1/log(Ss*500))*1.05



1.3*(0.16+0.0076)*(0.0076+1)/(0.16*0.45*50)




parameters <- c(#a = 0.026,  #0.11 # birth  rate  (1.34/year) Childs et al 2000
                #a_2 = 0,
                b = 0.0076,  # death  rate: Rabies Econ = 0.0076,  Childs et al 2000 = 0.016/week, Ruan et al = 0.026/week 
                gamma = ddp,   # density-dependence in mortality (0.0397km2/yr) Childs et al 2000
                sigma = 0.16, #inverse incubation  period  (1/sigma = 0.16) Rabies econ,  0.13  Childs et al 2000
                alpha = 1, # rabies-induced mortality (66.36/year) Childs et al 2000
                beta =  1.4*(0.16+0.0076)*(0.0076+1)/(0.16*0.45*Ss), # #transmission rate Ruan et al = 0.03, Rabies econ = 0.03, Childs et al = 33.25/ year 
                rho = (0.98), # 1- proportion  of  raccoon  that  develop  natural  immunity 
                
                ### Raccoon to human transmission
                
                b_h	<-  0.000326923,  #(Human_birth/52)/1000, #Human birth rate	0.000326923
                lambda_h <- 0, 	#Human loss of vaccination immunity rate	0
                m_h	<- (1/78)/52, #Human mortality rate	0.000274725
                v_h	<- 0.969, # (REQUIRES INPUT)Human vaccine efficacy	0.969
                alpha_h	<- 0.0, #Human prophylactic rate	0
                beta_dh <- 3E-05,  #Calculated for Omaha: 33/52/11250 kmsq/raccoon week
                P10 <- 0.95, 	#(REQUIRES INPUT)#PEP vaccination rate	25%
                
                mu_h <- 	1, # (REQUIRES INPUT)Inverse of average infective period, rabid human mortality rate	1.0000
                
                p_ExptoInf <-	0.015, #0.025,   #Rates of developing rabies from exposure (bites), per week
                p_ExptoNoInf <- 0.101) #0.097) #Rates of not developing rabies from exposure (bites), per week)  

1.5*((0.16+0.016)*(0.016+1.2))/(0.16*64)

 #(number of rabid raccoon contacting humans/week/km sq)


lotka_volterra <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    if (time %% 53 <= 12) #time %% 53 >= 12  & time %% 53 < 24  #time %% 53 <= 12
      a <- 0.12
    else
      a <- 0
    
    
    
    
    
    #dBeta <- 1.5*((sigma+b)*(b+alpha))/(sigma*X)
    
    #beta <- dBeta
    
   
    dX <- (a * (X + I)) - ((b + (beta * Y) + (gamma * N)) * X)  ##susceptible raccoonhosts
    dH1 <- (rho * beta * X * Y) - (b + sigma + (gamma * N)) * H1  #exposed hosts (i.e., infected but not infectious; H1)
    
    dH2  <- ((1- rho) * beta * X * Y) - (b + sigma + (gamma * N)) * H2  #hosts  exposed  that  eventually  develop  immunity
    
    dY <- sigma*H1 - (b + alpha + (gamma * N)) * Y  #rabid hosts
    
    dI <- sigma*H2 - (b + (gamma * N)) * I  #hosts that are immune
    
    dN = dX + dH1 + dH2 + dY + dI
    
    dSh <- (b_h*(Sh+Eh+Rh))+(lambda_h*Rh) -(m_h*Sh)-(v_h*alpha_h*Sh)-(beta_dh*Sh*Y) +(Eh*p_ExptoNoInf)  ##Susceptible human population 
    
    dEh <- (beta_dh*Sh*Y) -(m_h*Eh) -(Eh*p_ExptoInf*P10*v_h)-(Eh*p_ExptoInf*(1-P10*v_h)) -(Eh*p_ExptoNoInf)  ##Exposed human population per 
    
    #dExp <- (beta_dh*Sh*X) #((beta_dh*Sh*Y) -(m_h*Eh) -(Eh*p_ExptoInf*P10*v_h)-(Eh*p_ExptoInf*(1-P10*v_h)) -(Eh*p_ExptoNoInf))* 3.5
    
    dIh <- (Eh*p_ExptoInf*(1-P10*v_h))-(m_h*Ih)-(mu_h*Ih)  #Rabid human population 
    
    #dDh <- (Eh*p_ExptoInf*(1-P10*v_h))   #
    
    dRh <- (Eh*p_ExptoInf*P10*v_h)+(v_h*alpha_h*Sh)-(m_h*Rh)-(lambda_h*Rh)
    
    dNh = dSh+dEh+dIh+dRh
    
    
    return(list(c(dX, dH1, dH2, dY, dI, dN, dSh, dEh,  dIh, dRh, dNh)))
  })
}





initial_state <- c(X = Ss, H1 = 0, H2 = 0, Y = 0.002, I = 0, N = Ss+0.002,
                   Sh = Sh, Eh = 0,   Ih = 0,  Rh = 0, Nh = Sh)  # initial number of rabbits and foxes

# Set time points for which to solve the ODEs
times <- seq(0, 4000, by = 1)  # from 0 to 200, in increments of 0.1

# Solve the ODEs
output <- ode(y = initial_state, times = times, func = lotka_volterra, parms = parameters)

output <- as.data.frame(output)

output

write.csv(output, "C:/Users/mfl3/Downloads//periurban.csv")




output$Exp <- output$Sh * output$Y * 3E-05


intro <- 52*3
est <- 52*14
endm <- 52*40

52*15

750/52

as.data.frame(output) %>%
  ggplot(aes(x=times, y=Y * 500)) +
  geom_line(color = "red") + xlim(0, endm) +
  geom_vline(xintercept = intro, linetype=4) +
  geom_vline(xintercept = est, linetype=4) +
  labs(x = "Weeks", y = "Rabid raccoons")



as.data.frame(output) %>%
  ggplot(aes(x=times, y=(Exp* 500))) +
  geom_line() + xlim(0, 1040)




#sum(output[,"Dh"][1:200]*162)

sum(output[,"Y"][0:intro]*500)/(intro/52)  
sum(output[,"Exp"][0:intro]*500)/(intro/52)
sum(output[,"Ih"][0:intro]*500)/(intro/52)



sum(output[,"Y"][intro:est]*500)/((est-intro)/52)  # 260, 1040
sum(output[,"Exp"][intro:est]*500)/((est-intro)/52)
sum(output[,"Ih"][intro:est]*500)/((est-intro)/52)

sum(output[,"Y"][est:endm]*500)/((endm-est)/52)  # 260, 1040
sum(output[,"Exp"][est:endm]*500)/((endm-est)/52)
sum(output[,"Ih"][est:endm]*500)/((endm-est)/52)


45*3.5

as.data.frame(output) %>%
  ggplot(aes(x=times, y=Y * 162)) +
  geom_line() + xlim(1, 260) + ylim(0, 50)


as.data.frame(output) %>%
  ggplot(aes(x=times, y=(Eh* 162))) +
  geom_line() + xlim(1, 4000)



rac_deaths <- as.data.frame(output) %>%
  ggplot(aes(x=times, y=Y * 100)) +
  geom_line() + xlim(1, 2080) #+ ylim(0, 300)


hum_exp   <- as.data.frame(output) %>%
  ggplot(aes(x=times, y= Exp* 100)) +
  geom_line() + xlim(1, 2080)

hum_exp_rab   <- as.data.frame(output) %>%
  ggplot(aes(x=times, y=(Eh* 100))) +
  geom_line() + xlim(1, 2080)


ggarrange(rac_deaths, hum_exp, hum_exp_rab, ncol = 1, labels = c("Raccoon deaths", "Cumulative human exposures","Cumulative human exposures"))





######################################

0.0001491/4.17992E-05



(55+40+47+42+40+70)/6


(40+90+55+77+50+60)/6


(5+2+15+7+11+5)/6

RR_graph <- read_excel("C:/Users/mfl3/Downloads/RR_graph.xlsx")

HE_graph <- read_excel("C:/Users/mfl3/Downloads/RR_graph.xlsx", sheet = "Sheet2")


ggplot(RR_graph, aes(x=time, y=Y * 500, group = area)) +
  geom_line(aes(color = area), size = 1.5)+
  xlim(0,1500) + ylim(0,2500)+ geom_vline(xintercept = intro, linetype=4) +
  geom_vline(xintercept = est, linetype=4) +
  labs(x = "Weeks", y = "Rabid raccoons")+ 
  theme(legend.title = element_text(size = 1), 
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold"))

ggplot(HE_graph, aes(x=time, y=Exp * 500, group = area)) +
  geom_line(aes(color = area), size = 1.5)+
  xlim(0,1500) + #+ ylim(0,2500)+ 
  geom_vline(xintercept = intro, linetype=4) +
  geom_vline(xintercept = est, linetype=4) +
  labs(x = "Weeks", y = "Human exposures") +
  theme(legend.title = element_text(size = 1), 
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold"))


