
library(tidyverse)
library(ggforce)
SIRD <- function(t, x, parms){
  
  
  with(as.list(c(parms,x)),{
    dS <- - beta*S*I
    dI <- + beta*S*I - gamma*I - mu*I
    dR <- gamma*I   
    dD <- mu * I
    der <- c(dS, dI,dR, dD)
    list(der)  
  }) 
  
}  # end of function definition



library(deSolve)

### INITIALIZE PARAMETER SETTINGS

parms <- c(beta=1.5e-5, gamma=0.5e-1, mu = 1e-2)		# set the parameters of the model
inits <- c(S=10000, I=1, R=0, D=0)		# set the initial values
dt    <- seq(0,300,0.1)			# set the time points for evaluation



simulation <- as.data.frame(lsoda(inits, dt, SIRD, parms=parms)) # this way our set 'parms' will be used as default

my.palette <- c("blue", "darkgreen", "orange", "black")
simulation%>% 
  pivot_longer(-time) %>% 
  ggplot() +
  geom_line(aes(x = time, y = value, col = name), size = 1) +
  theme_bw() +
  ggtitle("Epidémie de Grovid-42") +
  labs(col = "Compartiment") +
  xlab("Temps") +
  ylab("Individus dans le compartiment")  +
  scale_color_manual(values = my.palette[c(4,3,2,1)])

ggsave("grpsvidnormal.png", width = 1920/80 , height = 1080/80, units = "cm")


simulation%>% 
  mutate(R = beta1*S/(parms[2])+parms[3]) %>% 
  ggplot() +
  geom_line(aes(x = time, y = R), size = 1) +
  theme_bw() +
  ggtitle("Epidémie de Grovid-42", "Evolution du Rt") +
  labs(col = "Compartiment") +
  xlab("Temps") +
  ylab("Individus dans le compartiment") +
  geom_hline(aes(yintercept = 1), lty = 2, col = "red")



# simu confinement --------------------------------------------------------

beta1 = 1.5e-5
beta2 = 0.3e-5
beta3 = 1.15e-5
t1 = seq(0,70,0.1)
t2 = seq(70.1,170,0.1)
t3 = seq(180.1,700,0.1)

parms <- c(beta=beta1, gamma=0.5e-1, mu = 1e-2)		# set the parameters of the model
inits <- c(S=10000, I=1, R=0, D=0)		# set the initial values
dt    <- t1			# set the time points for evaluation

# Calculate and print R_0 on the screen
N <- sum(inits)
R_0 <- with(as.list(parms),{beta*N/gamma})
print(paste("R_0 =",R_0),quote=FALSE)


simulation1 <- as.data.frame(lsoda(inits, dt, SIRD, parms=parms)) # this way our set 'parms' will be used as default

#sim2
parms <- c(beta=beta2, gamma=0.5e-1, mu = 1e-2)		# set the parameters of the model
inits <- c(S=simulation1$S[nrow(simulation1)], 
           I=simulation1$I[nrow(simulation1)],
           R=simulation1$R[nrow(simulation1)], 
           D=simulation1$D[nrow(simulation1)])		# set the initial values
dt    <- t2			# set the time points for evaluation

simulation2 <- as.data.frame(lsoda(inits, dt, SIRD, parms=parms)) # this way our set 'parms' will be used as default

#sim3 

parms <- c(beta=beta3, gamma=0.5e-1, mu = 1e-2)		# set the parameters of the model
inits <- c(S=simulation2$S[nrow(simulation2)], 
           I=simulation2$I[nrow(simulation2)],
           R=simulation2$R[nrow(simulation2)], 
           D=simulation2$D[nrow(simulation2)])		# set the initial values
dt    <- t3			# set the time points for evaluation

simulation3 <- as.data.frame(lsoda(inits, dt, SIRD, parms=parms)) # this way our set 'parms' will be used as default




my.palette <- c("blue", "darkgreen", "orange", "black")


rbind(simulation1, simulation2, simulation3)%>% 
  select(time, I, D) %>% 
  pivot_longer(-time) %>% 
  ggplot() +
  geom_line(aes(x = time, y = value, col = name), size=1) +
  theme_bw()  +
  ggtitle("Epidémie de Grovid-42","Confinement en début de phase + reprise") +
  labs(col = "Compartiment") +
  scale_color_manual(values = my.palette[c(4,3)])

ggsave("grovidgerevue.png", width = 1920/80 , height = 1080/80, units = "cm")

rbind(simulation1, simulation2, simulation3)%>% 
  pivot_longer(-time) %>% 
  ggplot() +
  geom_line(aes(x = time, y = value, col = name), size=1) +
  theme_bw()  +
  ggtitle("Epidémie de Grovid-42","Confinement en début de phase + reprise") +
  labs(col = "Compartiment")+
  scale_color_manual(values = my.palette[c(4,3,2,1)])

ggsave("grovidgere.png", width = 1920/80 , height = 1080/80, units = "cm")            


rbind(simulation1, simulation2, simulation3)%>% 
  select(time, D) %>% 
  pivot_longer(-time) %>% 
  ggplot() +
  geom_line(aes(x = time, y = value, col = name), size=1) +
  theme_bw()  +
  ggtitle("Epidémie de Grovid-42","Confinement en début de phase + reprise") +
  labs(col = "Compartiment")+
  scale_color_manual(values = my.palette[c(4,3,2,1)]) +
  facet_zoom( ylim = c(0,250), xlim = c(40,220), split = FALSE, horizontal = FALSE)

ggsave("grovidgeremortszoom.png", width = 1920/80 , height = 1080/80, units = "cm")            


simulation1$beta <- beta1
simulation2$beta <- beta2
simulation3$beta <- beta3

df <- rbind(simulation1, simulation2, simulation3)

df%>% 
  mutate(R = beta*S/(parms[2]+parms[3])) %>% 
  ggplot() +
  geom_line(aes(x = time, y = R), size = 1) +
  theme_bw() +
  ggtitle("Epidémie de Grovid-42", "Evolution du Rt") +
  labs(col = "Compartiment") +
  xlab("Temps") +
  ylab("Individus dans le compartiment") +
  geom_hline(aes(yintercept = 1), lty = 2, col = "red")


library(gganimate)

anim <- rbind(simulation1, simulation2, simulation3)%>% 
  select(time, I, D) %>% 
  pivot_longer(-time) %>% 
  ggplot() +
  geom_line(aes(x = time, y = value, col = name), size=1) +
  theme_bw()  +
  ggtitle("Epidémie de Grovid-42") +
  labs(col = "Compartiment") +
  scale_color_manual(values = my.palette[c(4,3)]) +
  transition_reveal(time) +
  view_follow()


animate(
  anim,
  width = 1920, height = 1080, renderer = av_renderer(file = "epid.mp4"),
  fps = 24, nframes = 24*40, res= 250) 
