library(deSolve, warn.conflicts=F, quietly=T)
library(cowplot, warn.conflicts=F, quietly=T)
library(ggplot2, warn.conflicts=F, quietly=T)
library(gridExtra, warn.conflicts=F, quietly=T)
library(reshape2, warn.conflicts=F, quietly=T)
library(viridis, warn.conflicts=F, quietly=T)
library(tidyverse, warn.conflicts=F, quietly=T)

###------------------ Temperature model

# Define the system of differential equations
mosquito_model_temp = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total human population
    N_h = S + E + I + R
    
    # Temperature function
    T = 2.7 * sin(-0.04 * t) + 29.8
    
    # Development rate functions
    kappa_ae = (1 / (2.7 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (2.7^2)))
    kappa_alb = (1 / (2.4 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (2.4^2)))
    
    # Biting rate functions
    br_ae = (1 / (1.1 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (1.1^2)))
    br_alb = (1 / (1.4 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (1.4^2)))
    
    # Mortality rate functions
    mu_Lae = (0.05 * (T - 25.5))^2 + 0.09
    mu_Lalb = (0.05 * (T - 25.5))^2 + 0.05
    mu_Aae = (0.05 * (T - 25.5))^2 + 0.015
    mu_Aalb = (0.05 * (T - 25.5))^2 + 0.01
    
    # Fecundity rate functions
    beta_ae = 70 * (1 / (3.5 * sqrt(2 * pi))) * exp(-((T - 27)^2) / (2 * (3.5^2)))
    beta_alb = 64 * (1 / (3.5 * sqrt(2 * pi))) * exp(-((T - 27)^2) / (2 * (3.5^2)))
    
    # Transmission probability function
    p_D = 4.8 * (1 / (2.4 * sqrt(2 * pi))) * exp(-((T - 30)^2) / (2 * (2.4^2)))
    
    # Larval dynamics
    dL_ae = beta_ae * L_ae * (1 - (L_ae + omega_alb * L_alb) / C) - kappa_ae * L_ae - mu_Lae * L_ae
    dL_alb = beta_alb * L_alb * (1 - (L_alb + omega_ae * L_ae) / C) - kappa_alb * L_alb - mu_Lalb * L_alb
    
    # Adult susceptible dynamics
    dX_ae = kappa_ae * L_ae * (1 - (X_ae + phi_alb * (X_alb + Y_alb)) / C_ad) - zeta_ae * X_ae * I / N_h - mu_Aae * X_ae
    dX_alb = kappa_alb * L_alb * (1 - (X_alb + phi_ae * (X_ae + Y_ae)) / C_ad) - zeta_alb * X_alb * I / N_h - mu_Aalb * X_alb
    
    # Adult infected dynamics
    dY_ae = zeta_ae * X_ae * I / N_h * (1 - (Y_ae + phi_alb * (X_alb + Y_alb)) / C_ad) - mu_Aae * Y_ae
    dY_alb = zeta_alb * X_alb * I / N_h * (1 - (Y_alb + phi_ae * (X_ae + Y_ae)) / C_ad) - mu_Aalb * Y_alb
    
    # Human SEIR dynamics
    lambda = p_D * (BR_ae * Y_ae + BR_alb * Y_alb) / N_h # Force of infection
    dS = - lambda * S #+ rho * R
    dE = lambda * S - theta_h * E
    dI = theta_h * E - gamma_h * I
    dR = gamma_h * I #- rho * R
    
    # Return derivatives
    list(c(dL_ae, dL_alb, dX_ae, dX_alb, dY_ae, dY_alb, dS, dE, dI, dR))
  })
}

# Parameters (updated with temperature-dependent functions)
parameters = c(
  omega_ae = 0.300, # Competition coefficient aegypti on albopictus
  omega_alb = 0.500,# Competition coefficient albopictus on aegypti
  C = 2000,         # Carrying capacity for larvae
  zeta_ae = 0.700 * 0.470, # Transmission coefficient for Aedes aegypti
  zeta_alb = 0.700 * 0.350, # Transmission coefficient for Aedes albopictus
  p_D = 0.700,      # Probability of Dengue transmission from mosquito to human
  BR_ae = 0.470,    # Biting rate for Aedes aegypti 
  BR_alb = 0.350,   # Biting rate for Aedes albopictus 
  theta_h = 1/10,   # Incubation rate of dengue in humans
  gamma_h = 1/7,    # Recovery rate of dengue in humans
  # rho = 1/90,     # Inmmunity rate
  phi_ae = 0.300,   # Competition coefficient for adult Aedes aegypti
  phi_alb = 0.500,  # Competition coefficient for adult Aedes albopictus
  C_ad = 10000      # Carrying capacity for adult mosquitoes
)

# Initial conditions (commented)
initial_state = c(
  L_ae = 100,  # Initial larvae population for Aedes aegypti
  L_alb = 100, # Initial larvae population for Aedes albopictus
  X_ae = 2000, # Initial adult susceptible population for Aedes aegypti
  X_alb = 2000,# Initial adult susceptible population for Aedes albopictus
  Y_ae = 0,    # Initial adult infected population for Aedes aegypti
  Y_alb = 0,   # Initial adult infected population for Aedes albopictus
  S = 999,     # Initial susceptible human population
  E = 1,       # Initial exposed human population
  I = 0,       # Initial infected human population
  R = 0        # Initial recovered human population
)

# Time vector
time = seq(0, 500, by = 1) # From 0 to 500 days in 1-day steps

# Solve the model
output = ode(y = initial_state, times = time, func = mosquito_model_temp,
             parms = parameters)

# Convert output to a data frame for visualization
output_df = as.data.frame(output)

output_df = output_df %>% select(time, S,E,I,R) %>% 
  pivot_longer(cols = "S":"R", names_to = "Pop", values_to = "N")

output_df$Pop = factor(output_df$Pop, levels = c("S", "E", "I", "R"))

# Plot results
ggplot(output_df, aes(x = time, y = N, color = Pop)) +
  geom_line(linewidth = 1) +
  labs(title = "Human SEIR Dynamics", x = "Time (days)", y = "Population") +
  theme_bw() +
  scale_color_manual(values = c("firebrick", "olivedrab", "grey68", "royalblue")) +
  guides(color = guide_legend(title = "Stage"))
