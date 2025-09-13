library(deSolve, warn.conflicts=F, quietly=T)
library(cowplot, warn.conflicts=F, quietly=T)
library(ggplot2, warn.conflicts=F, quietly=T)
library(gridExtra, warn.conflicts=F, quietly=T)
library(reshape2, warn.conflicts=F, quietly=T)
library(viridis, warn.conflicts=F, quietly=T)
library(tidyverse, warn.conflicts=F, quietly=T)

## General models

###------------------ No temperature model

# Define the system of differential equations
mosquito_model = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total human population
    N_h = S + E + I + R
    
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

# Parameters (updated with your research)
parameters = c(
  beta_ae = 0.340,  # Reproductive rate of Aedes aegypti
  beta_alb = 0.253, # Reproductive rate of Aedes albopictus
  omega_ae = 0.300, # Competition coefficient aegypti on albopictus
  omega_alb = 0.500,# Competition coefficient albopictus on aegypti
  C = 2000,         # Carrying capacity for larvae
  kappa_ae = 0.170, # Larval development rate of Aedes aegypti
  kappa_alb = 0.167,# Larval development rate of Aedes albopictus
  mu_Lae = 0.090,   # Larval mortality rate of Aedes aegypti
  mu_Lalb = 0.050,  # Larval mortality rate of Aedes albopictus
  mu_Aae = 1/10,    # Adult mortality rate of Aedes aegypti (1/Lv_ae)
  mu_Aalb = 1/12,   # Adult mortality rate of Aedes albopictus (1/Lv_alb)
  zeta_ae = 0.700 * 0.470, # Transmission coefficient for Aedes aegypti
  zeta_alb = 0.700 * 0.350, # Transmission coefficient for Aedes albopictus
  p_D = 0.700,      # Probability of Dengue transmission from mosquito to human
  BR_ae = 0.470,    # Biting rate for Aedes aegypti 
  BR_alb = 0.350,   # Biting rate for Aedes albopictus 
  theta_h = 1/10,   # Incubation rate of dengue in humans
  gamma_h = 1/7,    # Recovery rate of dengue in humans
  #rho = 1/90,      # Inmmunity rate
  phi_ae = 0.000,   # Competition coefficient for adult Aedes aegypti
  phi_alb = 0.000,  # Competition coefficient for adult Aedes albopictus
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
time = seq(0, 200, by = 1) # From 0 to 200 days in 1-day steps

# Solve the model
output = ode(y = initial_state, times = time, func = mosquito_model, parms = parameters)

# Convert output to a data frame for visualization
output_df = as.data.frame(output)

# Plot SEIR results
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = S, color = "Susceptible"), linewidth = 1) +
  geom_line(aes(y = E, color = "Exposed"), linewidth = 1) +
  geom_line(aes(y = I, color = "Infected"), linewidth = 1) +
  geom_line(aes(y = R, color = "Recovered"), linewidth = 1) +
  labs(title = "Human SEIR Dynamics", x = "Time (days)", y = "Population") +
  theme_linedraw() +
  scale_color_manual(values = c("firebrick", "olivedrab", "grey68", "royalblue")) +
  guides(color = guide_legend(title = "SEIR Compartments"))

# Plot results for larval populations (Aedes aegypti and Aedes albopictus)
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = L_ae, color = "Larvae Aedes aegypti"), linewidth = 1) +
  geom_line(aes(y = L_alb, color = "Larvae Aedes albopictus"), linewidth = 1) +
  labs(title = "Larval Populations", x = "Time (days)", y = "Larvae Population") +
  theme_linedraw() +
  scale_color_manual(values = c("purple", "orange")) +
  guides(color = guide_legend(title = "Larvae"))

# Plot results for adult populations (Aedes aegypti and Aedes albopictus)
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = X_ae, color = "Adults Aedes aegypti"), linewidth = 1) +
  geom_line(aes(y = X_alb, color = "Adults Aedes albopictus"), linewidth = 1) +
  labs(title = "Adult Populations", x = "Time (days)", y = "Adult Population") +
  theme_linedraw() +
  scale_color_manual(values = c("purple", "orange")) +
  guides(color = guide_legend(title = "Adults"))

# Plot results for infected adult populations (Aedes aegypti and Aedes albopictus)
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = Y_ae, color = "Infected Adults Aedes aegypti"), linewidth = 1) +
  geom_line(aes(y = Y_alb, color = "Infected Adults Aedes albopictus"), linewidth = 1) +
  labs(title = "Infected Adult Populations", x = "Time (days)", y = "Infected Adult Population") +
  theme_linedraw() +
  scale_color_manual(values = c("purple", "orange")) +
  guides(color = guide_legend(title = "Infected Adults"))


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
    dS = - lambda * S # + rho * R
    dE = lambda * S - theta_h * E
    dI = theta_h * E - gamma_h * I
    dR = gamma_h * I # - rho * R
    
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
  #rho = 1/90,      # Inmmunity rate
  phi_ae = 0.000,   # Competition coefficient for adult Aedes aegypti
  phi_alb = 0.000,  # Competition coefficient for adult Aedes albopictus
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

# Plot results
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = S, color = "Susceptible"), linewidth = 1) +
  geom_line(aes(y = E, color = "Exposed"), linewidth = 1) +
  geom_line(aes(y = I, color = "Infected"), linewidth = 1) +
  geom_line(aes(y = R, color = "Recovered"), linewidth = 1) +
  labs(title = "Human SEIR Dynamics", x = "Time (days)", y = "Population") +
  theme_linedraw() +
  scale_color_manual(values = c("firebrick", "olivedrab", "grey68", "royalblue")) +
  guides(color = guide_legend(title = "SEIR Compartments"))

# Plot results for larval populations (Aedes aegypti and Aedes albopictus)
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = L_ae, color = "Larvae Aedes aegypti"), linewidth = 1) +
  geom_line(aes(y = L_alb, color = "Larvae Aedes albopictus"), linewidth = 1) +
  labs(title = "Larval Populations", x = "Time (days)", y = "Larvae Population") +
  theme_linedraw() +
  scale_color_manual(values = c("purple", "orange")) +
  guides(color = guide_legend(title = "Larvae"))

# Plot results for adult populations (Aedes aegypti and Aedes albopictus)
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = X_ae, color = "Adults Aedes aegypti"), linewidth = 1) +
  geom_line(aes(y = X_alb, color = "Adults Aedes albopictus"), linewidth = 1) +
  labs(title = "Adult Populations", x = "Time (days)", y = "Adult Population") +
  theme_linedraw() +
  scale_color_manual(values = c("purple", "orange")) +
  guides(color = guide_legend(title = "Adults"))

# Plot results for infected adult populations (Aedes aegypti and Aedes albopictus)
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = Y_ae, color = "Infected Adults Aedes aegypti"), linewidth = 1) +
  geom_line(aes(y = Y_alb, color = "Infected Adults Aedes albopictus"), linewidth = 1) +
  labs(title = "Infected Adult Populations", x = "Time (days)", y = "Infected Adult Population") +
  theme_linedraw() +
  scale_color_manual(values = c("purple", "orange")) +
  guides(color = guide_legend(title = "Infected Adults"))
