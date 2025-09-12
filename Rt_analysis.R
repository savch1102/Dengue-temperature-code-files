library(deSolve, warn.conflicts=F, quietly=T)
library(cowplot, warn.conflicts=F, quietly=T)
library(ggplot2, warn.conflicts=F, quietly=T)
library(gridExtra, warn.conflicts=F, quietly=T)
library(reshape2, warn.conflicts=F, quietly=T)
library(viridis, warn.conflicts=F, quietly=T)
library(tidyverse, warn.conflicts=F, quietly=T)


##------------ No Temperature - *Ae. aegypti* greater abundance



# Define the system of differential equations
mosquito_model = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total human population
    N_h = S + E + I + R
    
    # Larval dynamics
    dL_ae = beta_ae * L_ae * (1 - (L_ae + omega_alb * L_alb) / C) - kappa_ae * L_ae - mu_Lae * L_ae
    dL_alb = beta_alb * L_alb * (1 - (L_alb + omega_ae * L_ae) / C) - kappa_alb * L_alb - mu_Lalb * L_alb
    
    # Adult susceptible dynamics
    dX_ae = kappa_ae * L_ae * (1 - (X_ae + phi_alb * (X_alb + Y_alb)) / C_ad) - zeta_ae * X_ae * I / N_h - mu_ae * X_ae
    dX_alb = kappa_alb * L_alb * (1 - (X_alb + phi_ae * (X_ae + Y_ae)) / C_ad) - zeta_alb * X_alb * I / N_h - mu_alb * X_alb
    
    # Adult infected dynamics
    dY_ae = zeta_ae * X_ae * I / N_h * (1 - (Y_ae + phi_alb * (X_alb + Y_alb)) / C_ad) - mu_ae * Y_ae
    dY_alb = zeta_alb * X_alb * I / N_h * (1 - (Y_alb + phi_ae * (X_ae + Y_ae)) / C_ad) - mu_alb * Y_alb
    
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
  C = 3000,         # Carrying capacity for larvae
  kappa_ae = 0.170, # Larval development rate of Aedes aegypti
  kappa_alb = 0.167,# Larval development rate of Aedes albopictus
  mu_Lae = 0.090,   # Larval mortality rate of Aedes aegypti
  mu_Lalb = 0.050,  # Larval mortality rate of Aedes albopictus
  mu_ae = 1/10,     # Adult mortality rate of Aedes aegypti (1/Lv_ae)
  mu_alb = 1/12,    # Adult mortality rate of Aedes albopictus (1/Lv_alb)
  zeta_ae = 0.700 * 0.470, # Transmission coefficient for Aedes aegypti
  zeta_alb = 0.700 * 0.350, # Transmission coefficient for Aedes albopictus
  p_D = 0.700,      # Probability of Dengue transmission from mosquito to human
  BR_ae = 0.470,    # Biting rate for Aedes aegypti 
  BR_alb = 0.350,   # Biting rate for Aedes albopictus
  theta_h = 1/10,   # Incubation rate of dengue in humans
  gamma_h = 1/7,    # Recovery rate of dengue in humans
  #rho = 1/90,      # Immunity rate
  phi_ae = 0.00,   # Competition coefficient for adult Aedes aegypti
  phi_alb = 0.00,  # Competition coefficient for adult Aedes albopictus
  C_ad = 10000      # Carrying capacity for adult mosquitoes
)

# Initial conditions (commented)
initial_state = c(
  L_ae = 200,  # Initial larvae population for Aedes aegypti
  L_alb = 50,  # Initial larvae population for Aedes albopictus
  X_ae = 3000, # Initial adult susceptible population for Aedes aegypti
  X_alb = 1000,# Initial adult susceptible population for Aedes albopictus
  Y_ae = 0,    # Initial adult infected population for Aedes aegypti
  Y_alb = 0,   # Initial adult infected population for Aedes albopictus
  S = 999,     # Initial susceptible human population
  E = 1,       # Initial exposed human population
  I = 0,       # Initial infected human population
  R = 0        # Initial recovered human population
)

# Time vector
time = seq(0, 200, by = 1) # From 0 to 200 days in 1-day steps

# Values for omega_ae and omega_alb
omega_ae = seq(0.1, 1.0, 0.05)       # Competition coefficient of Aedes aegypti over Aedes albopictus
omega_alb = rev(seq(0.1, 1.0, 0.05)) # Competition coefficient of Aedes albopictus over Aedes aegypti

Rt_func = function(omega_ae, omega_alb) {
  parameters["omega_ae"] = omega_ae
  parameters["omega_alb"] = omega_alb
  
  # Solve the model until equilibrium is reached
  output = ode(y = initial_state, times = time, func = mosquito_model, parms = parameters)
  
  output = ode(y = initial_state, times = time, func = mosquito_model, parms = parameters)
  
  output = as.data.frame(output)
  
  max_index = which.max(output$I)
  last_row = output[max_index-25, ]
  
  X_ae = last_row$X_ae
  X_alb = last_row$X_alb
  N_h = last_row$S + last_row$E + last_row$I + last_row$R
  mu_Aae = 1/10     # Adult mortality rate of Aedes aegypti (1/Lv_ae)
  mu_Aalb = 1/12    # Adult mortality rate of Aedes albopictus (1/Lv_alb)
  zeta_ae = 0.700 * 0.470 # Transmission coefficient for Aedes aegypti
  zeta_alb = 0.700 * 0.350 # Transmission coefficient for Aedes albopictus
  p_D = 0.700      # Probability of Dengue transmission from mosquito to human
  BR_ae = 0.470    # Biting rate for Aedes aegypti
  BR_alb = 0.350   # Biting rate for Aedes albopictus
  gamma_h = 1/7    # Recovery rate of dengue in humans
  
 Rt =  p_D * sqrt((BR_ae^2 * mu_Aalb * X_ae + BR_alb^2 * mu_Aae * X_alb) /
    (N_h * gamma_h * mu_Aae * mu_Aalb)) * (last_row$S/(last_row$S+last_row$E+last_row$I+last_row$R))
 return(Rt)
}

# Zeros matrix
Rt_matrix = matrix(0, nrow = length(omega_ae), ncol = length(omega_alb))

rownames(Rt_matrix) = omega_ae  # Set row names as omega_ae
colnames(Rt_matrix) = omega_alb # Set column names as omega_alb

for (i in seq_along(omega_ae)) {
  for (j in seq_along(omega_alb)) {
    Rt_matrix[i, j] = Rt_func(omega_ae[i], omega_alb[j])
  }
}

# Convert the matrix to a dataframe for ggplot
Rt_df = as.data.frame(as.table(Rt_matrix))
colnames(Rt_df) = c("omega_ae", "omega_alb", "Rt")

# Convert columns to numeric
Rt_df$omega_ae = as.numeric(as.character(Rt_df$omega_ae))
Rt_df$omega_alb = as.numeric(as.character(Rt_df$omega_alb))

# Create the heatmap with ggplot2
f11 = ggplot(Rt_df, aes(x = omega_ae, y = omega_alb, fill = Rt)) +
  geom_tile(color = "black", linewidth = 0.05) +
  scale_fill_viridis(option = "magma") +  # Use the viridis color palette
  labs(x = expression(omega[ae]), y = expression(omega[alb]), 
       fill = expression(R[t]), subtitle = "A) No temperature") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.y = element_text(angle = 90, vjust = 0.5),   # Ensures y-axis label is horizontal and aligned
    axis.text.y = element_text(hjust = 1),                   # Adjusts y-axis text alignment
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),                   # Reduce legend text size
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0),  # Subtitle bold and aligned left
    legend.key.size = unit(0.3, "cm")                      # Reduce legend bar size
  )


##------------ Temperature - *Ae. aegypti* greater abundance

# Define the system of differential equations
mosquito_model_temp = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total human population
    N_h = S + E + I + R
    
    # Temperature function
    T = 2.7 * sin(-0.04 * t) + 29.8
    
    # Development rates
    kappa_ae = (1 / (2.7 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (2.7^2)))
    kappa_alb = (1 / (2.4 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (2.4^2)))
    
    # Biting rates
    br_ae = (1 / (1.1 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (1.1^2)))
    br_alb = (1 / (1.4 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (1.4^2)))
    
    # Mortality rates
    mu_Lae = (0.05 * (T - 25.5))^2 + 0.09
    mu_Lalb = (0.05 * (T - 25.5))^2 + 0.05
    mu_Aae = (0.05 * (T - 25.5))^2 + 0.015
    mu_Aalb = (0.05 * (T - 25.5))^2 + 0.01
    
    # Fecundity rates
    beta_ae = 70 * (1 / (3.5 * sqrt(2 * pi))) * exp(-((T - 27)^2) / (2 * (3.5^2)))
    beta_alb = 64 * (1 / (3.5 * sqrt(2 * pi))) * exp(-((T - 27)^2) / (2 * (3.5^2)))
    
    # Transmission probability
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
    lambda = p_D * (br_ae * Y_ae + br_alb * Y_alb) / N_h # Force of infection
    dS = - lambda * S #+ rho * R
    dE = lambda * S - theta_h * E
    dI = theta_h * E - gamma_h * I
    dR = gamma_h * I #- rho * R
    
    list(c(dL_ae, dL_alb, dX_ae, dX_alb, dY_ae, dY_alb, dS, dE, dI, dR), p_D = p_D, br_ae = br_ae, br_alb = br_alb, N_h = N_h, mu_Aalb = mu_Aalb, mu_Aae = mu_Aae, zeta_ae = zeta_ae, zeta_alb = zeta_alb, gamma_h = gamma_h)
  })
}

# Parameters
parameters = c(
  C = 3000,
  zeta_ae = 0.700 * 0.470,
  zeta_alb = 0.700 * 0.350,
  p_D = 0.700,
  theta_h = 1/10,
  gamma_h = 1/7,
  #rho = 1/90,
  phi_ae = 0.00,
  phi_alb = 0.00,
  C_ad = 10000
)

# Initial conditions
initial_state = c(
  L_ae = 200,
  L_alb = 50,
  X_ae = 3000,
  X_alb = 1000,
  Y_ae = 0,
  Y_alb = 0,
  S = 999,
  E = 1,
  I = 0,
  R = 0
)

# Time vector
time = seq(0, 500, by = 1)

# Competition coefficients
omega_ae = seq(0.1, 1.0, 0.05)
omega_alb = rev(seq(0.1, 1.0, 0.05))

Rt_func = function(omega_ae, omega_alb) {
  parameters["omega_ae"] = omega_ae
  parameters["omega_alb"] = omega_alb
  
  # Solve the model until equilibrium is reached
  output = ode(y = initial_state, times = time, func = mosquito_model_temp, parms = parameters)
  
  output = as.data.frame(output)
  
  max_index = which.max(output$I)
  last_row = output[max_index-25, ]
  
  X_ae = last_row$X_ae
  X_alb = last_row$X_alb
  N_h = last_row$S + last_row$E + last_row$I + last_row$R
  p_D = last_row$p_D
  BR_ae = last_row$br_ae
  BR_alb = last_row$br_alb
  mu_Aalb = last_row$mu_Aalb
  mu_Aae = last_row$mu_Aae
  zeta_ae = last_row$zeta_ae
  zeta_alb = last_row$zeta_alb
  gamma_h = last_row$gamma_h
  
  Rt =  p_D * sqrt((BR_ae^2 * mu_Aalb * X_ae + BR_alb^2 * mu_Aae * X_alb) /
    (N_h * gamma_h * mu_Aae * mu_Aalb)) * (last_row$S/(last_row$S+last_row$E+last_row$I+last_row$R))
  return(Rt)
}

# Zeros matrix
Rt_matrix = matrix(0, nrow = length(omega_ae), ncol = length(omega_alb))

rownames(Rt_matrix) = omega_ae  # Set row names as omega_ae
colnames(Rt_matrix) = omega_alb # Set column names as omega_alb

for (i in seq_along(omega_ae)) {
  for (j in seq_along(omega_alb)) {
    Rt_matrix[i, j] = Rt_func(omega_ae[i], omega_alb[j])
  }
}

# Convert the matrix to a dataframe for ggplot
Rt_df = as.data.frame(as.table(Rt_matrix))
colnames(Rt_df) = c("omega_ae", "omega_alb", "Rt")

# Convert columns to numeric
Rt_df$omega_ae = as.numeric(as.character(Rt_df$omega_ae))
Rt_df$omega_alb = as.numeric(as.character(Rt_df$omega_alb))

# Create the heatmap with ggplot2
f12 = ggplot(Rt_df, aes(x = omega_ae, y = omega_alb, fill = Rt)) +
  geom_tile(color = "black", linewidth = 0.05) +
  scale_fill_viridis(option = "magma") +  # Use the viridis color palette
  labs(x = expression(omega[ae]), y = expression(omega[alb]), 
       fill = expression(R[t]), subtitle = "B) Temperature") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.y = element_text(angle = 90, vjust = 0.5),   # Ensures y-axis label is horizontal and aligned
    axis.text.y = element_text(hjust = 1),                   # Adjusts y-axis text alignment
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),                   # Reduce legend text size
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0),  # Subtitle bold and aligned left
    legend.key.size = unit(0.3, "cm")                      # Reduce legend bar size
  )




### Plot

plot_grid(f11, f12)
