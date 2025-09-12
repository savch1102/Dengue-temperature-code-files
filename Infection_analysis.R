library(deSolve, warn.conflicts=F, quietly=T)
library(cowplot, warn.conflicts=F, quietly=T)
library(ggplot2, warn.conflicts=F, quietly=T)
library(gridExtra, warn.conflicts=F, quietly=T)
library(reshape2, warn.conflicts=F, quietly=T)
library(viridis, warn.conflicts=F, quietly=T)
library(tidyverse, warn.conflicts=F, quietly=T)


### Infected - No temperature - Same abundance

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

# Values for omega_ae and omega_alb
omega_ae = seq(0.1, 1.0, 0.05)       # Competition coefficient of Aedes aegypti over Aedes albopictus
omega_alb = rev(seq(0.1, 1.0, 0.05)) # Competition coefficient of Aedes albopictus over Aedes aegypti

infection = function(omega_ae_val, omega_alb_val) {
  parameters["omega_ae"] = omega_ae_val
  parameters["omega_alb"] = omega_alb_val
  
  # Solve the model until equilibrium is reached
  output = ode(y = initial_state, times = time, func = mosquito_model, parms = parameters)
  
  output = as.data.frame(output)
  return(output)
}

# Zeros matrix
infection_matrix = matrix(0, nrow = length(omega_ae), ncol = length(omega_alb))

rownames(infection_matrix) = omega_ae # Set row names as omega_ae
colnames(infection_matrix) = omega_alb # Set column names as omega_alb

for (i in seq_along(omega_ae)) {
  for (j in seq_along(omega_alb)) {
    infection_matrix[i,j] = max(infection(omega_ae[i], omega_alb[j])$I)
  }
}

# Convert the matrix to a dataframe for ggplot
infection_df = as.data.frame(as.table(infection_matrix))
colnames(infection_df) = c("omega_ae", "omega_alb", "I")

# Convert columns to numeric
infection_df$omega_ae = as.numeric(as.character(infection_df$omega_ae))
infection_df$omega_alb = as.numeric(as.character(infection_df$omega_alb))

# Restore the original order of omega_alb (from 1 to 0) as a factor
#infection_df$omega_alb = factor(infection_df$omega_alb, levels = sort(unique(infection_df$omega_alb), decreasing = TRUE))

# Create the heatmap with ggplot2
f5 = ggplot(infection_df, aes(x = omega_ae, y = omega_alb, fill = I)) +
  geom_tile(color="black", linewidth = 0.05) +
  scale_fill_viridis(option = "viridis") +  # Use the viridis color palette
  labs(y = expression(omega[alb]~""), 
       x = expression(omega[ae]~""), 
       fill = "Max Infection", subtitle = "A) No Temperature | Same abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.y = element_text(angle = 90, vjust = 0.5),   # Ensures y-axis label is horizontal and aligned
    axis.text.y = element_text(hjust = 1),                   # Adjusts y-axis text alignment
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),                   # Reduce legend text size
    plot.subtitle = element_text(size = 14, face = "bold", hjust = 0),  # Subtitle bold and aligned left
    legend.key.size = unit(0.3, "cm")                      # Reduce legend bar size
  )



### Infected - No temperature - *Ae. aegypti* greater abundance



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
  phi_ae = 0.000,   # Competition coefficient for adult Aedes aegypti
  phi_alb = 0.000,  # Competition coefficient for adult Aedes albopictus
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

infection = function(omega_ae_val, omega_alb_val) {
  parameters["omega_ae"] = omega_ae_val
  parameters["omega_alb"] = omega_alb_val
  
  # Solve the model until equilibrium is reached
  output = ode(y = initial_state, times = time, func = mosquito_model, parms = parameters)
  
  output = as.data.frame(output)
  return(output)
}

# Zeros matrix
infection_matrix = matrix(0, nrow = length(omega_ae), ncol = length(omega_alb))

rownames(infection_matrix) = omega_ae # Set row names as omega_ae
colnames(infection_matrix) = omega_alb # Set column names as omega_alb

for (i in seq_along(omega_ae)) {
  for (j in seq_along(omega_alb)) {
    infection_matrix[i,j] = max(infection(omega_ae[i], omega_alb[j])$I)
  }
}

# Convert the matrix to a dataframe for ggplot
infection_df = as.data.frame(as.table(infection_matrix))
colnames(infection_df) = c("omega_ae", "omega_alb", "I")

# Convert columns to numeric
infection_df$omega_ae = as.numeric(as.character(infection_df$omega_ae))
infection_df$omega_alb = as.numeric(as.character(infection_df$omega_alb))

# Restore the original order of omega_alb (from 1 to 0) as a factor
#infection_df$omega_alb = factor(infection_df$omega_alb, levels = sort(unique(infection_df$omega_alb), decreasing = TRUE))

# Create the heatmap with ggplot2
f6 = ggplot(infection_df, aes(x = omega_ae, y = omega_alb, fill = I)) +
  geom_tile(color="black", linewidth = 0.05) +
  scale_fill_viridis(option = "viridis") +  # Use the viridis color palette
  labs(y = expression(omega[alb]~""), 
       x = expression(omega[ae]~""), 
       fill = "Max Infection", subtitle = "A) No Temperature") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.y = element_text(angle = 90, vjust = 0.5),   # Ensures y-axis label is horizontal and aligned
    axis.text.y = element_text(hjust = 1),                   # Adjusts y-axis text alignment
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),                   # Reduce legend text size
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0),  # Subtitle bold and aligned left
    legend.key.size = unit(0.3, "cm")                      # Reduce legend bar size
  )



### Infected - No temperature - *Ae. albopictus* greater abundance



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
  phi_ae = 0.000,   # Competition coefficient for adult Aedes aegypti
  phi_alb = 0.000,  # Competition coefficient for adult Aedes albopictus
  C_ad = 10000      # Carrying capacity for adult mosquitoes
)

# Initial conditions (commented)
initial_state = c(
  L_ae = 50,   # Initial larvae population for Aedes aegypti
  L_alb = 200, # Initial larvae population for Aedes albopictus
  X_ae = 1000, # Initial adult susceptible population for Aedes aegypti
  X_alb = 3000,# Initial adult susceptible population for Aedes albopictus
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

infection = function(omega_ae_val, omega_alb_val) {
  parameters["omega_ae"] = omega_ae_val
  parameters["omega_alb"] = omega_alb_val
  
  # Solve the model until equilibrium is reached
  output = ode(y = initial_state, times = time, func = mosquito_model, parms = parameters)
  
  output = as.data.frame(output)
  return(output)
}

# Zeros matrix
infection_matrix = matrix(0, nrow = length(omega_ae), ncol = length(omega_alb))

rownames(infection_matrix) = omega_ae # Set row names as omega_ae
colnames(infection_matrix) = omega_alb # Set column names as omega_alb

for (i in seq_along(omega_ae)) {
  for (j in seq_along(omega_alb)) {
    infection_matrix[i,j] = max(infection(omega_ae[i], omega_alb[j])$I)
  }
}

# Convert the matrix to a dataframe for ggplot
infection_df = as.data.frame(as.table(infection_matrix))
colnames(infection_df) = c("omega_ae", "omega_alb", "I")

# Convert columns to numeric
infection_df$omega_ae = as.numeric(as.character(infection_df$omega_ae))
infection_df$omega_alb = as.numeric(as.character(infection_df$omega_alb))

# Restore the original order of omega_alb (from 1 to 0) as a factor
#infection_df$omega_alb = factor(infection_df$omega_alb, levels = sort(unique(infection_df$omega_alb), decreasing = TRUE))

# Create the heatmap with ggplot2
f7 = ggplot(infection_df, aes(x = omega_ae, y = omega_alb, fill = I)) +
  geom_tile(color="black", linewidth = 0.05) +
  scale_fill_viridis(option = "viridis") +  # Use the viridis color palette
  labs(y = expression(omega[alb]~""), 
       x = expression(omega[ae]~""), 
       fill = "Max Infection", subtitle = expression(bold("C) No Temperature |  "* bolditalic("Ae. albopictus") *" > abundance"))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.y = element_text(angle = 90, vjust = 0.5),   # Ensures y-axis label is horizontal and aligned
    axis.text.y = element_text(hjust = 1),                   # Adjusts y-axis text alignment
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),                   # Reduce legend text size
    plot.subtitle = element_text(size = 14, face = "bold", hjust = 0),  # Subtitle bold and aligned left
    legend.key.size = unit(0.3, "cm")                      # Reduce legend bar size
  )



### Infected - Temperature - Same abundance



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
  C = 3000,         # Carrying capacity for larvae
  zeta_ae = 0.700 * 0.470, # Transmission coefficient for Aedes aegypti
  zeta_alb = 0.700 * 0.350, # Transmission coefficient for Aedes albopictus
  p_D = 0.700,      # Probability of Dengue transmission from mosquito to human
  BR_ae = 0.470,    # Biting rate for Aedes aegypti 
  BR_alb = 0.350,   # Biting rate for Aedes albopictus 
  theta_h = 1/10,   # Incubation rate of dengue in humans
  gamma_h = 1/7,    # Recovery rate of dengue in humans
  #rho = 1/90,      # Immunity rate
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
time = seq(0, 500, by = 1) # From 0 to 200 days in 1-day steps

# Values for omega_ae and omega_alb
omega_ae = seq(0.1, 1.0, 0.05)       # Competition coefficient of aegypti over albopictus
omega_alb = rev(seq(0.1, 1.0, 0.05)) # Competition coefficient of albopictus over aegypti

infection = function(omega_ae_val, omega_alb_val) {
  parameters["omega_ae"] = omega_ae_val
  parameters["omega_alb"] = omega_alb_val
  
  # Solve the model until equilibrium is reached
  output = ode(y = initial_state, times = time, func = mosquito_model_temp, parms = parameters)
  
  output = as.data.frame(output)
  return(output)
}

# Zeros matrix
infection_matrix = matrix(0, nrow = length(omega_ae), ncol = length(omega_alb))

rownames(infection_matrix) = omega_ae # Set row names as omega_ae
colnames(infection_matrix) = omega_alb # Set column names as omega_alb

for (i in seq_along(omega_ae)) {
  for (j in seq_along(omega_alb)) {
    infection_matrix[i,j] = max(infection(omega_ae[i], omega_alb[j])$I)
  }
}

# Convert the matrix to a data frame for ggplot
infection_df = as.data.frame(as.table(infection_matrix))
colnames(infection_df) = c("omega_ae", "omega_alb", "I")

# Convert columns to numeric
infection_df$omega_ae = as.numeric(as.character(infection_df$omega_ae))
infection_df$omega_alb = as.numeric(as.character(infection_df$omega_alb))

# Restore the original order of omega_alb (from 1 to 0) as a factor
#infection_df$omega_alb = factor(infection_df$omega_alb, levels = sort(unique(infection_df$omega_alb), decreasing = TRUE))

# Create the heatmap with ggplot2
f8 = ggplot(infection_df, aes(x = omega_ae, y = omega_alb, fill = I)) +
  geom_tile(color="black", linewidth = 0.05) +
  scale_fill_viridis(option = "viridis") +  # Use the viridis color palette
  labs(y = expression(omega[alb]~""), 
       x = expression(omega[ae]~""), 
       fill = "Max Infection", subtitle = "B) Temperature | Same abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.y = element_text(angle = 90, vjust = 0.5),   # Ensures y-axis label is horizontal and aligned
    axis.text.y = element_text(hjust = 1),                   # Adjusts y-axis text alignment
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),                   # Reduce legend text size
    plot.subtitle = element_text(size = 14, face = "bold", hjust = 0),  # Subtitle bold and aligned left
    legend.key.size = unit(0.3, "cm")                      # Reduce legend bar size
  )



### Infected - Temperature - *Ae. aegypti* greater abundance



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
  C = 3000,         # Carrying capacity for larvae
  zeta_ae = 0.700 * 0.470, # Transmission coefficient for Aedes aegypti
  zeta_alb = 0.700 * 0.350, # Transmission coefficient for Aedes albopictus
  p_D = 0.700,      # Probability of Dengue transmission from mosquito to human
  BR_ae = 0.470,    # Biting rate for Aedes aegypti 
  BR_alb = 0.350,   # Biting rate for Aedes albopictus 
  theta_h = 1/10,   # Incubation rate of dengue in humans
  gamma_h = 1/7,    # Recovery rate of dengue in humans
  #rho = 1/90,      # Immunity rate
  phi_ae = 0.000,   # Competition coefficient for adult Aedes aegypti
  phi_alb = 0.000,  # Competition coefficient for adult Aedes albopictus
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
time = seq(0, 500, by = 1) # From 0 to 200 days in 1-day steps

# Values for omega_ae and omega_alb
omega_ae = seq(0.1, 1.0, 0.05)       # Competition coefficient of aegypti over albopictus
omega_alb = rev(seq(0.1, 1.0, 0.05)) # Competition coefficient of albopictus over aegypti

infection = function(omega_ae_val, omega_alb_val) {
  parameters["omega_ae"] = omega_ae_val
  parameters["omega_alb"] = omega_alb_val
  
  # Solve the model until equilibrium is reached
  output = ode(y = initial_state, times = time, func = mosquito_model_temp, parms = parameters)
  
  output = as.data.frame(output)
  return(output)
}

# Zeros matrix
infection_matrix = matrix(0, nrow = length(omega_ae), ncol = length(omega_alb))

rownames(infection_matrix) = omega_ae # Set row names as omega_ae
colnames(infection_matrix) = omega_alb # Set column names as omega_alb

for (i in seq_along(omega_ae)) {
  for (j in seq_along(omega_alb)) {
    infection_matrix[i,j] = max(infection(omega_ae[i], omega_alb[j])$I)
  }
}

# Convert the matrix to a data frame for ggplot
infection_df = as.data.frame(as.table(infection_matrix))
colnames(infection_df) = c("omega_ae", "omega_alb", "I")

# Convert columns to numeric
infection_df$omega_ae = as.numeric(as.character(infection_df$omega_ae))
infection_df$omega_alb = as.numeric(as.character(infection_df$omega_alb))

# Restore the original order of omega_alb (from 1 to 0) as a factor
#infection_df$omega_alb = factor(infection_df$omega_alb, levels = sort(unique(infection_df$omega_alb), decreasing = TRUE))

# Create the heatmap with ggplot2
f9 = ggplot(infection_df, aes(x = omega_ae, y = omega_alb, fill = I)) +
  geom_tile(color="black", linewidth = 0.05) +
  scale_fill_viridis(option = "viridis") +  # Use the viridis color palette
  labs(y = expression(omega[alb]~""), 
       x = expression(omega[ae]~""), 
       fill = "Max Infection", subtitle = "B) Temperature") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.y = element_text(angle = 90, vjust = 0.5),   # Ensures y-axis label is horizontal and aligned
    axis.text.y = element_text(hjust = 1),                   # Adjusts y-axis text alignment
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),                   # Reduce legend text size
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0),  # Subtitle bold and aligned left
    legend.key.size = unit(0.3, "cm")                      # Reduce legend bar size
  )



### Infected - Temperature - *Ae. albopictus* greater abundance



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
  C = 3000,         # Carrying capacity for larvae
  zeta_ae = 0.700 * 0.470, # Transmission coefficient for Aedes aegypti
  zeta_alb = 0.700 * 0.350, # Transmission coefficient for Aedes albopictus
  p_D = 0.700,      # Probability of Dengue transmission from mosquito to human
  BR_ae = 0.470,    # Biting rate for Aedes aegypti 
  BR_alb = 0.350,   # Biting rate for Aedes albopictus 
  theta_h = 1/10,   # Incubation rate of dengue in humans
  gamma_h = 1/7,    # Recovery rate of dengue in humans
  #rho = 1/90,      # Immunity rate
  phi_ae = 0.000,   # Competition coefficient for adult Aedes aegypti
  phi_alb = 0.000,  # Competition coefficient for adult Aedes albopictus
  C_ad = 10000      # Carrying capacity for adult mosquitoes
)

# Initial conditions (commented)
initial_state = c(
  L_ae = 50,   # Initial larvae population for Aedes aegypti
  L_alb = 200, # Initial larvae population for Aedes albopictus
  X_ae = 1000, # Initial adult susceptible population for Aedes aegypti
  X_alb = 3000,# Initial adult susceptible population for Aedes albopictus
  Y_ae = 0,    # Initial adult infected population for Aedes aegypti
  Y_alb = 0,   # Initial adult infected population for Aedes albopictus
  S = 999,     # Initial susceptible human population
  E = 1,       # Initial exposed human population
  I = 0,       # Initial infected human population
  R = 0        # Initial recovered human population
)

# Time vector
time = seq(0, 500, by = 1) # From 0 to 200 days in 1-day steps

# Values for omega_ae and omega_alb
omega_ae = seq(0.1, 1.0, 0.05)       # Competition coefficient of aegypti over albopictus
omega_alb = rev(seq(0.1, 1.0, 0.05)) # Competition coefficient of albopictus over aegypti

infection = function(omega_ae_val, omega_alb_val) {
  parameters["omega_ae"] = omega_ae_val
  parameters["omega_alb"] = omega_alb_val
  
  # Solve the model until equilibrium is reached
  output = ode(y = initial_state, times = time, func = mosquito_model_temp, parms = parameters)
  
  output = as.data.frame(output)
  return(output)
}

# Zeros matrix
infection_matrix = matrix(0, nrow = length(omega_ae), ncol = length(omega_alb))

rownames(infection_matrix) = omega_ae # Set row names as omega_ae
colnames(infection_matrix) = omega_alb # Set column names as omega_alb

for (i in seq_along(omega_ae)) {
  for (j in seq_along(omega_alb)) {
    infection_matrix[i,j] = max(infection(omega_ae[i], omega_alb[j])$I)
  }
}

# Convert the matrix to a data frame for ggplot
infection_df = as.data.frame(as.table(infection_matrix))
colnames(infection_df) = c("omega_ae", "omega_alb", "I")

# Convert columns to numeric
infection_df$omega_ae = as.numeric(as.character(infection_df$omega_ae))
infection_df$omega_alb = as.numeric(as.character(infection_df$omega_alb))

# Restore the original order of omega_alb (from 1 to 0) as a factor
#infection_df$omega_alb = factor(infection_df$omega_alb, levels = sort(unique(infection_df$omega_alb), decreasing = TRUE))

# Create the heatmap with ggplot2
f10 = ggplot(infection_df, aes(x = omega_ae, y = omega_alb, fill = I)) +
  geom_tile(color="black", linewidth = 0.05) +
  scale_fill_viridis(option = "viridis") +  # Use the viridis color palette
  labs(y = expression(omega[alb]~""), 
       x = expression(omega[ae]~""), 
       fill = "Max Infection", subtitle = expression(bold("D) Temperature | "* bolditalic("Ae. albopictus") *" > abundance"))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.y = element_text(angle = 90, vjust = 0.5),   # Ensures y-axis label is horizontal and aligned
    axis.text.y = element_text(hjust = 1),                   # Adjusts y-axis text alignment
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),                   # Reduce legend text size
    plot.subtitle = element_text(size = 14, face = "bold", hjust = 0),  # Subtitle bold and aligned left
    legend.key.size = unit(0.3, "cm")                      # Reduce legend bar size
  )


### Plots

plot_grid(f6,f9,
          ncol = 2)

plot_grid(f5,f8, f7,f10,
          ncol = 2)
