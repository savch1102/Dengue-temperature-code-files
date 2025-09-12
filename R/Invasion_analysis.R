library(deSolve, warn.conflicts=F, quietly=T)
library(cowplot, warn.conflicts=F, quietly=T)
library(ggplot2, warn.conflicts=F, quietly=T)
library(gridExtra, warn.conflicts=F, quietly=T)
library(reshape2, warn.conflicts=F, quietly=T)
library(viridis, warn.conflicts=F, quietly=T)
library(tidyverse, warn.conflicts=F, quietly=T)

##------------------ No temperature

# Define the system of differential equations (mosquito dynamics)
mosquito_model = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Larval dynamics
    dL_ae = beta_ae * L_ae * (1 - (L_ae + omega_alb * L_alb) / C) - kappa_ae * L_ae - mu_Lae * L_ae  # Change in Aedes aegypti larvae
    dL_alb = beta_alb * L_alb * (1 - (L_alb + omega_ae * L_ae) / C) - kappa_alb * L_alb - mu_Lalb * L_alb  # Change in Aedes albopictus larvae
    
    # Adult susceptible dynamics
    dX_ae = kappa_ae * L_ae * (1 - (X_ae + phi_alb * X_alb) / C_ad) - mu_Aae * X_ae  # Change in susceptible Aedes aegypti adults
    dX_alb = kappa_alb * L_alb * (1 - (X_alb + phi_ae * X_ae) / C_ad) - mu_Aalb * X_alb  # Change in susceptible Aedes albopictus adults
    
    # Return derivatives (without infected mosquitoes)
    list(c(dL_ae, dL_alb, dX_ae, dX_alb))
  })
}

# Parameters (excluding infected mosquitoes)
parameters = c(
  beta_ae = 0.340,  # Reproductive rate of Aedes aegypti
  beta_alb = 0.253, # Reproductive rate of Aedes albopictus
  C = 5000,         # Larval carrying capacity
  kappa_ae = 0.170, # Larval development rate of Aedes aegypti
  kappa_alb = 0.167,# Larval development rate of Aedes albopictus
  mu_Lae = 0.090,   # Larval mortality rate of Aedes aegypti
  mu_Lalb = 0.050,  # Larval mortality rate of Aedes albopictus
  mu_Aae = 1/10,    # Adult mortality rate of Aedes aegypti (1/lifespan_ae)
  mu_Aalb = 1/12,   # Adult mortality rate of Aedes albopictus (1/lifespan_alb)
  phi_ae = 0.000,   # Competition coefficient for Aedes aegypti adults
  phi_alb = 0.000,  # Competition coefficient for Aedes albopictus adults
  C_ad = 10000      # Adult carrying capacity
)

# Initial conditions (without infected mosquitoes)
initial_state = c(
  L_ae = 1000,  # Initial population of Aedes aegypti larvae
  L_alb = 50,   # Initial population of Aedes albopictus larvae
  X_ae = 2000,  # Initial population of susceptible Aedes aegypti adults
  X_alb = 100   # Initial population of susceptible Aedes albopictus adults
)

# Time vector
time = seq(0, 200, by = 1) # From 0 to 200 days with 1-day steps

# Values for omega_ae and omega_alb (competition coefficients)
omega_ae = seq(0.1, 1.0, 0.05)       # Competition coefficient of aegypti on albopictus
omega_alb = rev(seq(0.1, 1.0, 0.05)) # Competition coefficient of albopictus on aegypti (reverse order)

# Function to calculate invasion exponent
calculate_invasion_exponent = function(omega_ae, omega_alb) {
  # Update parameters with current omega_ae and omega_alb values
  parameters["omega_ae"] = omega_ae
  parameters["omega_alb"] = omega_alb
  
  # Solve the system of differential equations
  out = ode(y = initial_state, times = time, func = mosquito_model, parms = parameters)
  
  # Calculate growth rate of invading population (albopictus)
  final_L_alb = tail(out[, "L_alb"], 1)  # Final larval population of albopictus
  initial_L_alb = initial_state["L_alb"] # Initial larval population of albopictus
  invasion_exponent = log(final_L_alb / initial_L_alb) / max(time)  # Invasion exponent calculation
  
  return(invasion_exponent)
}

# Create a matrix to store PIP (Pairwise Invasibility Plot) results
pip_matrix = matrix(NA, nrow = length(omega_ae), ncol = length(omega_alb))

# Fill matrix with invasion exponent values
for (i in 1:length(omega_ae)) {
  for (j in 1:length(omega_alb)) {
    pip_matrix[i, j] = calculate_invasion_exponent(omega_ae[i], omega_alb[j])
  }
}

# Create data frame for PIP plotting
pip_data = expand.grid(omega_ae = omega_ae, omega_alb = omega_alb)
pip_data$invasion_exponent = as.vector(pip_matrix)  # Flatten matrix into vector

# Plot PIP using ggplot2
p1 = ggplot(pip_data, aes(x = omega_ae, y = omega_alb, fill = invasion_exponent)) +
  geom_tile(color="black", linewidth = 0.1) +  # Create colored tiles with thin borders
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +  # Color gradient
  labs(subtitle = "A) No Temperature",
       y = expression(omega[alb]~"(Invader)"),  # y-axis label with mathematical notation
       x = expression(omega[ae]~"(Resident)"),  # x-axis label with mathematical notation
       fill = expression("Growth"~rho)) +  # Legend title
  theme_minimal() +  # Minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.y = element_text(angle = 90, vjust = 0.5),  # Ensure y-axis label is vertical
    axis.text.y = element_text(hjust = 1),  # Adjust y-axis text alignment
    axis.title = element_text(size = 11),  # Axis title size
    legend.title = element_text(size = 10),  # Legend title size
    legend.text = element_text(size = 8),  # Legend text size
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0),  # Subtitle formatting
    legend.key.size = unit(0.5, "cm")  # Legend color bar size
  )

##------------------ Temperature

# Define the system of differential equations with temperature dependence
mosquito_model_temp = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Temperature function (seasonal variation)
    T = 2.7 * sin(-0.04 * t) + 29.8  # Sinusoidal temperature pattern with mean ~30°C
    
    # Development rate functions (Gaussian curves centered at 33°C)
    kappa_ae = (1 / (2.7 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (2.7^2)))  # Aedes aegypti development rate
    kappa_alb = (1 / (2.4 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (2.4^2)))  # Aedes albopictus development rate
    
    # Biting rate functions (Gaussian curves centered at 33°C)
    br_ae = (1 / (1.1 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (1.1^2)))  # Aedes aegypti biting rate
    br_alb = (1 / (1.4 * sqrt(2 * pi))) * exp(-((T - 33)^2) / (2 * (1.4^2)))  # Aedes albopictus biting rate
    
    # Mortality rate functions (quadratic functions of temperature)
    mu_Lae = (0.05 * (T - 25.5))^2 + 0.09    # Aedes aegypti larval mortality
    mu_Lalb = (0.05 * (T - 25.5))^2 + 0.05   # Aedes albopictus larval mortality
    mu_Aae = (0.05 * (T - 25.5))^2 + 0.015   # Aedes aegypti adult mortality
    mu_Aalb = (0.05 * (T - 25.5))^2 + 0.01   # Aedes albopictus adult mortality
    
    # Fecundity rate functions (Gaussian curves centered at 27°C)
    beta_ae = 70 * (1 / (3.5 * sqrt(2 * pi))) * exp(-((T - 27)^2) / (2 * (3.5^2)))  # Aedes aegypti fecundity
    beta_alb = 64 * (1 / (3.5 * sqrt(2 * pi))) * exp(-((T - 27)^2) / (2 * (3.5^2)))  # Aedes albopictus fecundity
    
    # Larval dynamics
    dL_ae = beta_ae * L_ae * (1 - (L_ae + omega_alb * L_alb) / C) - kappa_ae * L_ae - mu_Lae * L_ae  # Change in Aedes aegypti larvae
    dL_alb = beta_alb * L_alb * (1 - (L_alb + omega_ae * L_ae) / C) - kappa_alb * L_alb - mu_Lalb * L_alb  # Change in Aedes albopictus larvae
    
    # Adult susceptible dynamics (no infected mosquitoes)
    dX_ae = kappa_ae * L_ae * (1 - (X_ae + phi_alb * X_alb) / C_ad) - mu_Aae * X_ae  # Change in susceptible Aedes aegypti adults
    dX_alb = kappa_alb * L_alb * (1 - (X_alb + phi_ae * X_ae) / C_ad) - mu_Aalb * X_alb  # Change in susceptible Aedes albopictus adults
    
    # Return derivatives
    list(c(dL_ae, dL_alb, dX_ae, dX_alb))
  })
}

# Parameters (temperature-independent parameters)
parameters = c(
  C = 5000,         # Larval carrying capacity
  phi_ae = 0.000,   # Competition coefficient for Aedes aegypti adults
  phi_alb = 0.000,  # Competition coefficient for Aedes albopictus adults
  C_ad = 10000      # Adult carrying capacity
)

# Initial conditions (no infected mosquitoes)
initial_state = c(
  L_ae = 1000,  # Initial population of Aedes aegypti larvae
  L_alb = 10,   # Initial population of Aedes albopictus larvae
  X_ae = 2000,  # Initial population of susceptible Aedes aegypti adults
  X_alb = 50    # Initial population of susceptible Aedes albopictus adults
)

# Time vector
time = seq(0, 1000, by = 1) # From 0 to 1000 days with 1-day steps

# Values for competition coefficients
omega_ae = seq(0.1, 1.0, 0.05)       # Competition coefficient of aegypti over albopictus (0.1 to 1.0 in 0.05 increments)
omega_alb = rev(seq(0.1, 1.0, 0.05)) # Competition coefficient of albopictus over aegypti (reverse order)

# Function to calculate invasion exponent
calculate_invasion_exponent = function(omega_ae, omega_alb) {
  # Update parameters with current omega values
  parameters["omega_ae"] = omega_ae
  parameters["omega_alb"] = omega_alb
  
  # Solve the system of differential equations
  out = ode(y = initial_state, times = time, func = mosquito_model_temp, parms = parameters)
  
  # Calculate growth rate of invading population (albopictus)
  final_L_alb = tail(out[, "L_alb"], 1)  # Final larval population of albopictus
  initial_L_alb = initial_state["L_alb"] # Initial larval population of albopictus
  invasion_exponent = log(final_L_alb / initial_L_alb) / max(time)  # Per-capita growth rate
  
  return(invasion_exponent)
}

# Create matrix to store Pairwise Invasibility Plot (PIP) results
pip_matrix = matrix(NA, nrow = length(omega_ae), ncol = length(omega_alb))

# Fill matrix with invasion exponent values
for (i in 1:length(omega_ae)) {
  for (j in 1:length(omega_alb)) {
    pip_matrix[i, j] = calculate_invasion_exponent(omega_ae[i], omega_alb[j])
  }
}

# Create data frame for plotting PIP
pip_data = expand.grid(omega_ae = omega_ae, omega_alb = omega_alb)
pip_data$invasion_exponent = as.vector(pip_matrix)  # Flatten matrix into vector

# Plot PIP using ggplot2
p2 = ggplot(pip_data, aes(x = omega_ae, y = omega_alb, fill = invasion_exponent)) +
  geom_tile(color="black", linewidth = 0.1) +  # Create heatmap tiles
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +  # Color scale
  labs(subtitle = "B) Temperature",
       y = expression(omega[alb]~"(Invader)"),  # y-axis label with math notation
       x = expression(omega[ae]~"(Resident)"),  # x-axis label with math notation
       fill = expression("Growth"~rho)) +  # Legend title
  theme_minimal() +  # Clean minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.y = element_text(angle = 90, vjust = 0.5),  # Vertical y-axis label
    axis.text.y = element_text(hjust = 1),  # Adjust y-axis text alignment
    axis.title = element_text(size = 11),  # Axis title size
    legend.title = element_text(size = 10),  # Legend title size
    legend.text = element_text(size = 8),  # Legend text size
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0),  # Subtitle formatting
    legend.key.size = unit(0.5, "cm")  # Legend color bar size
  )


plot_grid(p1,p2,
          ncol = 2)


