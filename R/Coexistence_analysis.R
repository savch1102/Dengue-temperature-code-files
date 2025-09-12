library(deSolve, warn.conflicts=F, quietly=T)
library(cowplot, warn.conflicts=F, quietly=T)
library(ggplot2, warn.conflicts=F, quietly=T)
library(gridExtra, warn.conflicts=F, quietly=T)
library(reshape2, warn.conflicts=F, quietly=T)
library(viridis, warn.conflicts=F, quietly=T)
library(tidyverse, warn.conflicts=F, quietly=T)

##----------------------- No Temperature

# Define the system of differential equations (mosquito dynamics only)
mosquito_model = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Larval dynamics
    dL_ae = beta_ae * L_ae * (1 - (L_ae + omega_alb * L_alb) / C) - kappa_ae * L_ae - mu_Lae * L_ae
    dL_alb = beta_alb * L_alb * (1 - (L_alb + omega_ae * L_ae) / C) - kappa_alb * L_alb - mu_Lalb * L_alb
    
    # Adult susceptible dynamics
    dX_ae = kappa_ae * L_ae * (1 - (X_ae + phi_alb * X_alb) / C_ad) - mu_Aae * X_ae
    dX_alb = kappa_alb * L_alb * (1 - (X_alb + phi_ae * X_ae) / C_ad) - mu_Aalb * X_alb
    
    # Return derivatives
    list(c(dL_ae, dL_alb, dX_ae, dX_alb))
  })
}

# Parameters
parameters = c(
  beta_ae = 0.340,  beta_alb = 0.253,  C = 5000,
  kappa_ae = 0.170, kappa_alb = 0.167,
  mu_Lae = 0.090,  mu_Lalb = 0.050,
  mu_Aae = 1/10,   mu_Aalb = 1/12,
  phi_ae = 0.000,  phi_alb = 0.000,
  C_ad = 10000
)

# Initial conditions
initial_state = c(L_ae = 1000, L_alb = 50, X_ae = 2000, X_alb = 100)

# Time vector
time = seq(0, 200, by = 1)

# Values for omega_ae and omega_alb
omega_ae_vals = seq(0.1, 1.0, 0.05)
omega_alb_vals = rev(seq(0.1, 1.0, 0.05))

# Dataframe to store larval abundances and percentages
results_df = data.frame(
  omega_ae = numeric(),
  omega_alb = numeric(),
  L_ae = numeric(),
  L_alb = numeric(),
  total = numeric(),
  perc_ae = numeric(),
  perc_alb = numeric()
)

# Run simulations and collect results
for (omega_ae_val in omega_ae_vals) {
  for (omega_alb_val in omega_alb_vals) {
    parameters["omega_ae"] = omega_ae_val
    parameters["omega_alb"] = omega_alb_val
    
    # Solve the model
    output = ode(y = initial_state, times = time, func = mosquito_model, parms = parameters)
    output_df = as.data.frame(output)
    
    # Calculate final larval abundances and percentages
    L_ae_final = output_df[nrow(output_df), "L_ae"]
    L_alb_final = output_df[nrow(output_df), "L_alb"]
    total = L_ae_final + L_alb_final
    perc_ae = ifelse(total > 0, L_ae_final / total * 100, 0)
    perc_alb = ifelse(total > 0, L_alb_final / total * 100, 0)
    
    # Append results to the dataframe
    results_df = rbind(results_df, data.frame(
      omega_ae = omega_ae_val,
      omega_alb = omega_alb_val,
      L_ae = L_ae_final,
      L_alb = L_alb_final,
      total = total,
      perc_ae = perc_ae,
      perc_alb = perc_alb
    ))
  }
}

# Melt the dataframe to long format for plotting
library(reshape2)
results_long = melt(results_df, id.vars = c("omega_ae", "omega_alb"), 
                    measure.vars = c("perc_ae", "perc_alb"),
                    variable.name = "Species", value.name = "Percentage")

results_long$Species = factor(results_long$Species, 
                              levels = c("perc_ae", "perc_alb"), 
                              labels = c("Ae. aegypti", "Ae. albopictus"))

# Create violin plot with jitter
f3 = ggplot(results_long, aes(x = Species, y = Percentage, fill = Species, color = Species)) +
  geom_violin(alpha = 0.4, trim = TRUE) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.7, color = "grey1") +
  stat_summary(fun = mean, geom = "point", size = 2)+
  stat_summary(fun = mean, geom = "text", aes(label = "Mean"), 
               hjust = -0.5, color = "black", size = 4) +
  scale_fill_manual(values = c("Ae. aegypti" = "white", "Ae. albopictus" = "black")) +
  scale_color_manual(values = c("Ae. aegypti" = "black", "Ae. albopictus" = "black")) +
  labs(x = "Species", y = "Percentage (%)", subtitle = "A) No Temperature | Larval %") +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "italic", size = 12),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(face = "italic", size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.subtitle = element_text(size = 14, face = "bold", hjust = 0))

results_long_notemp = results_long

##----------------------- Temperature

# Define the system of differential equations
mosquito_model_temp = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
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
    
    # Larval dynamics
    dL_ae = beta_ae * L_ae * (1 - (L_ae + omega_alb * L_alb) / C) - kappa_ae * L_ae - mu_Lae * L_ae
    dL_alb = beta_alb * L_alb * (1 - (L_alb + omega_ae * L_ae) / C) - kappa_alb * L_alb - mu_Lalb * L_alb
    
    # Adult susceptible dynamics (no infected mosquitoes)
    dX_ae = kappa_ae * L_ae * (1 - (X_ae + phi_alb * X_alb) / C_ad) - mu_Aae * X_ae
    dX_alb = kappa_alb * L_alb * (1 - (X_alb + phi_ae * X_ae) / C_ad) - mu_Aalb * X_alb
    
    # Return derivatives
    list(c(dL_ae, dL_alb, dX_ae, dX_alb))
  })
}

# Parameters (updated with temperature-dependent functions)
parameters = c(
  C = 5000,         # Carrying capacity for larvae
  phi_ae = 0.00,   # Competition coefficient for Aedes aegypti adults
  phi_alb = 0.00,  # Competition coefficient for Aedes albopictus adults
  C_ad = 10000      # Carrying capacity for adult mosquitoes
)

# Initial conditions (no infected mosquitoes)
initial_state = c(
  L_ae = 1000,  # Initial larvae population of Aedes aegypti
  L_alb = 50,   # Initial larvae population of Aedes albopictus
  X_ae = 2000,  # Initial susceptible adult population of Aedes aegypti
  X_alb = 100   # Initial susceptible adult population of Aedes albopictus
)

# Time vector
time = seq(0, 200, by = 1) # From 0 to 200 days with steps of 1 day

# Values for omega_ae and omega_alb
omega_ae_vals = seq(0.1, 1.0, 0.05)
omega_alb_vals = rev(seq(0.1, 1.0, 0.05))

# Dataframe to store larval abundances and percentages
results_df = data.frame(
  omega_ae = numeric(),
  omega_alb = numeric(),
  L_ae = numeric(),
  L_alb = numeric(),
  total = numeric(),
  perc_ae = numeric(),
  perc_alb = numeric()
)

# Run simulations and collect results
for (omega_ae_val in omega_ae_vals) {
  for (omega_alb_val in omega_alb_vals) {
    parameters["omega_ae"] = omega_ae_val
    parameters["omega_alb"] = omega_alb_val
    
    # Solve the model
    output = ode(y = initial_state, times = time, func = mosquito_model_temp, parms = parameters)
    output_df = as.data.frame(output)
    
    # Calculate final larval abundances and percentages
    L_ae_final = output_df[nrow(output_df), "L_ae"]
    L_alb_final = output_df[nrow(output_df), "L_alb"]
    total = L_ae_final + L_alb_final
    perc_ae = ifelse(total > 0, L_ae_final / total * 100, 0)
    perc_alb = ifelse(total > 0, L_alb_final / total * 100, 0)
    
    # Append results to the dataframe
    results_df = rbind(results_df, data.frame(
      omega_ae = omega_ae_val,
      omega_alb = omega_alb_val,
      L_ae = L_ae_final,
      L_alb = L_alb_final,
      total = total,
      perc_ae = perc_ae,
      perc_alb = perc_alb
    ))
  }
}

# Melt the dataframe to long format for plotting
library(reshape2)
results_long = melt(results_df, id.vars = c("omega_ae", "omega_alb"), 
                    measure.vars = c("perc_ae", "perc_alb"),
                    variable.name = "Species", value.name = "Percentage")

results_long$Species = factor(results_long$Species, 
                              levels = c("perc_ae", "perc_alb"), 
                              labels = c("Ae. aegypti", "Ae. albopictus"))


# Create violin plot with jitter
f4 = ggplot(results_long, aes(x = Species, y = Percentage, fill = Species, color = Species)) +
  geom_violin(alpha = 0.4, trim = TRUE) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.7, color = "grey1") +
  stat_summary(fun = mean, geom = "point", size = 2)+
  stat_summary(fun = mean, geom = "text", aes(label = "Mean"), 
               hjust = -0.45, color = "black", size = 4) +
  scale_fill_manual(values = c("Ae. aegypti" = "white", "Ae. albopictus" = "black")) +
  scale_color_manual(values = c("Ae. aegypti" = "black", "Ae. albopictus" = "black")) +
  labs(x = "Species", y = "Percentage (%)", subtitle = "B) Temperature | Larval %") +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "italic", size = 12),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(face = "italic", size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.subtitle = element_text(size = 14, face = "bold", hjust = 0))

results_long_temp = results_long

### Stats

results_long_temp$model = "Temp"
results_long_notemp$model = "NoTemp"

combined_data = rbind(results_long_temp, results_long_notemp)

# No temperature
wilcox.test(Percentage ~ Species, data = subset(results_long_notemp))

# Temperature
wilcox.test(Percentage ~ Species, data = subset(results_long_temp))

### Plots

a = plot_grid(f3, f4)
a
