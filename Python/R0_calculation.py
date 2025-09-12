import sympy as sp

sp.init_printing()

# Variables
S, E, I, x_ae, x_alb, y_ae, y_alb, N = sp.symbols('S E I x_ae x_alb y_ae y_alb N_h')

# Parameters
pd, br_ae, br_alb, zeta_ae, zeta_alb, theta, gamma, mu_Aae, mu_Aalb, K_ad= sp.symbols('p_D BR_ae BR_alb zeta_ae zeta_alb theta gamma mu_Aae mu_Aalb K_ad')

zeta_ae = pd * br_ae
zeta_alb = pd * br_alb

# --- Vector F(x): new infections ---
F = sp.Matrix([
    (pd * ((br_ae*y_ae) + (br_alb*y_alb)) * S / N),
    0,
    (zeta_ae * x_ae * I / N) * (1 - (y_ae / K_ad)), 
    (zeta_alb * x_alb * I / N )* (1 - (y_alb / K_ad))
])

# --- Vector V(x): transitions ---
V = sp.Matrix([
    theta * E,
    (-theta * E) + (gamma * I),
    mu_Aae * y_ae,
    mu_Aalb * y_alb
])

# --- Jacobians ---
infected_vars = [E, I, y_ae, y_alb]
Fx = F.jacobian(infected_vars)
Vx = V.jacobian(infected_vars)

# DFE
Fx_dfe = Fx.subs({E: 0, I: 0, y_ae: 0, y_alb: 0, S: N})
Vx_dfe = Vx.subs({E: 0, I: 0, y_ae: 0, y_alb: 0, S: N})

# --- NGM ---
K = Fx_dfe * Vx_dfe.inv()

# --- R0 ---
eigs = K.eigenvals()
print(len(eigs))
sp.print_latex(eigs)
sp.print_latex(K)

for ev in K.eigenvals().keys():
    sp.print_latex(sp.simplify(ev))
