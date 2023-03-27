import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Le système pour odeint
    rho_a_tot, rho_b_tot, rho_ab = rho
    f1 = rho_a_tot*(-1 + lambda_a*(1 - rho_a_tot) + lambda_a_delta*rho_a_tot*(1 - rho_a_tot))
    f2 = rho_b_tot*(-1 + lambda_b*(1 - rho_b_tot) + lambda_b*(epsilon_ab - 1)*(rho_a_tot - rho_ab))
    f3 = -2*rho_ab + epsilon_ab*lambda_b*(rho_a_tot - rho_ab)*rho_b_tot \
         + lambda_a*(rho_b_tot - rho_ab)*rho_a_tot \
         + lambda_a_delta*(rho_b_tot - rho_ab)*rho_a_tot**2
    return [f1, f2, f3]

# Les paramètres ( on a lambda qui varie pour tracer le graphe)
lambda_a_vals = np.linspace(0, 3, 100)
lambda_b = 0.8
lambda_a_delta = 0
epsilon_ab_vals = [1.5, 2, 4.5]

# Ca ça résoud
def solve_system(lambda_a, lambda_b, lambda_a_delta, epsilon_ab):
    rho0 = [0.1, 0.1, 0.1]
    t = np.linspace(0, 500, 10000)
    sol = odeint(f, rho0, t, args=(lambda_a, lambda_b, lambda_a_delta, epsilon_ab))
    rho_b_tot_asymptote = sol[-1, 1]
    return rho_b_tot_asymptote

# La on trace
for epsilon_ab in epsilon_ab_vals:
    rho_b_tot_asymptote_vals = [solve_system(lambda_a, lambda_b, lambda_a_delta, epsilon_ab) for lambda_a in lambda_a_vals]
    plt.plot(lambda_a_vals, rho_b_tot_asymptote_vals, label=f'epsilon_ab = {epsilon_ab}')


plt.xlabel('lambda_a')
plt.ylabel('rho_b_tot_asymptote')
plt.legend()


plt.show()

