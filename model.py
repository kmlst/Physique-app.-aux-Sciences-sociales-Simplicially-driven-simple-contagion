import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from tqdm import tqdm

# Le syst diffs
def f(rho, t, lambda_a, lambda_b, lambda_a_delta, epsilon_ab):
    rho_a_tot, rho_b_tot, rho_ab = rho
    f1 = rho_a_tot*(-1 + lambda_a*(1 - rho_a_tot) + lambda_a_delta*rho_a_tot*(1 - rho_a_tot))
    f2 = rho_b_tot*(-1 + lambda_b*(1 - rho_b_tot) + lambda_b*(epsilon_ab - 1)*(rho_a_tot - rho_ab))
    f3 = -2*rho_ab + epsilon_ab*lambda_b*(rho_a_tot - rho_ab)*rho_b_tot \
         + lambda_a*(rho_b_tot - rho_ab)*rho_a_tot \
         + lambda_a_delta*(rho_b_tot - rho_ab)*rho_a_tot**2
    return [f1, f2, f3]



lambda_b = 0.8
lambda_a_delta = 2.5

rho0 = [0.02, 0.02, 0.01]

# Le truc qui résoud
def solve_system(lambda_a, lambda_b, lambda_a_delta, epsilon_ab):
    t = np.linspace(0, 500, 10000)
    sol = odeint(f, rho0, t, args=(lambda_a, lambda_b, lambda_a_delta, epsilon_ab))
    rho_b_tot_asymptote = sol[-1, 1]
    return rho_b_tot_asymptote

# Différenets valeurs de eps
def graph_type1() :

    epsilon_ab_vals = [1.4, 1.5, 1.75]
    lambda_a_vals = np.linspace(0.5, 1.25, 100)


    for epsilon_ab in epsilon_ab_vals:
        rho_b_tot_asymptote_vals = [solve_system(lambda_a, lambda_b, lambda_a_delta, epsilon_ab) for lambda_a in lambda_a_vals]
        plt.plot(lambda_a_vals, rho_b_tot_asymptote_vals, label=f'epsilon_ab = {epsilon_ab}')

    plt.xlabel('lambda_a')
    plt.ylabel('rho_b_tot_asymptote')
    plt.legend()


    plt.show()


def plot_rho_lambda_epsilon(lambda_a_values, epsilon_ab_values):
    # Initialiser la grille de points
    lambda_a_grid, epsilon_ab_grid = np.meshgrid(lambda_a_values, epsilon_ab_values)
    rho_b_tot_asymptote_grid = np.zeros_like(lambda_a_grid)

    # Calculer la valeur de rho_b_tot_asymptote pour chaque point de la grille
    for i, epsilon_ab in tqdm(enumerate(epsilon_ab_values)):
        for j, lambda_a in enumerate(lambda_a_values):
            # Résolution du système d'équations différentielles
            rho_b_tot_asymptote  = solve_system(lambda_a, lambda_b, lambda_a_delta, epsilon_ab)

            # Asymptotes des densités

            rho_b_tot_asymptote_grid[i, j] = rho_b_tot_asymptote
        print(i)
    # Représenter graphiquement la valeur de rho_b_tot_asymptote pour chaque point de la grille en utilisant une carte de couleurs
    plt.pcolormesh(lambda_a_grid, epsilon_ab_grid, np.max(rho_b_tot_asymptote_grid)-rho_b_tot_asymptote_grid, cmap='Blues')
    plt.colorbar(label=r'$\rho_{B,tot}^{asymptote}$')
    plt.xlabel(r'$\lambda_A$')
    plt.ylabel(r'$\epsilon_{AB}$')
    plt.show()


def graph_type2() :

    lambda_a_values = np.linspace(0.5, 1.25, 100)
    epsilon_ab_values = np.linspace(1.3, 1.6, 100)
    plot_rho_lambda_epsilon(lambda_a_values, epsilon_ab_values)


graph_type2()


