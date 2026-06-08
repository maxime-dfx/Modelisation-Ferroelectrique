import matplotlib.pyplot as plt
import numpy as np 

def G(T, P, G0, Tc, alpha, beta, gamma, E):
    return G0 + 0.5 * alpha * (T - Tc) * P**2 + 0.25 * beta * P**4 + (1/6) * gamma * P**6 - E * P

# Paramètres thermodynamiques
T = 0
Tc = 80
alpha = 1
beta = 1
gamma = 0
E = 0
G0 = 0

P_array = np.linspace(-15, 15, 400)
G_array = G(T, P_array, G0, Tc, alpha, beta, gamma, E)

# Configuration et tracé
plt.figure(figsize=(8, 5))
plt.plot(P_array, G_array, label=f'T = {T} K (Phase Ferroélectrique)', color='b', linewidth=2)

# Ajout des axes (x=0 et y=0) pour la lisibilité
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')

# Esthétique du graphique
plt.xlabel('Polarisation $P$')
plt.ylabel('Énergie Libre $G(P)$')
plt.title('Modèle de Landau-Devonshire : Double puits de potentiel')
plt.grid(True, alpha=0.3)
plt.legend()

# Affichage
plt.show()