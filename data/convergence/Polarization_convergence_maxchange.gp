# --- Configuration de la sortie ---
set terminal pngcairo size 1000, 600 enhanced font 'Verdana,12'
set output '../data/convergence/Polarization_convergence_maxchange.png'

# --- Titres et Légendes ---
set title "Evolution de la convergence (Variation Max de la Polarisation)"
set xlabel "Nombre de pas de temps (Steps)"
set ylabel "Variation Max (Delta P)"

# --- Échelle logarithmique ---
# Très important car ta variation varie de 10^-6 à 10^-4
set logscale y

# --- Grille et style ---
set grid
set key right top

# --- Paramètres de tracé ---

# On ignore la première ligne (l'en-tête "step,max_variation")
plot "../data/convergence/Polarization_convergence_maxchange.dat" skip 1 using 1:2 with lines lw 2 lc rgb "blue" title "Delta P"

# --- Message de fin ---
print "Graphique généré : convergence_plot.png"