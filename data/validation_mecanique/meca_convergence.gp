# --- Configuration de la sortie ---
set terminal pngcairo size 1200,600 enhanced font 'Arial,12'
# Le fichier image sera enregistré dans le dossier data/validation
set output '../data/validation_mecanique/meca_validation.png'

# --- Format des données ---
set datafile separator "\t"

# --- Configuration de l'affichage en deux panneaux ---
set multiplot layout 1,2 title "Validation Mecanique - Paroi de Domaine 180°" font ",16" offset 0,-1

# ==========================================
# Graphique 1 : Profils de Polarisation
# ==========================================
set title "Comparaison Analytique vs Numérique" font ",14"
set xlabel "Position x"
set ylabel "Déplacement"
set grid
set key bottom right box

# Tracé avec chemins relatifs depuis "build/"
plot '../data/validation_mecanique/meca_convergence.dat' skip 1 using 1:2 with lines lw 2 lc rgb "black" title "Solution Analytique", \
     '../data/validation_mecanique/meca_convergence.dat' skip 1 using 1:3 with points pt 7 ps 0.8 lc rgb "red" title "Solution Numérique"

# ==========================================
# Graphique 2 : Erreur Locale
# ==========================================
set title "Erreur Spatiale Locale" font ",14"
set xlabel "Position x"
set ylabel "Différence (u_{exact} - u_{num})"
set grid
set autoscale y
# Sécurité mathématique contre un axe Y vide (si erreur parfaitement nulle)
set offsets 0, 0, 1e-4, 1e-4 
set key top right box

# Tracé avec chemin relatif depuis "build/"
plot '../data/validation_mecanique/meca_convergence.dat' skip 1 using 1:4 with lines lw 2 lc rgb "blue" title "Erreur absolue"

unset multiplot

# --- Message de fin dans le terminal ---
print "Graphique généré avec succès : ../data/validation_mecanique/meca_validation.png"