source("fonction_rk_mixte.R")

#Initial parameters
library(truncnorm)
library(animation)

# Setup parallel backend
library(foreach)

# booleen qui detecte l'environnement: machine locale ou cluster
is_cluster <- function() {
  cluster_vars <- c("PBS_JOBID", "SLURM_JOB_ID", "LSB_JOBID", "SGE_TASK_ID")
  any(sapply(cluster_vars, function(var) Sys.getenv(var) != ""))
}

####################
if (is_cluster()) { # we're on a cluster
  #doMPI initialization
  require(doMPI)
  cl <- startMPIcluster()
  registerDoMPI(cl)
} else {  #we're on a local machine assuming it is multicore
  require(doParallel)
  
  cl <- makeCluster(1) 
  registerDoParallel(cl)
}




##################### Paramètres
nb_ressources<-250 #nombre de points attribué aux cellules ressources

# --- Type 4 : Cellules immobiles ---
nombre_cellules_immobiles <- 25
mean_C_immobiles <- 100
borne_inf_im <- 85
borne_sup_im <- 120

points_initiaux_im <- rtruncnorm(nombre_cellules_immobiles, a=borne_inf_im, b=borne_sup_im, mean=mean_C_immobiles, sd=10)
nb_C_immobiles <- round(points_initiaux_im * (nombre_cellules_immobiles * mean_C_immobiles / sum(points_initiaux_im)))
ecart <- sum(nb_C_immobiles) - (nombre_cellules_immobiles * mean_C_immobiles)
indices <- sample(1:nombre_cellules_immobiles, abs(ecart))
nb_C_immobiles[indices] <- nb_C_immobiles[indices] - sign(ecart)


# --- Type 8 : Cellules immobiles type 8 ---
nombre_cellules_immobiles_type8 <- 25
mean_C_immobiles_type8 <- 100

points_initiaux_im_type8 <- rtruncnorm(nombre_cellules_immobiles_type8, a=borne_inf_im, b=borne_sup_im, mean=mean_C_immobiles_type8, sd=10)
nb_C_immobiles_type8 <- round(points_initiaux_im_type8 * (nombre_cellules_immobiles_type8 * mean_C_immobiles_type8 / sum(points_initiaux_im_type8)))
ecart <- sum(nb_C_immobiles_type8) - (nombre_cellules_immobiles_type8 * mean_C_immobiles_type8)
indices <- sample(1:nombre_cellules_immobiles_type8, abs(ecart))
nb_C_immobiles_type8[indices] <- nb_C_immobiles_type8[indices] - sign(ecart)


# --- Type 1 : Cellules mobiles ---
nombre_cellules_mobiles <- 3
mean_C_mobiles <- 100
borne_inf_mo <- 80
borne_sup_mo <- 120

points_initiaux_mo <- rtruncnorm(nombre_cellules_mobiles, a=borne_inf_mo, b=borne_sup_mo, mean=mean_C_mobiles, sd=10)
nb_C_mobiles <- round(points_initiaux_mo * (nombre_cellules_mobiles * mean_C_mobiles / sum(points_initiaux_mo)))
ecart <- sum(nb_C_mobiles) - (nombre_cellules_mobiles * mean_C_mobiles)
indices <- sample(1:nombre_cellules_mobiles, abs(ecart))
nb_C_mobiles[indices] <- nb_C_mobiles[indices] - sign(ecart)

# --- Type 7 : Cellules mobiles type 7 ---
nombre_cellules_mobiles_type7 <- 3
mean_C_mobiles_type7 <- 100

points_initiaux_mo_type7 <- rtruncnorm(nombre_cellules_mobiles_type7, a=borne_inf_mo, b=borne_sup_mo, mean=mean_C_mobiles_type7, sd=10)
nb_C_mobiles_type7 <- round(points_initiaux_mo_type7 * (nombre_cellules_mobiles_type7 * mean_C_mobiles_type7 / sum(points_initiaux_mo_type7)))
ecart <- sum(nb_C_mobiles_type7) - (nombre_cellules_mobiles_type7 * mean_C_mobiles_type7)
indices <- sample(1:nombre_cellules_mobiles_type7, abs(ecart))
nb_C_mobiles_type7[indices] <- nb_C_mobiles_type7[indices] - sign(ecart)


#paramètres de vie

cout_motilite<-1#cout supplémentaire à chaque itération pour les cellules mobiles
cout_motilite_type7<-5

cout_maintenance<-1 #cout commun au cellules mobiles et immobiles
cout_maintenance_type7<-3
cout_maintenance_type8<-3

cout_dormance<-1 #cout de la dormance (possible de mettre un cout de 0)
cout_dormance_type7<-1 
cout_dormance_type8<-1 

seuil_detection<-5 #seuil à partir duquel les cellules sont attirées par les ressources (attention avec un seuil trop bas les cellules ont du mal à s'orienter vers des ressources suffisement concentrées pour leur survie et avec un seuil trop haut les cellules mobiles ne sortent pas de leur dormance)

cout_reproduction<-55 #doit être supérieur au nombre de point cellule naissante, cout pour la reproduction d'une cellule
cout_reproduction_type7<-70
cout_reproduction_type8<-70

pts_c_naissante<-50 #nombre de points qui est attribué à la cellule issue de la reproduction de cellules vivantes ou immobiles

seuil_mort<-30 #seuil à partir duquel les points de la ressource sont trop faibles, elle meurt et devient une ressource
seuil_mort_type7<-30
seuil_mort_type8<-30

recyclage_seuil<-100/100 #doit être donné en %, il s'agit de la quantité de ressources effectivement utilisables par les autres cellules vivantes après la mort de la cellule
recyclage_seuil_type7<-100/100
recyclage_seuil_type8<-100/100

seuil_dormance<-60 #seuil à partir duquel les cellules entrent en dormance (les cellules mobiles ne peuvent plus bouger dans cette condition)
seuil_dormance_type7<-60 
seuil_dormance_type8<-60 

seuil_reproduction<-100 #seuil à partir duquel les cellules ont enmagasiné suffisement de points pour se reproduire
seuil_reproduction_type7<-250
seuil_reproduction_type8<-250

#caractéristiques de la diffusion et des ressources
vitesse_diffusion<-1 #nombre de points envoyés par la cellule ressource à ses voisins adjacents
seuil_diffusion<-10 #seuil à partir duquel une cellule ressource commence à diffuser, comme une cellule a 8 voisins vivants et qu'elle diffuse à tous ses voisins éviter de decendre en-dessous de 9 (à adapter avec la vitesse de diffusion)

# Paramètres d'ajout de ressources (attention à mettre des conditions logiques, l'aggrégation est le rayon dans lequel les cellules se trouvent on ne peut pas ajouter 10 cellules avec une agregation de 1)
ressource_value <- 250 # Points par ressource apporté
num_cells <- 8  # Nombre cellules sur lesquelles les points ressources sont apportés
aggregation <- 9  # Facteur d'agrégation des points de ressource (à 1 les cellules seront voisines)
frequency <- 1000#fréquence d'apport de nouvelles ressources, à chaque ite, toutes les 2 ite, toutes les 3 ite, 4... ect si pas d'apport mettre une frequence supérieure au nombre d'ite de la simu

#caractéristiques simulation
nombre_ite<-200
nombre_repetition<-1

# Fixed maximum values for color ranges
MAX_RESOURCE_POINTS <- 250
MAX_LIVING_CELL_POINTS <- 400
MAX_IMMOBILE_CELL_POINTS <- 400



# Create color palettes with fixed ranges
purple_palette <- colorRampPalette(c("lavender", "purple4"))(MAX_RESOURCE_POINTS + 1)
living_palette <- colorRampPalette(c("darkolivegreen1", "darkolivegreen4"))(MAX_LIVING_CELL_POINTS + 1)
living_palette_type7 <- colorRampPalette(c("#99FFCC", "#009966"))(MAX_LIVING_CELL_POINTS + 1)
immobile_palette <- colorRampPalette(c("#FF6666", "#CC0000"))(MAX_IMMOBILE_CELL_POINTS + 1)
immobile_palette_type8 <- colorRampPalette(c("#FFCC00", "#CC9933"))(MAX_IMMOBILE_CELL_POINTS + 1)




#Execution générale pour les G100 et P100
chemin_base <- "matrices"
fichiers_csv <- list.files(path = chemin_base, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
fichiers_cibles <- fichiers_csv[grepl("(G100|P100|G000|P000)", basename(fichiers_csv))]

cat("Fichiers trouvés :", fichiers_cibles, sep = "\n")


# Groupe les fichiers par dossier slice
fichiers_par_slice <- split(fichiers_cibles, dirname(fichiers_cibles))

# Simulation params
motility_costs <- c(1, 2, 3)
n_iter <- nombre_ite
n_repeats <- nombre_repetition
n_costs <- length(motility_costs)




# Résultats globaux sur toutes les slices
resultats_total <- foreach(dossier_slice = names(fichiers_par_slice), 
                           .combine = rbind,  # pour concaténer les data.frames
                           .packages = c("tools", "dplyr", "data.table", "animation")) %dopar% {
                             
                             fichiers <- fichiers_par_slice[[dossier_slice]]
                             nom_slice <- basename(dossier_slice)
                             
                             # Initialisation locale pour concaténer les résultats G100 + P100 dans une slice
                             resultats_slice <- data.frame()
                             
                             for (fichier in fichiers) {
                               type_matrice <- substr(fichier, nchar(fichier) - 7, nchar(fichier) - 4)
                               cat("\nTraitement du fichier :", fichier, "\n")  
                               
    # Lecture intelligente
    matrix_G <- read.csv2(fichier, sep = ",", header = TRUE)
    X_g <- as.matrix(matrix_G)
    if (nrow(X_g) != 100 || ncol(X_g) != 100) {
      matrix_G <- read.csv2(fichier, sep = ",", header = FALSE)
      X_g <- as.matrix(matrix_G)
      if (nrow(X_g) != 100 || ncol(X_g) != 100) {
        warning(paste("⚠️ Le fichier", fichier, "n'est pas une matrice 100x100. Ignoré."))
        next
      }
    }
    
    n <- 100
    # Initialiser la matrice des points de ressource avec 100 points pour chaque cellule ressource (X_g == 2)
    ressource_points <- matrix(0, nrow = n, ncol = n)  # Matrice de points pour les ressources
    ressource_points[X_g == 2] <- 250  # Donner 100 points aux cellules ressource
    plot_grid_with_shades <- function(X_g, ressource_points, living_cell_points, living_cell_points_type7, immobile_cell_points, immobile_cell_points_type8) {
      n <- nrow(X_g)  # Grid size
      
      # Choix de la palette selon type de matrice
      if (grepl("000", type_matrice)) {
        color_palette <- c("lightcyan1", purple_palette, living_palette, immobile_palette,
                           living_palette_type7, immobile_palette_type8, "grey", "white")
        has_white <- TRUE
      } else {
        color_palette <- c("lightcyan1", purple_palette, living_palette, immobile_palette,
                           living_palette_type7, immobile_palette_type8, "grey")  # SANS white
        has_white <- FALSE
      }
      # Longueurs des palettes
      len_purple <- length(purple_palette)
      len_living <- length(living_palette)
      len_living7 <- length(living_palette_type7)
      len_immobile <- length(immobile_palette)
      len_immobile8 <- length(immobile_palette_type8)
      
      # Offsets corrigés : commencer à 2 (car index 1 = lightcyan1)
      offset_purple <- 2
      offset_living <- offset_purple + len_purple
      offset_immobile <- offset_living + len_living
      offset_living7 <- offset_immobile + len_immobile
      offset_immobile8 <- offset_living7 + len_living7
      offset_grey <- offset_immobile8 + len_immobile8
      
      # Initialiser la matrice de couleurs
      color_indices <- matrix(1, nrow = n, ncol = n)  # valeur par défaut : "lightcyan1"
      
      # Type 2 : Ressources
      for (i in 1:n) {
        for (j in 1:n) {
          if (X_g[i, j] == 2 && ressource_points[i, j] > 0) {
            capped <- min(ressource_points[i, j], MAX_RESOURCE_POINTS)
            index <- ceiling(capped / MAX_RESOURCE_POINTS * (len_purple - 1)) + 1
            color_indices[i, j] <- offset_purple + index
          }
        }
      }
      
      # Type 1 : Cellules vivantes normales
      for (i in 1:n) {
        for (j in 1:n) {
          if (X_g[i, j] == 1) {
            capped <- min(living_cell_points[i, j], MAX_LIVING_CELL_POINTS)
            index <- ceiling(capped / MAX_LIVING_CELL_POINTS * (len_living - 1)) + 1
            color_indices[i, j] <- offset_living + index
          }
        }
      }
      
      # Type 4 : Cellules immobiles
      for (i in 1:n) {
        for (j in 1:n) {
          if (X_g[i, j] == 4) {
            capped <- min(immobile_cell_points[i, j], MAX_IMMOBILE_CELL_POINTS)
            index <- ceiling(capped / MAX_IMMOBILE_CELL_POINTS * (len_immobile - 1)) + 1
            color_indices[i, j] <- offset_immobile + index
          }
        }
      }
      
      # Type 7 : Cellules vivantes type 7
      for (i in 1:n) {
        for (j in 1:n) {
          if (X_g[i, j] == 7) {
            capped <- min(living_cell_points_type7[i, j], MAX_LIVING_CELL_POINTS)
            index <- ceiling(capped / MAX_LIVING_CELL_POINTS * (len_living7 - 1)) + 1
            color_indices[i, j] <- offset_living7 + index
          }
        }
      }
      
      # Type 8 : Cellules immobiles type 8
      for (i in 1:n) {
        for (j in 1:n) {
          if (X_g[i, j] == 8) {
            capped <- min(immobile_cell_points_type8[i, j], MAX_IMMOBILE_CELL_POINTS)
            index <- ceiling(capped / MAX_IMMOBILE_CELL_POINTS * (len_immobile8 - 1)) + 1
            color_indices[i, j] <- offset_immobile8 + index
          }
        }
      }
      
      # Apply grey color for immutable cells (X_g == 3)
      color_indices[X_g == 3] <- length(color_palette) - grepl("000", type_matrice)*1
      
      if (grepl("000", type_matrice)) {
        # Apply white color for cells with value 30 (preserved from second script)
        color_indices[X_g == 30] <- length(color_palette)
      }
 
      
      # Affichage
      par(mar = c(5, 5, 2, 2))
      image(1:n, 1:n, t(color_indices)[, n:1], col = color_palette, axes = FALSE, xlab = "", ylab = "")
      grid(nx = n, ny = n, col = "grey", lty = "solid")
      box()
    }
    
    
    # Fonction pour compter les cellules vivantes
    count_living_cells <- function(X_g) {
      return(sum(X_g == 1))
    }
    
    # Fonction pour compter les cellules immobiles
    count_immobile_cells <- function(X_g) {
      return(sum(X_g == 4))
    }
    
    
    # Fonction pour compter les cellules vivantes
    count_living_cells_type7 <- function(X_g) {
      return(sum(X_g == 7))
    }
    
    # Fonction pour compter les cellules immobiles
    count_immobile_cells_type8 <- function(X_g) {
      return(sum(X_g == 8))
    }
    
    # Fonction pour calculer le total des points de ressource
    count_total_ressource_points <- function(ressource_points) {
      return(sum(ressource_points))
    }
    
    # Fonction pour calculer le total des points des cellules vivantes
    count_total_living_cell_points <- function(living_cell_points) {
      return(sum(living_cell_points))
    }
    
    # Fonction pour calculer le total des points des cellules immobiles
    count_total_immobile_cell_points <- function(immobile_cell_points) {
      return(sum(immobile_cell_points))
    }
    
    # Fonction pour calculer le total des points des cellules vivantes
    count_total_living_cell_points_type7 <- function(living_cell_points_type7) {
      return(sum(living_cell_points_type7))
    }
    
    # Fonction pour calculer le total des points des cellules immobiles
    count_total_immobile_cell_points_type8 <- function(immobile_cell_points_type8) {
      return(sum(immobile_cell_points_type8))
    }
    # Add new parameter at the top with the motility costs you want to test
    motility_costs <- c(1)  # Three different motility costs to test
    
    # Modify the create_results_df function to include motility_cost column
    create_results_df <- function(n_iter, n_repeats, n_costs) {
      # Calculate total number of rows needed (iterations+1) * repeats * different costs
      total_rows <- (n_iter + 1) * n_repeats * n_costs
      
      results <- data.frame(
        motility_cost = numeric(total_rows),  # Add this new column
        repetition = integer(total_rows),
        iteration = integer(total_rows),
        
        # Rest of your columns remain the same
        living_cells = integer(total_rows),
        immobile_cells = integer(total_rows),
        living_cells_type7 = integer(total_rows),
        immobile_cells_type8 = integer(total_rows),
        
        total_ressource_points = numeric(total_rows),
        total_living_cell_points = numeric(total_rows),
        total_immobile_cell_points = numeric(total_rows),
        total_living_cell_points_type7 = numeric(total_rows),
        total_immobile_cell_points_type8 = numeric(total_rows),
        
        # Living Cells Resource Consumption
        living_cells_motility_consumption = numeric(total_rows),
        living_cells_maintenance_consumption = numeric(total_rows),
        living_cells_dormancy_consumption = numeric(total_rows),
        living_cells_reproduction_consumption = numeric(total_rows),
        living_cells_total_consumption = numeric(total_rows),
        cumulated_living_cells_motility_consumption = numeric(total_rows),
        cumulated_living_cells_maintenance_consumption = numeric(total_rows),
        cumulated_living_cells_dormancy_consumption = numeric(total_rows),
        cumulated_living_cells_reproduction_consumption = numeric(total_rows),
        cumulated_living_cells_total_consumption = numeric(total_rows),
        
        # Living Cells Resource Consumption type7
        living_cells_motility_consumption_type7 = numeric(total_rows),
        living_cells_maintenance_consumption_type7 = numeric(total_rows),
        living_cells_dormancy_consumption_type7 = numeric(total_rows),
        living_cells_reproduction_consumption_type7 = numeric(total_rows),
        living_cells_total_consumption_type7 = numeric(total_rows),
        cumulated_living_cells_motility_consumption_type7 = numeric(total_rows),
        cumulated_living_cells_maintenance_consumption_type7 = numeric(total_rows),
        cumulated_living_cells_dormancy_consumption_type7 = numeric(total_rows),
        cumulated_living_cells_reproduction_consumption_type7 = numeric(total_rows),
        cumulated_living_cells_total_consumption_type7 = numeric(total_rows),
        
        # immobile
        immobile_cells_maintenance_consumption = numeric(total_rows),
        immobile_cells_dormancy_consumption = numeric(total_rows),
        immobile_cells_reproduction_consumption = numeric(total_rows),
        immobile_cells_total_consumption = numeric(total_rows),
        cumulated_immobile_cells_maintenance_consumption = numeric(total_rows),
        cumulated_immobile_cells_dormancy_consumption = numeric(total_rows),
        cumulated_immobile_cells_reproduction_consumption = numeric(total_rows),
        cumulated_immobile_cells_total_consumption = numeric(total_rows),
        
        # immobile type8
        immobile_cells_maintenance_consumption_type8 = numeric(total_rows),
        immobile_cells_dormancy_consumption_type8 = numeric(total_rows),
        immobile_cells_reproduction_consumption_type8 = numeric(total_rows),
        immobile_cells_total_consumption_type8 = numeric(total_rows),
        cumulated_immobile_cells_maintenance_consumption_type8 = numeric(total_rows),
        cumulated_immobile_cells_dormancy_consumption_type8 = numeric(total_rows),
        cumulated_immobile_cells_reproduction_consumption_type8 = numeric(total_rows),
        cumulated_immobile_cells_total_consumption_type8 = numeric(total_rows),
        
        resources_uptake_living_cells = numeric(total_rows),
        resources_uptake_immobile_cells = numeric(total_rows),
        resources_uptake_living_cells_type7 = numeric(total_rows),
        resources_uptake_immobile_cells_type8 = numeric(total_rows),
        resources_diffused = numeric(total_rows),
        resources_added = numeric(total_rows),
        resources_from_dead_cells = numeric(total_rows),
        cumulated_resources_added = numeric(total_rows),
        cumulated_resources_from_dead_cells = numeric(total_rows),
        
        living_cells_reproduced = integer(total_rows),
        immobile_cells_reproduced = integer(total_rows),
        living_cells_reproduced_type7 = integer(total_rows),
        immobile_cells_reproduced_type8 = integer(total_rows),
        
        cells_died = integer(total_rows),
        
        warnings = character(total_rows),
        
        resource_cell_transfers = integer(total_rows),
        immobile_cell_transfers = integer(total_rows),
        living_cell_transfers = integer(total_rows),
        immobile_cell_transfers_type8 = integer(total_rows),
        living_cell_transfers_type7 = integer(total_rows),
        new_resource_cells = integer(total_rows)
      )
      return(results)
    }
    
    # Initialize the results with space for all motility costs
    n_iter <- nombre_ite
    n_repeats <- nombre_repetition
    n_costs <- length(motility_costs)
    results <- create_results_df(n_iter, n_repeats, n_costs)
    
    # Now we need to modify the main simulation loop to iterate through each motility cost
    overall_index <- 1  # Keep track of position in results dataframe
    
    
    for (cost_idx in 1:n_costs) {
      # Set the current motility cost
      cout_motilite <- motility_costs[cost_idx]
      
      # Repeat the simulation x10 for each motility cost
      for (repeat_idx in 1:n_repeats) {
        
        # Réinitialiser les matrices des points et des cellules vivantes pour chaque répétition
        X_g <- as.matrix(matrix_G)
        ressource_points <- matrix(0, nrow = n, ncol = n)
        ressource_points[X_g == 2] <- 100
        living_cell_points <- matrix(0, nrow = n, ncol = n)
        immobile_cell_points <- matrix(0, nrow = n, ncol = n)
        living_cell_points_type7 <- matrix(0, nrow = n, ncol = n)
        immobile_cell_points_type8 <- matrix(0, nrow = n, ncol = n)
        
        
        # Initialisation des cellules vivantes et des ressources
        initial_living_cells <- sample(which(X_g == 0), nombre_cellules_mobiles)
        X_g[initial_living_cells] <- 1
        living_cell_points[initial_living_cells] <- nb_C_mobiles
        initial_immobile_cells <- sample(which(X_g == 0), nombre_cellules_immobiles)
        X_g[initial_immobile_cells] <- 4
        immobile_cell_points[initial_immobile_cells] <- nb_C_immobiles
        
        initial_living_cells_type7 <- sample(which(X_g == 0), nombre_cellules_mobiles_type7)
        X_g[initial_living_cells_type7] <- 7
        living_cell_points_type7[initial_living_cells_type7] <- nb_C_mobiles_type7
        initial_immobile_cells_type8 <- sample(which(X_g == 0), nombre_cellules_immobiles_type8)
        X_g[initial_immobile_cells_type8] <- 8
        immobile_cell_points_type8[initial_immobile_cells_type8] <- nb_C_immobiles_type8
        
        # Initialize cumulated values for the current repetition
        cumulated_living_cells_motility_consumption <- 0
        cumulated_living_cells_maintenance_consumption <- 0
        cumulated_living_cells_dormancy_consumption <- 0
        cumulated_living_cells_reproduction_consumption <- 0
        cumulated_living_cells_total_consumption <- 0
        
        cumulated_immobile_cells_maintenance_consumption <- 0
        cumulated_immobile_cells_dormancy_consumption <- 0
        cumulated_immobile_cells_reproduction_consumption <- 0
        cumulated_immobile_cells_total_consumption <- 0
        
        # Initialize cumulated values for the current repetition
        cumulated_living_cells_motility_consumption_type7 <- 0
        cumulated_living_cells_maintenance_consumption_type7 <- 0
        cumulated_living_cells_dormancy_consumption_type7 <- 0
        cumulated_living_cells_reproduction_consumption_type7 <- 0
        cumulated_living_cells_total_consumption_type7 <- 0
        
        cumulated_immobile_cells_maintenance_consumption_type8 <- 0
        cumulated_immobile_cells_dormancy_consumption_type8 <- 0
        cumulated_immobile_cells_reproduction_consumption_type8 <- 0
        cumulated_immobile_cells_total_consumption_type8 <- 0
        
        cumulated_resources_added <- 0
        cumulated_resources_from_dead_cells <- 0
        
        # Record initial state (t0) in results
        results$motility_cost[overall_index] <- cout_motilite  # Store the motility cost
        results$motility_cost_type7[overall_index] <- cout_motilite_type7  # Store the motility cost
        results$repetition[overall_index] <- repeat_idx
        results$iteration[overall_index] <- 0  # t0 is iteration 0
        results$living_cells[overall_index] <- count_living_cells(X_g)
        results$immobile_cells[overall_index] <- count_immobile_cells(X_g)
        results$living_cells_type7[overall_index] <- count_living_cells_type7(X_g)
        results$immobile_cells_type8[overall_index] <- count_immobile_cells_type8(X_g)
        results$total_ressource_points[overall_index] <- count_total_ressource_points(ressource_points)
        results$total_living_cell_points[overall_index] <- count_total_living_cell_points(living_cell_points)
        results$total_living_cell_points_type7[overall_index] <- count_total_living_cell_points_type7(living_cell_points_type7)
        results$total_immobile_cell_points[overall_index] <- count_total_immobile_cell_points(immobile_cell_points)
        results$total_immobile_cell_points_type8[overall_index] <- count_total_immobile_cell_points_type8(immobile_cell_points_type8)
        
        overall_index <- overall_index + 1  # Increment the index
        
   
        # Sauvegarder l'animation pour chaque répétition et chaque coût de motilité
        gif_file <- paste0(nom_slice, type_matrice, "_rep", repeat_idx, ".gif")
        saveGIF({
          # Plot initial state (t0)
          plot_grid_with_shades(X_g, ressource_points, living_cell_points, immobile_cell_points,  living_cell_points_type7, immobile_cell_points_type8)
          title(main = paste("Coût Motilité:", cout_motilite, "Répétition", repeat_idx, "Itération 0 (État initial)"))
          
          for (ite in 1:n_iter) {
            # Exécuter la simulation pour une itération
            result <- iterate_diffusion_and_movement(X_g, ressource_points, living_cell_points, immobile_cell_points, living_cell_points_type7, immobile_cell_points_type8, frequency)
            # Mise à jour de la perte de ressources par cellule
            X_g <- result$X_g
            ressource_points <- result$ressource_points
            living_cell_points <- result$living_cell_points
            immobile_cell_points <- result$immobile_cell_points
            living_cell_points_type7 <- result$living_cell_points_type7
            immobile_cell_points_type8 <- result$immobile_cell_points_type8
            
            
            # Cumulate consumption and resource dynamics
            cumulated_living_cells_motility_consumption <- 
              cumulated_living_cells_motility_consumption + result$stats$resources_consumed_by_living$motility
            cumulated_living_cells_maintenance_consumption <- 
              cumulated_living_cells_maintenance_consumption + result$stats$resources_consumed_by_living$maintenance
            cumulated_living_cells_dormancy_consumption <- 
              cumulated_living_cells_dormancy_consumption + result$stats$resources_consumed_by_living$dormancy
            cumulated_living_cells_reproduction_consumption <- 
              cumulated_living_cells_reproduction_consumption + result$stats$resources_consumed_by_living$reproduction
            cumulated_living_cells_total_consumption <- 
              cumulated_living_cells_total_consumption + result$stats$resources_consumed_by_living$total_consumed
            
            cumulated_immobile_cells_maintenance_consumption <- 
              cumulated_immobile_cells_maintenance_consumption + result$stats$resources_consumed_by_immobile$maintenance
            cumulated_immobile_cells_dormancy_consumption <- 
              cumulated_immobile_cells_dormancy_consumption + result$stats$resources_consumed_by_immobile$dormancy
            cumulated_immobile_cells_reproduction_consumption <- 
              cumulated_immobile_cells_reproduction_consumption + result$stats$resources_consumed_by_immobile$reproduction
            cumulated_immobile_cells_total_consumption <- 
              cumulated_immobile_cells_total_consumption + result$stats$resources_consumed_by_immobile$total_consumed
            
            
            #type 7 et 8
            cumulated_living_cells_motility_consumption_type7 <- 
              cumulated_living_cells_motility_consumption_type7 + result$stats$resources_consumed_by_living_type7$motility
            cumulated_living_cells_maintenance_consumption_type7 <- 
              cumulated_living_cells_maintenance_consumption_type7 + result$stats$resources_consumed_by_living_type7$maintenance
            cumulated_living_cells_dormancy_consumption_type7 <- 
              cumulated_living_cells_dormancy_consumption_type7 + result$stats$resources_consumed_by_living_type7$dormancy
            cumulated_living_cells_reproduction_consumption_type7 <- 
              cumulated_living_cells_reproduction_consumption_type7 + result$stats$resources_consumed_by_living_type7$reproduction
            cumulated_living_cells_total_consumption_type7 <- 
              cumulated_living_cells_total_consumption_type7 + result$stats$resources_consumed_by_living_type7$total_consumed
            
            cumulated_immobile_cells_maintenance_consumption_type8 <- 
              cumulated_immobile_cells_maintenance_consumption_type8 + result$stats$resources_consumed_by_immobile_type8$maintenance
            cumulated_immobile_cells_dormancy_consumption_type8 <- 
              cumulated_immobile_cells_dormancy_consumption_type8 + result$stats$resources_consumed_by_immobile_type8$dormancy
            cumulated_immobile_cells_reproduction_consumption_type8 <- 
              cumulated_immobile_cells_reproduction_consumption_type8 + result$stats$resources_consumed_by_immobile_type8$reproduction
            cumulated_immobile_cells_total_consumption_type8 <- 
              cumulated_immobile_cells_total_consumption_type8 + result$stats$resources_consumed_by_immobile_type8$total_consumed
            
            
            
            cumulated_resources_added <- 
              cumulated_resources_added + result$stats$resources_added
            cumulated_resources_from_dead_cells <- 
              cumulated_resources_from_dead_cells + result$stats$resources_from_dead_cells
            
            # Basic Cell and Resource Metrics
            results$motility_cost[overall_index] <- cout_motilite  # Store the motility cost
            results$motility_cost_type7[overall_index] <- cout_motilite_type7
            
            results$repetition[overall_index] <- repeat_idx
            results$iteration[overall_index] <- ite
            results$living_cells[overall_index] <- count_living_cells(X_g)
            results$immobile_cells[overall_index] <- count_immobile_cells(X_g)
            results$living_cells_type7[overall_index] <- count_living_cells_type7(X_g)
            results$immobile_cells_type8[overall_index] <- count_immobile_cells_type8(X_g)
            results$total_ressource_points[overall_index] <- count_total_ressource_points(ressource_points)
            results$total_living_cell_points[overall_index] <- count_total_living_cell_points(living_cell_points)
            results$total_immobile_cell_points[overall_index] <- count_total_immobile_cell_points(immobile_cell_points)
            results$total_living_cell_points_type7[overall_index] <- count_total_living_cell_points_type7(living_cell_points_type7)
            results$total_immobile_cell_points_type8[overall_index] <- count_total_immobile_cell_points_type8(immobile_cell_points_type8)
            
            # Living Cells Resource Consumption
            results$living_cells_motility_consumption[overall_index] <- result$stats$resources_consumed_by_living$motility
            results$living_cells_maintenance_consumption[overall_index] <- result$stats$resources_consumed_by_living$maintenance
            results$living_cells_dormancy_consumption[overall_index] <- result$stats$resources_consumed_by_living$dormancy
            results$living_cells_reproduction_consumption[overall_index] <- result$stats$resources_consumed_by_living$reproduction
            results$living_cells_total_consumption[overall_index] <- result$stats$resources_consumed_by_living$total_consumed
            
            # Cumulated Living Cells Resource Consumption
            results$cumulated_living_cells_motility_consumption[overall_index] <- cumulated_living_cells_motility_consumption
            results$cumulated_living_cells_maintenance_consumption[overall_index] <- cumulated_living_cells_maintenance_consumption
            results$cumulated_living_cells_dormancy_consumption[overall_index] <- cumulated_living_cells_dormancy_consumption
            results$cumulated_living_cells_reproduction_consumption[overall_index] <- cumulated_living_cells_reproduction_consumption
            results$cumulated_living_cells_total_consumption[overall_index] <- cumulated_living_cells_total_consumption
            
            
            # Living Cells Resource Consumption_type7
            results$living_cells_motility_consumption_type7[overall_index] <- result$stats$resources_consumed_by_living_type7$motility
            results$living_cells_maintenance_consumption_type7[overall_index] <- result$stats$resources_consumed_by_living_type7$maintenance
            results$living_cells_dormancy_consumption_type7[overall_index] <- result$stats$resources_consumed_by_living_type7$dormancy
            results$living_cells_reproduction_consumption_type7[overall_index] <- result$stats$resources_consumed_by_living_type7$reproduction
            results$living_cells_total_consumption_type7[overall_index] <- result$stats$resources_consumed_by_living_type7$total_consumed
            
            # Cumulated Living Cells Resource Consumption_type7
            results$cumulated_living_cells_motility_consumption_type7[overall_index] <- cumulated_living_cells_motility_consumption_type7
            results$cumulated_living_cells_maintenance_consumption_type7[overall_index] <- cumulated_living_cells_maintenance_consumption_type7
            results$cumulated_living_cells_dormancy_consumption_type7[overall_index] <- cumulated_living_cells_dormancy_consumption_type7
            results$cumulated_living_cells_reproduction_consumption_type7[overall_index] <- cumulated_living_cells_reproduction_consumption_type7
            results$cumulated_living_cells_total_consumption_type7[overall_index] <- cumulated_living_cells_total_consumption_type7
            
            # Immobile Cells Resource Consumption
            results$immobile_cells_maintenance_consumption[overall_index] <- result$stats$resources_consumed_by_immobile$maintenance
            results$immobile_cells_dormancy_consumption[overall_index] <- result$stats$resources_consumed_by_immobile$dormancy
            results$immobile_cells_reproduction_consumption[overall_index] <- result$stats$resources_consumed_by_immobile$reproduction
            results$immobile_cells_total_consumption[overall_index] <- result$stats$resources_consumed_by_immobile$total_consumed
            
            # Cumulated Immobile Cells Resource Consumption
            results$cumulated_immobile_cells_maintenance_consumption[overall_index] <- cumulated_immobile_cells_maintenance_consumption
            results$cumulated_immobile_cells_dormancy_consumption[overall_index] <- cumulated_immobile_cells_dormancy_consumption
            results$cumulated_immobile_cells_reproduction_consumption[overall_index] <- cumulated_immobile_cells_reproduction_consumption
            results$cumulated_immobile_cells_total_consumption[overall_index] <- cumulated_immobile_cells_total_consumption
            
            # Immobile Cells Resource Consumption_type8
            results$immobile_cells_maintenance_consumption_type8[overall_index] <- result$stats$resources_consumed_by_immobile_type8$maintenance
            results$immobile_cells_dormancy_consumption_type8[overall_index] <- result$stats$resources_consumed_by_immobile_type8$dormancy
            results$immobile_cells_reproduction_consumption_type8[overall_index] <- result$stats$resources_consumed_by_immobile_type8$reproduction
            results$immobile_cells_total_consumption_type8[overall_index] <- result$stats$resources_consumed_by_immobile_type8$total_consumed
            
            # Cumulated Immobile Cells Resource Consumption
            results$cumulated_immobile_cells_maintenance_consumption_type8[overall_index] <- cumulated_immobile_cells_maintenance_consumption_type8
            results$cumulated_immobile_cells_dormancy_consumption_type8[overall_index] <- cumulated_immobile_cells_dormancy_consumption_type8
            results$cumulated_immobile_cells_reproduction_consumption_type8[overall_index] <- cumulated_immobile_cells_reproduction_consumption_type8
            results$cumulated_immobile_cells_total_consumption_type8[overall_index] <- cumulated_immobile_cells_total_consumption_type8
            
            # Resource Dynamics
            results$resources_uptake_living_cells[overall_index] <- result$stats$resources_uptake$living_cells
            results$resources_uptake_immobile_cells[overall_index] <- result$stats$resources_uptake$immobile_cells
            results$resources_uptake_living_cells_type7[overall_index] <- result$stats$resources_uptake$living_cells_type7
            results$resources_uptake_immobile_cells_type8[overall_index] <- result$stats$resources_uptake$immobile_cells_type8
            results$resources_diffused[overall_index] <- result$stats$resources_diffused
            results$resources_added[overall_index] <- result$stats$resources_added
            results$resources_from_dead_cells[overall_index] <- result$stats$resources_from_dead_cells
            
            # Cumulated Resource Dynamics
            results$cumulated_resources_added[overall_index] <- cumulated_resources_added
            results$cumulated_resources_from_dead_cells[overall_index] <- cumulated_resources_from_dead_cells
            
            # Cell Reproduction
            results$living_cells_reproduced[overall_index] <- result$stats$living_cells_reproduced
            results$immobile_cells_reproduced[overall_index] <- result$stats$immobile_cells_reproduced
            results$living_cells_reproduced_type7[overall_index] <- result$stats$living_cells_reproduced_type7
            results$immobile_cells_reproduced_type8[overall_index] <- result$stats$immobile_cells_reproduced_type8
            
            # Cells Lifecycle
            results$cells_died[overall_index] <- result$stats$cells_died
            
            # Warnings
            results$warnings[overall_index] <- paste(result$stats$warnings, collapse = "; ")
            
            overall_index <- overall_index + 1  # Increment the index
            
            # Afficher l'état actuel de la grille
            plot_grid_with_shades(X_g, ressource_points, living_cell_points, immobile_cell_points, living_cell_points_type7, immobile_cell_points_type8)
            title(main = paste("Coût Motilité:", cout_motilite, "Répétition", repeat_idx, "Itération", ite))
            
          }
          
          
        }, interval = 0.1, movie.name = gif_file)
      }
    }
    
    # Sauvegarder les résultats dans un fichier CSV
    nom_fichier <- paste0("simulation_results_", nom_slice, "_", type_matrice, "_all_costs.csv")
    write.csv(results, file = nom_fichier, row.names = FALSE)
    
    
  }
}

if (is_cluster()) {
  closeCluster(cl)
  mpi.quit()
} else {
  stopCluster(cl)
}
