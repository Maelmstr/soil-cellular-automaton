#Fonction pour calculer la distance entre deux points
calculate_distance <- function(i1, j1, i2, j2) {
  return(sqrt((i1 - i2)^2 + (j1 - j2)^2))
}

#Fonction pour obtenir les voisins d'une cellule
get_neighbors <- function(grid, i, j) {
  n <- nrow(grid)
  neighbors <- list()
  
  # Vérifier les 8 directions (haut, bas, gauche, droite, diagonales)
  directions <- list(
    c(-1, 0), c(1, 0), c(0, -1), c(0, 1),  # Cardinal directions
    c(-1, -1), c(-1, 1), c(1, -1), c(1, 1)  # Diagonal directions
  )
  
  for (dir in directions) {
    ni <- i + dir[1]
    nj <- j + dir[2]
    
    # Vérifier si le voisin est dans la grille
    if (ni >= 1 && ni <= n && nj >= 1 && nj <= n) {
      neighbors <- c(neighbors, list(c(ni, nj)))
    }
  }
  
  return(neighbors)
}

#Fonction pour ajouter des ressources aléatoires
add_random_ressources <- function(new_X_g, new_ressource_points, num_cells, ressource_value, aggregation, frequency) {
  # Récupérer les indices des cellules de type 0 ou 2
  available_cells <- which(new_X_g == 0 | new_X_g == 2, arr.ind = TRUE)
  
  # Calculer la somme totale des ressources avant ajout pour vérification
  total_resources_before <- sum(new_ressource_points)
  
  if (nrow(available_cells) >= 1) {
    # Si moins de num_cells disponibles, on en prend autant que possible
    num_cells_to_add <- min(nrow(available_cells), num_cells)
    
    # Sélectionner une cellule centrale aléatoire pour commencer l'agrégation
    central_index <- sample(1:nrow(available_cells), 1)
    central_cell <- available_cells[central_index, ]
    i_central <- as.integer(central_cell[1])
    j_central <- as.integer(central_cell[2])
    
    # Liste pour suivre les cellules déjà sélectionnées
    selected_cells <- matrix(0, nrow = nrow(new_X_g), ncol = ncol(new_X_g))
    
    # Ajouter des ressources à la cellule centrale
    new_ressource_points[i_central, j_central] <- new_ressource_points[i_central, j_central] + ressource_value
    new_X_g[i_central, j_central] <- 2
    selected_cells[i_central, j_central] <- 1
    
    # Fonction pour vérifier si une cellule est disponible
    is_available <- function(i, j) {
      return(i > 0 && i <= nrow(new_X_g) && j > 0 && j <= ncol(new_X_g) && 
               (new_X_g[i, j] == 0 || new_X_g[i, j] == 2) && selected_cells[i, j] == 0)
    }
    
    # Boucle pour ajouter des ressources aux cellules restantes
    for (k in 2:num_cells_to_add) {
      # Créer une liste de cellules dans la zone de voisinage
      neighbors <- list()
      for (i in max(1, i_central - aggregation):min(nrow(new_X_g), i_central + aggregation)) {
        for (j in max(1, j_central - aggregation):min(ncol(new_X_g), j_central + aggregation)) {
          if (is_available(i, j)) {
            neighbors <- c(neighbors, list(c(i, j)))
          }
        }
      }
      
      if (length(neighbors) > 0) {
        # Sélectionner une cellule valide au hasard parmi les voisins
        selected_index <- sample(1:length(neighbors), 1)
        selected_neighbor <- neighbors[[selected_index]]
        i <- as.integer(selected_neighbor[1]) 
        j <- as.integer(selected_neighbor[2])  
        
        # Ajouter des ressources à la cellule sélectionnée
        new_ressource_points[i, j] <- new_ressource_points[i, j] + ressource_value
        new_X_g[i, j] <- 2
        selected_cells[i, j] <- 1
      } else {
        # Si aucun voisin disponible dans le rayon, choisir une nouvelle cellule aléatoire
        remaining_cells <- which(new_X_g == 0 | new_X_g == 2, arr.ind = TRUE)
        remaining_cells <- remaining_cells[selected_cells[remaining_cells] == 0, ]
        
        if (nrow(remaining_cells) > 0) {
          random_index <- sample(1:nrow(remaining_cells), 1)
          random_cell <- remaining_cells[random_index, ]
          i <- as.integer(random_cell[1])
          j <- as.integer(random_cell[2])
          
          new_ressource_points[i, j] <- new_ressource_points[i, j] + ressource_value
          new_X_g[i, j] <- 2
          selected_cells[i, j] <- 1
        }
      }
    }
  } else {
    print("Aucune cellule de type 0 ou 2 n'est disponible pour ajouter des ressources.")
  }
  
  # Vérifier que les ressources ont été correctement ajoutées
  total_resources_after <- sum(new_ressource_points)
  expected_increase <- num_cells_to_add * ressource_value
  actual_increase <- total_resources_after - total_resources_before
  
  if (abs(actual_increase - expected_increase) > 0.001) {
    print(paste("AVERTISSEMENT: Augmentation inattendue des ressources. Attendu:", expected_increase, 
                "Réel:", actual_increase, "Différence:", actual_increase - expected_increase))
  }
  
  return(list(new_X_g = new_X_g, new_ressource_points = new_ressource_points))
}

#Fonction principale pour le mouvement, la diffusion

#Fonction pour trouver la ressource la plus proche
find_nearest_ressource <- function(X_g, ressource_points, i, j, seuil_detection) {
  n <- nrow(X_g)
  min_distance <- Inf
  nearest <- NULL
  
  for (x in 1:n) {
    for (y in 1:n) {
      if (X_g[x, y] == 2 && ressource_points[x, y] > seuil_detection) {
        distance <- sqrt((x-i)^2 + (y-j)^2)
        if (distance < min_distance) {
          min_distance <- distance
          nearest <- c(x, y)
        }
      }
    }
  }
  return(nearest)
}

#Fonction pour le suivi des stats
update_global_stats <- function(global_stats, stats, delta, status) {
  # Update cumulative consumption for living cells
  global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$motility <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$motility + 
    stats$resources_consumed_by_living$motility
  
  global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$maintenance <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$maintenance + 
    stats$resources_consumed_by_living$maintenance
  
  global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$dormancy <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$dormancy + 
    stats$resources_consumed_by_living$dormancy
  
  global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$reproduction <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$reproduction + 
    stats$resources_consumed_by_living$reproduction
  
  # Update cumulative consumption for immobile cells
  global_stats$total_resource_flow$cumulative_resources_consumed$immobile_cells$maintenance <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$immobile_cells$maintenance + 
    stats$resources_consumed_by_immobile$maintenance
  
  global_stats$total_resource_flow$cumulative_resources_consumed$immobile_cells$dormancy <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$immobile_cells$dormancy + 
    stats$resources_consumed_by_immobile$dormancy
  
  global_stats$total_resource_flow$cumulative_resources_consumed$immobile_cells$reproduction <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$immobile_cells$reproduction + 
    stats$resources_consumed_by_immobile$reproduction
  
  # Update total system consumption
  global_stats$total_resource_flow$cumulative_resources_consumed$total_system_consumption <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$living_cells$total_consumed + 
    global_stats$total_resource_flow$cumulative_resources_consumed$immobile_cells$total_consumed
  #Necromasse
  global_stats$total_resource_flow$cumulative_resources_consumed$necromass_losses <- 
    global_stats$total_resource_flow$cumulative_resources_consumed$necromass_losses + stats$necromass_losses
  
  # Update cumulative resource uptake
  global_stats$total_resource_flow$cumulative_resources_uptake$living_cells <- 
    global_stats$total_resource_flow$cumulative_resources_uptake$living_cells + 
    stats$resources_uptake$living_cells
  
  global_stats$total_resource_flow$cumulative_resources_uptake$immobile_cells <- 
    global_stats$total_resource_flow$cumulative_resources_uptake$immobile_cells + 
    stats$resources_uptake$immobile_cells
  
  global_stats$total_resource_flow$cumulative_resources_uptake$total_uptake <- 
    global_stats$total_resource_flow$cumulative_resources_uptake$total_uptake + 
    stats$resources_uptake$total_uptake
  
  # Update mass conservation tracking
  global_stats$total_resource_flow$mass_conservation_tracking$energy_deltas <- 
    c(global_stats$total_resource_flow$mass_conservation_tracking$energy_deltas, delta)
  
  global_stats$total_resource_flow$mass_conservation_tracking$conservation_status <- 
    c(global_stats$total_resource_flow$mass_conservation_tracking$conservation_status, status)
  
  if (!status) {
    global_stats$total_resource_flow$mass_conservation_tracking$conservation_violations <- 
      global_stats$total_resource_flow$mass_conservation_tracking$conservation_violations + 1
  }
  
  # Add iteration tracking
  global_stats$total_resource_flow$iterations_resource_tracking[[length(global_stats$total_resource_flow$iterations_resource_tracking) + 1]] <- stats
  
  return(global_stats)
}

#Fonction principale iterate diffusion et mouvement
iterate_diffusion_and_movement <- function(X_g, ressource_points, living_cell_points, 
                                           immobile_cell_points, living_cell_points_type7, immobile_cell_points_type8, frequency, custom_seuil_reproduction = 100, global_stats = NULL) {
  n <- nrow(X_g)
  loss_matrix <- matrix(0, n, n)
  
  # Use the custom reproduction threshold
  seuil_reproduction <- custom_seuil_reproduction
  
  # Initialize global statistics with comprehensive resource tracking
  if (is.null(global_stats)) {
    global_stats <- list(
      total_resource_flow = list(
        initial_grid_resources = sum(ressource_points),
        initial_cell_resources = sum(living_cell_points) + sum(immobile_cell_points_type8)+sum(living_cell_points_type7) + sum(immobile_cell_points),
        cumulative_resources_consumed = list(
          living_cells = list(
            motility = 0,
            maintenance = 0,
            dormancy = 0,
            reproduction = 0,
            total_consumed = 0
          ),
          living_cells_type7 = list(
            motility = 0,
            maintenance = 0,
            dormancy = 0,
            reproduction = 0,
            total_consumed = 0
          ),
          immobile_cells = list(
            maintenance = 0,
            dormancy = 0,
            reproduction = 0,
            total_consumed = 0
          ),
          immobile_cells_type8 = list(
            maintenance = 0,
            dormancy = 0,
            reproduction = 0,
            total_consumed = 0
          ),
          total_system_consumption = 0, 
          necromass_losses=0
        ),
        cumulative_resources_uptake = list(
          living_cells = 0,
          immobile_cells = 0,
          living_cells_type7=0,
          immobile_cells_type8 = 0,
          total_uptake = 0
        ),
        iterations_resource_tracking = list(),
        cumulated_resources_consumed = 0,
        mass_conservation_tracking = list(
          initial_total_resources = sum(ressource_points) + sum(living_cell_points) + sum(immobile_cell_points)+ sum(living_cell_points_type7) + sum(immobile_cell_points_type8),
          energy_deltas = numeric(0),
          conservation_status = logical(0),
          conservation_violations = 0,
          conservation_warnings = character(0)
        )
      )
    )
  }
  
  # Initialize iteration tracking with strict resource accounting
  stats <- list(
    resources_consumed_by_living = list(
      motility = 0,
      maintenance = 0,
      dormancy = 0,
      reproduction = 0,
      total_consumed = 0
    ),
    resources_consumed_by_immobile = list(
      maintenance = 0,
      dormancy = 0,
      reproduction = 0,
      total_consumed = 0
    ),
    resources_consumed_by_living_type7 = list(
      motility = 0,
      maintenance = 0,
      dormancy = 0,
      reproduction = 0,
      total_consumed = 0
    ),
    resources_consumed_by_immobile_type8 = list(
      maintenance = 0,
      dormancy = 0,
      reproduction = 0,
      total_consumed = 0
    ),
    resources_uptake = list(
      living_cells = 0,
      immobile_cells = 0,     
      living_cells_type7 = 0,
      immobile_cells_type8 = 0,
      total_uptake = 0
    ),
    resources_diffused = 0,
    cells_died = 0,
    resources_from_dead_cells = 0,
    resources_added = 0,
    living_cells_reproduced = 0,
    immobile_cells_reproduced = 0,
    living_cells_reproduced_type7 = 0,
    immobile_cells_reproduced_type8 = 0,
    total_living_points_before = sum(living_cell_points),
    total_immobile_points_before = sum(immobile_cell_points),
    total_living_points_before_type7 = sum(living_cell_points_type7),
    total_immobile_points_before_type8 = sum(immobile_cell_points_type8),
    total_resource_points_before = sum(ressource_points),
    warnings = character(0),
    new_resource_cells = 0,
    resource_cell_transfers = 0,
    immobile_cell_transfers = 0,
    living_cell_transfers = 0, 
    immobile_cell_transfers_type8 = 0,
    living_cell_transfers_type7 = 0, 
    necromass_losses=0
  )
  
  # Create deep copies of all grids to prevent reference issues
  new_X_g <- X_g
  new_ressource_points <- ressource_points
  new_living_cell_points <- living_cell_points
  new_immobile_cell_points <- immobile_cell_points
  new_living_cell_points_type7 <- living_cell_points_type7
  new_immobile_cell_points_type8 <- immobile_cell_points_type8
  
  # Calculate initial total resources for mass conservation check
  initial_total_resources <- sum(new_ressource_points) + 
    sum(new_living_cell_points) + 
    sum(new_immobile_cell_points)+ 
    sum(new_living_cell_points_type7) + 
    sum(new_immobile_cell_points_type8)
  
  # PHASE 1: Living Cells Maintenance - Strict resource tracking
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 1) {  # Living cell
        current_points <- new_living_cell_points[i, j]
        
        # Dormancy cost
        if (current_points >= seuil_mort && current_points <= seuil_dormance) {
          new_living_cell_points[i, j] <- current_points - cout_dormance
          stats$resources_consumed_by_living$dormancy <- stats$resources_consumed_by_living$dormancy + cout_dormance
        } 
        # Active maintenance and motility cost
        else if (current_points > seuil_dormance) {
          new_living_cell_points[i, j] <- current_points - (cout_maintenance + cout_motilite)
          stats$resources_consumed_by_living$maintenance <- stats$resources_consumed_by_living$maintenance + cout_maintenance
          stats$resources_consumed_by_living$motility <- stats$resources_consumed_by_living$motility + cout_motilite
        }
        
        # Verify no negative points after consumption
        if (new_living_cell_points[i, j] < 0) {
          stats$warnings <- c(stats$warnings, sprintf("Negative points after consumption at living cell [%d,%d]", i, j))
          new_living_cell_points[i, j] <- 0
        }
      }
    }
  }
  
  
  # PHASE 1 BIS: Living Cells Type 7 Maintenance - Strict resource tracking
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 7) {  # Living cell type 7
        current_points <- new_living_cell_points_type7[i, j]
        
        # Dormancy cost
        if (current_points >= seuil_mort_type7 && current_points <= seuil_dormance_type7) {
          new_living_cell_points_type7[i, j] <- current_points - cout_dormance_type7
          stats$resources_consumed_by_living_type7$dormancy <- stats$resources_consumed_by_living_type7$dormancy + cout_dormance_type7
        } 
        # Active maintenance and motility cost
        else if (current_points > seuil_dormance_type7) {
          new_living_cell_points_type7[i, j] <- current_points - (cout_maintenance_type7 + cout_motilite_type7)
          stats$resources_consumed_by_living_type7$maintenance <- stats$resources_consumed_by_living_type7$maintenance + cout_maintenance_type7
          stats$resources_consumed_by_living_type7$motility <- stats$resources_consumed_by_living_type7$motility+ cout_motilite_type7
        }
        
        # Vérification sur points négatifs
        if (new_living_cell_points_type7[i, j] < 0) {
          stats$warnings <- c(stats$warnings, sprintf("Negative points after consumption at type 7 living cell [%d,%d]", i, j))
          new_living_cell_points_type7[i, j] <- 0
        }
      }
    }
  }
  
  
  # PHASE 2: Immobile Cells Maintenance - Strict resource tracking
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 4) {  # Immobile cell
        current_points <- new_immobile_cell_points[i, j]
        
        # Maintenance cost for active cells
        if (current_points > seuil_dormance && current_points <= seuil_reproduction) {
          new_immobile_cell_points[i, j] <- current_points - cout_maintenance
          stats$resources_consumed_by_immobile$maintenance <- stats$resources_consumed_by_immobile$maintenance + cout_maintenance
        } 
        # Dormancy cost
        else if (current_points >= seuil_mort && current_points <= seuil_dormance) {
          new_immobile_cell_points[i, j] <- current_points - cout_dormance
          stats$resources_consumed_by_immobile$dormancy <- stats$resources_consumed_by_immobile$dormancy + cout_dormance
        }
        
        # Verify no negative points after consumption
        if (new_immobile_cell_points[i, j] < 0) {
          stats$warnings <- c(stats$warnings, sprintf("Negative points after consumption at immobile cell [%d,%d]", i, j))
          new_immobile_cell_points[i, j] <- 0
        }
      }
    }
  }
  
  # PHASE 2 BIS: Immobile Cells Maintenance - type 8
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 8) {  # Immobile cell type 8
        current_points <- new_immobile_cell_points_type8[i, j]
        
        # Maintenance cost for active immobile cells
        if (current_points > seuil_dormance_type8 && current_points <= seuil_reproduction_type8) {
          new_immobile_cell_points_type8[i, j] <- current_points - cout_maintenance_type8
          stats$resources_consumed_by_immobile_type8$maintenance <- stats$resources_consumed_by_immobile_type8$maintenance + cout_maintenance_type8
        } 
        # Dormancy cost
        else if (current_points >= seuil_mort_type8 && current_points <= seuil_dormance_type8) {
          new_immobile_cell_points_type8[i, j] <- current_points - cout_dormance_type8
          stats$resources_consumed_by_immobile_type8$dormancy<- stats$resources_consumed_by_immobile_type8$dormancy + cout_dormance_type8
        }
        
        # Vérification sur points négatifs
        if (new_immobile_cell_points_type8[i, j] < 0) {
          stats$warnings <- c(stats$warnings, sprintf("Negative points after consumption at immobile cell [%d,%d]", i, j))
          new_immobile_cell_points_type8[i, j] <- 0
        }
      }
    }
  }
  
  # PHASE 3: Living Cells Reproduction - Strict parent/child resource tracking
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 1) {  # Living cell
        current_points <- new_living_cell_points[i, j]
        
        if (current_points > seuil_reproduction) {
          neighbors <- get_neighbors(X_g, i, j)
          reproduced <- FALSE
          
          for (neighbor in neighbors) {
            ni <- neighbor[1]
            nj <- neighbor[2]
            
            # Only reproduce into empty or resource cells
            if (new_X_g[ni, nj] %in% c(0, 2)) {
              # Deduct reproduction cost from parent
              new_living_cell_points[i, j] <- current_points - cout_reproduction
              stats$resources_consumed_by_living$reproduction <- stats$resources_consumed_by_living$reproduction + cout_reproduction
              
              # Create new cell
              new_X_g[ni, nj] <- 1
              new_living_cell_points[ni, nj] <- pts_c_naissante
              
              stats$living_cells_reproduced <- stats$living_cells_reproduced + 1
              reproduced <- TRUE
              break
            }
          }
          
          if (!reproduced) {
            stats$warnings <- c(stats$warnings, 
                                sprintf("Living cell at [%d,%d] couldn't reproduce (points: %.1f)", i, j, current_points))
          }
        }
      }
    }
  }
  
  # PHASE 3 BIS: Living Cells Reproduction - Strict parent/child resource tracking
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 7) {  # Living cell type 7
        current_points <- new_living_cell_points_type7[i, j]
        
        if (current_points > seuil_reproduction_type7) {
          neighbors <- get_neighbors(X_g, i, j)
          reproduced <- FALSE
          
          for (neighbor in neighbors) {
            ni <- neighbor[1]
            nj <- neighbor[2]
            
            # Only reproduce into empty or resource cells
            if (new_X_g[ni, nj] %in% c(0, 2)) {
              
              # Reproduction successful
              new_living_cell_points_type7[i, j] <- current_points - cout_reproduction_type7
              stats$resources_consumed_by_living_type7$reproduction <- stats$resources_consumed_by_living_type7$reproduction + cout_reproduction_type7
              
              new_X_g[ni, nj] <- 7
              new_living_cell_points_type7[ni, nj] <- pts_c_naissante
              
              stats$living_cells_reproduced_type7 <- stats$living_cells_reproduced_type7 + 1
              reproduced <- TRUE
              break
            }
          }
          
          if (!reproduced) {
            stats$warnings <- c(stats$warnings, 
                                sprintf("Living cell at [%d,%d] couldn't reproduce (points: %.1f)", i, j, current_points))
          }
        }
      }
    }
  }
  
  
  # PHASE 4: Immobile Cells Reproduction - Strict parent/child resource tracking
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 4) {  # Immobile cell
        current_points <- new_immobile_cell_points[i, j]
        
        if (current_points > seuil_reproduction) {
          neighbors <- get_neighbors(X_g, i, j)
          reproduced <- FALSE
          
          for (neighbor in neighbors) {
            ni <- neighbor[1]
            nj <- neighbor[2]
            
            if (new_X_g[ni, nj] %in% c(0, 2)) {
              # Deduct reproduction cost from parent
              new_immobile_cell_points[i, j] <- current_points - cout_reproduction
              stats$resources_consumed_by_immobile$reproduction <- stats$resources_consumed_by_immobile$reproduction + cout_reproduction
              
              # Create new cell
              new_X_g[ni, nj] <- 4
              new_immobile_cell_points[ni, nj] <- pts_c_naissante
              
              stats$immobile_cells_reproduced <- stats$immobile_cells_reproduced + 1
              reproduced <- TRUE
              break
            }
          }
          
          if (!reproduced) {
            stats$warnings <- c(stats$warnings, 
                                sprintf("Immobile cell at [%d,%d] couldn't reproduce (points: %.1f)", i, j, current_points))
          }
        }
      }
    }
  }
  
  # PHASE 4 BIS: Immobile Cells Reproduction - type 8
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 8) {  # Immobile cell type 8
        current_points_ <- new_immobile_cell_points_type8[i, j]
        
        if (current_points_ > seuil_reproduction_type8) {
          neighbors <- get_neighbors(X_g, i, j)
          reproduced <- FALSE
          
          for (neighbor in neighbors) {
            ni <- neighbor[1]
            nj <- neighbor[2]
            
            if (new_X_g[ni, nj] %in% c(0, 2)) {
              # Deduct reproduction cost from parent
              new_immobile_cell_points_type8[i, j] <- current_points_ - cout_reproduction_type8
              stats$resources_consumed_by_immobile_type8$reproduction <- stats$resources_consumed_by_immobile_type8$reproduction + cout_reproduction_type8
              
              # Create new immobile cell (same type 8)
              new_X_g[ni, nj] <- 8
              new_immobile_cell_points_type8[ni, nj] <- pts_c_naissante
              
              stats$immobile_cells_reproduced_type8 <- stats$immobile_cells_reproduced_type8 + 1
              reproduced <- TRUE
              break
            }
          }
          
          if (!reproduced) {
            stats$warnings <- c(stats$warnings, 
                                sprintf("Immobile cell at [%d,%d] couldn't reproduce (points: %.1f)", i, j, current_points_))
          }
        }
      }
    }
  }
  
  
  
  
  # PHASE 5: Living Cell Movement and Resource Uptake - Precise point transfers
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 1) {  # Living cell
        current_points <- new_living_cell_points[i, j]
        nearest_ressource <- find_nearest_ressource(X_g, ressource_points, i, j, seuil_detection)
        
        if (!is.null(nearest_ressource)) {
          target_i <- nearest_ressource[1]
          target_j <- nearest_ressource[2]
          
          # Calculate movement direction
          delta_i <- sign(target_i - i)
          delta_j <- sign(target_j - j)
          new_i <- i + delta_i
          new_j <- j + delta_j
          
          # Validate new position
          if (new_i >= 1 && new_i <= n && new_j >= 1 && new_j <= n) {
            target_type <- new_X_g[new_i, new_j]
            
            # Handle movement to empty or resource cell
            if (target_type %in% c(0, 2)) {
              # Consume resources if moving to resource cell
              if (target_type == 2) {
                consumed <- new_ressource_points[new_i, new_j]
                current_points <- current_points + consumed
                stats$resources_uptake$living_cells <- stats$resources_uptake$living_cells + consumed
                stats$resources_uptake$total_uptake <- stats$resources_uptake$total_uptake + consumed
                new_ressource_points[new_i, new_j] <- 0
              }
              
              # Move the cell
              new_X_g[new_i, new_j] <- 1
              new_X_g[i, j] <- ifelse(target_type == 2, 2, 0)
              new_living_cell_points[new_i, new_j] <- current_points
              new_living_cell_points[i, j] <- 0
            }
          }
        }
      }
    }
  }
  
  
  
  # PHASE 5 BIS: Living Cell Movement and Resource Uptake Type7
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 7) {  # Living cell type 7
        current_points <- new_living_cell_points_type7[i, j]
        nearest_ressource <- find_nearest_ressource(X_g, ressource_points, i, j, seuil_detection)
        
        if (!is.null(nearest_ressource)) {
          target_i <- nearest_ressource[1]
          target_j <- nearest_ressource[2]
          
          # Calculate movement direction
          delta_i <- sign(target_i - i)
          delta_j <- sign(target_j - j)
          new_i <- i + delta_i
          new_j <- j + delta_j
          
          # Validate new position
          if (new_i >= 1 && new_i <= n && new_j >= 1 && new_j <= n) {
            target_type <- new_X_g[new_i, new_j]
            
            # Handle movement to empty or resource cell
            if (target_type %in% c(0, 2)) {
              # Consume resources if moving to resource cell
              if (target_type == 2) {
                consumed <- new_ressource_points[new_i, new_j]
                current_points <- current_points + consumed
                stats$resources_uptake$living_cells_type7 <- stats$resources_uptake$living_cells_type7 + consumed
                stats$resources_uptake$total_uptake <- stats$resources_uptake$total_uptake + consumed
                new_ressource_points[new_i, new_j] <- 0
              }
              
              # Move the cell
              new_X_g[new_i, new_j] <- 7  # cell type 7 moves
              new_X_g[i, j] <- 0
              new_living_cell_points_type7[new_i, new_j] <- current_points
              new_living_cell_points_type7[i, j] <- 0
            }
          }
        }
      }
    }
  }
  
  
  
  
  # PHASE 6: Cell Death and Resource Conversion - Exact point transfers
  for (i in 1:n) {
    for (j in 1:n) {
      # Living cell death
      if (X_g[i, j] == 1 && new_living_cell_points[i, j] < seuil_mort) {
        points_converted <- new_living_cell_points[i, j]
        new_X_g[i, j] <- 2
        new_ressource_points[i, j] <- points_converted*recyclage_seuil
        new_living_cell_points[i, j] <- 0
        
        stats$cells_died <- stats$cells_died + 1
        stats$resources_from_dead_cells <- stats$resources_from_dead_cells + points_converted * recyclage_seuil
        
        loss <- points_converted * (1 - recyclage_seuil)
        loss_matrix[i, j] <- loss_matrix[i, j] + loss  # ← ici on ajoute, on n'écrase plus
      }
      
      # Immobile cell death
      if (X_g[i, j] == 4 && new_immobile_cell_points[i, j] < seuil_mort) {
        points_converted <- new_immobile_cell_points[i, j]
        new_X_g[i, j] <- 2
        new_ressource_points[i, j] <- points_converted*recyclage_seuil
        new_immobile_cell_points[i, j] <- 0
        
        stats$cells_died <- stats$cells_died + 1
        stats$resources_from_dead_cells <- stats$resources_from_dead_cells + points_converted * recyclage_seuil
        
        
        loss <- points_converted * (1 - recyclage_seuil)
        loss_matrix[i, j] <- loss_matrix[i, j] + loss  # ← ici on ajoute, on n'écrase plus
      }
    }
  }
  
  
  
  # PHASE 6 BIS: Cell Death and Resource Conversion - Exact point transfers
  for (i in 1:n) {
    for (j in 1:n) {
      # Living cell death (type 7)
      if (X_g[i, j] == 7 && new_living_cell_points_type7[i, j] < seuil_mort_type7) {
        points_converted <- new_living_cell_points_type7[i, j]
        new_X_g[i, j] <- 2
        new_ressource_points[i, j] <- points_converted * recyclage_seuil_type7
        new_living_cell_points_type7[i, j] <- 0
        
        stats$cells_died <- stats$cells_died + 1
        stats$resources_from_dead_cells <- stats$resources_from_dead_cells + points_converted * recyclage_seuil_type7
        
        loss <- points_converted * (1 - recyclage_seuil_type7)
        loss_matrix[i, j] <- loss_matrix[i, j] + loss
      }
      
      # Immobile cell death (type 8)
      if (X_g[i, j] == 8 && new_immobile_cell_points_type8[i, j] < seuil_mort_type8) {
        points_converted <- new_immobile_cell_points_type8[i, j]
        new_X_g[i, j] <- 2
        new_ressource_points[i, j] <- points_converted * recyclage_seuil_type8
        new_immobile_cell_points_type8[i, j] <- 0
        
        stats$cells_died <- stats$cells_died + 1
        stats$resources_from_dead_cells <- stats$resources_from_dead_cells + points_converted * recyclage_seuil_type8
        
        loss <- points_converted * (1 - recyclage_seuil_type8)
        loss_matrix[i, j] <- loss_matrix[i, j] + loss
      }
    }
  }
  
  
  # PHASE 7: Resource Diffusion - Controlled point transfers
  for (i in 1:n) {
    for (j in 1:n) {
      if (X_g[i, j] == 2 && new_ressource_points[i, j] >= seuil_diffusion) {
        neighbors <- get_neighbors(X_g, i, j)
        
        for (neighbor in neighbors) {
          ni <- neighbor[1]
          nj <- neighbor[2]
          
          # Calculate maximum transferable points
          transfer_points <- min(vitesse_diffusion, new_ressource_points[i, j])
          if (transfer_points <= 0) next
          
          # Handle different neighbor types
          if (new_X_g[ni, nj] == 0) {
            # Diffusion to empty cell
            new_X_g[ni, nj] <- 2
            new_ressource_points[ni, nj] <- transfer_points
            new_ressource_points[i, j] <- new_ressource_points[i, j] - transfer_points
            
            stats$new_resource_cells <- stats$new_resource_cells + 1
            stats$resources_diffused <- stats$resources_diffused + transfer_points
          } 
          else if (new_X_g[ni, nj] == 2 && new_ressource_points[ni, nj] < new_ressource_points[i, j]) {
            # Transfer between resource cells
            new_ressource_points[ni, nj] <- new_ressource_points[ni, nj] + transfer_points
            new_ressource_points[i, j] <- new_ressource_points[i, j] - transfer_points
            
            stats$resource_cell_transfers <- stats$resource_cell_transfers + 1
            stats$resources_diffused <- stats$resources_diffused + transfer_points
          } 
          else if (new_X_g[ni, nj] == 4) {
            # Uptake by immobile cell
            new_immobile_cell_points[ni, nj] <- new_immobile_cell_points[ni, nj] + transfer_points
            new_ressource_points[i, j] <- new_ressource_points[i, j] - transfer_points
            
            stats$resources_uptake$immobile_cells <- stats$resources_uptake$immobile_cells + transfer_points
            stats$resources_uptake$total_uptake <- stats$resources_uptake$total_uptake + transfer_points
            stats$immobile_cell_transfers <- stats$immobile_cell_transfers + 1
            stats$resources_diffused <- stats$resources_diffused + transfer_points
          } 
          else if (new_X_g[ni, nj] == 8) {
            # Uptake by immobile cell type 8
            new_immobile_cell_points_type8[ni, nj] <- new_immobile_cell_points_type8[ni, nj] + transfer_points
            new_ressource_points[i, j] <- new_ressource_points[i, j] - transfer_points
            
            stats$resources_uptake$immobile_cells_type8 <- stats$resources_uptake$immobile_cells_type8 + transfer_points
            stats$resources_uptake$total_uptake <- stats$resources_uptake$total_uptake + transfer_points
            stats$immobile_cell_transfers_type8 <- stats$immobile_cell_transfers_type8 + 1
            stats$resources_diffused <- stats$resources_diffused + transfer_points
          }
          else if (new_X_g[ni, nj] == 1) {
            # Uptake by living cell
            new_living_cell_points[ni, nj] <- new_living_cell_points[ni, nj] + transfer_points
            new_ressource_points[i, j] <- new_ressource_points[i, j] - transfer_points
            
            stats$resources_uptake$living_cells <- stats$resources_uptake$living_cells + transfer_points
            stats$resources_uptake$total_uptake <- stats$resources_uptake$total_uptake + transfer_points
            stats$living_cell_transfers <- stats$living_cell_transfers + 1
            stats$resources_diffused <- stats$resources_diffused + transfer_points
          }
          else if (new_X_g[ni, nj] == 7) {
            # Uptake by living cell type 7
            new_living_cell_points_type7[ni, nj] <- new_living_cell_points_type7[ni, nj] + transfer_points
            new_ressource_points[i, j] <- new_ressource_points[i, j] - transfer_points
            
            stats$resources_uptake$living_cells <- stats$resources_uptake$living_cells + transfer_points
            stats$resources_uptake$total_uptake <- stats$resources_uptake$total_uptake + transfer_points
            stats$living_cell_transfers <- stats$living_cell_transfers + 1
            stats$resources_diffused <- stats$resources_diffused + transfer_points
          }
        }
      }
    }
  }
  # PHASE 8: Periodic Resource Addition - Track added points precisely
  if (ite %% frequency == 0) {
    added_points <- num_cells * ressource_value
    result <- add_random_ressources(new_X_g, new_ressource_points, num_cells, ressource_value, aggregation, frequency)
    new_X_g <- result$new_X_g
    new_ressource_points <- result$new_ressource_points
    
    stats$resources_added <- added_points
  }
  
  # Cleanup: Ensure cells match their grid positions
  new_living_cell_points[new_X_g != 1] <- 0
  new_immobile_cell_points[new_X_g != 4] <- 0
  new_living_cell_points_type7[new_X_g != 7] <- 0
  new_immobile_cell_points_type8[new_X_g != 8] <- 0
  
  # Calculate total consumption for this iteration
  stats$resources_consumed_by_living$total_consumed <- 
    stats$resources_consumed_by_living$motility +
    stats$resources_consumed_by_living$maintenance +
    stats$resources_consumed_by_living$dormancy +
    stats$resources_consumed_by_living$reproduction
  
  stats$resources_consumed_by_immobile$total_consumed <- 
    stats$resources_consumed_by_immobile$maintenance +
    stats$resources_consumed_by_immobile$dormancy +
    stats$resources_consumed_by_immobile$reproduction
  
  # Calculate total consumption for this iteration
  stats$resources_consumed_by_living_type7$total_consumed <- 
    stats$resources_consumed_by_living_type7$motility +
    stats$resources_consumed_by_living_type7$maintenance +
    stats$resources_consumed_by_living_type7$dormancy +
    stats$resources_consumed_by_living_type7$reproduction
  
  stats$resources_consumed_by_immobile_type8$total_consumed <- 
    stats$resources_consumed_by_immobile_type8$maintenance +
    stats$resources_consumed_by_immobile_type8$dormancy +
    stats$resources_consumed_by_immobile_type8$reproduction
  
  # Mass conservation verification with more detailed tracking
  # Calculate total consumption and additions
  total_living_cell_consumption <- 
    stats$resources_consumed_by_living$motility +
    stats$resources_consumed_by_living$maintenance +
    stats$resources_consumed_by_living$dormancy +
    stats$resources_consumed_by_living$reproduction
  
  total_immobile_cell_consumption <- 
    stats$resources_consumed_by_immobile$maintenance +
    stats$resources_consumed_by_immobile$dormancy +
    stats$resources_consumed_by_immobile$reproduction
  
  # Calculate total consumption and additions
  total_living_cell_consumption_type7 <- 
    stats$resources_consumed_by_living_type7$motility +
    stats$resources_consumed_by_living_type7$maintenance +
    stats$resources_consumed_by_living_type7$dormancy +
    stats$resources_consumed_by_living_type7$reproduction
  
  total_immobile_cell_consumption_type8 <- 
    stats$resources_consumed_by_immobile_type8$maintenance +
    stats$resources_consumed_by_immobile_type8$dormancy +
    stats$resources_consumed_by_immobile_type8$reproduction
  
  total_consumption <- total_living_cell_consumption + total_immobile_cell_consumption+total_immobile_cell_consumption_type8+total_living_cell_consumption_type7
  
  # Calculate current total resources
  current_total_resources <- sum(new_ressource_points) + 
    sum(new_living_cell_points) + 
    sum(new_living_cell_points_type7) +
    sum(new_immobile_cell_points) +
    sum(new_immobile_cell_points_type8)
  
  
  # Retrieve initial total resources from global stats
  initial_total_resources <- global_stats$total_resource_flow$mass_conservation_tracking$initial_total_resources
  
  # Calculate expected total resources
  expected_total_resources <- initial_total_resources - 
    total_consumption + 
    stats$resources_added
  
  # Calculate mass conservation delta with more precise comparison
  conservation_delta <- abs(current_total_resources - expected_total_resources)
  
  # Define a relative tolerance for mass conservation
  relative_tolerance <- 1e-10 * initial_total_resources
  conservation_status <- conservation_delta < relative_tolerance
  
  # Generate detailed mass conservation warning
  if (!conservation_status) {
    warning_message <- sprintf(
      "Mass Conservation Warning:\n  Initial Resources: %.4f\n  Current Resources: %.4f\n  Expected Resources: %.4f\n  Absolute Difference: %.4f\n  Relative Difference: %.6f%%",
      initial_total_resources,
      current_total_resources,
      expected_total_resources,
      conservation_delta,
      (conservation_delta / initial_total_resources) * 100
    )
    
    # Add warning to global stats
    global_stats$total_resource_flow$mass_conservation_tracking$conservation_warnings <- 
      c(global_stats$total_resource_flow$mass_conservation_tracking$conservation_warnings, 
        warning_message)
  }
  
  # Update global statistics with comprehensive tracking
  global_stats <- update_global_stats(global_stats, stats, conservation_delta, conservation_status)
  
  return(list(
    X_g = new_X_g, 
    ressource_points = new_ressource_points, 
    living_cell_points = new_living_cell_points, 
    immobile_cell_points = new_immobile_cell_points,
    living_cell_points_type7 = new_living_cell_points_type7, 
    immobile_cell_points_type8 = new_immobile_cell_points_type8,
    stats = stats,
    global_stats = global_stats,
    loss_matrix = loss_matrix
  ))
}
