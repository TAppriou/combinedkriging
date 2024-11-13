# Weights of the LOOCV method for the combination of two sub-models
weight_LOOCV_2comb <- function(e1,e2) {
  w_LOOCV <- (crossprod(e2) - crossprod(e1,e2)) / (crossprod(e2) - 2*crossprod(e1,e2) + crossprod(e1))
  if (w_LOOCV < 0) {w_LOOCV <- 0}
  if (w_LOOCV > 1) {w_LOOCV <- 1}
  return(w_LOOCV)
}


# Weights of the LOOCV method for a binary tree structure
weight_LOOCV_tree <- function(n_tree_levels, model) {

  n_nodes <- 2^n_tree_levels - 1 # Number of nodes in the tree
  eLOOv_nodes <- matrix(0,nrow=model$n_train,ncol=n_nodes)
  w_nodes <- rep(0,length=n_nodes) # Weights of each node in the combination
  w_tot <- rep(1,length=model$n_comb_max)
  ind <- 1

  # Compute the weights at each nodes
  for (level in 1:n_tree_levels) { # For each tree level
    for (p in 1:(2^(n_tree_levels-level))) { # For each node in the tree level

      # At the first level, no children
      if (level == 1) {
        eLOOv_nodes[,ind] <- model$Kinv_Y[,ind] / model$Kinv_diag[,ind]
      }

      # Parent eLOO is the weighted sum of two child eLOO
      if (level > 1) {
        child1 <- 2^n_tree_levels - (2*(2^n_tree_levels - ind) + 1)
        child2 <- 2^n_tree_levels - (2*(2^n_tree_levels - ind))
        eLOOv_nodes[,ind] <- w_nodes[child1]*eLOOv_nodes[,child1] + w_nodes[child2]*eLOOv_nodes[,child2]
      }

      # Compute the weights (by groups of 2) using the eLOO
      if (ind %% 2 == 0) {
        w <- weight_LOOCV_2comb(eLOOv_nodes[,ind-1],eLOOv_nodes[,ind])
        w_nodes[ind-1] <- w
        w_nodes[ind] <- 1-w
      }
      ind <- ind+1
    }
  }
  w_nodes[n_nodes] <- 1 #Weight of the finale combined model (fixed to 1 since no combination)

  return(w_nodes)
}
