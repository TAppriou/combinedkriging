compute_alpha_tree <- function(w_nodes, n_tree_levels, model) {

  n_nodes <- 2^n_tree_levels - 1 # Number of nodes in the tree
  alpha_nodes <- rep(0,length=n_nodes) # Weights of each node in the combination
  CholK_nodes <- array(0,dim=c(model$n_train,model$n_train,n_nodes))
  Kinv_diag_nodes <- matrix(0,nrow=model$n_train,ncol=n_nodes)
  ind <- 1

  # Compute the weights at each nodes
  for (level in 1:n_tree_levels) { # For each tree level
    for (p in 1:(2^(n_tree_levels-level))) { # For each node in the tree level

      # At the first level, no children
      if (level == 1) {
        CholK_nodes[,,ind] <- model$CholK[,,ind]
        Kinv_diag_nodes[,ind] <- model$Kinv_diag[,ind]
      }

      # Parent covariance is the weighted sum of the child covariances
      if (level > 1) {
        child1 <- 2^n_tree_levels - (2*(2^n_tree_levels - ind) + 1)
        child2 <- 2^n_tree_levels - (2*(2^n_tree_levels - ind))
        CholK_nodes[,,ind] <- chol.default(alpha_nodes[child1]*crossprod(CholK_nodes[,,child1]) + alpha_nodes[child2]*crossprod(CholK_nodes[,,child2]))
        Kinv_diag_nodes[,ind] <- diag(chol2inv(CholK_nodes[,,ind]))
      }

      # Compute the covariance weights (by group of 2) using the covariances and the combination weights
      if (ind %% 2 == 0) {
        a1 <- sum( 1/Kinv_diag_nodes[,ind-1] - diag(crossprod(CholK_nodes[,,ind-1] %*% chol2inv(CholK_nodes[,,ind])))/(Kinv_diag_nodes[,ind])^2 )
        a2 <- sum( 1/Kinv_diag_nodes[,ind] - diag(crossprod(CholK_nodes[,,ind] %*% chol2inv(CholK_nodes[,,ind-1])))/(Kinv_diag_nodes[,ind-1])^2 )
        gamma <- (a2*w_nodes[ind-1]^2 - sum(1/Kinv_diag_nodes[,ind]))  / (a2*w_nodes[ind-1]^2 - sum(1/Kinv_diag_nodes[,ind]) + a1*w_nodes[ind]^2 - sum(1/Kinv_diag_nodes[,ind-1]))
        alpha_nodes[ind-1] <- gamma^2
        alpha_nodes[ind] <- (1-gamma)^2
      }
      ind <- ind+1
    }
  }
  alpha_nodes[n_nodes] <- 1 #Weight of the final covariance (fixed to 1 since no combination)
  return(alpha_nodes)
}
