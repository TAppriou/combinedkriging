weight_total <- function(n_tree_levels, w) {

  w_tot <- rep(1,length=2^(n_tree_levels-1))

  # For each tree level
  for (level in 1:n_tree_levels) {
    for (p in 1:(2^(n_tree_levels-1))) {
      parent <- 2^n_tree_levels - (2^n_tree_levels - p)%/%(2^(level-1))
      w_tot[p] <- w_tot[p] * w[parent]
    }
  }

  return(w_tot)

}
