# nmtd
N-mixture models with time-to-detection data

# Installation
Run the R commands below: 
pkgs <- c("devtools", "copula", "hypergeo", "tictoc")
options(warn = -1)
for (i in pkgs){
  if (!require(i, quietly = TRUE, character.only = TRUE)){
    install.packages(i)
  }
}
devtools::install_github("david-borchers/nmtd",build = TRUE)
