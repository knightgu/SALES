## Asymmetric squared error loss

ercls <- function(r, tau) {
  abs(tau - (r < 0)) * r^2
} 
