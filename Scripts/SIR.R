# implements the basic SIR model, and plots simulation results

# FUNCTION DEFINITIONS

SIR <- function(t, x, parms){
  # Give names to the variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
 
  # Compte derivatives  
  dS <- - parms$beta*S*I
  dI <- + parms$beta*S*I - parms$r*I
  dR <- parms$r*I # Note: because S+I+R=constant, this equation could actually be omitted,
                    # and R at any time point could simply be calculated as N-S-I.
  # Prepare output
  derivs <- c(dS, dI, dR)
  list(derivs)  # the output must be returned    
}  # end of function definition

# MAIN PROGRAM

# Load R library for ordinary differential equation solvers
library(deSolve)
# If not already installed, type
# install.packages("deSolve")

# Parameters
parms <- as.list(c(beta=1e-3, r=1e-1))	# set the parameters of the model
inits <- c(S=499, I=1, R=0)	# set the initial values
dt    <- seq(0, 100, 0.1) # set the time points for evaluation

# Calculate and print R_0 on the screen
N <- sum(inits) # Total population size
R0 <- parms$beta*N/parms$r
print(R0)

# Simulate the model
sim <- as.data.frame(ode(inits, dt, SIR, parms=parms)) # this way our set 'parms' will be used as default

# Plot the output
par(las = 1) # (horizontal labels)
lwdcurves <- 3
plot(sim$time, sim$S, type="l", col="blue", ylim=c(0,sum(inits)), xlab="time", ylab="number of individuals", lwd=lwdcurves)
lines(sim$time, sim$I, col="red", lwd=lwdcurves)
lines(sim$time, sim$R, col="darkgreen", lwd=lwdcurves)

# Add a legend to the graph
legend(70, 400, legend=c("S","I","R"), col=c("blue", "red", "darkgreen"), lty=1, lwd=2)
#dev.off()
