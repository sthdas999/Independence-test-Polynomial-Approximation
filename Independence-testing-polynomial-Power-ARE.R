###############################################################
## Local asymptotic power of the proposed CVM test
###############################################################

set.seed(123)

alpha <- 0.05
B <- 100000
Delta <- 0.80
theta <- seq(0, 4, by = 0.5)

## Null distribution
T0 <- rnorm(B)^2

## Critical value
c.alpha <- quantile(T0, 1 - alpha)

## Local power
power <- numeric(length(theta))

for(i in seq_along(theta)){
  
  Z <- rnorm(B)
  
  T.theta <- (Z + theta[i] * sqrt(Delta))^2
  
  power[i] <- mean(T.theta > c.alpha)
}

results <- data.frame(theta, power)

print(results)

plot(theta,
     power,
     type = "b",
     pch = 19,
     lwd = 2,
     ylim = c(0, 1),
     xlab = expression(theta),
     ylab = "Estimated Local Power")

abline(h = 0.05, col = 2, lty = 2)


###############################################################
## Pitman Asymptotic Relative Efficiency
###############################################################

Delta.CVM <- 0.80
Delta0 <- 0.60

ARE <- Delta.CVM / Delta0

ARE

Delta.CVM <- c(0.60, 0.65, 0.70, 0.75, 0.80,
               0.84, 0.87, 0.89, 0.90)

Delta0 <- rep(0.60, length(Delta.CVM))

ARE <- Delta.CVM / Delta0

data.frame(theta, ARE)


###############################################################
## Monte Carlo Power Study
## Empirical CvM Statistic based on Polynomial Approximation
## Varying Polynomial Degree r and Difference Order d
###############################################################

set.seed(12345)


###############################################################
## Simulation Parameters
###############################################################

n.grid <- c(50, 100, 200, 500, 1000)

r.grid <- c(1, 3, 5, 8, 10)

d.grid <- c(1, 2, 5, 8, 10)

R <- 500
B <- 2000

alpha <- 0.05

theta <- 1.5


###############################################################
## Data Generation
###############################################################

generate.data <- function(n, theta) {
  
  Z <- runif(n, -1, 1)
  
  eps <- rnorm(n)
  
  ## Alternative hypothesis
  eps <- eps + theta * sin(pi * Z)
  
  ## Response
  Y <- sin(2 * pi * Z) + eps
  
  list(Y = Y, Z = Z)
}


###############################################################
## Polynomial Approximation
###############################################################

poly.map <- function(z, r) {
  
  out <- rep(0, length(z))
  
  for (k in 0:r) {
    
    out <- out +
      z^(2 * k + 1) / factorial(2 * k + 1)
  }
  
  return(out)
}


###############################################################
## Proposed Empirical CvM Statistic
###############################################################

compute.statistic <- function(Y, Z, r, d) {
  
  #############################################################
  ## Polynomial approximation
  #############################################################
  
  U <- poly.map(Z, r)
  
  
  #############################################################
  ## Finite difference
  #############################################################
  
  W <- diff(Y, differences = d)
  
  U <- U[(d + 1):length(U)]
  
  n1 <- length(U)
  
  
  #############################################################
  ## Marginal empirical CDFs
  #############################################################
  
  FnU <- rank(U) / (n1 + 1)
  
  FnW <- rank(W) / (n1 + 1)
  
  
  #############################################################
  ## Empirical Process
  #############################################################
  
  G <- numeric(n1)
  
  for (i in 1:n1) {
    
    Fn <- mean(U <= U[i] & W <= W[i])
    
    G[i] <- sqrt(n1) *
      (Fn - FnU[i] * FnW[i])
  }
  
  
  #############################################################
  ## Empirical CvM Statistic
  #############################################################
  
  T <- mean(G^2)
  
  return(T)
}


###############################################################
## Bootstrap Critical Value
###############################################################

bootstrap.cv <- function(n, r, d, B = 2000) {
  
  stat <- numeric(B)
  
  for (b in 1:B) {
    
    ## Generate data under H0
    dat <- generate.data(n, 0)
    
    stat[b] <- compute.statistic(
      dat$Y,
      dat$Z,
      r,
      d
    )
  }
  
  ## (1-alpha) bootstrap critical value
  return(as.numeric(
    quantile(stat, 1 - alpha)
  ))
}


###############################################################
## Monte Carlo Power Calculation
###############################################################

cat("Starting simulation...\n\n")

Power.Results <- data.frame()


for (n in n.grid) {
  
  cat("=======================================================\n")
  cat("Sample size n =", n, "\n")
  cat("=======================================================\n")
  
  
  #############################################################
  ## Bootstrap critical values
  #############################################################
  
  cat("Computing bootstrap critical values...\n")
  
  crit <- numeric(length(r.grid) * length(d.grid))
  
  counter <- 1
  
  for (d in d.grid) {
    
    cat("-------------------------------------------------------\n")
    cat("Difference order d =", d, "\n")
    cat("-------------------------------------------------------\n")
    
    for (r in r.grid) {
      
      cat("  r =", r, "\n")
      
      crit[counter] <- bootstrap.cv(
        n = n,
        r = r,
        d = d,
        B = B
      )
      
      counter <- counter + 1
    }
  }
  
  
  #############################################################
  ## Monte Carlo power
  #############################################################
  
  cat("Computing empirical power...\n")
  
  power <- numeric(length(r.grid) * length(d.grid))
  
  counter <- 1
  
  for (d in d.grid) {
    
    cat("-------------------------------------------------------\n")
    cat("Difference order d =", d, "\n")
    cat("-------------------------------------------------------\n")
    
    for (r in r.grid) {
      
      cat("  r =", r, "\n")
      
      reject <- numeric(R)
      
      for (rep in 1:R) {
        
        ## Generate data under alternative
        dat <- generate.data(n, theta)
        
        Tstat <- compute.statistic(
          dat$Y,
          dat$Z,
          r,
          d
        )
        
        reject[rep] <- (Tstat > crit[counter])
      }
      
      power[counter] <- mean(reject)
      
      counter <- counter + 1
    }
  }
  
  
  #############################################################
  ## Store results
  #############################################################
  
  temp <- data.frame(
    n = n,
    d = rep(d.grid, each = length(r.grid)),
    Degree = rep(r.grid, times = length(d.grid)),
    Critical.Value = round(crit, 4),
    Power = round(power, 3)
  )
  
  Power.Results <- rbind(
    Power.Results,
    temp
  )
  
  print(temp)
}


###############################################################
## Final Power Results
###############################################################

cat("\n=======================================================\n")
cat("FINAL POWER RESULTS\n")
cat("=======================================================\n")

print(Power.Results)


###############################################################
## Power Table: r as columns, separate d
###############################################################

Power.Table <- reshape(
  Power.Results[, c("n", "d", "Degree", "Power")],
  idvar = c("n", "d"),
  timevar = "Degree",
  direction = "wide"
)


###############################################################
## Rename columns
###############################################################

names(Power.Table) <- c(
  "n",
  "d",
  "r=1",
  "r=3",
  "r=5",
  "r=8",
  "r=10"
)


###############################################################
## Print final power table
###############################################################

cat("\n=======================================================\n")
cat("EMPIRICAL POWER TABLE\n")
cat("=======================================================\n")

print(Power.Table)


###############################################################
## Save Power Results
###############################################################

write.csv(
  Power.Results,
  "Power_Results_All_n_r_d.csv",
  row.names = FALSE
)

write.csv(
  Power.Table,
  "Power_Table_All_n_d.csv",
  row.names = FALSE
)


cat("\nPower results saved successfully.\n")


###############################################################
## Power Plot
## Power versus Polynomial Degree for each d
###############################################################

for (d.value in d.grid) {
  
  temp.plot <- subset(
    Power.Results,
    d == d.value
  )
  
  plot(
    r.grid,
    temp.plot$Power[temp.plot$n == max(n.grid)],
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "Polynomial Degree (r)",
    ylab = "Empirical Power",
    ylim = c(0, 1),
    main = paste(
      "Power versus Polynomial Degree, d =",
      d.value
    )
  )
  
  grid()
}


###############################################################
## Save Results
###############################################################

write.csv(
  Power.Table,
  "Power_vs_Polynomial_Degree_and_d.csv",
  row.names = FALSE
)


###############################################################
###############################################################
## Estimated Asymptotic Relative Efficiency
## Varying Polynomial Degree r and Difference Order d
###############################################################
###############################################################

R.ARE <- 1000


###############################################################
## Reference Statistic
###############################################################

reference.statistic <- function(Y, Z) {
  
  ## Classical empirical CvM statistic
  
  n1 <- length(Y)
  
  Fy <- rank(Y) / (n1 + 1)
  
  Fz <- rank(Z) / (n1 + 1)
  
  G <- numeric(n1)
  
  for (i in 1:n1) {
    
    Fn <- mean(Y <= Y[i] & Z <= Z[i])
    
    G[i] <- sqrt(n1) *
      (Fn - Fy[i] * Fz[i])
  }
  
  return(mean(G^2))
}


###############################################################
## Storage for ARE Results
###############################################################

ARE.Results <- data.frame()


###############################################################
## Estimated ARE for Different Sample Sizes,
## Polynomial Degrees and Difference Orders
###############################################################

for (n in n.grid) {
  
  cat("=======================================================\n")
  cat("Sample size n =", n, "\n")
  cat("=======================================================\n")
  
  
  for (d in d.grid) {
    
    cat("-------------------------------------------------------\n")
    cat("Difference order d =", d, "\n")
    cat("-------------------------------------------------------\n")
    
    ARE <- numeric(length(r.grid))
    
    
    for (j in seq_along(r.grid)) {
      
      r <- r.grid[j]
      
      cat("  Computing ARE for r =", r,
          "and d =", d, "\n")
      
      
      prop <- numeric(R.ARE)
      
      ref <- numeric(R.ARE)
      
      
      #########################################################
      ## Monte Carlo Replications
      #########################################################
      
      for (rep in 1:R.ARE) {
        
        dat <- generate.data(n, theta)
        
        
        ## Proposed statistic
        prop[rep] <- compute.statistic(
          dat$Y,
          dat$Z,
          r,
          d
        )
        
        
        ## Reference statistic
        ref[rep] <- reference.statistic(
          dat$Y,
          dat$Z
        )
      }
      
      
      #########################################################
      ## Variance-based ARE
      #########################################################
      
      var.prop <- var(prop)
      
      var.ref <- var(ref)
      
      ARE[j] <- var.ref / var.prop
    }
    
    
    ###########################################################
    ## Store results for current n and d
    ###########################################################
    
    temp.ARE <- data.frame(
      n = n,
      d = d,
      Degree = r.grid,
      ARE = round(ARE, 3)
    )
    
    ARE.Results <- rbind(
      ARE.Results,
      temp.ARE
    )
    
    print(temp.ARE)
  }
}


###############################################################
## Final ARE Results
###############################################################

cat("\n=======================================================\n")
cat("FINAL ESTIMATED ARE RESULTS\n")
cat("=======================================================\n")

print(ARE.Results)


###############################################################
## ARE Table: r as Columns, separate d
###############################################################

ARE.Table <- reshape(
  ARE.Results[, c("n", "d", "Degree", "ARE")],
  idvar = c("n", "d"),
  timevar = "Degree",
  direction = "wide"
)


###############################################################
## Rename columns
###############################################################

names(ARE.Table) <- c(
  "n",
  "d",
  "r=1",
  "r=3",
  "r=5",
  "r=8",
  "r=10"
)


###############################################################
## Print Final ARE Table
###############################################################

cat("\n=======================================================\n")
cat("ESTIMATED ARE TABLE\n")
cat("=======================================================\n")

print(ARE.Table)


###############################################################
## Save ARE Results
###############################################################

write.csv(
  ARE.Results,
  "ARE_Results_All_n_r_d.csv",
  row.names = FALSE
)

write.csv(
  ARE.Table,
  "ARE_Table_All_n_d.csv",
  row.names = FALSE
)


cat("\nARE results saved successfully.\n")