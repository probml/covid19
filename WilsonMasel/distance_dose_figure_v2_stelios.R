#install required packages if not already installed
if ("ggplot2" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggplot2")
  require(ggplot2)
} else{
  require(ggplot2)
}
if ("ggpubr" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggpubr")
  require(ggpubr)
} else{
  require(ggpubr)
}
if ("reshape2" %in% rownames(installed.packages()) == FALSE) {
  install.packages("reshape2")
  require(reshape2)
} else{
  require(reshape2)
}
if ("useful" %in% rownames(installed.packages()) == FALSE) {
  install.packages("useful")
  require(useful)
} else{
  require(useful)
}
if ("truncdist" %in% rownames(installed.packages()) == FALSE) {
  install.packages("truncdist")
  require(truncdist)
} else{
  require(truncdist)
}
if ("triangle" %in% rownames(installed.packages()) == FALSE) {
  install.packages("triangle")
  require(triangle)
} else{
  require(triangle)
}


library(tidyverse)

threshold.contact <- function(randomangle = "close", finer_grid = F) {
  #Gaussian plume approach
  points <- 100000
  
  #---- using spherical coordinates -----------
  
  #randomly sampling attenuation ranges
  threshold = 1
  
  distances <- c(0.1, 0.2, 0.5, 1, 1.0001, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
  
  if (finer_grid) {
    
    distances <- seq(0.1, 10, 0.1)
    
  }
  
  rho <-
    sample(distances, points, replace = TRUE) #get wide range of distances
  
  #randomly sampling phi and theta
  phi <- rep(NA, points)
  theta <- rep(NA, points)
  
  if (randomangle == "close") {
    phi[rho > threshold] <- pi / 2
    phi[rho <= threshold] <-
      rtriangle(
        n = length(phi[rho <= threshold]),
        a = pi / 4,
        b = 3 * pi / 4,
        c = pi / 2
      ) #runif(length(phi[rho>1]), 0, pi)
    
    theta[rho > threshold] <- 0
    theta[rho <= threshold] <-
      runif(length(theta[rho <= threshold]), 0, 2 * pi)
    
  } else{
    phi[rho <= threshold] <- pi / 2
    phi[rho > threshold] <-
      rtriangle(
        n = length(phi[rho > threshold]),
        a = pi / 4,
        b = 3 * pi / 4,
        c = pi / 2
      ) #runif(length(phi[rho>1]), 0, pi)
    
    theta[rho <= threshold] <- 0
    theta[rho > threshold] <-
      runif(length(theta[rho > threshold]), 0, 2 * pi)
    
  }
  
  #convert to Cartesian for plume equation
  x <- rho * sin(phi) * cos(theta)
  y <- rho * sin(phi) * sin(theta)
  z <- rep(NA, points)
  
  if (randomangle == "close") {
    z[rho <= threshold] <- rho[rho <= threshold] * cos(phi[rho <= threshold])
    z[rho > threshold] <- 0
  } else{
    z[rho > threshold] <- rho[rho > threshold] * cos(phi[rho > threshold])
    z[rho <= threshold] <- 0
  }
  
  
  #source of plume is at (0,0,0)
  
  #check on distances...
  #test<-4
  #sqrt(x[test]^2 + y[test]^2 + z[test]^2)
  #rho[test]
  
  #holding constant for weight setting
  C.emit <- 50 #in arbitrary units
  
  #exhalation rates
  X <- (rtrunc(
    points,
    "norm",
    a = 0,
    mean = 16.3,
    sd = 4.15
  ) / (24 * 60 * 60))
  
  Q <-
    C.emit * X #viral particles/m^3 x m^3/s exhalation rates (Exposure Factors Handbook)
  
  #U = (m^3 / s) / m^2 of cross sectional area --> m/s
  #A (surface area of mouth cross sectional) - large bite
  A <- runif(points, 23, 59) / (100 ^ 2) #convert from cm^2 to m^2
  U <- X / A
  
  
  #moderately stable
  I.y <- runif(points, 0.08, 0.25)
  I.z <- runif(points, 0.03, 0.07)
  
  omega.y <- rep(NA, points)
  omega.z <- rep(NA, points)
  
  omega.y[x > 0] <- I.y[x > 0] * x[x > 0]
  omega.z[x > 0] <- I.z[x > 0] * x[x > 0]
  
  #concentration at given x, y, z points
  C <- rep(NA, points)
  C[x > 0] <-
    (Q[x > 0] / U[x > 0]) * (1 / (2 * pi * omega.y[x > 0] * omega.z[x > 0] *
                                    1)) * exp(-y[x > 0] ^ 2 / (2 * omega.y[x > 0] ^ 2)) * exp(-z[x > 0] ^ 2 /
                                                                                                (2 * omega.z[x > 0] ^ 2))
  C[x <= 0] <- 0
  
  I <- (rtrunc(
    points,
    "norm",
    a = 0,
    mean = 16.3,
    sd = 4.15
  ) / (24 * 60))
  
  duration <- 30
  
  #viral particles/m^3 x m^3/min x min
  Dose <- C * I * duration
  
  k <- 3.50E-6
  
  #Infection risk
  risk <- 1 - exp(-Dose * k)
  
  frame <-
    data.frame(
      x = x,
      y = y,
      z = z,
      C = C,
      rho = rho,
      theta = theta,
      phi = phi,
      Dose = Dose,
      omega.y = omega.y,
      omega.z = omega.z,
      Q.U = Q / U,
      I = I,
      X = X,
      A = A,
      risk = risk
    )
  
  dose.mean <- rep(NA, length(distances))
  dose.sd <- rep(NA, length(distances))
  infection.mean <- rep(NA, length(distances))
  infection.sd <- rep(NA, length(distances))
  
  for (i in 1:length(distances)) {
    dose.mean[i] <- mean(frame$Dose[frame$rho == distances[i]])
    dose.sd[i] <- sd(frame$Dose[frame$rho == distances[i]])
    infection.mean[i] <- mean(frame$risk[frame$rho == distances[i]])
    infection.sd[i] <- sd(frame$risk[frame$rho == distances[i]])
  }
  
  frame.sum <-
    data.frame(
      mean = c(dose.mean, infection.mean),
      sd = c(dose.sd, infection.sd),
      distance = distances,
      param = c(rep("dose", length(distances)), rep("risk", length(distances)))
    )
  
  data.frame(
    dose = dose.mean,
    risk = infection.mean,
    distance = distances,
    threshold = threshold,
    randomangle = randomangle
  )
}



# Simulate ----------------------------------------------------------------

# Simulate
frame_coarse <- map_dfr(c("close", "far"), threshold.contact)
frame_fine <- map_dfr(c("close", "far"), threshold.contact, finer_grid = T)

distance_dose_figure_v2 <- 
  frame_coarse %>% 
  bind_rows(frame_fine, .id = "is_coarse") %>% 
  mutate(is_coarse = is_coarse == 1)

# Save
write_csv(distance_dose_figure_v2, "output/distance-dose-figure-v2.csv")


# Plot --------------------------------------------------------------------

# Original
distance_dose_figure_v2 %>%
  filter(is_coarse) %>% 
  ggplot() +
  geom_point(aes(
    x = distance,
    y = dose,
    color = "dose",
    group = randomangle,
    alpha = randomangle
  ),
  size = 5) +
  geom_line(aes(
    x = distance,
    y = dose,
    color = "dose",
    group = randomangle,
    alpha = randomangle
  )) +
  geom_point(
    aes(
      x = distance,
      y = risk * 10,
      color = "risk*5",
      group = randomangle,
      alpha = randomangle
    ),
    size = 5
  ) +
  geom_line(aes(
    x = distance,
    y = risk * 10,
    color = "risk*5",
    group = randomangle,
    alpha = randomangle
  )) +
  scale_alpha_discrete(range = c(0.3, 1), name = "") +
  scale_colour_manual(
    values = c("black", "red"),
    labels = c("Dose", "Probability of Infection"),
    name = ""
  ) + 
  scale_y_continuous(
    name = "Dose (arbitrary units)",
    sec.axis = sec_axis(trans =  ~ . / 10, name = "Probability of Infection"),
    trans = "log10"
  ) +
  scale_x_continuous(name = "Distance from Emitter (m)") +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.direction = "horizontal"
  ) + 
  guides(alpha = FALSE)

# Save
ggsave("output/distance-dose-figure-v2_original.jpg")


# Updated
distance_dose_figure_v2 %>% 
  mutate(risk = risk * 10) %>% 
  mutate(is_coarse = if_else(is_coarse, "Coarse", "Fine")) %>% 
  gather(key, value, dose, risk) %>% 
  ggplot(aes(
    x = distance, 
    y = value, 
    color = key, 
    fill = key,
    alpha = randomangle
  )) +
  geom_point() +
  geom_line() +
  facet_wrap(~ is_coarse, nrow = 2, scales = "free") +
  scale_x_continuous(name = "Distance from Emitter (m)") +
  scale_y_continuous(
    name = "Dose (arbitrary units)",
    sec.axis = sec_axis(trans =  ~ . / 10, name = "Probability of Infection"),
    trans = "log10"
  ) +
  scale_alpha_discrete(range = c(0.3, 1), name = "") +
  scale_colour_manual(
    values = c("black", "red"),
    labels = c("Dose", "Probability of Infection"),
    name = "",
    aesthetics = c("colour", "fill")
  ) + 
  theme_pubr() +
  guides(alpha = F)

# Save
ggsave("output/distance-dose-figure-v2_updated.jpg")
