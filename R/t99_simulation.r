# # Load necessary package
# set.seed(123) # For reproducibility

# # Parameters
# time_steps <- 100  # Number of time steps
# initial_population <- 100  # Starting number of individuals
# initial_biomass <- 10  # Initial biomass per individual
# temperature_amplitude <- 5  # Amplitude of temperature oscillation
# mean_temperature <- 20  # Mean temperature
# biomass_decrease_rate <- 0.1  # Biomass decrease rate per unit temperature increase
# reproduction_rate <- 0.2  # Reproduction increase rate per unit temperature increase

# # Initialize vectors
# temperature <- numeric(time_steps)
# population_list <- vector("list", time_steps)  # List to store individual biomasses

# # Initial values
# temperature[1] <- mean_temperature
# population_list[[1]] <- rep(initial_biomass, initial_population)  # Initial biomasses

# # Simulation
# for (t in 2:time_steps) {
#   # Oscillating temperature
#   temperature[t] <- mean_temperature + temperature_amplitude * sin(2 * pi * t / 20)
  
#   # Calculate new biomasses for each individual
#   previous_population <- population_list[[t - 1]]
#   # new_biomasses <- previous_population * (1 - biomass_decrease_rate * (temperature[t] - mean_temperature))
#   new_biomasses <- pmax(new_biomasses, 0.1)  # Ensure positive biomass
  
#   # Calculate reproduction: Add new individuals proportional to reproduction rate
#   num_new_individuals <- round(length(previous_population) * reproduction_rate * (temperature[t] - mean_temperature))
#   if (num_new_individuals > 0) {
#     new_individuals <- rep(0.1, num_new_individuals)  # New individuals start with a small biomass
#     new_biomasses <- c(new_biomasses, new_individuals)
#   }
  
#   # Store updated population
#   population_list[[t]] <- new_biomasses
# }

# # Example output
# str(population_list)  # Check the structure of the population list

# # Plot the total population biomass over time
# total_population_biomass <- sapply(population_list, sum)
# plot(1:time_steps, total_population_biomass, type = "l", col = "red", ylab = "Total Biomass", xlab = "Time")
