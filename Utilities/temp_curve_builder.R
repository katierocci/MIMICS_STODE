generate_temp_curve <- function(winter_low, summer_high, mean_temp) {
  # Define the number of days in a year
  days_in_year <- 365
  
  # Generate a sequence for each day of the year
  day_of_year <- 1:days_in_year
  
  # Calculate the amplitude based on the given high and low temperatures
  amplitude <- (summer_high - winter_low) / 2
  
  # Adjust the phase to start with winter (assuming day 1 is Jan 1st)
  # Phase shift to align the sine wave with the calendar year
  # Assuming day 81 (around March 22) is the first day of spring
  phase_shift <- -2 * pi * (81 / days_in_year)
  
  # Create the temperature vector
  temperature <- mean_temp + amplitude * sin(2 * pi * day_of_year / days_in_year + phase_shift)
  
  # Make a plot
  plot(temperature)
  
  return(temperature)
}

