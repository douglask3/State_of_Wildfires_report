find_levels <- function(z, percentiles = seq(20, 80, by = 20)) {

    if (any(z < 0)) symetry = T else symetry = F
    # Find the corresponding quantiles
    quantiles = quantile(abs(z), probs = percentiles / 100)

    scaling_factor = 10^(ceiling(log10((quantiles))))  # Adjust based on the magnitude of your data

    # Round the scaled quantiles to the nearest nice looking numbers
    levels = round(quantiles / scaling_factor, 1) * scaling_factor
    if (symetry) levels = c(-rev(levels), levels)

    return(levels)
}

