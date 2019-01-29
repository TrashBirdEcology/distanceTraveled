#' @title Calculate the `Distance Travelled` by a multivariable system
#' @param dataIn A data frame containing columns `value`, `time`, and any number of state variables (the order of all columns is irrelevant).
#' @param group.vars Character vector of state variable(s) (e.g., species) to be used as grouping variables. Distances will be calculated indiviudally over time, according to these grouping variables. Default = 'variable'
#' @param value The values over which distance will be calculated (e.g., species abundances).
#' @import tidyr
#' @param derivs Logical (Default = T), calculate derivatives of distance traveled?
#' @export
#' @examples
#' time <- c(1:20)
#' value <- rnorm(n=20)
#' variable <- rep(c('a','b'), 20)
#' dist <- distance_traveled(data.frame(time, variable, value))
#' plot( dist$time, dist$s) # Cumulative distance traveled
#' plot( dist$time, dist$dsdt)  # Velocity = First derivative of s (Ds/t or dsdt or s')
#' plot( dist$time, dist$d2sdt2)  # Acceleration = Second derivative of s (Ds/t or dsdt or s')
distance_travelled <-
    function(dataIn,
             group.vars = c('variable'),
             derivs = T,
             print = T) {
        groupingSyms <- rlang::syms(group.vars)

        distances <-
            dataIn %>%
            ungroup() %>%  # precautionary
            # Group data by
            group_by(!!!groupingSyms) %>%
            # Sort data for distance calculations
            arrange(time) %>%
            # Calculate the Euclidean distance for each variable between two time points
            # Calculates first difference
            mutate(dx = value - lag(value)) %>%
            # Remove the first dx for each variable because it == NA
            na.omit(dx) %>%
            # Remove the grouping factor
            ungroup() %>%
            group_by(time) %>%
            # Sum the dx across all species at each time point
            summarize(ds = sqrt(sum(dx ^ 2))) %>%
            # Calculate cumulative sum of ds
            mutate(s = cumsum(ds))

        if (derivs == T) {
            distances <- distances %>%
                mutate(dsdt = ((s - lag(s)) / (time - lag(time))),
                       d2sdt2 = ((dsdt - lag(dsdt)) / (time - lag(time)))) %>%
                ungroup()

        } # Calculate the derivatives of the cumulative sum (s)

        if (print == T) {
            head(distances)
        } # Print df to screen

        return(distances)

    }
