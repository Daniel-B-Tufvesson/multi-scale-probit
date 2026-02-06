source("R/internal.R")
source("R/util.R")

# Create a data structure for fitting the multi-scale probit model. This data structure
# defines the predictors (X) and the responses (Y). It contains several datasets which 
# the MSP model will be jointly fitted to. Each dataset needs to contain all the
# predictors, and its corresponding ordinal response variable.
# 
# Arguments:
# predictorNames: A vector of predictor variable names.
# responseNames: A vector of response variable names corresponding to each dataset.
# datasets: A list of data frames, each containing the predictors and response variable 
#           by columns.
# seed: An optional random seed for reproducibility. Only used if data is synthetic.
# addLevelNamePredicate: If TRUE, adds a predicate to level names (default: TRUE). Useful 
#                        for distinguishing levels when levels from different scales share 
#                        the same names.
#
# Returns:
# An object of class 'mspm_data' containing the data for fitting the MSP model.
mspm_data <- function(
    predictorNames,
    responseNames,
    datasets,
    ...,
    seed = NULL,
    addLevelNamePredicate = TRUE
) {

    # Check that there is one response name for each dataset.
    if (length(responseNames) != length(datasets)) {
        stop("Number of response names must match number of datasets. Response names: ",
             length(responseNames), ", datasets: ", length(datasets))
    }

    # Check that the predictors are in all the data sets.
    for (data in datasets) {
        missing_preds <- setdiff(predictorNames, colnames(data))
        if (length(missing_preds) > 0) {
            stop(
                paste0(
                    "Missing predictor(s) in dataset: ",
                    paste(missing_preds, collapse = ", ")
                )
            )
        }
    }

    # Check that each response variable is in its corresponding dataset.
    for (i in seq_along(datasets)) {
        data <- datasets[[i]]
        resp <- responseNames[i]
        if (!(resp %in% colnames(data))) {
            stop(
                paste0(
                    "Response variable '",
                    resp,
                    "' not found in dataset ",
                    i,
                    "."
                )
            )
        }
    }

    # Extract predictor and response data.
    Xlist <- list()
    ylist <- list()
    nlevels = c()
    for (i in seq_along(datasets)) {
        data <- datasets[[i]]
        resp <- responseNames[i]
        # Ensure predictors are a matrix, not a data frame
        Xlist[[i]] <- as.matrix(data[, predictorNames])
        # Convert ordinal string responses to integers, preserving order
        ylist[[i]] <- as.integer(factor(data[[resp]], levels = unique(as.character(sort(data[[resp]])))))
        # Get levels for response variable.
        nlevels[i] <- length(levels(as.factor(data[[resp]])))
    }

    # Extract level names.
    levelNames <- list()
    for (i in seq_along(datasets)) {
        data <- datasets[[i]]
        resp <- responseNames[i]
        levels_i <- levels(as.factor(data[[resp]]))
        if (addLevelNamePredicate) {
            levels_i <- paste0("D", i, "_", levels_i)
        }
        levelNames[[i]] <- levels_i
    }

    new_mspm_data(
        predictorNames = predictorNames,
        responseNames = responseNames,
        Xlist = Xlist,
        ylist = ylist,
        levelNames = levelNames,
        nlevels = nlevels,
        ntargets = length(datasets),
        seed = seed,
        call = match.call()
    )
}


# Split the data into training sets and test sets. 
#
# Arguments:
# data    : An mspm_data object.
# prop    : Proportion of data to include in the training set (default: 0.8).
# seed    : Random seed for reproducibility (default: NULL). 
# 
# Returns:
# A list with two elements: 'train' and 'test', each an mspm_data object.
split_data.mspm_data <- function(
    data,
    ...,
    prop = 0.8,
    seed = NULL
) {
    # Set seed for reproducibility.
    if (!is.null(seed)) {
        set.seed(seed)
    }


    # Create train/test splits for each dataset.
    Xtrain = list()
    ytrain = list()
    Xtest = list()
    ytest = list()
    for (i in seq_along(data$Xlist)) {
        X <- data$Xlist[[i]]
        y <- data$ylist[[i]]
        nobs <- nrow(X)
        train_indices <- sample(1:nobs, size = floor(prop * nobs), replace = FALSE)
        test_indices <- setdiff(1:nobs, train_indices)

        Xtrain[[i]] <- X[train_indices, , drop = FALSE]
        ytrain[[i]] <- y[train_indices]
        Xtest[[i]] <- X[test_indices, , drop = FALSE]
        ytest[[i]] <- y[test_indices]
    }

    train <- new_mspm_data(
        predictorNames = predictorNames(data),
        responseNames = responseNames(data),
        Xlist = Xtrain,
        ylist = ytrain,
        levelNames = levelNames(data),
        nlevels = nlevels(data),
        ntargets = ntargets(data),
        seed = data$seed,
        call = match.call()
    )

    test <- new_mspm_data(
        predictorNames = predictorNames(data),
        responseNames = responseNames(data),
        Xlist = Xtest,
        ylist = ytest,
        levelNames = levelNames(data),
        nlevels = nlevels(data),
        ntargets = ntargets(data),
        seed = data$seed,
        call = match.call()
    )

    # Return as mspm_data objects.
    list(train = train, test = test)
}


# Generate synthetic experiment data.
#
# Arguments:
#   nobs    : Number of observations to generate for the dataset. Note, there needs to be 
#             at least one observation per level in each target.
#   ncov    : Number of covariates (features) in the dataset.
#   data_sd : Standard deviation for generated covariate values (default: 1).
#   beta    : Numeric vector of regression coefficients (optional; generated if NULL).
#   gammas  : List of threshold values for each group (optional; generated if NULL).
#   ngamma  : Vector specifying number of gamma thresholds per group (default: c(1, 3, 3)).
#   seed    : Random seed for reproducibility (optional; generated if NULL).
#   maxProbitTries : Maximum number of attempts to generate valid probit data (default: 10000).
#
# Returns:
# An mspm_data object containing the generated synthetic data.
generate_synthetic_data <- function(
    nobs,
    ncov,
    ...,
    data_sd = 1,
    beta = NULL,
    gammas = NULL,
    ngamma = c(1, 3, 3),
    seed = NULL,
    maxProbitTries = 10000
) {

    # Check that nobs is sufficient for the number of levels.
    min_required_obs <- sum(ngamma+1)
    if (nobs < min_required_obs) {
        stop(paste0(
            "Number of observations (nobs = ", nobs, 
            ") is less than the minimum required (", min_required_obs, 
            ") to accommodate all levels in each target."
        ))
    }

    # Set seed.
    if (is.null(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }
    set.seed(seed)

    # If beta is not provided, generate random beta coefficients
    if (is.null(beta)) {
        beta <- rnorm(ncov, 0, 1)
    }
    # If gammas is not provided, generate random gamma thresholds for each group
    if (is.null(gammas)) {
        gammas <- list()
        for (i in 1:length(ngamma)) {
            gammas[[i]] <- sort(rnorm(ngamma[i], 0, 5))
        }
    }
    # Create covariate names
    X.names <- paste("X", 1:ncov, sep = "")
    data <- list()
    # For each gamma group, generate data and apply roprobit
    for (i in 1:length(gammas)) {
        # Generate covariate matrix X
        # X <- as.matrix(sapply(1:ncov, function(i) {
        #     return(rnorm(nobs, sd = data_sd))
        # }))
        X <- matrix(rnorm(nobs * ncov, sd = data_sd), nrow = nobs, ncol = ncov)
        colnames(X) <- X.names
        
        # Apply roprobit to training and test sets
        tmp_data <- roprobit(X, beta, gammas[[i]])
        # If roprobit fails, regenerate gammas and try again
        tries <- 0
        while (!is.list(tmp_data)) {
            if (tries > maxProbitTries) {
                stop("Exceeded maximum number of tries to generate valid probit data.")
            }
            tries <- tries + 1

            gammas[[i]] <- sort(rnorm(ngamma[i], 0, 5))
            tmp_data <- roprobit(X, beta, gammas[[i]])
        }

        # Convert list to data frame.
        df <- as.data.frame(tmp_data$X)
        df$Y <- tmp_data$Y
        colnames(df) <- c(X.names, "Y")

        # Store results
        data[[i]] <- df
    }

    return(mspm_data(
        predictorNames = X.names,
        responseNames = rep("Y", length(gammas)),
        datasets = data,
        seed = seed
    ))
}
