

run_all_data_tests <- function() {
    .test_generate_data()
    .test_split_data()
    .test_level_names()
    .test_minimal_data()
    .test_accessors()
    .test_split_replication()
}


# Test case 1: Generate two sets of synthetic data with the same seed and verify they 
# are identical.
.test_generate_data <- function() {
    data1 <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )
    data2 <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )

    # Check that the generated datasets are identical.
    if (!all.equal(data1$Xlist[[1]], data2$Xlist[[1]])) {
        stop("Test failed: Predictor matrices are not identical.")
    }
    if (!all.equal(data1$ylist[[1]], data2$ylist[[1]])) {
        stop("Test failed: Response vectors are not identical.")
    }

    cat("Test passed: Generated datasets are identical.\n")
}


# Test case 2: Generate synthetic data, split it test and train, and verify the split.
.test_split_data <- function() {
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )

    split <- split_data(
        data = data,
        prop = 0.7,
        seed = 42
    )

    train <- split$train
    test <- split$test

    for (i in 1:ntargets(data)) {
        # Check size of train and test sets.
        expected_train_size <- floor(0.7 * nrow(data$Xlist[[i]]))
        expected_test_size <- nrow(data$Xlist[[i]]) - expected_train_size
        if (nrow(train$Xlist[[i]]) != expected_train_size) {
            stop(paste("Test failed: Training set size for dataset", i, "is incorrect.\n"))
        }
        if (nrow(test$Xlist[[i]]) != expected_test_size) {
            stop(paste("Test failed: Test set size for dataset", i, "is incorrect.\n"))
        }

        # Check that all split data is in the original set.
        combined_X <- rbind(train$Xlist[[i]], test$Xlist[[i]])
        combined_y <- c(train$ylist[[i]], test$ylist[[i]])
        # Sort rows before comparison to allow for unordered rows
        sort_matrix <- function(m) m[do.call(order, as.data.frame(m)), , drop=FALSE]
        orig_X_sorted <- sort_matrix(data$Xlist[[i]])
        combined_X_sorted <- sort_matrix(combined_X)
        if (!isTRUE(all.equal(orig_X_sorted, combined_X_sorted))) {
            stop(paste("Test failed: Predictor data for dataset", i, "does not match original after split (ignoring row order).\n"))
        }
        # Sort y before comparison (order by value)
        if (!isTRUE(all.equal(sort(data$ylist[[i]]), sort(combined_y)))) {
            stop(paste("Test failed: Response data for dataset", i, "does not match original after split (ignoring order).\n"))
        }
    }

    cat("Test passed: Data split correctly into training and test sets.\n")
}

# Test case 3: Generate synthetic data and verify that the level names are correct.
.test_level_names <- function() {

    # Generate synthetic data.
    ngamma <- c(2, 4)
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = ngamma,
        seed = 42
    )

    # Check level names for each dataset.
    levelNames <- levelNames(data)
    for (i in 1:ntargets(data)) {
        expected_level_names <- paste0("D", i, "_", 0:ngamma[i])
        if (!isTRUE(all.equal(levelNames[[i]], expected_level_names))) {
            stop(paste("Test failed: Level names for dataset", i, "are incorrect. Expected:", 
                       paste(expected_level_names, collapse = ", "), 
                       "Got:", 
                       paste(levelNames[[i]], collapse = ", ")))
        }
    }

    cat("Test passed: Level names are correct.\n")
}

# Test case 4: Generate minimal synthetic data (2 obs, 1 covariate, 1 gamma) and check structure
.test_minimal_data <- function() {
    data <- generate_synthetic_data(
        nobs = 2,
        ncov = 1,
        ngamma = c(1),
        seed = 2026
    )
    # Check that Xlist and ylist have correct dimensions
    if (nrow(data$Xlist[[1]]) != 2 || ncol(data$Xlist[[1]]) != 1) {
        stop("Test failed: Predictor matrix dimensions are incorrect for minimal data.")
    }
    if (length(data$ylist[[1]]) != 2) {
        stop("Test failed: Response vector length is incorrect for minimal data.")
    }

    # Check number of levels.
    if (nlevels(data)[1] != 2) {
        stop("Test failed: Number of levels is incorrect for minimal data.")
    }

    # Check levels.
    expected_level_names <- c("D1_0", "D1_1")
    if (!isTRUE(all.equal(levelNames(data)[[1]], expected_level_names))) {
        stop("Test failed: Level names are incorrect for minimal data.")
    }

    cat("Test passed: Minimal data (2 obs, 1 covariate, 1 gamma) generated correctly.\n")
}


# Test case 5: Test if the accessor functions work.
.test_accessors <- function() {
    # Create a mock mspm_data object.
    data <- new_mspm_data(
        predictorNames = c("X1", "X2"),
        responseNames = c("Y1", "Y2", "Y3"),
        Xlist = list(matrix(rnorm(100), nrow=50, ncol=2)),
        ylist = list(sample(1:3, 50, replace=TRUE),
                     sample(1:4, 50, replace=TRUE),
                     sample(1:2, 50, replace=TRUE)),
        levelNames = list(c("A", "B", "C"),
                          c("D", "E", "F", "G"),
                          c("H", "I")),
        nlevels = c(3, 4, 2),
        ntargets = 3,
        seed = 42,
        call = match.call()
    )

    # Test ntargets accessor.
    if (ntargets(data) != 3) {
        stop("Test failed: ntargets accessor returned incorrect value. Got: ", 
             ntargets(data), ", expected: 3")
    }

    # Test nlevels accessor.
    if (!isTRUE(all.equal(nlevels(data), c(3, 4, 2)))) {
        stop("Test failed: nlevels accessor returned incorrect value.")
    }

    # Test predictorNames accessor.
    if (!isTRUE(all.equal(predictorNames(data), c("X1", "X2")))) {
        stop("Test failed: predictorNames accessor returned incorrect value.")
    }

    # Test responseNames accessor.
    if (!isTRUE(all.equal(responseNames(data), c("Y1", "Y2", "Y3")))) {
        stop("Test failed: responseNames accessor returned incorrect value.")
    }

    # Test levelNames accessor.
    expected_level_names <- list(c("A", "B", "C"),
                                 c("D", "E", "F", "G"),
                                 c("H", "I"))
    if (!isTRUE(all.equal(levelNames(data), expected_level_names))) {
        stop("Test failed: levelNames accessor returned incorrect value.")
    }

    cat("Test passed: Accessor functions work correctly.\n")
}


# Test case 5: Test that two splits with the same seed produce the same results.
.test_split_replication <- function() {
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )

    split1 <- split_data(
        data = data,
        prop = 0.7,
        seed = 42
    )

    split2 <- split_data(
        data = data,
        prop = 0.7,
        seed = 42
    )

    # Check that the two splits are identical.
    for (i in 1:ntargets(data)) {
        if (!isTRUE(all.equal(split1$train$Xlist[[i]], split2$train$Xlist[[i]]))) {
            stop(paste("Test failed: Training predictor matrices for dataset", i, "are not identical across splits."))
        }
        if (!isTRUE(all.equal(split1$train$ylist[[i]], split2$train$ylist[[i]]))) {
            stop(paste("Test failed: Training response vectors for dataset", i, "are not identical across splits."))
        }
        if (!isTRUE(all.equal(split1$test$Xlist[[i]], split2$test$Xlist[[i]]))) {
            stop(paste("Test failed: Test predictor matrices for dataset", i, "are not identical across splits."))
        }
        if (!isTRUE(all.equal(split1$test$ylist[[i]], split2$test$ylist[[i]]))) {
            stop(paste("Test failed: Test response vectors for dataset", i, "are not identical across splits."))
        }
    }
    cat("Test passed: Data split is reproducible with the same seed.\n")

}
