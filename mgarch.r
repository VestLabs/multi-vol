packages <- c("Dict", "dendextend", "expm", "MCS", "mvtnorm", "parallel", "progress", "purrr", "rugarch", "rmgarch", "scoringRules", "stats")
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Loading data
returns <- c('multi-vol/data/returns_50_2023-03-04.csv', 'multi-vol/data/returns_50_2023-03-08.csv', 'multi-vol/data/returns_50_2023-03-12.csv')
prices <- c('multi-vol/data/prices_50_2023-03-04.csv', 'multi-vol/data/prices_50_2023-03-08.csv', 'multi-vol/data/prices_50_2023-03-12.csv')

# Concatenate returns and prices
get_rs_ps <- function(returns, prices) {
    rs <- list(length(returns))
    ps <- list(length(prices))
    for (i in 1:length(returns)) {
        print(returns[[i]])
        return_df <- read.csv(returns[[i]])
        price_df <- read.csv(prices[[i]])
        rs[[i]] <- return_df[, 2:ncol(return_df)]
        ps[[i]] <- price_df[nrow(price_df), 2:ncol(price_df)] |> unlist()
    }
    return(list(rs = rs, ps = ps))
}
rs_ps <- get_rs_ps(returns, prices)
returns_list <- rs_ps$rs
prices_list <- rs_ps$ps

# Multivariate GARCH models

unigarch <- function(data, p = 1, q = 1, n_sample = 100) {
    start_time <- Sys.time()
    num_assets <- ncol(data)
    spec <- ugarchspec(
        variance.model = list(model = "sGARCH", garchOrder = c(p, q)), 
        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)
    )
    mspec <- multispec(replicate(num_assets, spec))
    fit <- multifit(
        mspec, data, out.sample = n_sample - 1,
        fit.control = list(scale = 1), solver.control = list(n.restarts = 3)
    )
    forecast <- multiforecast(fit, n.ahead = 1, n.roll = n_sample - 1)
    rcovs <- list(length(n_sample))
    for (i in 1:n_sample) {
        rcov <- matrix(0, num_assets, num_assets)
        for (j in 1:num_assets) {
            rcov[j,j] <- sigma(forecast)[1,i,j]^2
        }
        rcovs[[i]] <- rcov
    }
    end_time <- Sys.time()
    result <- Dict$new(rcovs = rcovs, time = end_time - start_time)
    return(result)
}

dccgarch <- function(data, p = 1, q = 1, n_sample = 100) {
    start_time <- Sys.time()
    num_assets <- ncol(data)
    spec <- ugarchspec(
        variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)
    )
    mspec <- multispec(replicate(num_assets, spec))
    dspec = dccspec(uspec = mspec, dccOrder = c(1, 1), distribution = "mvnorm")

    # Parallelize fitting DCC
    num_cores <- detectCores()
    cl <- makeCluster(num_cores)
    fit <- dccfit(
        dspec, data, out.sample = n_sample - 1, cluster = cl,
        fit.control = list(scale = TRUE), solver.control = list(n.restarts = 3)
    )
    forecast <- dccforecast(fit, n.ahead = 1, n.roll = n_sample - 1)
    stopCluster(cl)
    
    raw_rcors <- rcor(forecast)
    raw_rcovs <- rcov(forecast)
    rcors <- list(length(n_sample))
    rcovs <- list(length(n_sample))
    for (i in 1:n_sample) {
        rcors[[i]] <- drop(raw_rcors[[i]][,,1])
        rcovs[[i]] <- drop(raw_rcovs[[i]][,,1])
    }
    end_time <- Sys.time()
    result <- Dict$new(
        rcors = rcors, rcovs = rcovs, fit = fit, time = end_time - start_time
    )
    return(result)
}

# Returns number of principal components to keep based on MP
n_comp_MP <- function(data) {
    eigenvals <- eigen(cor(data))$values # This is sorted in descending order

    # Calculate MP threshold from Marchenko-Pastur distribution
    lambda_plus <- (1 + sqrt(ncol(data)/nrow(data)))^2
    return(max(c(2, sum(eigenvals > lambda_plus))))
}

gogarch <- function(data, p = 1, q = 1, n_sample = 100, mp = FALSE) {
    start_time <- Sys.time()
    gspec <- gogarchspec(
        mean.model=list(model = "constant", include.mean=FALSE),
        variance.model=list(model = "sGARCH", garchOrder = c(p, q)), 
        distribution.model = "mvnorm",
        ica = "fastica",
    )
    n_comp <- if (mp) n_comp_MP(data) else ncol(data)
    fit <- gogarchfit(
        gspec, data, out.sample = n_sample - 1, n.comp = n_comp,
        fit.control = list(scale = TRUE), solver.control = list(n.restarts = 3)
    )

    forecast <- gogarchforecast(fit, n.ahead = 1, n.roll = n_sample - 1)
    
    raw_rcors <- rcor(forecast)
    raw_rcovs <- rcov(forecast)
    rcors <- list(length(n_sample))
    rcovs <- list(length(n_sample))
    for (i in 1:n_sample) {
        rcors[[i]] <- drop(raw_rcors[[i]][,,1])
        rcovs[[i]] <- drop(raw_rcovs[[i]][,,1])
    }
    end_time <- Sys.time()
    result <- Dict$new(
        rcors = rcors, rcovs = rcovs, fit = fit, time = end_time - start_time
    )
    return(result)
}

# Returns clustering based on statistical distance
get_clusters <- function(data, n_clusters = 3) {
    num_assets <- ncol(data)
    cor_matrix <- cor(data)
    dist_matrix <- as.dist(sqrt((1 - cor_matrix)/2))
    tree <- hclust(dist_matrix, method = "ward.D2")
    clusters <- cutree(tree, k = n_clusters)
    while (!all(sapply(1:n_clusters, \(x) (sum(clusters == x) > 1)))) {
        n_clusters <- n_clusters - 1
        clusters <- cutree(tree, k = n_clusters)
    }  

    # Plot dendrogram
    # plot(tree, cex = 0.6, hang = -1)
    # rect.hclust(tree, k = 3, cluster = clusters)

    return(clusters)
}

hgarch_helper <- function(data, n_clusters = 3, mp = FALSE, round = FALSE) {
    # Split data into clusters
    clusters <- get_clusters(data, n_clusters)

    # Apply GO-GARCH within clusters and recover factors
    num_assets <- ncol(data)
    big_rcov <- matrix(0, num_assets, num_assets)
    mixing_matrices <- list()
    factors <- list()
    factor_indices <- list()

    for (i in 1:n_clusters) {
        cluster_data <- data[, clusters == i]
        cluster_idx <- which(clusters == i)

        cluster_gogarch <- gogarch(cluster_data, n_sample=1, mp = mp)
        
        big_rcov[cluster_idx, cluster_idx] <- cluster_gogarch$get("rcovs")[[1]]
        
        gogarch_fit <- cluster_gogarch$get("fit")
        mixing_matrices[[i]] <- rmgarch::as.matrix(gogarch_fit, which = "A")

        unmixing_matrix <- rmgarch::as.matrix(gogarch_fit, which = "W")
        factors[[i]] <- as.matrix(cluster_data) %*% t(unmixing_matrix)
        factor_indices[[i]] <- ifelse(i == 1, 1, factor_indices[[i - 1]] + ncol(factors[[i - 1]]))
    }

    # Apply DCC to factors
    factors_dccgarch <- dccgarch(do.call(cbind, factors), n_sample=1)
    factors_rcov <- factors_dccgarch$get("rcovs")[[1]]
    for (i in 1:n_clusters) {
        cluster_i_idx <- which(clusters == i)
        for (j in (i + 1):n_clusters) {
            if (j > n_clusters) {
                break
            }
            cluster_j_idx <- which(clusters == j)
            cov_ij <- factors_rcov[factor_indices[[i]]:(factor_indices[[i]] + ncol(factors[[i]]) - 1), factor_indices[[j]]:(factor_indices[[j]] + ncol(factors[[j]]) - 1)]

            A <- mixing_matrices[[i]]
            B <- mixing_matrices[[j]]

            # insert inter-sector covariance matrix into big covariance matrix
            cov_ij <- A %*% cov_ij %*% t(B)
            big_rcov[cluster_i_idx, cluster_j_idx] <- cov_ij
            big_rcov[cluster_j_idx, cluster_i_idx] <- t(cov_ij)
        }
    }

    psd_condition <- all(eigen(big_rcov)$values >= -10^(-8))
    if (!psd_condition) {
        print("Not PSD :(")
    }

    # Round inter-factor correlation coefficient to enforce PSD condition
    if (!psd_condition && round) {
        factors_rcor <- factors_dccgarch$get("rcors")[[1]]
        total_factors <- sapply(factors, \(x) ncol(x)) |> sum()
        for (i in 1:n_clusters) {
            cluster_i_idx <- which(clusters == i)
            for (j in (i + 1):n_clusters) {
                if (j > n_clusters) {
                    break
                }
                cluster_j_idx <- which(clusters == j)
                cor_ij <- factors_rcor[factor_indices[[i]]:(factor_indices[[i]] + ncol(factors[[i]]) - 1), factor_indices[[j]]:(factor_indices[[j]] + ncol(factors[[j]]) - 1)]
                cov_ij <- factors_rcov[factor_indices[[i]]:(factor_indices[[i]] + ncol(factors[[i]]) - 1), factor_indices[[j]]:(factor_indices[[j]] + ncol(factors[[j]]) - 1)]
                rounded <- FALSE
                for (k in 1:ncol(factors[[i]])) {
                    for (l in 1:ncol(factors[[j]])) {
                        max_corr <- sqrt(1 / ((total_factors - ncol(factors[[i]])) * (total_factors - ncol(factors[[j]]))))
                        if (abs(cor_ij[k, l]) > max_corr) {
                            new_cov <- cov_ij[k, l] * max_corr / abs(cor_ij[k, l])
                            cov_ij[k, l] <- new_cov
                            cov_ij[l, k] <- new_cov
                            rounded <- TRUE
                        }
                    }
                }
                if (rounded) {
                    A <- mixing_matrices[[i]]
                    B <- mixing_matrices[[j]]
                    cov_ij <- A %*% cov_ij %*% t(B)
                    big_rcov[cluster_i_idx, cluster_j_idx] <- cov_ij
                    big_rcov[cluster_j_idx, cluster_i_idx] <- t(cov_ij)
                }
            }
        }
    }
    
    return(big_rcov)
}

hgarch <- function(data, n_clusters = 3, n_sample=100, mp = FALSE, round = FALSE) {
    start_time <- Sys.time()
    rcovs <- list()
    pb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta", total = n_sample)
    for (i in 1:n_sample) {
        subdata <- data[i:(nrow(data)-(n_sample - i)),]
        rcovs[[i]] <- hgarch_helper(subdata, n_clusters=n_clusters, mp = mp, round = round)
        pb$tick()
    }
    end_time <- Sys.time()
    result <- Dict$new(rcovs = rcovs, time = end_time - start_time)
    return(result)
}

# Scoring rules
calculate_score <- function(data, n_sample, model, scoring_rule, T, steps_ahead) {
    data_test <- data[(T + 1):(T + n_sample + steps_ahead - 1), ]
    model_fit <- model(data[1:(T + n_sample - 1), ], n_sample = n_sample)
    forecasts <- model_fit$get("rcovs")
    scores <- replicate(n_sample, 0)
    for (i in 1:n_sample) {
        scores[i] <- scoring_rule(data_test[i:(i + steps_ahead - 1), ], forecasts[[i]])
    }
    return(scores)
}

energy_score <- function(data, model, T, n_sample, steps_ahead) {
    scoring_rule <- function(observations, forecast_cov_matrix, steps_ahead) {
        samples <- t(rmvnorm(n = 500, sigma = steps_ahead * forecast_cov_matrix))
        return(es_sample(y = as.numeric(observations[steps_ahead, ]), dat = samples))
    }
    return(calculate_score(data, n_sample, model, partial(scoring_rule, steps_ahead = steps_ahead), T, steps_ahead))
}

frobenius_loss <- function(data, model, T, n_sample, steps_ahead) {
    scoring_rule <- function(observations, forecast_cov_matrix, steps_ahead) {
        true_cov <- matrix(0, ncol(observations), ncol(observations))
        for (i in 1:steps_ahead) {
            true_cov <- true_cov + (t(as.matrix(observations[i, ])) %*% as.matrix(observations[i, ]))
        }
        return(norm(steps_ahead * forecast_cov_matrix - true_cov, type = "F"))
    }
    return(calculate_score(data, n_sample, model, partial(scoring_rule, steps_ahead = steps_ahead), T, steps_ahead))}

sse_loss <- function(data, model, T, n_sample, steps_ahead) {
    scoring_rule <- function(observations, forecast_cov_matrix, steps_ahead) {
        true_cov <- matrix(0, ncol(observations), ncol(observations))
        for (i in 1:steps_ahead) {
            true_cov <- true_cov + (t(as.matrix(observations[i, ])) %*% as.matrix(observations[i, ]))
        }
        return(sum((steps_ahead * forecast_cov_matrix - true_cov)^2))
    }
    return(calculate_score(data, n_sample, model, partial(scoring_rule, steps_ahead = steps_ahead), T, steps_ahead))}

mse_loss <- function(data, model, T, n_sample, steps_ahead) {
    scoring_rule <- function(observations, forecast_cov_matrix, steps_ahead) {
        true_cov <- matrix(0, ncol(observations), ncol(observations))
        for (i in 1:steps_ahead) {
            true_cov <- true_cov + (t(as.matrix(observations[i, ])) %*% as.matrix(observations[i, ]))
        }
        return(mean((steps_ahead * forecast_cov_matrix - true_cov)^2))
    }
    return(calculate_score(data, n_sample, model, partial(scoring_rule, steps_ahead = steps_ahead), T, steps_ahead))}

# MCS
calculate_mcs <- function(data, models, scoring_rule, T = 864, n_sample = 100, steps_ahead = 1, alpha = 0.1) {
    num_models <- length(models)
    loss_matrix <- matrix(0, nrow = n_sample, ncol = num_models)

    for (j in 1:num_models) {
        loss_matrix[,j] <- scoring_rule(data, models[[j]], T, n_sample, steps_ahead)
        print(sprintf("model %f: %f", j, mean(loss_matrix[, j])))
    }
    
    # Parallelize running mcs
    num_cores <- detectCores()
    cl <- makeCluster(num_cores)
    mcs <- MCSprocedure(Loss = loss_matrix, alpha = alpha, stat = "Tmax", cl = cl)
    stopCluster(cl)

    return(mcs)
}

# Example model evaluation
models <- list(
    unigarch,
    partial(dccgarch),
    gogarch,
    partial(gogarch, mp = TRUE),
    partial(hgarch, mp = TRUE, n_clusters = 2, round = TRUE)
)
mcs <- calculate_mcs(
    returns_list[[1]], models, energy_score,
    T = 864, n_sample = 100, steps_ahead = 12
)

# EVaR approximation
evar_multinorm_approx <- function(prices, mu, sigma, alpha, imbalances, cum_entry, new_sizes, lp = 0) {
    price_mu <- prices * mu 
    price_sigma <- diag(prices) %*% sigma %*% diag(prices)
    risk <- t(as.matrix(price_mu)) %*% as.matrix(imbalances + new_sizes) + sqrt(-2 * log(alpha) * t(as.matrix(imbalances + new_sizes)) %*% price_sigma %*% as.matrix(imbalances + new_sizes)) - (lp+cum_entry + t(as.matrix(prices)) %*% as.matrix(new_sizes))
    return(risk)
}

evar_premium <- function(prices, mu, sigma, alpha, imbalances, cum_entry, new_sizes, lp = 0) {
    num_assets <- length(prices)
    premium <- evar_multinorm_approx(prices, mu, sigma, alpha, imbalances, cum_entry, new_sizes, lp) - evar_multinorm_approx(prices, mu, sigma, alpha, imbalances, cum_entry, replicate(num_assets, 0), lp)
    return(premium)
}

evar_premium_single <- function(prices, mu, sigma, alpha, imbalances, cum_entry, index, new_size, lp = 0, absolute=FALSE) {
    num_assets <- length(prices)
    new_sizes <- replicate(num_assets, 0)
    new_sizes[index] <- new_size
    premium <- evar_premium(prices, mu, sigma, alpha, imbalances, cum_entry, new_sizes, lp)
    if (absolute) {
        return(premium)
    } else {
        return(ifelse(new_size != 0, prices[index] + premium / new_size, NaN))
    }
}

# Example pricing
sample_returns <- read.csv('multi-vol/data/returns_50_2023-03-20.csv')
sample_returns <- sample_returns[, 2:ncol(sample_returns)]
sample_prices <- read.csv('multi-vol/data/prices_50_2023-03-20.csv')
sample_prices <- sample_prices[nrow(sample_prices), 2:ncol(sample_prices)] |> unlist()

num_assets <- ncol(sample_returns)
mu <- replicate(num_assets, 1)
# Calculate risk as the potential loss against traders in the next hour
approx_sigma <- function(data, model) {
    num_assets <- ncol(data)
    rcov <- model(data[(nrow(data) - 863):nrow(data), ], n_sample = 1)$get("rcovs")[[1]]
    return(expm(rcov) - diag(num_assets))
}
uni_sigma <- approx_sigma(sample_returns, unigarch)
dcc_sigma <- approx_sigma(sample_returns, dccgarch)
go_sigma <- approx_sigma(sample_returns, partial(gogarch, mp = TRUE))
h_sigma <- approx_sigma(sample_returns, partial(hgarch, mp = TRUE, round = TRUE, n_clusters = 2))
alpha <- 0.005
btc_index <- 12
eth_index <- 21
sizes <- seq(-100, 100, by = 1)

model_sigma <- Dict$new(uni = uni_sigma, dcc = dcc_sigma, go = go_sigma, h = h_sigma)
imbalance <- replicate(num_assets, 0)
random_imbalance <- function(p) {
    random_imbalance <- replicate(num_assets, 0)
    for (i in 1:num_assets) {
        random_imbalance[i] <- (2 * rbinom(1, 1, p) - 1) * rexp(1, 1) * 1000 / sample_prices[i]
    }
    random_imbalance <- unlist(random_imbalance)
    return(random_imbalance)
}

mark_prices <- data.frame(sizes = sizes, index = sample_prices[eth_index])
imbalances <- Dict$new(
    eth_long_heavy = replace(imbalance, eth_index, 20), 
    eth_short_heavy = replace(imbalance, eth_index, -20),
    btc_long_heavy = replace(imbalance, btc_index, 2), 
    btc_short_heavy = replace(imbalance, btc_index, -2),
    random_long_heavy = random_imbalance(0.8),
    random_short_heavy = random_imbalance(0.2)
)

for (imb_id in imbalances$keys) {
    cum_entry <- imbalances$get(imb_id) %*% sample_prices
    for (model in model_sigma$keys) {
        id <- sprintf("%s_%s", model, imb_id)
        mark_prices[id] <- sapply(sizes, \(x) evar_premium_single(sample_prices, mu, model_sigma$get(model), alpha, imbalances$get(imb_id), cum_entry, eth_index, x))
    }
}

