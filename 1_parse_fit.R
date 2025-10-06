library(PARSE)
source("utility.r")
set.seed(1234)


# Data load
source("0 data loading.R")

# Parse Fit
K = 4
parse.fit = parse(K = K, lambda = c(0.1, 0.3, 0.5), y = data_norm,  cores = 8, pca_adjust = 0.5)
saveRDS(parse.fit, paste0(c("islet/parse_islet_K",K,".rds"),collapse = ""))


# NULL Data generation
data_norm_null = constructNull_Guassian(mat = t(data_norm))
saveRDS(data_norm_null, "islet/data_norm_null.rds")


# PARSE NULL fit
parse.fit.null = parse(K = K, lambda = parse.fit$lambda.best, y = data_norm_null, cores = 8, pca_adjust = 0.5)
saveRDS(parse.fit.null, paste0(c("islet/parse_islet_null_K",K,".rds"),collapse = ""))


