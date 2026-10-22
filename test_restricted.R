require(tidyverse)

res <- c()
for (model in 1:3) {
  file_name <- paste0("res_restrict_", model, ".txt")
  vx <- read_tsv(file_name)
  res <- rbind(res ,vx)
}




res$model <- as.factor(res$model)
df <- res %>%
  dplyr::select(-c("extinct", "sym_high", "sym_low", "allo", "jiggle", "water_rate", "model", "crown_age", "tree_height"))

to_remove <- c()
for (i in 1:ncol(df)) {
  av <- which(is.na(df[, i]))
  if (length(av) > 0) {
    cat(colnames(df)[i], "\n")
    to_remove <- c(to_remove, av)
  }
}
to_remove <- unique(to_remove)

if (length(to_remove)) df <- df[-to_remove, ]

vz <- prcomp(scale(df))

require(ggfortify)

if (length(to_remove)) {
  res2 <- res[-to_remove, ]
} else {
  res2 <- res
}

autoplot(vz,
         data = res2,
         colour = 'model',
  #   loadings = TRUE,
  #   loadings.label = TRUE,
     frame = TRUE,
     frame.type = 'norm') +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(colour = "") +
  labs(fill = "")

library(factoextra)
df2 <- scale(df)
data.pca <- princomp(df2)
summary(data.pca)
fviz_eig(data.pca, addlabels = TRUE)
fviz_pca_var(data.pca, col.var = "black")
fviz_cos2(data.pca, choice = "var", axes = 1:4)
fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

training <- res %>%
  dplyr::select(-c("extinct", "sym_high", "sym_low", "allo", "jiggle", "water_rate", "crown_age", "tree_height"))

training <- training %>%
  mutate_at(2:69, scale)

require(MASS)
linear <- MASS::lda(model~., training)
linear
wine.lda.values <- predict(linear)
ldahist(wine.lda.values$x[,2], g = model)
plot(wine.lda.values$x[,1], wine.lda.values$x[,2])

newdata <- data.frame(type = training[,1], lda = wine.lda.values$x)
library(ggplot2)
ggplot(newdata) + geom_point(aes(lda.LD1, lda.LD2, colour = model), size = 2.5)

to_keep <- 1:length(res[1, ])
to_remove <- c(1:6) # parameters
to_remove <- c(to_remove, which(colnames(res) == "crown_age"))
to_remove <- c(to_remove, which(colnames(res) == "tree_height"))

to_keep <- to_keep[-to_remove]
for_abc <- res[, to_keep]

for_abc <- subset(for_abc, for_abc$model != 4)

forest <- abcrf::abcrf(model ~ ., data = for_abc, ntree = 1000, paral = TRUE, lda = FALSE)
forest
abcrf::variableImpPlot(forest)

test <- ape::rphylo(n = 10, 1, 0)
brtstats <- treestats::calc_brts_stats(test)

to_keep <- c(names(brtstats), "model")

for_abc2 <- res[, to_keep]
forest2 <- abcrf::abcrf(model ~ ., data = for_abc2, ntree = 1000, paral = TRUE)
abcrf::variableImpPlot(forest2)
forest2
