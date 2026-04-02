# install.packages(c("gradientForest","ggplot2"))  # if needed
library(gradientForest)
library(ggplot2)

df <- read.csv("resistance_imputed_data.csv", check.names = FALSE)

env_cols <- c(
  "bio_17","bio_10","bio_01","bio_04","bio_15","proportion_agricultural_1000m","shannon_idx_1000m","proportion_natural_1000m","PPO inhibitor (Group 14)", "GLYPHOSATE", "ALS inhibitor (Group 2)"
)
# snp_cols <- grep("^AgCMH_", colnames(df), value = TRUE)
# snp_cols <- c(
#     "EPSPSPro106Ser_6_17480274"
# )
snp_cols <- c(
     "EPSPSPro106Ser_6_17480274", "ALSAla205Asp_8_10192948","ALSAsp376Glu_8_10193474","ALSTrp574Leu_8_10194067","ALSSer653Asn_8_10194301","PPO210Del_8_10446786"
)
if (length(env_cols) == 0) stop("No environmental columns found in gf_input.csv using the specified names.")
if (length(snp_cols) == 0) stop("No DroughtAnc_ SNP columns found in gf_input.csv.")

# Create mapping for column names (original -> valid R names)
name_mapping <- setNames(make.names(env_cols), env_cols)
# Rename columns in dataframe for gradientForest
df_renamed <- df
colnames(df_renamed)[match(env_cols, colnames(df_renamed))] <- name_mapping
env_cols_renamed <- name_mapping

# Create semantic names for display in plots
better_names <- c(
  "bio_17" = "Precipitation of Driest Quarter (BIO17)",
  "bio_10" = "Mean Temperature of Warmest Quarter (BIO10)",
  "bio_01" = "Annual Mean Temperature (BIO1)",
  "bio_04" = "Temperature Seasonality (BIO4)",
  "bio_15" = "Precipitation Seasonality (BIO15)",
  "proportion_agricultural_1000m" = "Proportion Agricultural Land Use (1000m radius)",
  "shannon_idx_1000m" = "Shannon Diversity Index (1000m radius)",
  "proportion_natural_1000m" = "Proportion of Natural Land Use (1000m radius)",
  "PPO inhibitor (Group 14)" = "PPO Inhibitors applied (county-level kg)",
  "GLYPHOSATE" = "Glyphosate applied(county-level kg)",
  "ALS inhibitor (Group 2)" = "ALS Inhibitors applied (county-level kg)"
)
# Map to semantic names in the same order as env_cols
env_cols_semantic <- better_names[env_cols]
            
used <- c(env_cols_renamed, snp_cols)
df_renamed <- df_renamed[complete.cases(df_renamed[, used]), ]
for (c in snp_cols) {
  df_renamed[[c]] <- pmin(pmax(df_renamed[[c]], 0), 2)
}

# compute max rf level
maxLevel <- log2(0.368 * nrow(df_renamed) / 2)

set.seed(42)
gf <- gradientForest(
  data = df_renamed[, used],
  predictor.vars = env_cols_renamed,
  response.vars  = snp_cols,
  ntree = 500,
  maxLevel = maxLevel,
  corr.threshold = 0.50,
  trace = TRUE
)

# Get number of environmental variables
n_env <- length(env_cols)

# Calculate appropriate plot layout (rows x cols)
n_cols <- 4  # 4 columns
n_rows <- ceiling(n_env / n_cols)  # Calculate rows needed

png("resistance_all_gf_cumulative_importance_all_vars.png", 
    width = 1600, height = 400 * n_rows, res = 130)

# Set up graphics parameters with space for title
par(mfrow = c(n_rows, n_cols),
    mgp = c(1.5, 0.5, 0),
    mar = c(3, 4.5, 0.1, 0.5),         # Margins: bottom, left, top, right
    omi = c(0.3, 0.5, 1.0, 0.3))        # Outer margins: bottom, left, top, right (increased top for title)

plot(gf, 
     plot.type = "Cumulative.Importance",
     par.args = list(
       mfrow = c(n_rows, n_cols),  
       mgp = c(1.5, 0.5, 0),
       mar = c(3, 4.5, 0.1, 0.5),         # Margins: bottom, left, top, right
       omi = c(0.3, 0.5, 1.0, 0.3)  # Outer margins: bottom, left, top, right (increased top for title)
     ),

     plot.args = list(
       show.species = FALSE,        # don't show individual SNP curves
       show.overall = TRUE,
       common.scale = FALSE,
       cex.axis = 0.6,
       cex.lab = 0.7,
       line.ylab = 1.9
     ),

     imp.vars = env_cols_renamed,
     imp.vars.names = env_cols_semantic 
)

# Add main title using mtext for more reliable placement
mtext("Cumulative Importance of Environmental Predictors on All Resistance SNPs", 
      side = 3, 
      outer = TRUE, 
      cex = 1.2,
      line = 3,
      font = 2)

dev.off()

png("agri_gf_performance.png", 
    width = 1600, height = 400 * n_rows, res = 130)

par(omi = c(0.3, 0.3, 0.8, 0.3))  # Outer margins: bottom, left, top, right (increased top for title)

plot(gf, plot.type = "P", show.names = F, horizontal = F, cex.axis = 1, cex.labels = 0.7, line = 2.5)

# Add main title
title(main = "Gradient Forest Model Performance", 
      outer = TRUE, 
      cex.main = 1.2,
      line = 2)

dev.off()

env_only <- df_renamed[, env_cols_renamed, drop = FALSE]
env_gf <- predict(gf, env_only)   # gradientForest transform
pc <- prcomp(env_gf, center = TRUE, scale. = TRUE)

ord <- data.frame(
  PC1 = pc$x[,1],
  PC2 = pc$x[,2]
)

png("agri_gf_turnover_ordination.png", width = 1000, height = 800, res = 130)
print(
  ggplot(ord, aes(PC1, PC2)) +
    geom_point(alpha = 0.8) +
    labs(title = "Genomic turnover (GF-transformed env → PCA)",
         subtitle = "Each point is a population (seq_names); distances reflect predicted turnover",
         x = "PC1", y = "PC2") +
    theme_minimal(base_size = 13)
)
dev.off()