## ============================================================
## 0. Setup
## ============================================================

library(tidyverse)
library(ggplot2)

gex <- readRDS("data/raw/pt16_emory_GEX_immune_FULL_v2.rds")

## 구조 확인
names(gex)
str(gex, max.level = 1)

## 필요한 객체 꺼내기
X_cat  <- gex$X_cat
Y      <- gex$Y
meta   <- gex$meta
M_expr <- gex$M_expr

## 분석용 data frame 만들기
eda_df <- meta %>%
  as.data.frame() %>%
  mutate(
    X_cat = factor(X_cat, levels = c("inside", "adjacent", "far")),
    Y = as.numeric(Y)
  )

## 좌표 변수명 확인
names(eda_df)

## ============================================================
## 0-1. Coordinate variable setting
## ============================================================

## 아래 두 줄은 실제 meta 변수명에 맞게 수정
## 예: coord_x <- "pxl_col_in_fullres"
##     coord_y <- "pxl_row_in_fullres"

coord_x <- "imagecol"
coord_y <- "imagerow"

eda_df <- eda_df %>%
  mutate(
    x_coord = .data[[coord_x]],
    y_coord = .data[[coord_y]]
  )

## ============================================================
## 1. Matrix-level weighted network summary
## ============================================================
M <- M_expr
weight_vec <- as.numeric(M)


##### PCA #####
## PCA
M_scaled <- scale(M_expr)

pca_fit <- prcomp(M_scaled, center = FALSE, scale. = FALSE)

## ============================================================
## Final PCA axis orientation
##
## PC1_R: stromal / ECM / humoral immune direction is positive
## PC2_R: antigen-presentation / adaptive immune direction is positive
## PC3  : mitochondrial / metabolic / plasma-cell direction is positive
## ============================================================

pca_df <- eda_df %>%
  mutate(
    PC1_original = pca_fit$x[, 1],
    PC2_original = pca_fit$x[, 2],
    PC3_original = pca_fit$x[, 3],
    
    PC1_R = -PC1_original,
    PC2_R = -PC2_original,
    PC3   =  PC3_original
  )



## Scree plot
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)

scree_df <- tibble(
  PC = seq_along(var_explained),
  var_explained = var_explained,
  cum_var_explained = cumsum(var_explained)
)

ggplot(scree_df[1:30, ], aes(x = PC, y = var_explained)) +
  geom_col() +
  labs(
    title = "Scree plot",
    x = "Principal component",
    y = "Proportion of variance explained"
  ) +
  theme_bw(base_size = 13)
## ============================================================
## Final loading orientation
## ============================================================

loading_df <- as.data.frame(pca_fit$rotation[, 1:3]) %>%
  rownames_to_column("gene") %>%
  rename(
    PC1_original = PC1,
    PC2_original = PC2,
    PC3_original = PC3
  ) %>%
  mutate(
    PC1_R = -PC1_original,
    PC2_R = -PC2_original,
    PC3   =  PC3_original
  )

## PC1_R
pc1_pos <- loading_df %>%
  arrange(desc(PC1_R)) %>%
  slice_head(n = 150) %>%
  pull(gene)

pc1_neg <- loading_df %>%
  arrange(PC1_R) %>%
  slice_head(n = 150) %>%
  pull(gene)

## PC2_R
pc2_pos <- loading_df %>%
  arrange(desc(PC2_R)) %>%
  slice_head(n = 100) %>%
  pull(gene)

pc2_neg <- loading_df %>%
  arrange(PC2_R) %>%
  slice_head(n = 100) %>%
  pull(gene)

## PC3
pc3_pos <- loading_df %>%
  arrange(desc(PC3)) %>%
  slice_head(n = 100) %>%
  pull(gene)

pc3_neg <- loading_df %>%
  arrange(PC3) %>%
  slice_head(n = 100) %>%
  pull(gene)

cat("PC1_R positive:\n")
cat(pc1_pos, sep = " ")

cat("\n\nPC1_R negative:\n")
cat(pc1_neg, sep = " ")

cat("\n\nPC2_R positive:\n")
cat(pc2_pos, sep = " ")

cat("\n\nPC2_R negative:\n")
cat(pc2_neg, sep = " ")

cat("\n\nPC3 positive:\n")
cat(pc3_pos, sep = " ")

cat("\n\nPC3 negative:\n")
cat(pc3_neg, sep = " ")
## ============================================================
## Spatial score plot function
## ============================================================
library(dplyr)
library(ggplot2)
library(sf)
library(RANN)

#--------------------------------------------------
# 1. plotting 좌표 만들기
#    scale_y_reverse()를 쓰지 않고 y 자체를 뒤집음
#--------------------------------------------------

plot_df <- pca_df |>
  mutate(
    x_plot = x_coord,
    y_plot = -y_coord
  )


#--------------------------------------------------
# 2. spot 사이의 실제 최근접 거리 계산
#--------------------------------------------------

xy_mat <- as.matrix(plot_df[, c("x_plot", "y_plot")])

nn_result <- RANN::nn2(
  data = xy_mat,
  query = xy_mat,
  k = 2
)

# 첫 번째는 자기 자신, 두 번째가 가장 가까운 spot
nearest_distance <- median(nn_result$nn.dists[, 2])

nearest_distance

#--------------------------------------------------
# 3. 각 spot을 sf point로 변환하고 buffer 생성
#--------------------------------------------------

spot_sf <- st_as_sf(
  plot_df,
  coords = c("x_plot", "y_plot"),
  remove = FALSE,
  crs = NA
)

# spot 간격 절반보다 조금 크게 설정해 인접 spot들이 연결되도록 함
buffer_radius <- nearest_distance * 0.58

spot_buffer_sf <- st_buffer(
  spot_sf,
  dist = buffer_radius
)


#--------------------------------------------------
# 4. X_cat별 spot 영역을 합쳐 category polygon 생성
#--------------------------------------------------

category_polygon_sf <- spot_buffer_sf |>
  group_by(X_cat) |>
  summarise(
    geometry = st_union(geometry),
    .groups = "drop"
  ) |>
  st_make_valid()

plot_pc_spatial <- function(
    data,
    pc_var,
    variance_explained,
    pc_label,
    boundary_sf
) {
  
  x_breaks <- pretty(data$x_plot, n = 5)
  y_breaks_original <- pretty(abs(data$y_plot), n = 5)
  
  ggplot() +
    geom_point(
      data = data,
      aes(
        x = x_plot,
        y = y_plot,
        color = .data[[pc_var]]
      ),
      size = 2.2,
      alpha = 0.9
    ) +
    geom_sf(
      data = boundary_sf,
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.35
    ) +
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = pc_label
    ) +
    scale_x_continuous(
      breaks = x_breaks
    ) +
    scale_y_continuous(
      breaks = -y_breaks_original,
      labels = function(z) abs(z)
    ) +
    coord_sf(
      datum = NA,
      expand = FALSE
    ) +
    labs(
      title = paste0(
        pc_label,
        " spatial score map\n(",
        sprintf("%.2f", 100 * variance_explained),
        "% variance explained)"
      ),
      x = "x coordinate",
      y = "y coordinate"
    ) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)

p_pc1 <- plot_pc_spatial(
  data = plot_df,
  pc_var = "PC1_R",
  variance_explained = var_explained[1],
  pc_label = "PC1_R",
  boundary_sf = category_polygon_sf
)

p_pc2 <- plot_pc_spatial(
  data = plot_df,
  pc_var = "PC2_R",
  variance_explained = var_explained[2],
  pc_label = "PC2_R",
  boundary_sf = category_polygon_sf
)

p_pc3 <- plot_pc_spatial(
  data = plot_df,
  pc_var = "PC3",
  variance_explained = var_explained[3],
  pc_label = "PC3",
  boundary_sf = category_polygon_sf
)

p_pc1
p_pc2
p_pc3


## ------------------------------------------------------------
## Helper function: linear product-of-coefficients decomposition
## ------------------------------------------------------------
## This is useful because it directly shows:
## X -> PC effect, PC -> Y effect, NIE = alpha * beta.
## For multiple PCs, NIE is sum_k alpha_k * beta_k.

decompose_linear <- function(data, exposure, outcome, mediators, covariates = NULL) {
  
  ## formulas
  cov_part <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }
  
  ## total effect model: Y ~ X (+ C)
  te_formula <- as.formula(
    paste(
      outcome, "~", exposure,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  
  fit_te <- lm(te_formula, data = data)
  TE <- coef(fit_te)[exposure]
  
  ## mediator models: M_k ~ X (+ C)
  alpha_tbl <- map_dfr(mediators, function(med) {
    m_formula <- as.formula(
      paste(
        med, "~", exposure,
        if (!is.null(cov_part)) paste("+", cov_part) else ""
      )
    )
    
    fit_m <- lm(m_formula, data = data)
    
    tibble(
      mediator = med,
      alpha_X_to_M = coef(fit_m)[exposure]
    )
  })
  
  ## outcome model: Y ~ X + M_1 + ... + M_K (+ C)
  y_formula <- as.formula(
    paste(
      outcome, "~", exposure, "+", paste(mediators, collapse = " + "),
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  
  fit_y <- lm(y_formula, data = data)
  
  beta_tbl <- tibble(
    mediator = mediators,
    beta_M_to_Y = coef(fit_y)[mediators]
  )
  
  path_tbl <- alpha_tbl %>%
    left_join(beta_tbl, by = "mediator") %>%
    mutate(
      indirect_component = alpha_X_to_M * beta_M_to_Y
    )
  
  NIE <- sum(path_tbl$indirect_component)
  NDE <- coef(fit_y)[exposure]
  PM  <- NIE / TE
  
  summary_tbl <- tibble(
    mediators = paste(mediators, collapse = " + "),
    covariates = ifelse(is.null(covariates) || length(covariates) == 0,
                        "none",
                        paste(covariates, collapse = " + ")),
    TE = TE,
    NDE = NDE,
    NIE = NIE,
    PM = PM
  )
  
  list(
    summary = summary_tbl,
    path = path_tbl,
    fit_total = fit_te,
    fit_outcome = fit_y
  )
}


###### full ######
## ============================================================
## Full 3-level exposure data
## Baseline/reference = inside
## X_adjacent = 1 if adjacent, 0 otherwise
## X_far      = 1 if far,      0 otherwise
## ============================================================

## ============================================================
## Analysis dataset using final oriented PCs
## ============================================================

dat_full_allpc <- pca_df %>%
  mutate(
    X_cat = factor(
      X_cat,
      levels = c("inside", "adjacent", "far")
    ),
    X_adjacent = as.integer(X_cat == "adjacent"),
    X_far      = as.integer(X_cat == "far")
  ) %>%
  select(
    Y,
    X_cat,
    X_adjacent,
    X_far,
    PC1_R,
    PC2_R,
    PC3,
    x_coord,
    y_coord
  ) %>%
  drop_na()

table(dat_full_allpc$X_cat)

## ============================================================
## Helper function: decomposition for multi-level categorical exposure
## with dummy-coded exposure variables
## ============================================================

decompose_linear_multix <- function(data,
                                    exposures,
                                    outcome,
                                    mediators,
                                    covariates = NULL) {
  
  cov_part <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }
  
  exposure_part <- paste(exposures, collapse = " + ")
  mediator_part <- paste(mediators, collapse = " + ")
  
  ## ------------------------------------------------------------
  ## Total effect model:
  ## Y ~ X_adjacent + X_far (+ C)
  ## ------------------------------------------------------------
  te_formula <- as.formula(
    paste(
      outcome, "~", exposure_part,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  
  fit_te <- lm(te_formula, data = data)
  te_coef <- coef(fit_te)[exposures]
  
  ## ------------------------------------------------------------
  ## Mediator models:
  ## M_k ~ X_adjacent + X_far (+ C)
  ## ------------------------------------------------------------
  alpha_tbl <- map_dfr(mediators, function(med) {
    
    m_formula <- as.formula(
      paste(
        med, "~", exposure_part,
        if (!is.null(cov_part)) paste("+", cov_part) else ""
      )
    )
    
    fit_m <- lm(m_formula, data = data)
    
    tibble(
      mediator = med,
      exposure = exposures,
      alpha_X_to_M = coef(fit_m)[exposures]
    )
  })
  
  ## ------------------------------------------------------------
  ## Outcome model:
  ## Y ~ X_adjacent + X_far + PC1 + PC2 + PC3 (+ C)
  ## ------------------------------------------------------------
  y_formula <- as.formula(
    paste(
      outcome, "~", exposure_part, "+", mediator_part,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  
  fit_y <- lm(y_formula, data = data)
  
  beta_tbl <- tibble(
    mediator = mediators,
    beta_M_to_Y = coef(fit_y)[mediators]
  )
  
  path_tbl <- alpha_tbl %>%
    left_join(beta_tbl, by = "mediator") %>%
    mutate(
      indirect_component = alpha_X_to_M * beta_M_to_Y
    )
  
  ## ------------------------------------------------------------
  ## Summarize by exposure contrast
  ## ------------------------------------------------------------
  summary_tbl <- path_tbl %>%
    group_by(exposure) %>%
    summarise(
      NIE = sum(indirect_component),
      .groups = "drop"
    ) %>%
    mutate(
      TE = te_coef[exposure],
      NDE = coef(fit_y)[exposure],
      PM = NIE / TE,
      mediators = paste(mediators, collapse = " + "),
      covariates = ifelse(
        is.null(covariates) || length(covariates) == 0,
        "none",
        paste(covariates, collapse = " + ")
      )
    ) %>%
    select(exposure, mediators, covariates, TE, NDE, NIE, PM)
  
  list(
    summary = summary_tbl,
    path = path_tbl,
    fit_total = fit_te,
    fit_outcome = fit_y
  )
}



## ============================================================
## Full 3-level exposure model: PC1-PC3, with coordinate controls
## ============================================================
res_full_pc123_coord <- decompose_linear_multix(
  data = dat_full_allpc,
  exposures = c("X_adjacent", "X_far"),
  outcome = "Y",
  mediators = c("PC1_R", "PC2_R", "PC3"),
  covariates = c("x_coord", "y_coord")
)

res_full_pc123_coord$summary
res_full_pc123_coord$path
summary(res_full_pc123_coord$fit_outcome)
## ============================================================
## Full 3-level exposure model: PC1-PC3, no controls
## ============================================================
res_full_pc123_nocov <- decompose_linear_multix(
  data = dat_full_allpc,
  exposures = c("X_adjacent", "X_far"),
  outcome = "Y",
  mediators = c("PC1_R", "PC2_R", "PC3"),
  covariates = NULL
)

res_full_pc123_nocov$summary
res_full_pc123_nocov$path
## ============================================================
## Full 3-level path contribution table
## ============================================================

full3_path_contribution_table <- bind_rows(
  res_full_pc123_nocov$path %>%
    mutate(model = "PC1-PC3, no controls"),
  res_full_pc123_coord$path %>%
    mutate(model = "PC1-PC3, coordinate controls")
) %>%
  mutate(
    contrast = case_when(
      exposure == "X_adjacent" ~ "adjacent vs inside",
      exposure == "X_far" ~ "far vs inside",
      TRUE ~ exposure
    ),
    across(c(alpha_X_to_M, beta_M_to_Y, indirect_component), ~ round(.x, 5))
  ) %>%
  select(contrast, model, mediator, alpha_X_to_M, beta_M_to_Y, indirect_component)

full3_path_contribution_table


##### category and covariates #####
## Ensure reference category is inside
dat_full_allpc$X_cat <- factor(dat_full_allpc$X_cat, levels = c("inside", "adjacent", "far"))

## Design matrix for exposure + coordinates
X_design <- model.matrix(~ X_cat + x_coord + y_coord, data = dat_full_allpc)

## Check dimensions and rank
dim(X_design)
qr(X_design)$rank
ncol(X_design)

## Full rank?
qr(X_design)$rank == ncol(X_design)
colnames(X_design)
fit_rank_check <- lm(
  rnorm(nrow(dat_full_allpc)) ~ X_cat + x_coord + y_coord,
  data = dat_full_allpc
)

alias(fit_rank_check)

## Remove intercept before scaling, or keep as-is for raw design condition
kappa(X_design)

## Scale continuous columns for more interpretable condition number
X_scaled <- X_design
continuous_cols <- c("x_coord", "y_coord")

X_scaled[, continuous_cols] <- scale(X_scaled[, continuous_cols])

kappa(X_scaled)

# install.packages("car")
library(car)

fit_vif <- lm(
  rnorm(nrow(dat_full_allpc)) ~ X_cat + x_coord + y_coord,
  data = dat_full_allpc
)

vif(fit_vif)

## Residualize each PC on exposure + coordinates
resid_PC1 <- resid(lm(PC1_R ~ X_cat + x_coord + y_coord, data = dat_full_allpc))
resid_PC2 <- resid(lm(PC2_R ~ X_cat + x_coord + y_coord, data = dat_full_allpc))
resid_PC3 <- resid(lm(PC3 ~ X_cat + x_coord + y_coord, data = dat_full_allpc))

## Conditional/residual correlation
cor(cbind(resid_PC1, resid_PC2, resid_PC3))


## Residualize each PC on exposure + coordinates
resid_PC1 <- resid(lm(PC1_R ~ X_cat , data = dat_full_allpc))
resid_PC2 <- resid(lm(PC2_R ~ X_cat , data = dat_full_allpc))
resid_PC3 <- resid(lm(PC3 ~ X_cat , data = dat_full_allpc))
cor(cbind(resid_PC1, resid_PC2, resid_PC3))

##### Bootstrap #####
set.seed(123)

B <- 1000
n <- nrow(dat_full_allpc)

boot_summary <- vector("list", B)
boot_path <- vector("list", B)

for (b in 1:B) {
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  dat_b <- dat_full_allpc[idx, ]
  
  fit_b <- decompose_linear_multix(
    data = dat_b,
    exposures = c("X_adjacent", "X_far"),
    outcome = "Y",
    mediators = c("PC1_R", "PC2_R", "PC3"),
    covariates = c("x_coord", "y_coord")
  )
  
  boot_summary[[b]] <- fit_b$summary
  boot_path[[b]] <- fit_b$path
}
library(dplyr)
library(purrr)

boot_summary_df <- bind_rows(boot_summary, .id = "boot_id")
boot_path_df <- bind_rows(boot_path, .id = "boot_id")

summary_ci <- boot_summary_df %>%
  group_by(exposure) %>%
  summarise(
    TE_mean = mean(TE, na.rm = TRUE),
    TE_lwr = quantile(TE, 0.025, na.rm = TRUE),
    TE_upr = quantile(TE, 0.975, na.rm = TRUE),
    
    NDE_mean = mean(NDE, na.rm = TRUE),
    NDE_lwr = quantile(NDE, 0.025, na.rm = TRUE),
    NDE_upr = quantile(NDE, 0.975, na.rm = TRUE),
    
    NIE_mean = mean(NIE, na.rm = TRUE),
    NIE_lwr = quantile(NIE, 0.025, na.rm = TRUE),
    NIE_upr = quantile(NIE, 0.975, na.rm = TRUE),
    
    PM_mean = mean(PM, na.rm = TRUE),
    PM_lwr = quantile(PM, 0.025, na.rm = TRUE),
    PM_upr = quantile(PM, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

path_ci <- boot_path_df %>%
  group_by(exposure, mediator) %>%
  summarise(
    alpha_mean = mean(alpha_X_to_M, na.rm = TRUE),
    alpha_lwr = quantile(alpha_X_to_M, 0.025, na.rm = TRUE),
    alpha_upr = quantile(alpha_X_to_M, 0.975, na.rm = TRUE),
    
    beta_mean = mean(beta_M_to_Y, na.rm = TRUE),
    beta_lwr = quantile(beta_M_to_Y, 0.025, na.rm = TRUE),
    beta_upr = quantile(beta_M_to_Y, 0.975, na.rm = TRUE),
    
    indirect_mean = mean(indirect_component, na.rm = TRUE),
    indirect_lwr = quantile(indirect_component, 0.025, na.rm = TRUE),
    indirect_upr = quantile(indirect_component, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

summary_ci
path_ci


##### pc3 #####
# 1. far 내부 PC3 하위 10%
pc3_cutoff <- pca_df |>
  filter(X_cat == "far") |>
  summarise(cutoff = quantile(PC3, 0.10, na.rm = TRUE)) |>
  pull(cutoff)

# 2. 오른쪽 아래 국소 cluster
pc3_target_cluster <- pca_df |>
  filter(
    X_cat == "far",
    PC3 <= pc3_cutoff,
    x_coord >= 7600,
    y_coord >= 7800
  )

# 3. 선택된 위치 확인
ggplot(pca_df, aes(x = x_coord, y = y_coord)) +
  geom_point(color = "grey85", size = 2) +
  geom_point(
    data = pc3_target_cluster,
    color = "black",
    size = 2.7
  ) +
  scale_y_reverse() +
  coord_fixed() +
  theme_bw()

pc3_far_cutoff_10 <- pca_df |>
  filter(X_cat == "far") |>
  summarise(cutoff = quantile(PC3, 0.10, na.rm = TRUE)) |>
  pull(cutoff)

pc3_far_low <- pca_df |>
  filter(
    X_cat == "far",
    PC3 <= pc3_far_cutoff_10
  )
pc3_negative_genes <- c(
  "CXCL9", "CXCL10", "CCL5",
  "GZMK", "GNLY",
  "ISG15", "IFITM1",
  "TAP1", "NLRC5", "B2M",
  "HLA-DRA", "CD74"
)

pc3_target_cluster <- pc3_far_low |>
  filter(
    x_coord >= 7600,
    y_coord >= 7800
  )
ggplot(pca_df, aes(x = x_coord, y = y_coord)) +
  geom_point(color = "grey85", size = 2) +
  geom_point(
    data = pc3_target_cluster,
    color = "black",
    size = 2.5
  ) +
  scale_y_reverse() +
  coord_fixed() +
  theme_bw() +
  labs(
    title = "Far spots with low PC3 scores",
    x = "x coordinate",
    y = "y coordinate"
  )


M_expr
gene_df <- as.data.frame(M_expr) |>
  tibble::rownames_to_column("spot_id")
pca_df2 <- pca_df |>
  mutate(
    spot_id = rownames(pca_fit$x)
  ) |>
  left_join(
    gene_df,
    by = "spot_id"
  )
comparison_df <- pca_df2 |>
  mutate(
    pc3_cluster = if_else(
      X_cat == "far" &
        x_coord >= 7600 &
        x_coord <= 8500 &
        y_coord >= 7800 &
        y_coord <= 9300 &
        PC3 <= pc3_far_cutoff_10,
      "PC3-negative local cluster",
      "Other far spots"
    )
  ) |>
  filter(X_cat == "far")
comparison_df |>
  group_by(pc3_cluster) |>
  summarise(
    across(
      all_of(pc3_negative_genes),
      ~ mean(.x, na.rm = TRUE)
    )
  )
comparison_df |>
  group_by(pc3_cluster) |>
  summarise(
    n = n(),
    mean_immune_score = mean(immune_score, na.rm = TRUE),
    median_immune_score = median(immune_score, na.rm = TRUE),
    sd_immune_score = sd(immune_score, na.rm = TRUE)
  )
comparison_df |>
  count(pc3_cluster)
