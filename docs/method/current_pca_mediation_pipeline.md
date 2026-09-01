# Current PCA Mediation Pipeline: Statistical Analysis Specification

이 문서는 `scripts/current/run_pca_mediation.R`, `scripts/current/pcaaxis_0722.R`, `R/pca_mediation_pipeline.R`에 구현된 current frequentist PCA mediation pipeline을 기준으로 작성한 statistical analysis specification이다. `R/sensitivity_v2.R`은 main analysis에 source되지 않으므로 sensitivity analysis의 존재 여부와 pipeline boundary를 확인하는 참고 자료로만 사용하였다.

## 1. Analysis objective

현재 분석의 목표는 spatial transcriptomics spot 단위에서 tumor proximity가 gene-expression-derived PCA mediator를 통해 immune infiltration score와 어떻게 연관되는지를 linear mediation decomposition 형태로 요약하는 것이다.

- Exposure `X`: tumor proximity category, `X_cat`.
- Mediators `M`: `M_expr`에서 PCA로 구성한 `PC1_R`, `PC2_R`, `PC3`.
- Outcome `Y`: immune infiltration score. 코드에서는 `gex$Y`를 numeric으로 사용하며, 현재 RDS에서는 `meta$immune_score`와 동일하다.
- Spatial coordinates: `imagecol`, `imagerow`를 각각 `x_coord`, `y_coord`로 복사하여 main coordinate-adjusted mediation models의 covariates로 사용한다. 좌표는 PCA 자체에는 사용되지 않는다. Spatial PC plot에서는 `y_plot = -y_coord`로 뒤집어 시각화 좌표를 만든다.

현재 분석은 엄밀한 counterfactual natural direct/indirect effect를 별도로 추정하는 절차가 아니라, 동일한 spot-level data에서 선형회귀를 적합한 뒤 `alpha x beta` 곱과 그 합으로 indirect component를 계산하는 frequentist linear product-of-coefficients decomposition이다. 코드에서는 요약량 이름으로 `TE`, `NDE`, `NIE`, `PM`을 사용한다.

## 2. Analysis sample and variables

사용 데이터 파일은 `data/raw/pt16_emory_GEX_immune_FULL_v2.rds`이다. 이 RDS는 current run script에서 `rds_path`로 지정된다.

RDS object는 current code가 요구하는 `X_cat`, `Y`, `meta`, `M_expr`를 포함해야 한다. 현재 확인된 RDS에는 추가로 `X_design`, `params`, `M_expr_raw`, `resid_info`가 포함되어 있다.

분석에 포함되는 spot은 `build_full_analysis_data()`에서 다음 변수들을 선택한 뒤 `drop_na()`를 통과한 spot이다.

- `Y`
- `X_cat`
- `X_adjacent`
- `X_far`
- `PC1_R`
- `PC2_R`
- `PC3`
- `x_coord`
- `y_coord`

현재 RDS 기준으로 이 변수들에는 결측이 없으므로 1,507개 spot이 모두 포함된다. Category count는 `inside = 609`, `adjacent = 276`, `far = 622`이다. Current pipeline 자체는 추가 spot filtering을 다시 수행하지 않는다. RDS `params`에는 upstream construction parameter로 `MIN_SPOT_UMI = 500`이 기록되어 있으나, current main pipeline 코드에서 이 필터를 재적용하지 않는다.

Exposure 변수는 `X_cat`이며 factor level은 `inside`, `adjacent`, `far` 순서로 고정된다. Reference category는 `inside`이다. Mediation model에서는 다음 two dummy variables를 사용한다.

\[
A_i = I(X_{cat,i} = \text{adjacent}), \quad F_i = I(X_{cat,i} = \text{far})
\]

따라서 `X_adjacent`는 adjacent vs inside contrast, `X_far`는 far vs inside contrast를 나타낸다.

Outcome 변수는 `Y`이다. Spatial coordinate 변수는 `meta$imagecol`과 `meta$imagerow`이며 각각 `x_coord`, `y_coord`로 사용된다.

Gene-expression mediator source는 `gex$M_expr`이다. Current pipeline은 `M_expr`의 모든 column을 PCA에 사용하며, pipeline 안에서 유전자 선택을 새로 수행하지 않는다. 현재 RDS 기준 `M_expr`는 1,507 spot x 1,985 gene matrix이다.

RDS metadata에서 확인되는 gene preprocessing 정보는 다음과 같다.

- `params$M_GENE_COUNT = 2000`.
- `params$EXCLUDE_DEFINING_GENES = TRUE`.
- `params$DEFINING_GENES_REMOVED`에는 `MLANA`, `CD3D`, `CD3E`, `CD8A`, `CD4`, `TRAC`, `TRBC2`, `CD19`, `MS4A1`, `CD79A`, `CD79B`, `CD14`, `CD163`, `LYZ`, `CSF1R`가 기록되어 있다.
- `params$M_genes` length는 1,985이다.
- `resid_info$vars_regressed = "percent.mt"`.
- `resid_info$dropped_genes`는 empty character vector이다.
- `resid_info$note`에는 `M_expr`가 SCT log-normalized expression에서 technical covariate residualization을 거친 matrix이며, `log_total_umi`는 X-confounding 때문에 제외되었다고 기록되어 있다.

다만 current main pipeline 코드만으로는 upstream에서 2,000개 gene이 어떤 ranking rule로 선택되었는지 확인할 수 없다. Current analysis specification에서 확정할 수 있는 것은 PCA 입력이 이미 RDS에 저장된 `M_expr` 1,985개 gene이라는 점이다.

Missing value 처리는 mediation analysis dataset 생성 시 `drop_na()`로 수행된다. PCA 단계 전에는 별도 missing-value filtering 또는 imputation이 없다. 현재 RDS에서는 `X_cat`, `Y`, `meta`, `M_expr`에 결측이 없는 것으로 확인되었다.

## 3. PCA mediator construction

PCA 입력 행렬은 spot x gene matrix `gex$M_expr`이다.

먼저 gene별 standardization을 수행한다.

\[
Z_{ig} = \frac{M_{ig} - \bar{M}_{g}}{s_g}
\]

이는 R의 `scale(M_expr)`에 해당하므로 gene column별 centering과 scaling이 모두 수행된다. 이후 `prcomp(M_scaled, center = FALSE, scale. = FALSE)`를 실행한다. 즉 `prcomp()` 내부에서는 추가 centering 또는 추가 scaling을 하지 않는다.

PCA score는 `pca_fit$x`에서 추출된다. Current mediation pipeline은 첫 세 개 principal components만 mediator로 사용한다. Explained variance는 모든 PC에 대해 다음과 같이 계산된다.

\[
\text{PVE}_k = \frac{sdev_k^2}{\sum_{\ell} sdev_\ell^2}
\]

Loading은 `pca_fit$rotation[, 1:3]`에서 첫 세 PC에 대해 추출한다.

PCA axis orientation은 코드에서 명시적으로 다음과 같이 바뀐다.

\[
PC1\_R_i = -PC1\_{original,i}
\]

\[
PC2\_R_i = -PC2\_{original,i}
\]

\[
PC3_i = PC3\_{original,i}
\]

Loading도 동일한 부호 convention으로 변환한다.

\[
loading(PC1\_R)_g = -loading(PC1\_{original})_g
\]

\[
loading(PC2\_R)_g = -loading(PC2\_{original})_g
\]

\[
loading(PC3)_g = loading(PC3\_{original})_g
\]

`pcaaxis_0722.R`의 주석에 따르면 final orientation label은 다음과 같다.

- `PC1_R`: stromal / ECM / humoral immune direction이 positive.
- `PC2_R`: antigen-presentation / adaptive immune direction이 positive.
- `PC3`: mitochondrial / metabolic / plasma-cell direction이 positive.

`PC1_R`, `PC2_R`의 `_R`은 original PCA score와 loading의 sign reversal을 의미한다. PCA score를 만든 뒤 mediation analysis 전에 PC score를 다시 z-score로 표준화하는 단계는 없다. 따라서 mediator regression과 outcome regression에서 사용하는 PC score scale은 `prcomp()` score scale이며, 각 PC score의 variance는 해당 PC eigenvalue에 대응한다.

Sign reversal은 downstream mediation analysis에 직접 반영된다. 즉 alpha coefficient는 reversed score 기준의 `X -> PC_R` association이고, beta coefficient도 reversed score 기준의 `PC_R -> Y` partial association이다. Original PC axis와 비교하면 `PC1_R`, `PC2_R`에서는 alpha와 beta가 각각 부호가 반전되지만, 같은 mediator에 대한 product `alpha x beta`는 수학적으로 sign reversal에 대해 invariant하다. 다만 mediator direction의 biological interpretation은 reversed axis label을 기준으로 해야 한다.

Conditional PCA, Bayesian PCA, exposure-adjusted PCA, coordinate-adjusted PCA는 current main pipeline에 구현되어 있지 않다.

## 4. Mediator models: X to M

Main coordinate-adjusted analysis에서 각 mediator \(j \in \{1,2,3\}\)에 대해 별도의 linear regression을 적합한다. 여기서 \(M_{i1}=PC1\_R_i\), \(M_{i2}=PC2\_R_i\), \(M_{i3}=PC3_i\)이다.

\[
M_{ij} =
\alpha_{0j}
+ \alpha_{Aj} A_i
+ \alpha_{Fj} F_i
+ \delta_{xj} x_i
+ \delta_{yj} y_i
+ \varepsilon^M_{ij}
\]

여기서 \(x_i = x\_coord_i\), \(y_i = y\_coord_i\)이다.

Exposure coding은 `X_adjacent`와 `X_far` two dummy variables이며, `inside` spot은 \(A_i=0, F_i=0\)인 reference group이다. Current main model의 covariates는 `x_coord + y_coord`뿐이다. 추가 biological 또는 technical covariate는 current mediation regression에 포함되지 않는다.

Mediator model은 mediator별 separate regression이다. `PC1_R`, `PC2_R`, `PC3`를 multivariate response로 함께 적합하는 joint mediator model은 아니다.

\(\alpha_{Aj}\)는 coordinate-adjusted adjacent vs inside contrast가 mediator \(j\) score와 갖는 association이다. 즉 같은 coordinates를 조건으로, adjacent spot의 mediator score가 inside spot 대비 얼마나 다른지를 나타내는 linear model coefficient이다.

\(\alpha_{Fj}\)는 coordinate-adjusted far vs inside contrast가 mediator \(j\) score와 갖는 association이다.

코드의 `alpha_X_to_M`는 각 mediator model에서 `X_adjacent` 또는 `X_far` coefficient를 추출한 값이며, 각각 \(\alpha_{Aj}\), \(\alpha_{Fj}\)에 해당한다.

Current pipeline은 no-control version도 보조 결과로 계산한다. 이 경우 식은 spatial covariate 없이 \(M_{ij} = \alpha_{0j} + \alpha_{Aj} A_i + \alpha_{Fj} F_i + \varepsilon^M_{ij}\)이다. Main analysis specification은 coordinate-adjusted version이다.

## 5. Outcome model: M to Y

Main coordinate-adjusted outcome regression은 두 exposure dummy, 세 PC mediators, spatial coordinates를 동시에 포함하는 single linear model이다.

\[
Y_i =
\theta_0
+ c'_A A_i
+ c'_F F_i
+ \sum_{j=1}^{3} \beta_j M_{ij}
+ \eta_x x_i
+ \eta_y y_i
+ \varepsilon^Y_i
\]

Exposure는 outcome model에 직접 포함된다. 따라서 \(c'_A\)와 \(c'_F\)는 mediators와 coordinates를 조건으로 한 adjacent vs inside 및 far vs inside direct exposure coefficient이다.

세 PC mediators는 동시에 outcome model에 들어간다. 따라서 \(\beta_j\)는 다른 PC mediators, exposure dummies, spatial coordinates를 모두 조건으로 했을 때 mediator \(j\) score 1-unit 증가와 `Y`의 association이다.

Current multiple mediator model은 outcome side에서 `PC1_R + PC2_R + PC3`를 같은 regression에 동시에 넣고, mediator side에서는 각 PC에 대해 별도 regression을 적합한 뒤 같은 outcome-model beta를 mediator-specific product에 사용하는 방식으로 구현된다.

## 6. Mediator-specific indirect components

각 exposure contrast와 mediator \(j\)에 대한 mediator-specific indirect component는 `alpha_X_to_M * beta_M_to_Y`로 계산된다.

Adjacent vs inside contrast:

\[
IC_{Aj} = \alpha_{Aj} \beta_j
\]

Far vs inside contrast:

\[
IC_{Fj} = \alpha_{Fj} \beta_j
\]

여기서 \(\alpha_{Aj}\), \(\alpha_{Fj}\)는 section 4의 mediator-specific exposure coefficients이고, \(\beta_j\)는 section 5의 common multiple-mediator outcome model에서 추출한 mediator coefficient이다.

`PC1_R`와 `PC2_R`는 original PC axis에서 sign-reversed된 mediator이다. Current code는 reversed score를 mediator model과 outcome model 양쪽에 모두 사용하므로, mediator-specific indirect product는 original axis로 표현해도 동일한 numerical value를 갖는다. 그러나 \(\alpha\)와 \(\beta\) 각각의 부호 해석은 `PC1_R`, `PC2_R`의 positive biological direction을 기준으로 해야 한다.

## 7. Total mediation decomposition

Current code는 exposure contrast별로 다음 quantities를 계산한다.

Total effect model:

\[
Y_i =
\tau_0
+ \tau_A A_i
+ \tau_F F_i
+ \kappa_x x_i
+ \kappa_y y_i
+ \varepsilon^{TE}_i
\]

`TE`는 이 model에서 추출한 exposure coefficient이다.

\[
TE_A = \tau_A, \quad TE_F = \tau_F
\]

`NDE` 또는 reported direct effect는 total effect와 total indirect effect의 차이로 정의된다. Section 5의 outcome model에서 mediator들을 함께 조정한 exposure coefficient는 `outcome_direct_coef`로 별도로 보존되지만, 보고되는 `NDE`에는 사용하지 않는다.

\[
NDE_A = TE_A - NIE_A, \quad NDE_F = TE_F - NIE_F
\]

`NIE` 또는 total indirect effect는 mediator-specific indirect components의 합으로 계산된다.

\[
NIE_A = \sum_{j=1}^{3} \alpha_{Aj}\beta_j
\]

\[
NIE_F = \sum_{j=1}^{3} \alpha_{Fj}\beta_j
\]

`PM`은 proportion mediated로 계산된다.

\[
PM_A = \frac{NIE_A}{TE_A}, \quad PM_F = \frac{NIE_F}{TE_F}
\]

분모는 각 exposure contrast의 `TE`이다. 코드에는 `TE = 0` 또는 sign-changing decomposition에 대한 별도 예외 처리나 truncation이 없다.

`TE`, `NDE`, `NIE`, `PM`은 `X_adjacent`와 `X_far`에 대해 각각 별도로 계산된다. 따라서 `adjacent vs inside`와 `far vs inside`는 같은 model 안에서 동시에 parameterized된 두 contrast이다.

코드는 `TE`를 total effect model에서 따로 추정하고, `NIE`를 \(\sum_j \alpha_j \beta_j\)로 계산한 뒤 `NDE = TE - NIE`로 계산한다. 따라서 현재 구현에서는 `TE = NDE + NIE`가 정의상 성립한다. Outcome model의 exposure coefficient는 `outcome_direct_coef`라는 내부 진단량으로 별도 보존된다. 동일 sample, 동일 covariates, additive linear models, no exposure-mediator interaction의 product-of-coefficients decomposition에서는 이 관계가 선형회귀 분해로 해석된다. 이 문서의 `NDE`와 `NIE` 용어는 코드의 naming을 유지한 것이며, 별도 counterfactual natural effect estimator가 구현되어 있다는 뜻은 아니다.

## 8. Bootstrap uncertainty

Bootstrap은 `dat_full_allpc`의 row, 즉 spot을 resampling unit으로 사용한다. Current run script에서는 `bootstrap_B = 1000`, `bootstrap_seed = 123`이다.

각 replicate \(b\)에서 다음을 수행한다.

1. `set.seed(123)` 이후, 전체 spot 수 \(n\)과 같은 크기로 row index를 `sample(seq_len(n), size = n, replace = TRUE)`로 추출한다.
2. Resampled dataset `dat_b`를 만든다.
3. `dat_b`에서 coordinate-adjusted `decompose_linear_multix()`를 다시 실행한다.
4. Total effect model, mediator models, outcome model을 모두 다시 적합한다.
5. Replicate별 `summary`와 `path`를 저장한다.

PCA 자체는 bootstrap replicate 안에서 다시 적합하지 않는다. Bootstrap은 original full data에서 한 번 계산된 `PC1_R`, `PC2_R`, `PC3` score를 fixed variables로 두고 row resampling만 수행한다. 따라서 PCA loading, PCA explained variance, gene scaling, PC projection, upstream gene selection과 preprocessing의 uncertainty는 bootstrap에 포함되지 않는다.

Confidence interval은 bootstrap replicate distribution의 empirical percentile interval이다. 각 exposure별 `TE`, `NDE`, `NIE`, `PM`에 대해 2.5%와 97.5% quantile을 계산한다. 각 exposure-mediator pair별 `alpha_X_to_M`, `beta_M_to_Y`, `indirect_component`에 대해서도 2.5%와 97.5% quantile을 계산한다. 코드에서는 bootstrap mean도 함께 계산하지만, percentile interval 외의 BCa, normal approximation, studentized interval은 구현되어 있지 않다.

현재 bootstrap이 포함하는 uncertainty는 fixed PC scores와 fixed selected analysis dataset을 조건으로 한 spot-level row resampling variability이다. 현재 bootstrap이 포함하지 않는 uncertainty는 PCA fitting/projection uncertainty, gene selection uncertainty, upstream residualization uncertainty, spatially correlated residual 또는 neighborhood dependence를 명시적으로 반영하는 block/bootstrap design uncertainty이다.

## 9. Outputs

Current pipeline의 주요 결과는 file로 자동 저장되지 않고, `run_pca_mediation_pipeline()`의 returned list와 interactive plot objects로 생성된다.

생성되는 주요 analysis outputs는 다음과 같다.

- PCA-scaled expression matrix `M_scaled`.
- PCA fit object, PC scores, PC loadings.
- Explained variance와 scree data.
- Final oriented `PC1_R`, `PC2_R`, `PC3` scores.
- Final oriented loading table과 top loading genes.
- Spatial PC score figures for `PC1_R`, `PC2_R`, `PC3` when `make_plots = TRUE`.
- Coordinate-adjusted `X -> M` estimates.
- Coordinate-adjusted `M -> Y` estimates from the multiple mediator outcome model.
- Mediator-specific indirect components for `adjacent vs inside` and `far vs inside`.
- `TE`, `NDE`, `NIE`, `PM` summaries for both exposure contrasts.
- No-control decomposition results as auxiliary comparison.
- Bootstrap percentile confidence intervals for summary decomposition quantities and path-level quantities.
- Design diagnostics, including design rank, condition number, VIF, and residual PC correlations.

Pipeline에는 PC3 local cluster exploratory summary도 포함되어 있으나, 이는 main mediation decomposition의 구성 요소가 아니라 PC3 해석을 보조하는 diagnostic/exploratory output이다.

## 10. Current statistical assumptions

### Modeling assumptions

Current regression decomposition은 다음 modeling assumptions에 의존한다.

- Relationships among exposure dummies, PC mediators, coordinates, and outcome are adequately represented by additive linear models.
- Exposure-mediator interaction terms are absent from the implemented outcome model.
- Spatial coordinate adjustment is linear in `x_coord` and `y_coord`.
- PC mediators are treated as fixed observed variables after full-sample PCA construction.
- The same complete-case analysis dataset is used across total effect, mediator, outcome, and bootstrap models.
- OLS coefficients are interpreted as conditional linear associations. Bootstrap CIs use row resampling rather than analytic normal-theory standard errors.

### Causal identification assumptions

Interpreting the product-of-coefficients decomposition causally as mediation would additionally require assumptions not verified by the code, including:

- No relevant unmeasured confounding of tumor proximity and immune score conditional on included covariates.
- No relevant unmeasured confounding of tumor proximity and PCA mediators conditional on included covariates.
- No relevant unmeasured confounding of PCA mediators and immune score conditional on exposure and included covariates.
- No mediator-outcome confounder affected by exposure, unless otherwise handled outside the current code.
- The spatial coordinates `x_coord` and `y_coord` are sufficient for the intended spatial confounding adjustment, if causal interpretation is claimed.

The current code does not test these causal assumptions. It implements the specified linear regression decomposition.

## 11. Current limitations

The following limitations are directly visible from the current implementation.

- PCA is unconditional with respect to exposure and coordinates. `M_expr` is scaled gene-wise and passed to `prcomp()` without `X_cat`, `x_coord`, or `y_coord`.
- \(K = 3\) mediators is analytically fixed in the current pipeline. The mediation models hard-code `PC1_R`, `PC2_R`, `PC3`.
- PCA projection uncertainty is not included in bootstrap because bootstrap resamples fixed PC score rows and does not refit PCA.
- Mediator selection uncertainty is not included. The selected mediator set `PC1_R`, `PC2_R`, `PC3` is fixed before bootstrapping.
- Spatial dependence is handled only through linear adjustment for `x_coord` and `y_coord` in the main models. The bootstrap is ordinary row bootstrap, not spatial block bootstrap.
- Main pipeline does not include mediator-outcome unmeasured-confounding sensitivity analysis. `R/sensitivity_v2.R` contains sensitivity helper functions, but it is not sourced or run by `scripts/current/run_pca_mediation.R`.
- Current pipeline does not include exposure-mediator interaction terms in the outcome model.
- Upstream selection of the original 2,000 genes is not implemented in the current pipeline files. The current code uses the RDS-provided `M_expr`; exact upstream ranking rule is 현재 코드 기준 확인 필요.

## 12. Boundary of the current pipeline

현재 main pipeline에 구현된 것은 frequentist PCA mediation pipeline이다. It consists of full-sample PCA mediator construction, sign-oriented PC1 to PC3 scores, linear product-of-coefficients decomposition for `adjacent vs inside` and `far vs inside`, coordinate-adjusted main models, no-control auxiliary comparison, diagnostics, plots, and ordinary row bootstrap percentile intervals conditional on fixed PC scores.

아직 current main pipeline에 포함되지 않은 것은 다음이다.

- Bayesian mediation.
- Updated sensitivity analysis as part of the main run.
- Alternative K analysis.
- Full PCA uncertainty propagation.
- Spatial block bootstrap.
- Conditional PCA.
- Shiny interface.

These items are outside the current main analysis boundary and are not specified here as implemented methods.
