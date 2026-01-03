## dominance_and_residuals.R
## Reproduces reviewer-requested dominance analysis and residual diagnostics.
##
## Linear model:
##   ResponseTime ~ HRPE_Density + PromLength_Scaled
##
## Outputs (saved in the current working directory):
##   - FigS_dominance_barplot.pdf / .png (600 dpi)
##   - FigS_residual_diagnostics.pdf / .png (600 dpi)   [NO point numbers]
##   - FigS_residuals_vs_fitted.pdf / .png (600 dpi)    [gene names shown]
##   - dominance_results.csv
##   - residual_influence_table.csv
##
## Requirements:
##   install.packages(c("tidyverse","dominanceanalysis"))
##
## Run:
##   source("code/analysis/dominance_and_residuals.R")

## ============================================================
## Dominance analysis + publication-quality figures (PDF + PNG)
## + residual diagnostics (reviewer requested)
##
## Linear model used for dominance analysis:
##   ResponseTime ~ HRPE_Density + PromLength_Scaled
##
## Outputs:
##   - FigS_dominance_barplot.pdf / .png (600 dpi)
##   - FigS_residual_diagnostics.pdf / .png (600 dpi)   [NO point numbers]
##   - FigS_residuals_vs_fitted.pdf / .png (600 dpi)    [gene names shown]
##   - dominance_results.csv
##   - residual_influence_table.csv
## ============================================================

library(tidyverse)
library(dominanceanalysis)

## ------------------------------------------------------------
## 1) Data (exactly as MATLAB)
## ------------------------------------------------------------
geneNames <- c("ADH1","HB1","HRE2","HUP7","PCO1","HRA1","LBD41","PDC1","SUS4")

ResponseTime <- c(68.54,94.79,57.79,83.81,47.62,13.88,14.57,73.83,64.31)

geneNames_hrpe <- c("SUS4","HB1","PCO1","HUP7","HRE2","HRA1","LBD41","PDC1","ADH1")
hrpe_number    <- c(3,1,3,5,3,8,10,4,3)
prom_length_bp <- c(2348,518,1212,1614,2402,1052,2262,2632,867)

idx <- match(geneNames, geneNames_hrpe)
if (any(is.na(idx))) stop("Gene name mismatch detected.")

HRPE_Density <- hrpe_number[idx] / prom_length_bp[idx]
PromLength_Scaled <- as.numeric(scale(prom_length_bp[idx]))  # z-score

dat <- tibble(
  Gene = geneNames,
  ResponseTime = ResponseTime,
  HRPE_Density = HRPE_Density,
  PromLength_Scaled = PromLength_Scaled
)

cat("\n=== Analysis dataset ===\n")
print(dat)

## ------------------------------------------------------------
## 2) Linear model (dominance model)
## ------------------------------------------------------------
lm_fit <- lm(ResponseTime ~ HRPE_Density + PromLength_Scaled, data = dat)

cat("\n=== Linear model summary ===\n")
print(summary(lm_fit))
cat("\nLM R²:", summary(lm_fit)$r.squared, "\n")

## ------------------------------------------------------------
## 3) Dominance analysis
## ------------------------------------------------------------
da <- dominanceAnalysis(lm_fit)
ac <- averageContribution(da)[["r2"]]

dom_df <- tibble(
  Predictor = names(ac),
  GeneralDominance_R2 = as.numeric(ac),
  RelativeImportance_pct = 100 * GeneralDominance_R2 / sum(GeneralDominance_R2)
) %>%
  mutate(
    Predictor_label = case_when(
      Predictor == "HRPE_Density" ~ "HRPE density",
      Predictor == "PromLength_Scaled" ~ "Promoter length",
      TRUE ~ Predictor
    )
  )

cat("\n=== Dominance results ===\n")
print(dom_df)

## ------------------------------------------------------------
## 4) Dominance bar plot (no grid; publication style)
## ------------------------------------------------------------
p_dom <- ggplot(dom_df,
                aes(x = reorder(Predictor_label, -RelativeImportance_pct),
                    y = RelativeImportance_pct)) +
  geom_col(width = 0.60, colour = "black", linewidth = 0.40) +
  geom_text(aes(label = sprintf("%.1f%%", RelativeImportance_pct)),
            vjust = -0.50, size = 4.0) +
  labs(
    x = NULL,
    y = "Relative importance (% of explained variance)",
    title = "Dominance analysis of response time"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5),
    axis.line = element_blank(),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text.x  = element_text(face = "plain"),
    plot.title   = element_text(face = "plain", hjust = 0.5),
    plot.margin  = margin(t = 8, r = 8, b = 8, l = 20)
  ) +
  ylim(0, max(dom_df$RelativeImportance_pct) * 1.25)

print(p_dom)

ggsave("FigS_dominance_barplot.pdf", p_dom, width = 7.6, height = 4.4)
ggsave("FigS_dominance_barplot.png", p_dom, width = 7.6, height = 4.4, dpi = 600)

## ------------------------------------------------------------
## 5) Residual diagnostics (4-panel) — NO numbers on points
## ------------------------------------------------------------
make_resid_diag <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(
    mfrow = c(2, 2),
    mar = c(5.0, 6.0, 3.2, 2.0),
    cex.lab = 1.1,
    cex.axis = 1.0,
    cex.main = 1.0,
    xpd = FALSE
  )
  
  plot(lm_fit, which = 1, id.n = 0)  # Residuals vs Fitted
  plot(lm_fit, which = 2, id.n = 0)  # Normal Q-Q
  plot(lm_fit, which = 3, id.n = 0)  # Scale-Location
  plot(lm_fit, which = 5, id.n = 0)  # Residuals vs Leverage
}

pdf("FigS_residual_diagnostics.pdf", width = 8.5, height = 8.5)
make_resid_diag()
dev.off()

png("FigS_residual_diagnostics.png",
    width = 8.5, height = 8.5, units = "in", res = 600)
make_resid_diag()
dev.off()

## ------------------------------------------------------------
## 6) Residuals vs fitted values (gene names shown)
## ------------------------------------------------------------
res <- resid(lm_fit)
fit <- fitted(lm_fit)

y_pad <- 0.12 * diff(range(res))
ylim_use <- c(min(res) - y_pad, max(res) + y_pad)

y_threshold <- quantile(res, 0.85)
pos_vec <- ifelse(res >= y_threshold, 1, 3)  # 1=below, 3=above

make_resid_label_plot <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(
    mar = c(5.0, 6.0, 3.2, 2.0),
    cex.lab = 1.1,
    cex.axis = 1.0,
    cex.main = 1.0,
    xpd = FALSE
  )
  
  plot(fit, res,
       pch = 19,
       xlab = "Fitted values",
       ylab = "Residuals",
       main = "Residuals vs fitted values",
       ylim = ylim_use)
  
  abline(h = 0, lty = 2)
  
  par(xpd = NA)
  text(fit, res, labels = dat$Gene, pos = pos_vec, cex = 0.9, offset = 0.4)
  par(xpd = FALSE)
}

pdf("FigS_residuals_vs_fitted.pdf", width = 7.6, height = 5.4)
make_resid_label_plot()
dev.off()

png("FigS_residuals_vs_fitted.png",
    width = 7.6, height = 5.4, units = "in", res = 600)
make_resid_label_plot()
dev.off()

## ------------------------------------------------------------
## 7) Influence table + export
## ------------------------------------------------------------
infl_df <- tibble(
  Gene = dat$Gene,
  Fitted = fit,
  Residual = res,
  StudentizedResidual = rstudent(lm_fit),
  Leverage = hatvalues(lm_fit),
  CooksD = cooks.distance(lm_fit)
) %>%
  arrange(desc(CooksD))

cat("\n=== Influence table ===\n")
print(infl_df)

cat("\nMax Cook's distance:",
    round(max(infl_df$CooksD), 3),
    "at gene",
    infl_df$Gene[which.max(infl_df$CooksD)], "\n")

write.csv(dom_df,  "dominance_results.csv",        row.names = FALSE)
write.csv(infl_df, "residual_influence_table.csv", row.names = FALSE)

cat("\nSaved files (current working directory):\n",
    "- FigS_dominance_barplot.pdf / .png\n",
    "- FigS_residual_diagnostics.pdf / .png\n",
    "- FigS_residuals_vs_fitted.pdf / .png\n",
    "- dominance_results.csv\n",
    "- residual_influence_table.csv\n",
    "\nTip: run getwd() in R to see the folder path.\n")
