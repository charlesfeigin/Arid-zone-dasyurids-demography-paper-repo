library(geodata)
library(usdm)
library(ENMeval)
library(maxnet)
library(pastclim)

# Set up pastclim and download the Beyer et al. (2020) climate and vegetation data
dir.create("pastclim_data", showWarnings = FALSE)
set_data_path(path_to_nc = "pastclim_data", ask = FALSE)
env_vars <- c("bio01", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11",
              "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19",
              "npp", "lai")
download_dataset(dataset = "Beyer2020", bio_variables = env_vars)

aus <- ext(112, 155, -45, -10)
get_clim <- function(time_bp) {
  x <- region_slice(time_bp = time_bp, bio_variables = env_vars, dataset = "Beyer2020")
  crop(x, aus)
}

# Load present-day climate and vegetation
bio <- get_clim(0)

# Load and clean occurrence records
d <- read.csv("nridei_records-2026-09-29.csv")
d <- d[!is.na(d$decimalLatitude) & !is.na(d$decimalLongitude), ]
d <- d[!duplicated(d[, c("decimalLatitude", "decimalLongitude")]), ]
d <- d[is.na(d$coordinateUncertaintyInMeters) | d$coordinateUncertaintyInMeters <= 5000, ]


# Load the IUCN range map shapefile for N. ridei
# UPDATE THIS PATH to point to your downloaded IUCN shapefile
range_map <- vect("redlist_species_data_24f88d4b-270b-4071-afbd-05f70a0b59df/data_0.shp")

# IUCN range map with a buffer (terra buffer widths are in metres)
range_map_buffered <- buffer(range_map, width = 200 * 1000)

# Convert occurrence data frame to a spatial points object
occ_pts <- vect(d, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

# Extract polygon attributes for each point
# Points outside the buffered polygon will have NA for all polygon attributes
ext <- extract(range_map_buffered, occ_pts)

# Filter the original data frame to keep only points inside the buffered range map
# We check if at least one polygon attribute column is not NA
inside_range <- rowSums(!is.na(ext[, -1])) > 0
d <- d[inside_range, ]


# Thin to one record per climate grid cell
cell <- cellFromXY(bio, as.matrix(d[, c("decimalLongitude", "decimalLatitude")]))
d <- d[!is.na(cell) & !duplicated(cell), ]
occ <- d[, c("decimalLongitude", "decimalLatitude")]

# Remove variables correlated above 0.7
vals <- na.omit(extract(bio, occ, ID = FALSE))
vc <- vifcor(vals, th = 0.7)
vars <- setdiff(names(bio), vc@excluded)
bio <- bio[[vars]]



# Background: all present-day land cells in the study area (continental Australia)
bg <- as.data.frame(bio, xy = TRUE, na.rm = TRUE)[, c("x", "y")]
colnames(bg) <- c("decimalLongitude", "decimalLatitude")

# Fit MaxEnt models with spatial block cross-validation
e <- ENMevaluate(occs = occ, envs = bio, bg = bg,
                 algorithm = "maxnet", partitions = "block",
                 tune.args = list(fc = c("L", "LQ", "LQH"), rm = 1:4))

# Select the best model (lowest AICc)
results <- eval.results(e)
best <- results[which(results$delta.AICc == 0), ]
mod <- eval.models(e)[[as.character(best$tune.args)]]

# Permutation importance: drop in AUC when each variable is shuffled (10 shuffles each)
p_env <- extract(bio, occ, ID = FALSE)
b_env <- extract(bio, bg, ID = FALSE)
env_all <- rbind(p_env, b_env)
pa <- c(rep(1, nrow(p_env)), rep(0, nrow(b_env)))
keep <- complete.cases(env_all)
env_all <- env_all[keep, ]
pa <- pa[keep]

auc_fun <- function(pred, pa) {
  r <- rank(pred)
  n1 <- sum(pa == 1)
  n0 <- sum(pa == 0)
  (sum(r[pa == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

base_auc <- auc_fun(as.numeric(predict(mod, env_all, type = "cloglog")), pa)
set.seed(1)
perm_imp <- sapply(vars, function(v) mean(replicate(10, {
  shuffled <- env_all
  shuffled[[v]] <- sample(shuffled[[v]])
  base_auc - auc_fun(as.numeric(predict(mod, shuffled, type = "cloglog")), pa)
})))

var_importance <- data.frame(variable = vars, auc_drop = perm_imp,
                             percent = 100 * pmax(perm_imp, 0) / sum(pmax(perm_imp, 0)))
var_importance[order(-var_importance$auc_drop), ]




# Present-day suitability
suit_cur <- predict(bio, mod, type = "cloglog", na.rm = TRUE)

# MSMC2 sample location
msmc_sample <- data.frame(decimalLongitude = 122, decimalLatitude = -24.45)

# Check results
nrow(occ)
vars
best[, c("tune.args", "auc.val.avg", "AICc")]

# Time points for the SDM maps, every 20 ky from 0 to 100 kya
times <- seq(0, 100000, by = 20000)

# MSMC2 Ne vs mean climatic suitability at the MSMC2 sample site, per MSMC2 time step (midpoints up to 100 kya)
mu <- 5.95e-9
gen <- 1
msmc <- read.table("whole_genome_msmc2_wongai_i50_p1-2_25-1_1-2_1-3.final.txt", header = TRUE)
msmc$start_years <- msmc$left_time_boundary / mu * gen
msmc$end_years <- msmc$right_time_boundary / mu * gen
msmc$mid_years <- (msmc$start_years + msmc$end_years) / 2
msmc$Ne <- (1 / msmc$lambda) / (2 * mu)
msmc_s <- msmc[msmc$mid_years <= 100000, ]

# Suitability at the sample site for every climate time step these MSMC2 steps span
clim_times_s <- get_time_bp_steps(dataset = "Beyer2020")
clim_times_s <- clim_times_s[-clim_times_s <= max(msmc_s$end_years)]
clim_suit_s <- sapply(clim_times_s, function(t)
  as.numeric(predict(mod, extract(get_clim(t)[[vars]], msmc_sample, ID = FALSE), type = "cloglog")))

# Mean suitability over each MSMC2 step (nearest climate time step if none fall inside)
ne_vs_suit <- data.frame(
  kya = msmc_s$mid_years / 1000,
  Ne = msmc_s$Ne,
  suitability = sapply(seq_len(nrow(msmc_s)), function(i) {
    inside <- -clim_times_s >= msmc_s$start_years[i] & -clim_times_s < msmc_s$end_years[i]
    if (any(inside)) mean(clim_suit_s[inside]) else clim_suit_s[which.min(abs(-clim_times_s - msmc_s$mid_years[i]))]
  })
)
ne_vs_suit

# Spearman rank correlation of Ne and mean suitability across MSMC2 steps
cor.test(ne_vs_suit$Ne, ne_vs_suit$suitability, method = "spearman")





# Write each plot as both PNG and SVG
# Sizes are in inches; PNGs are written at 250 pixels per inch
save_png_svg <- function(name, width_in, height_in, plot_fun) {
  png(paste0(name, ".png"), width = width_in * 250, height = height_in * 250, res = 250)
  plot_fun()
  dev.off()
  svg(paste0(name, ".svg"), width = width_in, height = height_in)
  plot_fun()
  dev.off()
}

# 1) Present-day SDM with occurrence records (black) and the MSMC2 sample (red)
save_png_svg("nridei_ausbg_sdm_present_records", 8, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(3, 3, 2, 5))
  plot(suit_cur, main = expression("Present-day SDM suitability (" * italic("N. ridei") * ")"), range = c(0, 1))
  points(occ, pch = 20, cex = 1.0, col = "black")
  points(msmc_sample, pch = 20, col = "red", cex = 1.5)
})

# 2) SDM at each time point, 100 to 0 kya in 20 ky steps, as 2 x 3 panels, with the MSMC2 sample (red)
# Each map is masked to present-day land so all panels show modern Australia's outline
suit_times <- lapply(rev(times), function(t) mask(predict(get_clim(-t)[[vars]], mod, type = "cloglog", na.rm = TRUE), suit_cur))
save_png_svg("nridei_ausbg_sdm_timepoints", 14.4, 8.8, function() {
  par(mfrow = c(2, 3), mar = c(3, 3, 2, 5))
  for (i in seq_along(suit_times)) {
    t <- rev(times)[i]
    label <- if (t == 0) "Present" else paste0(t / 1000, " ka")
    plot(suit_times[[i]], main = label, range = c(0, 1))
    points(msmc_sample, pch = 20, col = "red", cex = 1.5)
  }
})

# 3) Ne vs mean suitability at the MSMC2 sample site per MSMC2 step
save_png_svg("nridei_ausbg_ne_vs_suitability", 9.6, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 2, 5), cex.lab = 1.25, cex.axis = 1.2)
  plot(ne_vs_suit$kya, ne_vs_suit$Ne, type = "n",
       xlim = c(0, 100), ylim = c(0, max(ne_vs_suit$Ne) * 1.1),
       xlab = "Thousand years before present", ylab = "Effective population size (Ne)",
       main = expression(bold("Ne and SDM suitability per MSMC2 time step (" * bolditalic("N. ridei") * ")")))
  lines(ne_vs_suit$kya, ne_vs_suit$Ne, type = "o", pch = 19, col = "darkgreen")
  par(new = TRUE)
  plot(ne_vs_suit$kya, ne_vs_suit$suitability, type = "o", pch = 17, col = "darkorange",
       xlim = c(0, 100), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext("Mean climatic suitability at MSMC2 sample site", side = 4, line = 3, cex = 1.25)
  legend("topright", legend = c(expression(italic("N. ridei") ~ "Ne"), "SDM suitability"),
         pch = c(19, 17), lty = c(1, 1),
         col = c("darkgreen", "darkorange"), pt.cex = c(1, 1), cex = 0.95)
})

# 4) MSMC2 Ne through time with the 100 bootstrap replicates, for N. ridei and S. murina
# type = "s" draws each Ne value from its left time boundary to the next boundary
read_msmc <- function(f) {
  x <- read.table(f, header = TRUE)
  data.frame(years = x$left_time_boundary / mu * gen, Ne = (1 / x$lambda) / (2 * mu))
}
nr_main <- read_msmc("whole_genome_msmc2_wongai_i50_p1-2_25-1_1-2_1-3.final.txt")
nr_boot <- lapply(1:100, function(i)
  read_msmc(paste0("wongai/wongai_i50_p1-2_25-1_1-2_1-3.bootstrap_", i, ".final.txt")))
sm_main <- read_msmc("whole_genome_msmc2_slendertailed_i50_p1-2_25-1_1-2_1-3.final.txt")
sm_boot <- lapply(1:100, function(i)
  read_msmc(paste0("slendertailed/slendertailed_i50_p1-2_25-1_1-2_1-3.bootstrap_", i, ".final.txt")))

x_ticks <- seq(0, 240000, by = 20000)
y_ticks <- seq(0, 1e6, by = 2e5)

save_png_svg("nridei_ausbg_smurina_msmc2_bootstraps", 9.6, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(5, 6, 2, 2))
  plot(nr_main$years, nr_main$Ne, type = "n",
       xlim = c(0, 250000), ylim = c(0, 1e6), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
       xlab = "Years before present", ylab = "",
       main = "Demographic history inferred with MSMC2")
  axis(1, at = x_ticks, labels = formatC(x_ticks, format = "d", big.mark = ","), cex.axis = 0.75, gap.axis = -1)
  axis(2, at = y_ticks, labels = formatC(y_ticks, format = "d", big.mark = ","), las = 1)
  mtext("Effective population size (Ne)", side = 2, line = 4.5)
  for (b in nr_boot) lines(b$years, b$Ne, type = "s", col = adjustcolor("darkgreen", alpha.f = 0.15), lwd = 0.5)
  for (b in sm_boot) lines(b$years, b$Ne, type = "s", col = adjustcolor("darkviolet", alpha.f = 0.15), lwd = 0.5)
  lines(nr_main$years, nr_main$Ne, type = "s", col = "darkgreen", lwd = 1.75)
  lines(sm_main$years, sm_main$Ne, type = "s", col = "darkviolet", lwd = 1.75)
  abline(v = c(1000, 100000), lty = 2, col = adjustcolor("black", alpha.f = 0.75))
})

library(writexl)

# Publication summary of the N. ridei SDM
sdm_summary <- data.frame(
  Measure = c("Feature classes", "Regularization multiplier", "Selection method", "Candidate models",
              "Delta AICc of selected model", "Number of parameters",
              "Validation AUC (mean)", "Validation AUC (SD)",
              "Validation CBI (mean)", "Validation CBI (SD)",
              "Training AUC", "AUC difference, training minus validation (mean)", "AUC difference (SD)",
              "Occurrence records (after cleaning, filtering and thinning)", "Background points",
              "Background", "Partitioning", "Number of blocks", "Correlation filter"),
  Value = c(as.character(best$fc), as.character(best$rm), "Lowest AICc", nrow(results),
            round(best$delta.AICc, 2), best$ncoef,
            round(best$auc.val.avg, 3), round(best$auc.val.sd, 3),
            round(best$cbi.val.avg, 3), round(best$cbi.val.sd, 3),
            round(best$auc.train, 3), round(best$auc.diff.avg, 3), round(best$auc.diff.sd, 3),
            nrow(occ), nrow(bg),
            "All present-day land cells, continental Australia", "Spatial block",
            length(unique(e@occs.grp)), "Pairwise |r| > 0.7 (vifcor)")
)



# MSMC2 Ne through time for N. ridei and S. murina under two mutation rates:
# devil (5.95e-09, existing colours) and koala (6.12e-09, alternative shades)
mu_devil <- 5.95e-9
mu_koala <- 6.12e-9

# Read an MSMC2 output file, scaling time and Ne with a given mutation rate
read_msmc_mu <- function(f, mu_rate) {
  x <- read.table(f, header = TRUE)
  data.frame(years = x$left_time_boundary / mu_rate * gen, Ne = (1 / x$lambda) / (2 * mu_rate))
}
nr_devil <- read_msmc_mu("whole_genome_msmc2_wongai_i50_p1-2_25-1_1-2_1-3.final.txt", mu_devil)
nr_koala <- read_msmc_mu("whole_genome_msmc2_wongai_i50_p1-2_25-1_1-2_1-3.final.txt", mu_koala)
sm_devil <- read_msmc_mu("whole_genome_msmc2_slendertailed_i50_p1-2_25-1_1-2_1-3.final.txt", mu_devil)
sm_koala <- read_msmc_mu("whole_genome_msmc2_slendertailed_i50_p1-2_25-1_1-2_1-3.final.txt", mu_koala)

x_ticks <- seq(0, 240000, by = 20000)
y_ticks <- seq(0, 1e6, by = 2e5)

save_png_svg("nridei_ausbg_smurina_msmc2_mutation_rates", 9.6, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(5, 6, 2, 2))
  plot(nr_devil$years, nr_devil$Ne, type = "n",
       xlim = c(0, 250000), ylim = c(0, 1e6), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
       xlab = "Years before present", ylab = "",
       main = "Demographic history inferred with MSMC2")
  axis(1, at = x_ticks, labels = formatC(x_ticks, format = "d", big.mark = ","), cex.axis = 0.75, gap.axis = -1)
  axis(2, at = y_ticks, labels = formatC(y_ticks, format = "d", big.mark = ","), las = 1)
  mtext("Effective population size (Ne)", side = 2, line = 4.5)
  lines(sm_devil$years, sm_devil$Ne, type = "s", col = "darkviolet", lwd = 1.75)
  lines(sm_koala$years, sm_koala$Ne, type = "s", col = "#4B0082", lwd = 1.75)
  lines(nr_devil$years, nr_devil$Ne, type = "s", col = "darkgreen", lwd = 1.75)
  lines(nr_koala$years, nr_koala$Ne, type = "s", col = "darkolivegreen4", lwd = 1.75)
  abline(v = c(1000, 100000), lty = 2, col = adjustcolor("black", alpha.f = 0.75))
  legend("topright",
         legend = c("S. murina, devil \u03bc", "S. murina, koala \u03bc",
                    "N. ridei, devil \u03bc", "N. ridei, koala \u03bc"),
         col = c("darkviolet", "#4B0082", "darkgreen", "darkolivegreen4"),
         lty = 1, lwd = 1.75, cex = 0.8, bg = "white")
  legend("bottomright",
         legend = c(paste0("devil \u03bc = ", sprintf("%.2e", mu_devil)),
                    paste0("koala \u03bc = ", sprintf("%.2e", mu_koala))),
         cex = 0.8, bg = "white")
})



# MSMC2 Ne through time for N. ridei and S. murina under two generation times (1 and 1.25 years),
# with the devil mutation rate (5.95e-09) fixed: 1 year in the existing colours, 1.25 years in alternative shades
mu_devil <- 5.95e-9
gen_short <- 1
gen_long <- 1.25

# Read an MSMC2 output file, scaling time and Ne with a given mutation rate and generation time
read_msmc_gen <- function(f, mu_rate, gen_time) {
  x <- read.table(f, header = TRUE)
  data.frame(years = x$left_time_boundary / mu_rate * gen_time, Ne = (1 / x$lambda) / (2 * mu_rate))
}
nr_gen_short <- read_msmc_gen("whole_genome_msmc2_wongai_i50_p1-2_25-1_1-2_1-3.final.txt", mu_devil, gen_short)
nr_gen_long <- read_msmc_gen("whole_genome_msmc2_wongai_i50_p1-2_25-1_1-2_1-3.final.txt", mu_devil, gen_long)
sm_gen_short <- read_msmc_gen("whole_genome_msmc2_slendertailed_i50_p1-2_25-1_1-2_1-3.final.txt", mu_devil, gen_short)
sm_gen_long <- read_msmc_gen("whole_genome_msmc2_slendertailed_i50_p1-2_25-1_1-2_1-3.final.txt", mu_devil, gen_long)

x_ticks <- seq(0, 240000, by = 20000)
y_ticks <- seq(0, 1e6, by = 2e5)

save_png_svg("nridei_ausbg_smurina_msmc2_generation_times", 9.6, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(5, 6, 2, 2))
  plot(nr_gen_short$years, nr_gen_short$Ne, type = "n",
       xlim = c(0, 250000), ylim = c(0, 1e6), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
       xlab = "Years before present", ylab = "",
       main = "Demographic history inferred with MSMC2")
  axis(1, at = x_ticks, labels = formatC(x_ticks, format = "d", big.mark = ","), cex.axis = 0.75, gap.axis = -1)
  axis(2, at = y_ticks, labels = formatC(y_ticks, format = "d", big.mark = ","), las = 1)
  mtext("Effective population size (Ne)", side = 2, line = 4.5)
  lines(sm_gen_short$years, sm_gen_short$Ne, type = "s", col = "darkviolet", lwd = 1.75)
  lines(sm_gen_long$years, sm_gen_long$Ne, type = "s", col = "#4B0082", lwd = 1.75)
  lines(nr_gen_short$years, nr_gen_short$Ne, type = "s", col = "darkgreen", lwd = 1.75)
  lines(nr_gen_long$years, nr_gen_long$Ne, type = "s", col = "darkolivegreen4", lwd = 1.75)
  abline(v = c(1000, 100000), lty = 2, col = adjustcolor("black", alpha.f = 0.75))
  legend("topright",
         legend = c("S. murina, 1 yr generation", "S. murina, 1.25 yr generation",
                    "N. ridei, 1 yr generation", "N. ridei, 1.25 yr generation"),
         col = c("darkviolet", "#4B0082", "darkgreen", "darkolivegreen4"),
         lty = 1, lwd = 1.75, cex = 0.8, bg = "white")
  legend("bottomright",
         legend = paste0("devil \u03bc = ", sprintf("%.2e", mu_devil)),
         cex = 0.8, bg = "white")
})




# Retained variables and their permutation importance (%)
sdm_variables <- var_importance[order(-var_importance$percent), c("variable", "percent")]
sdm_variables$percent <- round(sdm_variables$percent, 1)
names(sdm_variables) <- c("Variable", "Permutation importance (%)")

# Write to Excel (sheet 1: model summary, sheet 2: variables)
write_xlsx(list(Summary = sdm_summary, Variables = sdm_variables), "nridei_sdm_publication_tables.xlsx")



library(writexl)

# Beyer et al. (2020) variable descriptions from pastclim
beyer_vars <- get_vars_for_dataset(dataset = "Beyer2020", details = TRUE)

# Retained variables with permutation importance and their Beyer description
variable_table <- var_importance[order(-var_importance$percent), c("variable", "percent")]
variable_table$percent <- round(variable_table$percent, 1)
variable_table$description <- beyer_vars$long_name[match(variable_table$variable, beyer_vars$variable)]
names(variable_table) <- c("Variable", "Permutation importance (%)", "Description")

# N. ridei SDM: write to Excel
write_xlsx(list(Variables = variable_table), "nridei_ausbg_sdm_variable_table.xlsx")







