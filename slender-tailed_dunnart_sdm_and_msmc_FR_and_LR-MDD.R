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
d <- read.csv("smurina_records-2026-09-29.csv")
d <- d[!is.na(d$decimalLatitude) & !is.na(d$decimalLongitude), ]
d <- d[!duplicated(d[, c("decimalLatitude", "decimalLongitude")]), ]
d <- d[is.na(d$coordinateUncertaintyInMeters) | d$coordinateUncertaintyInMeters <= 5000, ]


# Load the IUCN range map shapefile for S. murina
# UPDATE THIS PATH to point to your downloaded IUCN shapefile
range_map <- vect("redlist_species_data_fed6d4c4-8e08-47c0-bc43-2e0174de1d5f/data_0.shp")

# IUCN range map with a buffer (terra buffer widths are in metres)
range_map_buffered <- buffer(range_map, width = 200 * 1000)

# Convert occurrence data frame to a spatial points object
occ_pts <- vect(d, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

# Extract polygon attributes for each point
# Points outside the buffered polygon will have NA for all polygon attributes
ext <- extract(range_map_buffered, occ_pts)

# Full range: keep only points inside the buffered range map
# We check if at least one polygon attribute column is not NA
inside_range <- rowSums(!is.na(ext[, -1])) > 0
d_fullrange <- d[inside_range, ]


# Murray Darling Depression bioregion (IBRA 7), merged from its subregions
ibra <- vect("IBRA7_subregions/ibra7_subregions.shp")
mdd <- project(aggregate(ibra[ibra$REG_NAME_7 == "Murray Darling Depression", ]), "EPSG:4326")

# Murray Darling Depression with a buffer (terra buffer widths are in metres)
mdd_buffered <- buffer(mdd, width = 200 * 1000)

# Murray Darling: keep only full range records that are also inside the buffered Murray Darling Depression
d_pts <- vect(d_fullrange, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
d_murraydarling <- d_fullrange[is.related(d_pts, mdd_buffered, "intersects"), ]


# Thin to one record per climate grid cell
cell <- cellFromXY(bio, as.matrix(d_fullrange[, c("decimalLongitude", "decimalLatitude")]))
d_fullrange <- d_fullrange[!is.na(cell) & !duplicated(cell), ]
occ_fullrange <- d_fullrange[, c("decimalLongitude", "decimalLatitude")]

cell <- cellFromXY(bio, as.matrix(d_murraydarling[, c("decimalLongitude", "decimalLatitude")]))
d_murraydarling <- d_murraydarling[!is.na(cell) & !duplicated(cell), ]
occ_murraydarling <- d_murraydarling[, c("decimalLongitude", "decimalLatitude")]

# Remove variables correlated above 0.7
vals <- na.omit(extract(bio, occ_fullrange, ID = FALSE))
vc <- vifcor(vals, th = 0.7)
vars_fullrange <- setdiff(names(bio), vc@excluded)
bio_fullrange <- bio[[vars_fullrange]]

vals <- na.omit(extract(bio, occ_murraydarling, ID = FALSE))
vc <- vifcor(vals, th = 0.7)
vars_murraydarling <- setdiff(names(bio), vc@excluded)
bio_murraydarling <- bio[[vars_murraydarling]]



# Background: all present-day land cells in the study area (continental Australia)
bg_fullrange <- as.data.frame(bio_fullrange, xy = TRUE, na.rm = TRUE)[, c("x", "y")]
colnames(bg_fullrange) <- c("decimalLongitude", "decimalLatitude")

bg_murraydarling <- as.data.frame(bio_murraydarling, xy = TRUE, na.rm = TRUE)[, c("x", "y")]
colnames(bg_murraydarling) <- c("decimalLongitude", "decimalLatitude")

# Fit MaxEnt models with spatial block cross-validation
e_fullrange <- ENMevaluate(occs = occ_fullrange, envs = bio_fullrange, bg = bg_fullrange,
                           algorithm = "maxnet", partitions = "block",
                           tune.args = list(fc = c("L", "LQ", "LQH"), rm = 1:4))

e_murraydarling <- ENMevaluate(occs = occ_murraydarling, envs = bio_murraydarling, bg = bg_murraydarling,
                               algorithm = "maxnet", partitions = "block",
                               tune.args = list(fc = c("L", "LQ", "LQH"), rm = 1:4))

# Select the best model (lowest AICc)
results_fullrange <- eval.results(e_fullrange)
best_fullrange <- results_fullrange[which(results_fullrange$delta.AICc == 0), ]
mod_fullrange <- eval.models(e_fullrange)[[as.character(best_fullrange$tune.args)]]

results_murraydarling <- eval.results(e_murraydarling)
best_murraydarling <- results_murraydarling[which(results_murraydarling$delta.AICc == 0), ]
mod_murraydarling <- eval.models(e_murraydarling)[[as.character(best_murraydarling$tune.args)]]

# Permutation importance: drop in AUC when each variable is shuffled (10 shuffles each)
auc_fun <- function(pred, pa) {
  r <- rank(pred)
  n1 <- sum(pa == 1)
  n0 <- sum(pa == 0)
  (sum(r[pa == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Full range
p_env <- extract(bio_fullrange, occ_fullrange, ID = FALSE)
b_env <- extract(bio_fullrange, bg_fullrange, ID = FALSE)
env_all <- rbind(p_env, b_env)
pa <- c(rep(1, nrow(p_env)), rep(0, nrow(b_env)))
keep <- complete.cases(env_all)
env_all <- env_all[keep, ]
pa <- pa[keep]

base_auc <- auc_fun(as.numeric(predict(mod_fullrange, env_all, type = "cloglog")), pa)
set.seed(1)
perm_imp <- sapply(vars_fullrange, function(v) mean(replicate(10, {
  shuffled <- env_all
  shuffled[[v]] <- sample(shuffled[[v]])
  base_auc - auc_fun(as.numeric(predict(mod_fullrange, shuffled, type = "cloglog")), pa)
})))

var_importance_fullrange <- data.frame(variable = vars_fullrange, auc_drop = perm_imp,
                                       percent = 100 * pmax(perm_imp, 0) / sum(pmax(perm_imp, 0)))
var_importance_fullrange[order(-var_importance_fullrange$auc_drop), ]

# Murray Darling
p_env <- extract(bio_murraydarling, occ_murraydarling, ID = FALSE)
b_env <- extract(bio_murraydarling, bg_murraydarling, ID = FALSE)
env_all <- rbind(p_env, b_env)
pa <- c(rep(1, nrow(p_env)), rep(0, nrow(b_env)))
keep <- complete.cases(env_all)
env_all <- env_all[keep, ]
pa <- pa[keep]

base_auc <- auc_fun(as.numeric(predict(mod_murraydarling, env_all, type = "cloglog")), pa)
set.seed(1)
perm_imp <- sapply(vars_murraydarling, function(v) mean(replicate(10, {
  shuffled <- env_all
  shuffled[[v]] <- sample(shuffled[[v]])
  base_auc - auc_fun(as.numeric(predict(mod_murraydarling, shuffled, type = "cloglog")), pa)
})))

var_importance_murraydarling <- data.frame(variable = vars_murraydarling, auc_drop = perm_imp,
                                           percent = 100 * pmax(perm_imp, 0) / sum(pmax(perm_imp, 0)))
var_importance_murraydarling[order(-var_importance_murraydarling$auc_drop), ]




# Present-day suitability
suit_cur_fullrange <- predict(bio_fullrange, mod_fullrange, type = "cloglog", na.rm = TRUE)
suit_cur_murraydarling <- predict(bio_murraydarling, mod_murraydarling, type = "cloglog", na.rm = TRUE)

# MSMC2 sample location
msmc_sample <- data.frame(decimalLongitude = 140.511497, decimalLatitude = -33.633467)

# Check results
nrow(occ_fullrange)
vars_fullrange
best_fullrange[, c("tune.args", "auc.val.avg", "AICc")]

nrow(occ_murraydarling)
vars_murraydarling
best_murraydarling[, c("tune.args", "auc.val.avg", "AICc")]

# Time points for the SDM maps, every 20 ky from 0 to 100 kya
times <- seq(0, 100000, by = 20000)

# MSMC2 Ne vs mean climatic suitability at the MSMC2 sample site, per MSMC2 time step (midpoints up to 100 kya)
mu <- 5.95e-9
gen <- 1
msmc <- read.table("whole_genome_msmc2_slendertailed_i50_p1-2_25-1_1-2_1-3.final.txt", header = TRUE)
msmc$start_years <- msmc$left_time_boundary / mu * gen
msmc$end_years <- msmc$right_time_boundary / mu * gen
msmc$mid_years <- (msmc$start_years + msmc$end_years) / 2
msmc$Ne <- (1 / msmc$lambda) / (2 * mu)
msmc_s <- msmc[msmc$mid_years <= 100000, ]

# Climate and vegetation at the sample site for every climate time step these MSMC2 steps span
clim_times_s <- get_time_bp_steps(dataset = "Beyer2020")
clim_times_s <- clim_times_s[-clim_times_s <= max(msmc_s$end_years)]
site_clim_s <- do.call(rbind, lapply(clim_times_s, function(t)
  extract(get_clim(t), msmc_sample, ID = FALSE)))

# Suitability at the sample site for every climate time step these MSMC2 steps span
clim_suit_s_fullrange <- as.numeric(predict(mod_fullrange, site_clim_s[, vars_fullrange, drop = FALSE], type = "cloglog"))
clim_suit_s_murraydarling <- as.numeric(predict(mod_murraydarling, site_clim_s[, vars_murraydarling, drop = FALSE], type = "cloglog"))

# Mean suitability over each MSMC2 step (nearest climate time step if none fall inside)
ne_vs_suit_fullrange <- data.frame(
  kya = msmc_s$mid_years / 1000,
  Ne = msmc_s$Ne,
  suitability = sapply(seq_len(nrow(msmc_s)), function(i) {
    inside <- -clim_times_s >= msmc_s$start_years[i] & -clim_times_s < msmc_s$end_years[i]
    if (any(inside)) mean(clim_suit_s_fullrange[inside]) else clim_suit_s_fullrange[which.min(abs(-clim_times_s - msmc_s$mid_years[i]))]
  })
)
ne_vs_suit_fullrange

ne_vs_suit_murraydarling <- data.frame(
  kya = msmc_s$mid_years / 1000,
  Ne = msmc_s$Ne,
  suitability = sapply(seq_len(nrow(msmc_s)), function(i) {
    inside <- -clim_times_s >= msmc_s$start_years[i] & -clim_times_s < msmc_s$end_years[i]
    if (any(inside)) mean(clim_suit_s_murraydarling[inside]) else clim_suit_s_murraydarling[which.min(abs(-clim_times_s - msmc_s$mid_years[i]))]
  })
)
ne_vs_suit_murraydarling

# Spearman rank correlation of Ne and mean suitability across MSMC2 steps
cor.test(ne_vs_suit_fullrange$Ne, ne_vs_suit_fullrange$suitability, method = "spearman")
cor.test(ne_vs_suit_murraydarling$Ne, ne_vs_suit_murraydarling$suitability, method = "spearman")





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
save_png_svg("smurina_fullrange_sdm_present_records", 8, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(3, 3, 2, 5))
  plot(suit_cur_fullrange, main = expression("Present-day SDM suitability (" * italic("S. murina") * ", full range)"), range = c(0, 1))
  points(occ_fullrange, pch = 20, cex = 1.0, col = "black")
  points(msmc_sample, pch = 20, col = "red", cex = 1.5)
})

save_png_svg("smurina_murraydarling_sdm_present_records", 8, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(3, 3, 2, 5))
  plot(suit_cur_murraydarling, main = expression("Present-day SDM suitability (" * italic("S. murina") * ", local range)"), range = c(0, 1))
  points(occ_murraydarling, pch = 20, cex = 1.0, col = "black")
  points(msmc_sample, pch = 20, col = "red", cex = 1.5)
})

# 2) SDM at each time point, 100 to 0 kya in 20 ky steps, as 2 x 3 panels, with the MSMC2 sample (red)
# Each map is masked to present-day land so all panels show modern Australia's outline
suit_times_fullrange <- lapply(rev(times), function(t) mask(predict(get_clim(-t)[[vars_fullrange]], mod_fullrange, type = "cloglog", na.rm = TRUE), suit_cur_fullrange))
save_png_svg("smurina_fullrange_sdm_timepoints", 14.4, 8.8, function() {
  par(mfrow = c(2, 3), mar = c(3, 3, 2, 5))
  for (i in seq_along(suit_times_fullrange)) {
    t <- rev(times)[i]
    label <- if (t == 0) "Present" else paste0(t / 1000, " ka")
    plot(suit_times_fullrange[[i]], main = label, range = c(0, 1))
    points(msmc_sample, pch = 20, col = "red", cex = 1.5)
  }
})

suit_times_murraydarling <- lapply(rev(times), function(t) mask(predict(get_clim(-t)[[vars_murraydarling]], mod_murraydarling, type = "cloglog", na.rm = TRUE), suit_cur_murraydarling))
save_png_svg("smurina_murraydarling_sdm_timepoints", 14.4, 8.8, function() {
  par(mfrow = c(2, 3), mar = c(3, 3, 2, 5))
  for (i in seq_along(suit_times_murraydarling)) {
    t <- rev(times)[i]
    label <- if (t == 0) "Present" else paste0(t / 1000, " ka")
    plot(suit_times_murraydarling[[i]], main = label, range = c(0, 1))
    points(msmc_sample, pch = 20, col = "red", cex = 1.5)
  }
})

# 3) Ne (purple) vs mean suitability (dark orange) at the MSMC2 sample site per MSMC2 step, full range SDM
save_png_svg("smurina_fullrange_ne_vs_suitability", 9.6, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 2, 5), cex.lab = 1.25, cex.axis = 1.2)
  plot(ne_vs_suit_fullrange$kya, ne_vs_suit_fullrange$Ne, type = "n",
       xlim = c(0, 100), ylim = c(0, max(ne_vs_suit_fullrange$Ne) * 1.1),
       xlab = "Thousand years before present", ylab = "Effective population size (Ne)",
       main = expression(bold("Ne and SDM suitability per MSMC2 time step (" * bolditalic("S. murina") * ", full range)")))
  lines(ne_vs_suit_fullrange$kya, ne_vs_suit_fullrange$Ne, type = "o", pch = 19, col = "darkviolet")
  par(new = TRUE)
  plot(ne_vs_suit_fullrange$kya, ne_vs_suit_fullrange$suitability, type = "o", pch = 17, col = "darkorange",
       xlim = c(0, 100), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext("Mean climatic suitability at MSMC2 sample site", side = 4, line = 3, cex = 1.25)
  legend("topright", legend = c(expression(italic("S. murina") ~ "Ne"), "SDM suitability (full range)"),
         pch = c(19, 17), lty = c(1, 1),
         col = c("darkviolet", "darkorange"), pt.cex = c(1, 1), cex = 0.95)
})

# 3) Ne (purple) vs mean suitability (dark orange) at the MSMC2 sample site per MSMC2 step, local range SDM
save_png_svg("smurina_murraydarling_ne_vs_suitability", 9.6, 7.2, function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 2, 5), cex.lab = 1.25, cex.axis = 1.2)
  plot(ne_vs_suit_murraydarling$kya, ne_vs_suit_murraydarling$Ne, type = "n",
       xlim = c(0, 100), ylim = c(0, max(ne_vs_suit_murraydarling$Ne) * 1.1),
       xlab = "Thousand years before present", ylab = "Effective population size (Ne)",
       main = expression(bold("Ne and SDM suitability per MSMC2 time step (" * bolditalic("S. murina") * ", local range)")))
  lines(ne_vs_suit_murraydarling$kya, ne_vs_suit_murraydarling$Ne, type = "o", pch = 19, col = "darkviolet")
  par(new = TRUE)
  plot(ne_vs_suit_murraydarling$kya, ne_vs_suit_murraydarling$suitability, type = "o", pch = 17, col = "darkorange",
       xlim = c(0, 100), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext("Mean climatic suitability at MSMC2 sample site", side = 4, line = 3, cex = 1.25)
  legend("topright", legend = c(expression(italic("S. murina") ~ "Ne"), "SDM suitability (local range)"),
         pch = c(19, 17), lty = c(1, 1),
         col = c("darkviolet", "darkorange"), pt.cex = c(1, 1), cex = 0.95)
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

# Publication summary of one SDM: model summary values and retained variables with permutation importance
sdm_publication_tables <- function(e, results, best, occ, bg, var_importance) {
  summary_table <- data.frame(
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
  variables_table <- var_importance[order(-var_importance$percent), c("variable", "percent")]
  variables_table$percent <- round(variables_table$percent, 1)
  names(variables_table) <- c("Variable", "Permutation importance (%)")
  list(Summary = summary_table, Variables = variables_table)
}

# S. murina, full range SDM: write to Excel (sheet 1: model summary, sheet 2: variables)
sdm_fullrange <- sdm_publication_tables(e_fullrange, results_fullrange, best_fullrange,
                                        occ_fullrange, bg_fullrange, var_importance_fullrange)
write_xlsx(sdm_fullrange, "smurina_fullrange_sdm_publication_tables.xlsx")

# S. murina, local range SDM: write to Excel (sheet 1: model summary, sheet 2: variables)
sdm_murraydarling <- sdm_publication_tables(e_murraydarling, results_murraydarling, best_murraydarling,
                                            occ_murraydarling, bg_murraydarling, var_importance_murraydarling)
write_xlsx(sdm_murraydarling, "smurina_murraydarling_sdm_publication_tables.xlsx")



library(writexl)

# Beyer et al. (2020) variable descriptions from pastclim
beyer_vars <- get_vars_for_dataset(dataset = "Beyer2020", details = TRUE)

# Full range: retained variables with permutation importance and their Beyer description
variable_table_fullrange <- var_importance_fullrange[order(-var_importance_fullrange$percent), c("variable", "percent")]
variable_table_fullrange$percent <- round(variable_table_fullrange$percent, 1)
variable_table_fullrange$description <- beyer_vars$long_name[match(variable_table_fullrange$variable, beyer_vars$variable)]
names(variable_table_fullrange) <- c("Variable", "Permutation importance (%)", "Description")

# S. murina, full range SDM: write to Excel
write_xlsx(list(Variables = variable_table_fullrange), "smurina_fullrange_sdm_variable_table.xlsx")

# Local range: retained variables with permutation importance and their Beyer description
variable_table_murraydarling <- var_importance_murraydarling[order(-var_importance_murraydarling$percent), c("variable", "percent")]
variable_table_murraydarling$percent <- round(variable_table_murraydarling$percent, 1)
variable_table_murraydarling$description <- beyer_vars$long_name[match(variable_table_murraydarling$variable, beyer_vars$variable)]
names(variable_table_murraydarling) <- c("Variable", "Permutation importance (%)", "Description")

# S. murina, local range SDM: write to Excel
write_xlsx(list(Variables = variable_table_murraydarling), "smurina_murraydarling_sdm_variable_table.xlsx")


