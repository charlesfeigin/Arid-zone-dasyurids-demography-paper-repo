# MSMC2 Ne vs mean palaeo-to-modern river discharge ratio (Hesse et al. 2018), per MSMC2 time step
mu <- 5.95e-9
gen <- 1
msmc <- read.table("whole_genome_msmc2_slendertailed_i50_p1-2_25-1_1-2_1-3.final.txt", header = TRUE)
msmc$start_years <- msmc$left_time_boundary / mu * gen
msmc$end_years <- msmc$right_time_boundary / mu * gen
msmc$mid_years <- (msmc$start_years + msmc$end_years) / 2
msmc$Ne <- (1 / msmc$lambda) / (2 * mu)
msmc <- msmc[msmc$mid_years <= 100000, ]

# River channels with both an age and a discharge ratio
h <- read.csv("hesse2018_table3.csv")
h <- h[!is.na(h$median_age_ka) & !is.na(h$Qbf_ratio), ]

# Merge rows that are the same dated channel (same river and same age) into one averaged row
h <- aggregate(Qbf_ratio ~ river + median_age_ka, data = h, FUN = mean)
h_years <- h$median_age_ka * 1000

# Average discharge ratio of the channels dated within each MSMC2 step
msmc$n_channels <- sapply(seq_len(nrow(msmc)), function(i)
  sum(h_years >= msmc$start_years[i] & h_years < msmc$end_years[i]))
msmc$mean_Qbf_ratio <- sapply(seq_len(nrow(msmc)), function(i) {
  inside <- h_years >= msmc$start_years[i] & h_years < msmc$end_years[i]
  if (any(inside)) mean(h$Qbf_ratio[inside]) else NA
})

ne_vs_flow <- data.frame(start_kya = msmc$start_years / 1000,
                         end_kya = msmc$end_years / 1000,
                         mid_kya = msmc$mid_years / 1000,
                         Ne = msmc$Ne,
                         n_channels = msmc$n_channels,
                         mean_Qbf_ratio = msmc$mean_Qbf_ratio)
ne_vs_flow

# MSMC2 steps with at least one dated channel
flow_ok <- ne_vs_flow[!is.na(ne_vs_flow$mean_Qbf_ratio), ]

# Spearman rank correlation of Ne and discharge ratio across MSMC2 steps
cor.test(flow_ok$Ne, flow_ok$mean_Qbf_ratio, method = "spearman")


# Plot Ne (left axis) and mean discharge ratio (right axis, log scale)
plot_ne_vs_flow <- function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 2, 5), cex.lab = 1.25, cex.axis = 1.2)
  plot(ne_vs_flow$mid_kya, ne_vs_flow$Ne, type = "o", pch = 19, col = "purple",
       xlim = c(0, 100), ylim = c(0, max(ne_vs_flow$Ne) * 1.1),
       xlab = "Thousand years before present", ylab = "Effective population size (Ne)",
       main = expression(bold("Ne and river discharge ratio per MSMC2 time step (" * bolditalic("S. murina") * ")")))
  par(new = TRUE)
  plot(ne_vs_flow$mid_kya, ne_vs_flow$mean_Qbf_ratio, type = "o", pch = 17, col = "steelblue",
       xlim = c(0, 100), log = "y", axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext("Mean palaeo/modern bankfull discharge ratio (log scale)", side = 4, line = 3, cex = 1.25)
  legend("topright", legend = c(expression(italic("S. murina") ~ "Ne"), "River discharge ratio"),
         pch = c(19, 17), lty = 1, col = c("purple", "steelblue"), cex = 0.95)
}
plot_ne_vs_flow()

# Write the plot as both PNG and SVG
# Sizes are in inches; the PNG is written at 250 pixels per inch
save_png_svg <- function(name, width_in, height_in, plot_fun) {
  png(paste0(name, ".png"), width = width_in * 250, height = height_in * 250, res = 250)
  plot_fun()
  dev.off()
  svg(paste0(name, ".svg"), width = width_in, height = height_in)
  plot_fun()
  dev.off()
}
save_png_svg("smurina_ne_vs_river_discharge", 9.6, 7.2, plot_ne_vs_flow)