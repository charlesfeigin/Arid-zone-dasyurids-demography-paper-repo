library(terra)
library(geodata)

# Rivers with channels in the river discharge analysis (Hesse et al. 2018)
rivers <- c("Murray", "Darling", "Murrumbidgee", "Lachlan", "Macquarie", "Gwydir", "Namoi", "Goulburn")

# Area covering these rivers and the Murray Darling Depression
study_area <- ext(138, 152.5, -38.5, -28.5)

# Murray Darling Depression bioregion (IBRA 7), merged from its subregions
ibra <- vect("IBRA7_subregions/ibra7_subregions.shp")
mdd <- aggregate(ibra[ibra$REG_NAME_7 == "Murray Darling Depression", ])

# Read every line in the study area from Surface Hydrology Lines (National) by location
# (uses the .gdb spatial index), once, and save it; later runs read the saved file
lines_file <- "mdb_hydrolines.gpkg"
if (!file.exists(lines_file)) {
  gdb <- list.dirs(getwd(), recursive = TRUE)
  gdb <- gdb[grepl("SurfaceHydrologyLinesNational\\.gdb$", gdb)]
  gdb <- gdb[sapply(gdb, function(p) any(grepl("\\.gdbtable$", list.files(p))))][1]
  message("Reading lines in the study area from ", gdb)
  hl <- vect(gdb, layer = "HydroLines", extent = study_area)
  message("Read ", nrow(hl), " lines; saving to ", lines_file)
  writeVector(hl[, c("NAME", "FEATURETYPE", "HIERARCHY")], lines_file)
}
hl <- vect(lines_file)
mdd <- project(mdd, crs(hl))
message(nrow(hl), " lines loaded")

# Start and end point of every line
g <- geom(hl)
first_pt <- g[!duplicated(g[, "geom"]), c("x", "y")]
last_pt <- g[!duplicated(g[, "geom"], fromLast = TRUE), c("x", "y")]
start_key <- paste(round(first_pt[, "x"], 7), round(first_pt[, "y"], 7))
end_key <- paste(round(last_pt[, "x"], 7), round(last_pt[, "y"], 7))

# Flow direction test: where Murray River lines near Mildura meet, one line's end should meet the next line's start
murray_test <- which(hl$NAME == "MURRAY RIVER" & first_pt[, "x"] > 141.5 & first_pt[, "x"] < 143 &
                     first_pt[, "y"] > -34.6 & first_pt[, "y"] < -33.9)
end_to_start <- sum(end_key[murray_test] %in% start_key[murray_test])
message("Flow direction test: ", end_to_start, " of ", length(murray_test),
        " Murray lines near Mildura end where another starts")
if (end_to_start / length(murray_test) < 0.9) stop("Lines are not consistently drawn in the direction of flow")

# Junction test: the Darling's downstream end should be exactly where a Murray line starts
darling <- which(hl$NAME == "DARLING RIVER")
darling_ends <- end_key[darling][!(end_key[darling] %in% start_key[darling])]
darling_joins_murray <- any(darling_ends %in% start_key[hl$NAME == "MURRAY RIVER"])
message("Junction test: Darling's downstream end meets the start of a Murray line: ", darling_joins_murray)
if (!darling_joins_murray) stop("Rivers are not split at junctions, so tracing from line ends would miss joins")

# Rivers by name in the data (excluding the NSW Goulburn, a Hunter tributary)
river_names <- c(paste(toupper(rivers), "RIVER"), paste("RIVER", toupper(rivers)))
seeds <- which(hl$NAME %in% river_names)
nsw_goulburn <- hl$NAME[seeds] == "GOULBURN RIVER" & first_pt[seeds, "y"] > -35.5
seeds <- seeds[!nsw_goulburn]
message(length(seeds), " lines named as the rivers")

# Follow downstream: from each line, the lines that start where it ends, until nothing further is found
downstream <- split(seq_along(start_key), start_key)
traced <- rep(FALSE, nrow(hl))
traced[seeds] <- TRUE
current <- seeds
while (length(current) > 0) {
  nxt <- unique(unlist(downstream[end_key[current]], use.names = FALSE))
  nxt <- nxt[!traced[nxt]]
  traced[nxt] <- TRUE
  current <- nxt
}
message(sum(traced), " lines in the rivers and everything downstream of them (",
        sum(traced) - length(seeds), " added by following downstream)")

# Lines whose start point is inside the MDD (checked only for lines within the MDD's bounding box)
mdd_box <- ext(mdd)
in_box <- which(first_pt[, "x"] >= xmin(mdd_box) & first_pt[, "x"] <= xmax(mdd_box) &
                first_pt[, "y"] >= ymin(mdd_box) & first_pt[, "y"] <= ymax(mdd_box))
in_mdd <- rep(FALSE, nrow(hl))
in_mdd[in_box] <- relate(vect(first_pt[in_box, , drop = FALSE], crs = crs(hl)), mdd, "intersects")[, 1]

# Waterways inside the MDD that flow into the rivers or anything downstream of them:
# follow upstream (lines that end where a connected line starts), keeping only lines inside the MDD
upstream <- split(seq_along(end_key), end_key)
connected <- traced
current <- which(traced)
while (length(current) > 0) {
  nxt <- unique(unlist(upstream[start_key[current]], use.names = FALSE))
  nxt <- nxt[!connected[nxt] & in_mdd[nxt]]
  connected[nxt] <- TRUE
  current <- nxt
}
message(sum(connected) - sum(traced), " MDD waterways added that flow into the rivers or their downstream network")
river_network <- hl[connected]

# Label position for each river: the named line closest to the middle of that river's named lines
river_labels <- do.call(rbind, lapply(rivers, function(r) {
  idx <- seeds[hl$NAME[seeds] %in% c(paste(toupper(r), "RIVER"), paste("RIVER", toupper(r)))]
  cents <- crds(centroids(hl[idx]))
  mid <- colMeans(cents)
  closest <- which.min((cents[, 1] - mid[1])^2 + (cents[, 2] - mid[2])^2)
  data.frame(river = r, x = cents[closest, 1], y = cents[closest, 2])
}))
river_labels <- vect(river_labels, geom = c("x", "y"), crs = crs(hl))

# MSMC2 sample location (S. murina)
msmc_sample <- vect(data.frame(x = 140.511497, y = -33.633467), geom = c("x", "y"), crs = "EPSG:4326")

# Australia outline: GADM country boundary (level 0), downloaded once into gadm_data
dir.create("gadm_data", showWarnings = FALSE)
aus_outline <- project(gadm("AUS", level = 0, path = "gadm_data"), crs(hl))

# Map extent zoomed out 2x around the centre of the study area
centre <- c(mean(c(xmin(study_area), xmax(study_area))), mean(c(ymin(study_area), ymax(study_area))))
half_w <- (xmax(study_area) - xmin(study_area))
half_h <- (ymax(study_area) - ymin(study_area))
map_area <- ext(centre[1] - half_w, centre[1] + half_w, centre[2] - half_h, centre[2] + half_h)

# Map: land (grey), MDD outline (black), Hesse et al. (2018) rivers, everything downstream of them
# and the MDD waterways connected to them (blue), each river labelled, and the MSMC2 sample (red)
# The plot area is set up with base R (aspect ratio corrected for latitude), so its coordinates
# match the plot frame and the legend sits in the frame's top-right corner
plot_river_map <- function() {
  par(mfrow = c(1, 1), mar = c(3, 3, 2.5, 1), cex.axis = 1.2, cex.main = 1.45)
  plot.new()
  plot.window(xlim = c(xmin(map_area), xmax(map_area)), ylim = c(ymin(map_area), ymax(map_area)),
              asp = 1 / cos(mean(c(ymin(map_area), ymax(map_area))) * pi / 180), xaxs = "i", yaxs = "i")
  plot(aus_outline, add = TRUE, col = "grey75", border = "grey40")
  plot(mdd, add = TRUE, border = "black", lwd = 2)
  plot(river_network, add = TRUE, col = "blue", lwd = 1.2)
  text(river_labels, labels = river_labels$river, font = 2, cex = 1.05, halo = TRUE)
  points(msmc_sample, pch = 19, col = "red", cex = 1.2)
  box()
  axis(1)
  axis(2)
  title(main = "Rivers with dated palaeochannels flowing into the Murray Darling Depression")
  legend("topright",
         legend = c("Rivers (Hesse et al. 2018) and connected waterways", "Murray Darling Depression",
                    expression(italic("S. murina") ~ "MSMC2 sample")),
         col = c("blue", "black", "red"), lty = c(1, 1, NA), lwd = c(1.2, 2, NA), pch = c(NA, NA, 19),
         cex = 1.0, bg = "white")
}
plot_river_map()

# Write the map as both PNG and SVG
# Sizes are in inches; the PNG is written at 250 pixels per inch
save_png_svg <- function(name, width_in, height_in, plot_fun) {
  png(paste0(name, ".png"), width = width_in * 250, height = height_in * 250, res = 250)
  plot_fun()
  dev.off()
  svg(paste0(name, ".svg"), width = width_in, height = height_in)
  plot_fun()
  dev.off()
}
save_png_svg("smurina_rivers_mdd_map", 10, 8, plot_river_map)