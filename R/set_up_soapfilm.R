# Install and load some needed packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, sf, rmapshaper, ggplot2)

# Read in shapefiles
shp <- list.files("Geodata", pattern = ".shp$", full.names = TRUE)
nms <- tools::file_path_sans_ext(basename(shp))
shps <- lapply(shp, st_read, quiet = TRUE)
names(shps) <- nms

# Extract common column names
keep_cols <- Reduce(intersect, lapply(shps, names))
shps <- lapply(shps, function(i) select(i, one_of(keep_cols)))
shps <- do.call(rbind, shps)

# Dissolve them where they overlap to create study area (sa) polygon
sa <- shps %>% group_by(Elevation) %>% summarize() %>% st_zm()

# Correct geometry
sa <- st_buffer(sa, 0)
st_is_valid(sa)

# Simplify the boundaries or mgcv will revolt
# This is trial and error, so you may need to come back and reduce further
sa5 <- ms_simplify(sa, keep = 0.05)

# Plot them for comparison
p_sa <- ggplot() + geom_sf(data = sa[1], fill = "lightblue") + ggtitle("Original shapes")
p_sa5 <- ggplot() + geom_sf(data = sa5[1], fill = "lightblue") + ggtitle("5% of vertices retained")
gridExtra::grid.arrange(p_sa, p_sa5, nrow = 1)

# Put boundary in the format (list of lists) expected by mgcv
sa5_df <- st_coordinates(sa5) %>% as.data.frame()
pieces <- unique(sa5_df$L1)
bound <- lapply(pieces, function(i) {
  tmp <- filter(sa5_df, L1 == i)
  list(x = tmp$X, y = tmp$Y)
})

# Define boundary
source("R/soap.R")
knots <- make_knots(sa5, cellsize_xy = 5000)

# I suspect you'll want to add some knots manually to the river systems
# They're so small, relative to the study area, that knots rarely fall within them...
# Maybe add them every `cellsize_xy` meters along the channel away from the larger water bodies?
soap_check(bound, knots)

# Then you can pass these objects to the soap film smoother in your GAM
# For example:

# mygam <- gam(ndetects ~ s(x, y, k = 20, bs = "so", xt = list(bnd = bound)),
#              data = mydata, family = binomial, weights = navail, 
#              method = "REML", knots = knots)
