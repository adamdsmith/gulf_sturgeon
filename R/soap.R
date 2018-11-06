# Generate potential knot locations
make_knots <- function(sf_poly, cellsize_xy = 5000, 
                       bnd_buffer = 0,
                       random_start = TRUE, seed = NULL) {
  stopifnot(is.numeric(cellsize_xy) || is.integer(cellsize_xy))
  if (length(cellsize_xy) == 1) cellsize_xy <- rep(cellsize_xy, 2)
  
  # Calculate nrows/ncols
  xdiff <- diff(st_bbox(sf_poly)[c("xmin", "xmax")])
  ydiff <- diff(st_bbox(sf_poly)[c("ymin", "ymax")])
  nr <- floor(ydiff/cellsize_xy[2])
  nc <- floor(xdiff/cellsize_xy[1])
  
  if (random_start) {
    if (!is.null(seed)) set.seed(seed)
    shift <- runif(2, -cellsize_xy, cellsize_xy)
  } else shift <- c(0, 0)
  
  # Make grid of points based on extent of sf_poly and cellsize_xy
  knots <- st_make_grid(sf_poly, n = c(nc, nr), what = "centers",
                        offset = st_bbox(sf_poly)[1:2] + shift) %>% st_sf()

  if (bnd_buffer > 0)
    sf_poly <- st_buffer(sf_poly, -bnd_buffer)
  
  # Filter to within study area
  knots <- st_join(knots, sf_poly, left = FALSE)
  
  knots <- st_coordinates(knots)[, 1:2] %>%
    as.data.frame(knots)
  names(knots) <- c("x", "y")
  knots
}

# Check soap film surface
# Modified slightly from David Miller
# https://github.com/dill/soap_checker
soap_check <- function(bnd, knots=NULL, data=NULL, plot=TRUE,
                       tol=sqrt(.Machine$double.eps)){
  
  if (!requireNamespace("sp", quietly = TRUE)) install.packages("sp", quiet = TRUE)
  if (!requireNamespace("rgeos", quietly = TRUE)) install.packages("sp", quiet = TRUE)
  if (!requireNamespace("mgcv", quietly = TRUE)) install.packages("sp", quiet = TRUE)
  ## check that the boundary makes sense
  # check that boundary is a list
  stopifnot(is.list(bnd))
  
  # check that the boundary (or boundary part) have x and y elements
  lapply(bnd, function(x) stopifnot(c("x","y") %in% names(x)))
  # each boundary part must have at least 4 elements!
  lapply(bnd, function(x) stopifnot(length(x$x)>3, length(x$y)>3))
  
  # check that the boundary loops are actually loops
  check_ends <- function(x, tol){
    all.equal(c(x$x[1], x$y[1]),
              c(x$x[length(x$y)], x$y[length(x$y)]),tolerance=tol)
  }
  end_check <- unlist(lapply(bnd, check_ends, tol=tol))
  end_check_logical <- is.character(end_check)
  if(any(end_check_logical)){
    stop(paste("Boundary loop(s)",which(end_check_logical),
               "don't have identical start & end points",collapse=" "))
  }
  
  islands <- FALSE
  # check for intersections
  if(length(bnd)>1){
    inds <- combn(1:length(bnd),2)
    
    ## make the bnds into polys here
    make_bnd_poly <- function(bnd){
      
      bnd$x[length(bnd$x)] <-  bnd$x[1]
      bnd$y[length(bnd$y)] <-  bnd$y[1]
      sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(bnd)),ID=1)))
    }
    bnd_poly <- lapply(bnd, make_bnd_poly)
    
    # function to see if two polygons intersect
    intersects <- function(this.ind, bnd){
      poly1 <- bnd[[this.ind[1]]]
      poly2 <- bnd[[this.ind[2]]]
      
      rgeos::gIntersects(poly1, poly2)
    }
    # apply over all the combinations
    inter <- apply(inds, 2, intersects, bnd=bnd_poly)
    
    if(any(inter)){
      # get the index for the prospective "outer" loop
      outer_ind <- which.max(unlist(lapply(bnd_poly, rgeos::gArea)))
      outer_bnd <- bnd_poly[[outer_ind]]
      
      other_bnd <- bnd_poly
      other_bnd[[outer_ind]] <- NULL
      
      # is everything else inside that?
      islands <- unlist(lapply(other_bnd, rgeos::gWithin, spgeom2=outer_bnd))
      
      if(!all(islands)){
        stop(paste("Polygon parts",
                   paste0(apply(inds[,inter, drop=FALSE], 2, paste0,
                                collapse=" and "),
                          collapse=", "), "intersect"))
      }
      
      islands <- all(islands)
    }
  }
  
  ## plot what the boundary is
  # highlighting the area to be modelled
  if(plot){
    # colourblind-safe colours from colorbrewer2 "qualitative" map
    red <- "#d95f02"
    # if the boundary is only 1 part, plotting is rather easier
    if(!islands){
      plot(bnd[[1]], type="l", main="Red indicates soap film surface", asp=1)
      lapply(bnd, polygon, col=red)
    }else{
      outer_bnd <- bnd[[outer_ind]]
      other_bnd <- bnd
      other_bnd[[outer_ind]] <- NULL
      plot(outer_bnd, type="n", main="Red indicates soap film surface", asp=1)
      # plot the outer loop
      polygon(outer_bnd, col=red)
      # plot the other polygons on top in white
      lapply(other_bnd, polygon, col="white")
    }
  }
  
  # function to check if points are inside the boundary
  point_check <- function(bnd, x, y, type){
    
    if(length(bnd)>1){
      # inSide doesn't deal with edge points very well
      # but does handle multiple rings better
      inout <- mgcv::inSide(bnd, x, y)
    }else{
      # use sp::point.in.polygon
      # see ?point.in.polygon for returned codes, 1 is inside
      pip <- function(bnd, x, y){
        sp::point.in.polygon(x, y, bnd$x, bnd$y)==1
      }
      # apply over the parts of the polygon
      inout <- pip(bnd[[1]], x, y)
    }
    if(!all(inout)){
      stop(paste(type, paste(which(!inout),collapse=", "),
                 "are outside the boundary."))
    }
  }
  
  ## check the knots
  if(!is.null(knots)){
    # check that the points have x and y elements
    stopifnot(c("x","y") %in% names(knots))
    point_check(bnd, knots$x, knots$y, "Knots")
    if(plot) points(knots, col="#1b9e77", pch=19)
  }
  
  ## check the data
  if(!is.null(data)){
    # check that the points have x and y elements
    stopifnot(c("x","y") %in% names(data))
    point_check(bnd, data$x, data$y, "Data points")
    if(plot) points(data, col="#7570b3", pch=19)
  }
  
  return(TRUE)
}