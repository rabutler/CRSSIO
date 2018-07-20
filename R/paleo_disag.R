
#' Spatial and temporal diaggregattion of paleo flow data
#' 
#' The default parameter values are setup to perform the typical disaggregation
#' for CRSS, based on scaling the Upper Basin inflows, and performing a 
#' direct resampling of Lower Basin tributaries, i.e., no scale factor is 
#' applied to the LB tributaries. Sites 1-20 are scaled (`sf_sites = 1:20`); 
#' therefore the remaining sites (21:29) are not scaled.
#' 
#' @param x The annual paleo data to disaagregate.
#' @param ann_flw Observed annual flow data used for picking analog year.
#' @param mon_flw Intervening monthly natural flow. Used for spatially and 
#'   temporaly disaggregating paleo data based on the index year.
#' @param nsite The number of sites to diaggregate the data too. 
#' @param sf_sites The site numbers (indeces), that will scale the index year's
#'   volume based on the annual flow being disaggregated. The remaining sites
#'   will select the index year directly. See **Details**.
#' @param nsim Number of times to repeat the space/time disaggregation.
#' @param ofolder Optional. If specified, the disaggregated flow data and the 
#'   selected index years are saved to this folder as csv files.
#' @param index_years Optional. If specified, these index years will be used 
#'   instead of selecting years based on weights and sampling. 
#' @param k_weights If `NULL`, parameters are set based on definitions in Nowak
#'   et al. (2010). Users may force `k` and the `weights` by specifiying this 
#'   argument. It should be a list with two named entries: `k` and `weights`.
#' 
#' @author Ken Nowak
#' 
#' @references Nowak, K., Prairie, J., Rajagopalan, B., Lall, U. (2010).
#'   A nonparametric stochastic approach for multisite disaggregation of
#'   annual to daily streamflow. *Water Resources Research.*
#' 
#' @examples 
#' \dontrun{
#' # read in annual synthetic mon_flw for disag
#' x <- matrix(scan("data-raw/Meko.txt"), ncol = 2, byrow = TRUE) 
#' # intervening natural flow mon_flw - monthly CY text file
#' mon_flw <- as.matrix(read.table(
#'   "tests/dp/NFfullbasinWY0608intervening.txt", 
#'   sep = "\t"
#' ))
#' 
#' # observed annual flow for picking analog disag yr
#' ann_flw <- as.matrix(read.table("tests/dp/LFWYTotal.txt"))
#' zz <- paleo_disagg(x, ann_flw, mon_flw, 29, 1)
#' }
#' 
#' @export
paleo_disagg <- function(x, 
                         ann_flw, 
                         mon_flw, 
                         nsite = 29, 
                         sf_sites = 1:20,
                         nsim = 1,
                         ofolder = NULL, 
                         index_years = NULL,
                         k_weights = NULL)
{
  n_paleo_yrs <- nrow(x) # 1244 for meko; length of each simulation (yrs)
  
  # how many yrs of observed mon_flw
  n_obs_yrs <- nrow(mon_flw)/12 
  
  if (n_obs_yrs != nrow(ann_flw))
    stop(
      "`ann_flw` and `mon_flw` must have the same number of years.", 
      call. = FALSE
    )
  
  if (nsite != ncol(mon_flw))
    stop("`mon_flow` needs to have `nsite` columns.", call. = FALSE)
  
  if (!is.null(index_years)) {
    if (ncol(index_years) != nsim) 
      stop(
        "`index_years` must be specified for all simulations.\n",
        " So, `nsim` must equal the number of columns in `index_years`.",
        call. = FALSE
      )
    
    if (nrow(index_years) != n_paleo_yrs)
      stop(
        "`index_years` must be specified for every year in the paleo record.",
        call. = FALSE
      )
  }
  
  if (max(sf_sites) > nsite) {
    stop(
      "max(`sf_sites`), must be <= the number of sites (`site`).", 
      call. = FALSE
    )
  }
  
  if (!all(1:nsite %in% sf_sites)) {
    # set ind_sites to the remaining sites
    # ind_sites are selected directly
    ind_sites <- 1:nsite
    ind_sites <- ind_sites[!(1:nsite %in% sf_sites)]
    message(
      "Sites ", toString(ind_sites), "\n",
      "will be selected directly from the index years, i.e., not scaled."
    )
  } else {
    ind_sites <- NULL
  }
  
  # matrix for observed values - row(yr), col (month), index (site)
  dat_a <- array(data = NA, dim=c(n_obs_yrs, 12, nsite)) 

  # matrix for disag values - row(yr), col (month), index1 (site), index2(sim #)
  disag <- array(data = NA, dim = c(n_paleo_yrs, 12, nsite, nsim))
  
  # matrix for recording yr index for disag (optional)
  index_mat <- matrix(ncol = nsim, nrow = n_paleo_yrs) 
  
  # matrix for recording scale factors used (optional)
  sf_mat <- matrix(ncol = nsim, nrow = n_paleo_yrs)
  
  # this loop moves observed monthly mon_flw from 2d matrix to 3d array
  mgn <- length(dat_a[1,1,])
  
  for (j in seq_len(mgn)) {
  
    s <- 1
    e <- 12

  	for (i in seq_len(n_obs_yrs)) {
  
  	  dat_a[i, , j] <- mon_flw[s:e, j]
  
  	  s <- s + 12
  	  e <- e + 12
  	}
  }
  
  # temporary matrix for storing time only disag
  temp <- array(data = NA, dim = c(n_paleo_yrs, 12, nsim)) 
  
  # knn parameters (k) and weights
  if (is.null(k_weights)) {
    tmp <- knn_params(n_obs_yrs) #number of neighbors
    k <- tmp$k
    weights <- tmp$weights
  } else {
    k <- k_weights$k
    weights <- k_weights$weights
  }
  
  # loop through the number of simulations ---------------------
  for(j in seq_len(nsim)){
  
    # this picks the 1st year for disag based only on the annual flow
  	
  	Flow <- x[1, 2]
  	
  	# select the index year and scaling factor
  	index_atts <- get_index_sf(Flow, ann_flw, k, weights, index_years, j, 1)
  	
  	temp[1, , j] <- dat_a[index_atts$pos, , mgn] * index_atts$SF	
  	index_mat[1, j] <- index_atts$N
  	sf_mat[1, j] <- index_atts$SF
    disag[1, , sf_sites, j] <- dat_a[index_atts$pos, , sf_sites] * index_atts$SF
  	disag[1, , ind_sites, j] <- dat_a[index_atts$pos, , ind_sites]
  
    # now that one year has been disaggregated, the remaining years in the 
  	# trace use annual flow and also december of last yr (CY)
  	# *** I don't think this comment is true; not sure how it's any differnt
  	# than the first selection 
  	for(h in 2:n_paleo_yrs){
  
  		Flow <- x[h, 2]
  		# *** delete this next row? it does nothing...
  		# D <- 2:n_obs_yrs
  	
  		# select the index year and scaling factor
  		index_atts <- get_index_sf(Flow, ann_flw, k, weights, index_years, j, h)
      
  		index_mat[h, j] <- index_atts$N
  		sf_mat[h, j] <- index_atts$SF
  		
      disag[h, , sf_sites, j] <- dat_a[index_atts$pos, , sf_sites] * 
        index_atts$SF
  		disag[h, , ind_sites, j] <- dat_a[index_atts$pos, , ind_sites]
  	}
  }
  			
  # output to "flat" file
  
  disag_out <- lapply(seq_len(nsim), function(ii) {
    do.call(
      cbind, 
      lapply(seq_len(nsite), function(jj) as.vector(t(disag[,,jj,ii])))
    )
  })
  
  if (!is.null(ofolder)) {
    lapply(seq_len(nsim), function(ii) 
      utils::write.csv(
        disag_out[[ii]], 
        file = file.path(ofolder, paste0("paleo_disagg_", ii, ".csv"))
      )
    )
    utils::write.csv(index_mat, file = file.path(ofolder, "index_years.csv"))
  }

  invisible(list(paleo_disagg = disag_out, index_years = index_mat))
}


#' @param j The simulation number
#' 
get_index_sf <- function(flow, ann_flw, k, weights, index_years, j, h)
{
  # Selects a year to be "nearest neighbor"
  if (is.null(index_years)) {
    
    D <- abs(ann_flw[, 2] - flow)
    
    # combines difference and corresponding year into one matrix
    Delta <- cbind(ann_flw[, 1], D) 
    
    # reorders the delta matrix based on distances
    Delta_sort <- cbind(Delta[, 1][order(Delta[, 2])], sort(Delta[, 2])) 
    
    # selects the "k-nearest-neighbors" from Delta_sort 
    kmatrix <- Delta_sort[1:k, 1:2, drop = FALSE] 
    
    if (k != 1) {
      N <- sample(kmatrix[, 1, drop = FALSE], 1, replace = TRUE, prob = weights)
    } else {
      N <- kmatrix[, 1, drop = FALSE]
    }
  } else {
    N <- index_years[h, j]
  }
  
  # index for selected yr
  pos <- N - (ann_flw[1, 1] - 1) 
  
  # scaling factor to apply for disag
  SF <- flow/(ann_flw[pos, 2]) 

  list(pos = pos, SF = SF, N = N)
}
