
#' Spatial and temporal diaggregattion of paleo flow data
#' 
#' This is a direct resampling of LB tribs, i.e., no scale factor is applied to
#' them.
#' 
#' @param x The annual paleo data to disaagregate.
#' @param ann_flw Observed annual flow data used for picking analog year.
#' @param mon_flw Intervening monthly natural flow. Used for spatially and 
#'   temporaly disaggregating paleo data based on the index year.
#' @param nsite The number of sites to diaggregate the data too. 
#' @param nsim Number of times to repeat the space/time disaggregation.
#' 
#' @examples 
#' \dontrun{
#' # read in annual synthetic mon_flw for disag
#' x <- matrix(scan("data-raw/Meko.txt"), ncol = 2, byrow = TRUE) 
#' # intervening natural flow mon_flw - monthly CY text file
#' mon_flw <- matrix(
#'   scan("data-raw/CRB_CY_MON_IV_29_06_08.txt"), 
#'   ncol = 29, 
#'   byrow = TRUE
#' )
#' 
#' # observed annual flow for picking analog disag yr
#' ann_flw <- matrix(scan("data-raw/LF_06_08.txt"), ncol = 2, byrow = TRUE)
#' zz <- paleo_disagg(x, ann_flw, mon_flw, 29, 1)
#' }
#' 
#' @export
paleo_disagg <- function(x, ann_flw, mon_flw, nsite, nsim, ofolder = NULL)
{
  n_paleo_yrs <- nrow(x) # 1244 for meko; length of each simulation (yrs)
  
  # how many yrs of observed mon_flw
  n_obs_yrs <- nrow(mon_flw)/12 
  
  if (n_obs_yrs != nrow(ann_flw))
    stop("`ann_flw` and `mon_pttrn` must have the same number of years.")
  
  if (nsite != ncol(mon_flw))
    stop("`mon_flow` needs to have `nsite` columns.")
  
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
  
  for(j in 1:mgn){
  
    s <- 1
    e <- 12

  	for(i in 1:n_obs_yrs){
  
  	  dat_a[i,,j] <- mon_flw[s:e,j]
  
  	  s <- s + 12
  	  e <- e + 12
  	}
  }
  
  # temporary matrix for storing time only disag
  temp <- array(data = NA, dim = c(n_paleo_yrs, 12, nsim)) 
  
  for(j in 1:nsim){
  
    #this picks the 1st year for disag based only on the annual flow
  	k <- sqrt(n_obs_yrs) #number of neighbors
  	
  	Flow <- x[1, 2]
  	
  	D <- abs(ann_flw[, 2] - Flow)
  	
  	# combines difference and corresponding year into one matrix
  	Delta <- cbind(ann_flw[, 1], D) 
  	
  	# reorders the delta matrix based on distances
  	Delta_sort <- cbind(Delta[, 1][order(Delta[, 2])], sort(Delta[, 2])) 
  	
  	# selects the "k-nearest-neighbors" from Delta_sort 
  	kmatrix <- Delta_sort[1:k, 1:2] 
  	
  	# defines matrix for weights
  	weight <- matrix(nrow = k, ncol = 1) 
   	
  	# ranks distances for purpose of generating weights	
  	rnk <- rank(kmatrix[, 2]) 
  		
  		for(i in 1:k){
  		  # fills weighting matrix
  			weight[i, 1] <- 1/(rnk[i]) 
  	
  		}
  
  	z <- sum(weight) # sums weights 
  
  	# divides weights by sum of weights so cumulative probability = 1
  	weights <- weight/z	
  	
  	# Selects a year to be "nearest neighbor"
  	N <- sample(kmatrix[, 1], 1, replace = TRUE, prob=weights) 
  
  	# index for selected yr
  	pos <- N - (ann_flw[1, 1] - 1) 
    
  	# scaling factor to apply for disag
  	SF <- Flow/(ann_flw[pos, 2]) 
  	temp[1, , j] <- dat_a[pos, , mgn] * SF	
  	index_mat[1, j] <- N
  	sf_mat[1, j] <- SF
    disag[1, , 1:20, j] <- dat_a[pos, , 1:20]*SF
  	disag[1, , 21:29, j] <- dat_a[pos, , 21:29]
  
    # now that one year has been disaggregated, the remaining years in the 
  	# trace use annual flow and also december of last yr (CY)
  	for(h in 2:n_paleo_yrs){
  
  		k <- sqrt(n_obs_yrs) #number of neighbors
  
  		Flow <- x[h, 2]
  		D <- 2:n_obs_yrs
  	
  		# annual as the only selection criteria
      D <- abs(ann_flw[,2] - Flow)
      
      # combines difference and corresponding year into one matrix 
      # these use just m.a.f
  		Delta <- cbind(ann_flw[,1],D) 
  
  		# reorders the delta matrix based on distances
  		Delta_sort <- cbind(Delta[,1][order(Delta[,2])], sort(Delta[,2])) 
  		
  		# selects the "k-nearest-neighbors" from Delta_sort
  		kmatrix <- Delta_sort[1:k,1:2]  
  	  weight <- matrix(nrow=k, ncol=1) # defines matrix for weights
   		
  	  # ranks distances for purpose of generating weights
  		rnk <- rank(kmatrix[,2]) 
  		
  		for(i in 1:k){
  	
  			weight[i,1] <- 1/(rnk[i]) #fills weighting matrix
  	
  		}
  
  		z <- sum(weight) # sums weights 
  	
  		#divides weights by sum of weights so cumulative probability = 1
  		weights <- weight/z	
  		
  		# Selects a year to be "nearest neighbor"
  		N <- sample(kmatrix[, 1], 1, replace = TRUE, prob = weights) 
  		pos <- N - (ann_flw[1, 1] - 1) # index for selected yr
  		SF <- Flow/(ann_flw[pos, 2]) # scaling factor to apply for disag
      
  		index_mat[h, j] <- N
  		sf_mat[h, j] <- SF
  		
  		#disag[h, , ,j]=dat_a[pos, ,]*SF
      disag[h, , 1:20, j] <- dat_a[pos, , 1:20]*SF
  		disag[h, , 21:29, j] <- dat_a[pos, , 21:29]
  	}
  }
  			
  # output to "flat" file
  
  #disagmat=matrix(ncol=nsite, nrow=(n_paleo_yrs*12*nsim))
  
  disag_out <- lapply(seq_len(nsim), function(ii) {
    do.call(
      cbind, 
      lapply(seq_len(nsite), function(jj) as.vector(t(disag[,,jj,ii])))
    )
  })
  
  # for(i in 1:nsite){
  #   p <- 1
  # 
  #   for(k in 1:nsim){
  #     for(j in 1:n_paleo_yrs){
  #       disagmat[p:(p+11),i] <- disag[j,,i,k]
  #       p <- p + 12
  #     }
  #   }
  # }
  
  if (!is.null(ofolder)) {
    lapply(seq_len(nsim), function(ii) 
      write.csv(
        disag_out[[ii]], 
        file = file.path(ofolder, paste0("paleo_disagg_", ii, ".csv"))
      )
    )
    write.csv(index_mat, file = file.path(ofolder, "index_years.csv"))
  }

  invisible(list(paleo_disagg = disag_out, index_yrs = index_mat))
}
