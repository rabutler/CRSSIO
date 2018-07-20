context("check paleo space time disagg code")

# compare to previous results ----------------------------

x <- matrix(scan("../dp/Meko.txt", quiet = TRUE), ncol = 2, byrow = TRUE) 
# intervening natural flow mon_flw - monthly WY text file
mon_flw <- as.matrix(read.table(
  "../dp/NFfullbasinWY0608intervening.txt", 
  sep = "\t"
))

# observed annual flow for picking analog disag yr
ann_flw <- as.matrix(read.table("../dp/LFWYTotal.txt"))

# ** can I check old CRSS packages to see direct paleo files when these data 
# were used. Is the mon_flw file correct for the 1906-2008 natural flows?

# ** this contains weird numbers for the Grand Canyon reach
# zz <- as.matrix(read.table("../dp/MatrixSimDataCRBwithObsLB_DP.txt"))
zz <- as.matrix(read.csv(
  "../dp/MatrixSimDataCRBwithObsLB_DP_rab20180620.csv"
))
 
index_yrs <- matrix(scan("../dp/indexpick.txt", quiet = TRUE), ncol = 1)

# paleo_disagg <- function(x, 
#                          ann_flw, 
#                          mon_flw, 
#                          nsite = 29, 
#                          sf_sites = 1:20,
#                          nsim = 1,
#                          ofolder = NULL, 
#                          index_years = NULL,
#                          k_weights = NULL)

test_that("disagg matches previous code's results", {
  expect_equivalent(
    tmp <- paleo_disagg(
      x, 
      ann_flw = ann_flw, 
      mon_flw = mon_flw, 
      index_years = index_yrs)$paleo_disagg[[1]],
    zz,
    tolerance = 0.00001
  )
  expect_equivalent(round(tmp, 0), round(zz, 0))
  expect_equal(range(tmp - zz), c(0, 0))
})

# compare random selection -----------------------------

orig_index <- as.matrix(read.csv("../dp/index_years_rseed408.csv"))
dimnames(orig_index) <- NULL
set.seed(403) # this was the first entry of .Random.seed when implementing this

test_that("current random selection matches original random selection", {
  expect_equal(paleo_disagg(x, ann_flw, mon_flw)$index_years, orig_index)
  set.seed(403)
  expect_equal(knn_get_index(x, ann_flw), orig_index)
})

# ***** still need to make function much safer to the format of incoming data,
# i.e., which input need years associated with them, and which don't, matrices, 
# vs. vectors, etc.
# test for k = 1 and weights = 1
# Should also consider round to nearest AF, but what are the effects of that on 
# matching the inut Lees Ferry value

# should error if index_years and k_weights are specified by user
