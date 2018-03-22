context("check space time disagg code")

# compare to previous results ----------------------------

x <- matrix(scan("../dp/Meko.txt"), ncol = 2, byrow = TRUE) 
# intervening natural flow mon_flw - monthly WY text file
mon_flw <- as.matrix(read.table(
  "../dp/NFfullbasinWY0608intervening.txt", 
  sep = "\t"
))

# observed annual flow for picking analog disag yr
ann_flw <- as.matrix(read.table("../dp/LFWYTotal.txt"))

zz <- as.matrix(read.table("../dp/MatrixSimDataCRBwithObsLB_DP.txt"))
 
index_yrs <- matrix(scan("../dp/indexpick.txt"), ncol = 1)

test_that("disagg matches previous code's results", {
  expect_equivalent(
    paleo_disagg(x, ann_flw, mon_flw, 29, 1, index_years = index_yrs)$paleo_disagg[[1]],
    zz,
    tolerance = 1
  )
})
