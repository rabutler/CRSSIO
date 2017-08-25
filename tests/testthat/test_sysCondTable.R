library(RWDataPlyr)
library(dplyr)

context('Check system condition table creation')

slotAggList <- RWDataPlyr::createSlotAggList(CRSSIO::sysCondSALMatrix())
scenFolder <- 'ISM1988_2014,2007Dems,IG,Most'
scenName <- 'scen1'
scenPath <- system.file('extdata','Scenario/',package = 'RWDataPlyr')
sysData <- RWDataPlyr::getDataForAllScens(scenFolder, scenName, slotAggList,
                                          scenPath, 'tmp.feather', TRUE)
on.exit(file.remove("tmp.feather"))

yrs <- 2018:2022
sysCondTable <- createSysCondTable(sysData, yrs)

test_that('object dimensions and attributes are correct', {
  expect_equal(length(sysCondTable), 2)
  expect_equal(names(sysCondTable), c("fullTable", "limitedTable"))
  expect_equal(dim(sysCondTable$fullTable), c(length(CRSSIO:::slotNames())+4,length(yrs)))
  expect_equal(dim(sysCondTable$limitedTable), c(length(CRSSIO:::slotNames())+1,length(yrs)))
})

s2 <- sysData %>%
  mutate(Scenario = "scen2") %>%
  bind_rows(sysData)

test_that("warnings and errors are as expected", {
  expect_warning(createSysCondTable(s2, yrs),
                 'There are 2 Scenarios in the data. Please note, these scenarios will be averaged together when creating the system conditions table.'
  )
})

expVals <- matrix(c(
  rep(0, 5),
  rep(0, 5),
  rep(0, 5),
  c(100, 100, 25, 25, 50),
  c(100, 100, 0, 0, 50),
  c(0, 0, 25, 25, 0), 
  rep(0, 5),
  c(0, 0, 75, 75, 25),
  rep(0, 5),
  c(0, 0, 75, 75, 25),
  rep(0, 4), 25,
  rep(0, 5), 
  rep(0, 5),
  rep(0, 4), 25,
  c(0, 50, 75, 75, 100),
  c(0, 50, 75, 25, 50),
  c(0, 0, 0, 50, 50),
  rep(0, 5),
  rep(0, 5),
  rep(0, 5),
  c(100, 50, 25, 25, 0)),
  ncol = 5, 
  byrow = T
)



# thinkt hat checking against the orig values, and then selecting a few values
# to recompute should test it enough
# check UEB Total, Shortage 1, MER 7.48 and normal year
r1 <- apply(
  rdfSlotToMatrix(RWDataPlyr::sysRdf, "SummaryOutputData.UpperBalancingAbove823") +
  rdfSlotToMatrix(RWDataPlyr::sysRdf, "SummaryOutputData.UpperBalancingAt823") +
  rdfSlotToMatrix(RWDataPlyr::sysRdf, "SummaryOutputData.UpperBalancingBelow823"),
  1,
  mean
) * 100
r2 <- apply(
  rdfSlotToMatrix(RWDataPlyr::sysRdf, "SummaryOutputData.LBShortageStep1"),
  1,
  mean
) * 100
r3 <- apply(
  rdfSlotToMatrix(RWDataPlyr::sysRdf, "SummaryOutputData.MidElevationReleaseAt748"),
  1,
  mean
) * 100
r4 <- apply(
  rdfSlotToMatrix(RWDataPlyr::sysRdf, "SummaryOutputData.LBNormalCondition"),
  1,
  mean
) * 100
test_that("computations of chances are correct", {
  expect_equivalent(expVals, sysCondTable$fullTable)
  expect_equivalent(sysCondTable$fullTable[4,], r1)
  expect_equivalent(sysCondTable$fullTable[16,], r2)
  expect_equivalent(sysCondTable$fullTable[10,], r3)
  expect_equivalent(sysCondTable$fullTable[21,], r4)
})

test_that("rows sum together correctly", {
  expect_equal(sysCondTable$fullTable[1,], apply(sysCondTable$fullTable[2:3,], 2, sum))
  expect_equal(sysCondTable$fullTable[4,], apply(sysCondTable$fullTable[5:7,], 2, sum))
  expect_equal(sysCondTable$fullTable[8,], apply(sysCondTable$fullTable[9:10,], 2, sum))
  expect_equal(sysCondTable$fullTable[11,], apply(sysCondTable$fullTable[12:14,], 2, sum))
  expect_equal(sysCondTable$fullTable[15,], apply(sysCondTable$fullTable[16:18,], 2, sum))
  expect_equivalent(rep(100, 5), apply(sysCondTable$fullTable[c(1,4,8,11),], 2, sum))
  expect_equivalent(rep(100, 5), apply(sysCondTable$fullTable[c(15,19,21),], 2, sum))
  expect_true(all(sysCondTable$fullTable[20,] <= sysCondTable$fullTable[19,]))
})