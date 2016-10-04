library(CRSSIO)
context('Check NF file trimming function.')

dir.create('tmp')
dir.create('tmp/trace1')
dir.create('tmp/trace2')
dir.create('tmp/trace3')
dir.create('tmp/trace4')
p1 <- '..' # for automated tests
#p1 <- 'tests'
rr <- floor(runif(4,1,30)) # get 4 random nodes
message(cat('4 random nodes are:',rr))
rr <- CRSSNFInputNames()[rr]

f1 <- file.path('tmp/trace1',rr[1])
f2 <- file.path('tmp/trace2',rr[2])
f3 <- file.path('tmp/trace3',rr[3])
f4 <- file.path('tmp/trace4',rr[4])

f1.o <- file.path(p1,'trace1',rr[1])
f2.o <- file.path(p1,'trace2',rr[2])
f3.o <- file.path(p1,'trace3',rr[3])
f4.o <- file.path(p1,'trace4',rr[4])

# need to copy nf files since these functions overwrite existing files
file.copy(f1.o, f1)
file.copy(f2.o, f2)
file.copy(f3.o, f3)
file.copy(f4.o, f4)

zz <- trimSingleFile(f1, 2017, 2020)
zz <- trimSingleFile(f2, 2018, 2020)
zz <- trimSingleFile(f3, 2019, 2021)
zz <- trimSingleFile(f4, 2017, 2021)

readNF <- function(ff) as.matrix(read.table(ff, sep = '\t', skip = 2))

test_that('length of trimmed data is correct', {
  expect_equal(length(readNF(f1)),48)
  expect_equal(length(readNF(f2)),36)
  expect_equal(length(readNF(f3)),36)
  expect_equal(length(readNF(f4)),60)
})

test_that('single files match subset of untrimmed data', {
  expect_equal(readNF(f1)[,], readNF(f1.o)[1:48])
  expect_equal(readNF(f2)[,],readNF(f2.o)[13:48])
  expect_equal(readNF(f3)[,], readNF(f3.o)[25:60])
  expect_equal(readNF(f4)[,], readNF(f4.o)[,])
})

readHeader <- function(ff) scan(ff, what = 'char', nlines = 2, sep = '\t', quiet = T)

test_that('start date is correctly set in trimmed data', {
  expect_equal(readHeader(f1)[1],"start_date: 2017-1-31 24:00")
  expect_equal(readHeader(f2)[1],"start_date: 2018-1-31 24:00")
  expect_equal(readHeader(f3)[1],"start_date: 2019-1-31 24:00")
  expect_equal(readHeader(f4)[1],"start_date: 2017-1-31 24:00")
})

unlink('tmp',recursive = T) # delete the files that were created in the tests

dir.create('tmp')
file.copy(file.path(p1,'trace1'),'tmp', recursive = T)
file.copy(file.path(p1,'trace2'),'tmp', recursive = T)

test_that('all files are trimmed', {
  expect_equal(trimCCNFFiles(2018,2020,'tmp',2),58)
})

zz <- trimCCNFFiles(2018,2020,'tmp',2)

test_that('files are all trimmed as expected', {
  expect_equal(length(readNF(f1)), 36)
  expect_equal(length(readNF(f2)), 36)
})

unlink('tmp',recursive = T) # delete the files that were created in the tests
