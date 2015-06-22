context("mutation feature vector format reading and estimation")


test_that("read mutation feature vector format data and estimation (type: indepenent, numBases: 5, trDir: TRUE)", {
  inputFile <- system.file("extdata/Hoang_MFVF.ind.txt", package="pmsignature")
  G <- readMFVFile(inputFile, numBases = 5, trDir = TRUE, type = "independent")
  Param <- getPMSignature(G, K = 3)
  for (k in 1:3) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 6))
  }
})


test_that("read mutation feature vector format data and estimation (type: full, numBases: 5, trDir: TRUE)", {
  inputFile <- system.file("extdata/Hoang_MFVF.full.txt", package="pmsignature")
  G <- readMFVFile(inputFile, numBases = 5, trDir = TRUE, type = "full")
  Param <- getPMSignature(G, K = 3)
  for (k in 1:3) {
    expect_equal(sum(getSignatureValue(Param, k)), 1)
  }
})