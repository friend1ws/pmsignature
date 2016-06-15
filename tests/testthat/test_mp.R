context("mutation position format reading and estimation")

test_that("read mutation position format data and estimation (type: indepenent, numBases: 3, trDir: FALSE)", {
  
  inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
  G <- readMPFile(inputFile, numBases = 3, trDir = FALSE)
  BG_prob <- readBGFile(G)
  
  Param <- getPMSignature(G, K = 3)
  for (k in 1:3) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 3))
  }
  
  Param <- getPMSignature(G, K = 3, BG = BG_prob)
  for (k in 1:2) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 3))
  }
  
})

test_that("read mutation position format data and estimation (type: indepenent, numBases: 3, trDir: TRUE)", {
  
  inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
  G <- readMPFile(inputFile, numBases = 3, trDir = TRUE)
  BG_prob <- readBGFile(G)
  
  Param <- getPMSignature(G, K = 3)
  for (k in 1:3) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 4))
  }
  
  Param <- getPMSignature(G, K = 3, BG = BG_prob)
  for (k in 1:2) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 4))
  }
  
})

test_that("read mutation position format data and estimation (type: indepenent, numBases: 5, trDir: FALSE)", {
  
  inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
  G <- readMPFile(inputFile, numBases = 5, trDir = FALSE)
  BG_prob <- readBGFile(G)
  
  Param <- getPMSignature(G, K = 3)
  for (k in 1:3) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 5))
  }
  
  Param <- getPMSignature(G, K = 3, BG = BG_prob)
  for (k in 1:2) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 5))
  }
  
})

test_that("read mutation position format data and estimation (type: indepenent, numBases: 5, trDir: TRUE)", {
  
  inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
  G <- readMPFile(inputFile, numBases = 5, trDir = TRUE)
  BG_prob <- readBGFile(G)
  
  Param <- getPMSignature(G, K = 3)
  for (k in 1:3) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 6))
  }
  
  Param <- getPMSignature(G, K = 3, BG = BG_prob)
  for (k in 1:2) {
    expect_equal(rowSums(getSignatureValue(Param, k)), rep(1, 6))
  }
  
})
