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


test_that("read mutation position format data and estimation (type: full, numBases: 5, trDir: TRUE)", {
  
  inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
  G <- readMPFile(inputFile, numBases = 5, type = "full", trDir = TRUE)
  BG_prob <- readBGFile(G)
  
  Param <- getPMSignature(G, K = 3)
  for (k in 1:3) {
    expect_equal(sum(getSignatureValue(Param, k)), 1)
  }
  
  Param <- getPMSignature(G, K = 3, BG = BG_prob)
  for (k in 1:2) {
    expect_equal(sum(getSignatureValue(Param, k)), 1)
  }
  
})


test_that("read mutation position format data and estimation for hg18 (type: independent, numBases: 5, trDir: TRUE)", {
  
  inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.hg18.txt.gz", package="pmsignature")
  G <- readMPFile(inputFile, numBases = 5, trDir = TRUE, 
                  bs_genome = BSgenome.Hsapiens.UCSC.hg18::BSgenome.Hsapiens.UCSC.hg18,
                  txdb_transcript = TxDb.Hsapiens.UCSC.hg18.knownGene::TxDb.Hsapiens.UCSC.hg18.knownGene)
  
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
