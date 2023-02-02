test_that("sequencebox", {
    ref <- "GTCCG"
    expect_equal(align(ref, c("TCC", "GTC")), c("-TCC-", "GTC--"))
})
