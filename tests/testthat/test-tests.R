test_that("sequencebox", {
    expect_equal(align(c("GTCCG"), c("TCC", "GTCC")), c("-TCC-", "GTCC-"))
})
