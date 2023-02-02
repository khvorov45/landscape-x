test_that("sequencebox", {
    ref <- "GTCCG"
    expect_equal(align(ref, c("TCC", "GTCC")), c("-TCC-", "GTCC-"))
})
