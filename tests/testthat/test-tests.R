test_that("align", {
    result <- align(c("GTCCG"), c("TCC", "GTCC", "AGTCCG"))
    expect_equal(result$sequences, c("-TCC-", "GTCC-", "AGTCCG"))
    expect_equal(result$references, c("GTCCG", "GTCCG", "-GTCCG"))
})
