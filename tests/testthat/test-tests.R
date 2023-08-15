test_that("align", {
    {
        result <- align_sequences(c("GTCCG"), c("TCC", "GTCC", "AGTCCG"))
        expect_equal(result$sequences,  c("-TCC-", "GTCC-", "AGTCCG"))
        expect_equal(result$references, c("GTCCG", "GTCCG", "-GTCCG"))
    }

    {
        result <- align_sequences(c("GTCCG", "ABCDE"), c("TCC", "AABC"))
        expect_equal(result$sequences,  c("-TCC-", "AABC--"))
        expect_equal(result$references, c("GTCCG", "-ABCDE"))
    }

    {
        result <- align_sequences("GTCCG", c("TCC", "GTCC", "AGTC", "CCGA", "GTACCG"), "common")
        expect_equal(result$sequences,  c("--T-CC--", "-GT-CC--", "AGT-C---", "----CCGA", "-GTACCG-"))
        expect_equal(result$references, c("-GT-CCG-"))
    }
})
