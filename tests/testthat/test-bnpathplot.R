test_that("bnpathplot: not produce errors in variou parameters", {
    data("exampleEaRes");data("exampleGeneExp");
    R <- 5
    expect_error( bnpathplot(exampleEaRes, exampleGeneExp, R = R), NA)
    expect_error( bnpathplot(exampleEaRes, exampleGeneExp, R = R, edgeLink=TRUE), NA)
    expect_error( bnpathplot(exampleEaRes, exampleGeneExp, disc = T, R = R), NA)

    # Custom
    expect_error( bnpathplotCustom(exampleEaRes, exampleGeneExp,
    R = R, glowEdgeNum=3, hub=3), NA)
})
