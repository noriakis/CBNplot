test_that("bnpathplot: not produce errors in variou parameters", {
    data("exampleEaRes");data("exampleGeneExp");
    R <- 5
    expect_error( bnpathplot(exampleEaRes, exampleGeneExp, R = R, strThresh=0), NA)
    expect_error( bnpathplot(exampleEaRes, exampleGeneExp, R = R, edgeLink=TRUE, strThresh=0), NA)
    expect_error( bnpathplot(exampleEaRes, exampleGeneExp, disc = T, R = R, strThresh=0), NA)

    # Custom
    expect_error( bnpathplotCustom(exampleEaRes, exampleGeneExp,
    R = R, glowEdgeNum=3, hub=3, strThresh=0), NA)
})
