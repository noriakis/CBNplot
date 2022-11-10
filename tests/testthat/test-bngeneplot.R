test_that("bngeneplot: not produce errors in basic usages of variou parameters", {
    data("exampleEaRes");data("exampleGeneExp");
    R <- 5
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, hub=2), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, convertSymbol=FALSE), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, returnNet=T), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, onlyDf=T), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, showDir=T), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, strType="ms"), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1:2, R = R), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, strengthPlot=TRUE), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, delZeroDegree=FALSE), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, disc=T), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, compareRef=T), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R,
        sizeDep=T, dep=depmap::depmap_crispr()), NA)

    ## Custom
    expect_error( bngeneplotCustom(exampleEaRes, exampleGeneExp,
        pathNum = 2, R = R, glowEdgeNum=3, hub=3), NA)
    
})
