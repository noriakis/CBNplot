test_that("bngeneplot: not produce errors in basic usages of variou parameters", {
    data("exampleEaRes");data("exampleGeneExp");
    R <- 5
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, edgeLink=TRUE, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, hub=2, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, convertSymbol=FALSE, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, returnNet=T, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, onlyDf=T, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, showDir=T, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, strType="ms", strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, strengthPlot=TRUE, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, delZeroDegree=FALSE, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, disc=T, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R, compareRef=T, strThresh=0), NA)
    expect_error( bngeneplot(exampleEaRes, exampleGeneExp, pathNum = 1, R = R,
        sizeDep=T, dep=depmap::depmap_crispr(), strThresh=0), NA)

    ## Custom
    expect_error( bngeneplotCustom(exampleEaRes, exampleGeneExp,
        pathNum = 2, R = R, glowEdgeNum=3, hub=3, strThresh=0), NA)
    
})
