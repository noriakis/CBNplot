test_that("utilities: not produce errors in variou parameters", {
    data("exampleEaRes");data("exampleGeneExp");
    library(bnlearn)
    data("exampleEaRes")
    data("exampleGeneExp")
    net <- bngeneplot(exampleEaRes, exampleGeneExp, R=5,
                   pathNum=1, returnNet=TRUE)
    fitted <- bn.fit(net$av, net$df)
    res <- queryCpDistLw(fitted, candidate="ERCC4", evidence="ERCC2",
                          levels=c(0.1, 0.5, 0.8), n=50, alpha=TRUE)
    res <- queryCpDistLw(fitted, candidate="ERCC4", evidence="ERCC2",
                          levels=c(0.1, 0.5, 0.8), n=50, alpha=FALSE)
    res <- queryCpDistLw(fitted, candidate="ERCC4", evidence="ERCC2",
                          levels=c(0.1, 0.5, 0.8), n=50, point=TRUE)
    res <- queryCpDistLw(fitted, candidate="ERCC4", evidence="ERCC2",
                          levels=c(0.1, 0.5, 0.8), n=50, point=FALSE)
})
