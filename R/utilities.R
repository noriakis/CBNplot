#' queryCpDistLs
#'
#' produce a plot of bnlearn::cpdist by performing bnlearn::cpdist
#' on specified node, evidence and level.
#'
#' @param fitted bn.fit object
#' @param candidate name of node
#' @param evidences the evidences
#' @param discPalette palette to be used for plotting if the event is discrete 
#' @param ... other parameters passed to bnlearn cpdist
#'
#' @importFrom dplyr mutate
#' @return list of dataframe containing raw values
#' @examples
#' library(bnlearn)
#' data("exampleEaRes")
#' data("exampleGeneExp")
#' net <- bngeneplot(exampleEaRes, exampleGeneExp,
#'                   pathNum=1, returnNet=TRUE)
#' fitted <- bn.fit(net$av, net$df)
#' res <- queryCpDistLs(fitted, candidate="ERCC4",
#'               evidences=c("ERCC2<0.1","ERCC2>0.5","ERCC2>0.8"), n=500)
#' @export
#'
queryCpDistLs <- function (fitted, candidate,
                            evidences, discPalette="Set2",...) {
    evRenamed <- c()
    for (e in evidences){
        if (grepl("<=",e,fixed = TRUE)) {
            tmp <- unlist(strsplit(e,"<="));
            newe <- paste0("`",tmp[1],"`<= ",tmp[2]);
            evRenamed <- c(evRenamed, newe)}
        else if (grepl(">=",e,fixed = TRUE)) {
            tmp <- unlist(strsplit(e,">="));
            newe <- paste0("`",tmp[1],"`>= ",tmp[2]);
            evRenamed <- c(evRenamed, newe)}
        else if (grepl("<",e,fixed = TRUE)) {
            tmp <- unlist(strsplit(e,"<"));
            newe <- paste0("`",tmp[1],"`< ",tmp[2]);
            evRenamed <- c(evRenamed, newe)}
        else if (grepl(">",e,fixed = TRUE)) {
            tmp <- unlist(strsplit(e,">"));
            newe <- paste0("`",tmp[1],"`> ",tmp[2]);
            evRenamed <- c(evRenamed, newe)}
        else if (grepl("==",e,fixed = TRUE)) {
            tmp <- unlist(strsplit(e,"=="));
            newe <- paste0("`",tmp[1],"`== ",tmp[2]);
            evRenamed <- c(evRenamed, newe)}
    }

    ## If event is discrete
    if ("prob" %in% names(fitted[[candidate]])) {
        conc <- c()
        for (e in evRenamed){
            tmp <- eval(parse(text=paste0("cpdist(fitted, nodes='",
                candidate,"', evidence=",e,", method='ls', ...)")))
            conc <- rbind(conc, data.frame(tmp, rep(e, nrow(tmp))))
        }
        colnames(conc) <- c(candidate, "level")
        returnList <- list()
        returnList[["df"]] <- conc
        groupedDf <- conc %>%
                    group_by(.data$level, !!!syms(candidate)) %>%
                    summarise(n = n()) %>%
                    mutate(freq = n / sum(n))
        groupedPlot <- groupedDf %>% 
                    ggplot(aes_(~level,~freq,
                        fill=candidate))+
                    geom_bar(stat="identity",position='dodge')+
                    theme_bw()+scale_fill_brewer(palette=discPalette)

        returnList[["grouped"]] <- groupedDf
        returnList[["plot"]] <- groupedPlot
        return(returnList)
    }

    ## Otherwise
    conc <- c()
    for (e in evRenamed){
        tmp <- eval(parse(text=paste0("cpdist(fitted, nodes='",
            candidate,"', evidence=",e,", method='ls', ...)")))
        conc <- rbind(conc, data.frame(tmp, rep(e, nrow(tmp))))
    }
    colnames(conc) <- c(candidate, "level")
    
    returnList <- list()
    returnList[["df"]] <- conc

    if (length(candidate)==1){
        disMean <- conc %>% group_by(.data$level) %>%
                    summarise(mean=mean(get(candidate)), n=n())
        disMean$label <- paste0(disMean$level,
            " (mean=", round(disMean$mean,3), ")")
        lb <- c()
        for (nm in unique(conc$level)) {
            lab <- as.character(disMean %>%
                filter(.data$level==nm) %>% select(.data$label))
            num <- as.numeric(disMean %>%
                filter(.data$level==nm) %>% select(n))
            lb <- c(lb, rep(lab, num))
        }
        conc$label <- lb
        p <- conc %>%
            ggplot(aes_(y=~level, x=candidate,
                color=~label, fill=~label)) +
            ggdist::stat_dotsinterval()+
            theme_bw() + xlab(candidate)
        returnList[["plot"]] <- p
    }
    return(returnList)
}


#' queryCpDistLw
#'
#' produce a plot of bnlearn::cpdist by performing bnlearn::cpdist
#' on specified node, evidence and level.
#'
#'
#' @param fitted bn.fit object
#' @param candidate name of node
#' @param evidence evidence variable name
#' @param levels level to be listed
#' @param alpha whether to reflect the weights by alpha (TRUE) or color (FALSE)
#' @param point geom_point the weighted mean
#' @param pointSize point size for geom_point
#' @param ... other parameters passed to bnlearn cpdist
#' @importFrom dplyr summarise
#' @importFrom ggdist geom_dots
#' @importFrom ggdist scale_fill_ramp_continuous
#' @importFrom bnlearn cpdist
#' @return list of dataframe containing raw values
#' @examples
#' library(bnlearn)
#' data("exampleEaRes")
#' data("exampleGeneExp")
#' net <- bngeneplot(exampleEaRes, exampleGeneExp,
#'                   pathNum=1, returnNet=TRUE)
#' fitted <- bn.fit(net$av, net$df)
#' res <- queryCpDistLw(fitted, candidate="ERCC4", evidence="ERCC2",
#'                      levels=c(0.1, 0.5, 0.8), n=500)
#' @export
#'
queryCpDistLw <- function (fitted, candidate, evidence, levels, point=FALSE,
                            pointSize=5, alpha=TRUE, ...) {

    if ("prob" %in% names(fitted[[candidate]])) {
        stop("the cpdist(lw) for the categorical 
            event node is currently not supported.")}

    conc <- c()
    if(length(evidence)!=1){stop("evidence must be one node")}
    for (l in levels){
        ll <- list()
        ll[[evidence]] <- l
        tmp <- cpdist(fitted, nodes=c(candidate),
            evidence=ll, method="lw", ...)
        conc <- rbind(conc, data.frame(tmp,
            attributes(tmp)$weights, rep(l, nrow(tmp))))
    }
    
    colnames(conc) <- c(candidate, "weights", "level")
    conc$level <- factor(conc$level)
    wm <- conc %>% group_by(.data$level) %>%
            dplyr::summarise_at(candidate, ~ weighted.mean(., weights))
    returnList <- list()
    returnList[["df"]] <- conc
    returnList[["wm"]] <- wm
    if (length(candidate)==1){
        ## Points are colored | set alpha by their weights
        if (alpha){
            p <- conc %>% ggplot()+
                ggdist::geom_dots(aes_(y=~level,
                    x=candidate, fill=~level, alpha=~weights), color=NA)+
                    theme_bw()+xlab(candidate)+ylab(evidence)
            if (point) {p <- p + geom_point(wm,
                mapping=aes_(x=candidate, y=~level, fill=~level),
                shape=21, size=pointSize) }
        } else {
            p <- conc %>% ggplot()+
                ggdist::geom_dots(aes_(y=~level,
                    x=candidate, fill=~weights), color=NA)+
                    scale_fill_ramp_continuous()+
                    theme_bw()+xlab(candidate)+ylab(evidence)
            if (point) {p <- p + geom_point(wm,
                mapping=aes_(x=candidate, y=~level), color="tomato",
                size=pointSize) }
        }
        p 
        returnList[["plot"]] <- p
    }
    return(returnList)
}


#' inferMS
#'
#' multiscale bootstrap-based inference of Bayesian network
#'
#' @param data data.frame to perform inference
#' @param algo structure learning method used in boot.strength()
#' @param algorithm.args parameters to pass to
#'                       bnlearn structure learnng function
#' @param R the number of bootstrap
#' @param cl cluster object from parallel::makeCluster()
#' @param r vector for size of each bootstrap replicate
#' @importFrom utils packageVersion
#' @importFrom bnlearn inclusion.threshold
#' @importFrom pvclust msfit
#' @importFrom purrr reduce
#' @return object of class bn.strength
#'
#'
inferMS <- function(data, algo, algorithm.args, R, cl=NULL,
                    r=seq(0.5, 1.5, 0.1))
{
    #stop("Currently cannot be used until the release of bnlearn 4.7")
    #threshold <- utils::getFromNamespace("threshold", "bnlearn")
    nList <- list()
    for (s in r){
        sb <- boot.strength(data=data, cluster=cl, algorithm=algo, R=R,
                            m=as.integer(nrow(data)*s),
                            algorithm.args=algorithm.args)
        nList[[paste0("s",s)]] <- sb
    }
    md <- nList %>% reduce(dplyr::left_join, by = c("from","to"))

    stcol <- c("from", "to", colnames(md)[grepl("strength",colnames(md))])
    dicol <- c("from", "to", colnames(md)[grepl("direc",colnames(md))])

    st <- md[,stcol]
    ms <- c()
    for (i in seq_len(nrow(st))){
        tmpbp <- as.numeric(st[i,][3:length(st[i,])])
        ft <- msfit(tmpbp, r, R)
        ms <- c(ms, as.numeric(1-pnorm(ft$coef[1]-ft$coef[2])))
    }
    st$MS <- ms

    di <- md[,dicol]
    dis <- c()
    for (i in seq_len(nrow(di))){
        tmpbp <- as.numeric(di[i,][3:length(di[i,])])
        ft <- msfit(tmpbp, r, R)
        dis <- c(dis, as.numeric(1-pnorm(ft$coef[1]-ft$coef[2])))
    }
    di$DIR <- dis
    nn <- unique(c(di$from, di$to))
    res <- data.frame(di$from, di$to, st$MS, di$DIR)
    colnames(res) <- c("from","to","strength","direction")

    res <- structure(res, method = "bootstrap", threshold = 0,
                    nodes = nn,
                    class = c("bn.strength", class(res)))

    ## Using bnlearn version 4.7.20210803 or
    ## above for function inclusion.threshold().
    if (packageVersion("bnlearn")=="4.7"){
        attributes(res)$threshold <- inclusion.threshold(res)
    }
    return(res)
}

#' discDF
#' 
#' Discretize by arules::discretize, by the k-means clustering
#'
#' @noRd
discDF <- function(df, tr=NULL, remainCont=NULL) {
    if (is.null(tr)) {
        discMeth <- list()
        for ( i in names(df)){
            if (i %in% remainCont){
                discMeth[[i]] <- list(method="none")
            } else {
                discMeth[[i]] <- list(method="cluster")
            }
        }
        } else {
            discMeth <- tr
        }
    dData <- arules::discretizeDF(df=df, methods=discMeth)
    return(dData)
}


# #' chooseEdgeDir (deprecated)
# #' 
# #' Recursively perform choose.direction() from bnlearn on undirected arcs.
# #'
# #' @noRd
# chooseEdgeDir <- function(av, pcs, scoreType, debug=FALSE) {
#     if (dim(undirected.arcs(av))[1]!=0){
#         message("found undirected arc(s)")
#         undir <- undirected.arcs(av)
#         for (n in seq_len(dim(undir)[1])) {
#             av <- choose.direction(av, c(undir[,1][n], undir[,2][n]),
#                                pcs, criterion = scoreType, debug = debug)
#         }
#     }
#     return(av)
# }

#' returnReactomeIntersection
#'
#' @noRd
returnReactomeIntersection <- function(res, g) {
    url <- "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
    bfc <- BiocFileCache::BiocFileCache()
    path <- BiocFileCache::bfcrpath(bfc, url)
    pathRel <- read.csv(path, sep="\t", header=FALSE)
    pathRelG <- graph_from_edgelist(as.matrix(pathRel), directed = FALSE)
    incPath <-  V(pathRelG)[names(V(pathRelG)) %in% res$ID]
    incPathG <- igraph::subgraph(pathRelG, incPath)
    incPathG <- set.vertex.attribute(incPathG, "name",
        value=res[names(V(incPathG)), "Description"])
    undirG <- graph_from_edgelist(as_edgelist(g), directed = FALSE)
    ovlG <- as_edgelist(igraph::intersection(incPathG, undirG))
    ovlEdge <- rbind(cbind(ovlG[,1], ovlG[,2]), cbind(ovlG[,2], ovlG[,1]))
    ovlEMat <- c()
    for (i in seq_len(dim(ovlEdge)[1])){
        ovlEMat <- c(ovlEMat, paste(ovlEdge[i,], collapse="_"))
    }
    edgeLabel <- c()
    eList <- as_edgelist(g)
    for (i in seq_len(dim(eList)[1])){
        if (paste(eList[i,], collapse="_") %in% ovlEMat){
            edgeLabel <- c(edgeLabel, "solid")
        } else {
            edgeLabel <- c(edgeLabel, "dotted")
        }
    }
    return(edgeLabel)
}


#' obtainPath
#' 
#' obtain the analysis results including the queried gene symbol
#' 
#' @param res enrichment analysis result
#' @param geneSymbol the candidate gene
#'
#' @importFrom purrr map
#' @return subset of enrichment results
#' @examples
#' data("exampleEaRes")
#' obtainPath(res = exampleEaRes, geneSymbol="ERCC7")
#' @export

obtainPath <- function(res, geneSymbol) {
    # if (res@keytype=="kegg"){tot <- "ENTREZID"} else {tot <- res@keytype}
    # candidateGene <- clusterProfiler::bitr(geneSymbol, fromType = "SYMBOL",
    # toType = tot, org.Hs.eg.db)[tot]
    res@result <- res@result[unlist(map(map(res@result$geneID,
                            function (x) unlist(strsplit(x, "/"))),
                            function(x) geneSymbol %in% x)),]
    return(res)
}

#' compareBNs
#' 
#' Take the list of networks and returns the F-measures
#' 
#' @param listOfNets list of networks
#' @importFrom utils combn
#' @return F-measures of each combination of network
#' @examples
#' data("exampleEaRes");data("exampleGeneExp")
#' net1 <- bngeneplot(results = exampleEaRes,
#'         exp = exampleGeneExp, pathNum = 1, R = 10, returnNet=TRUE)
#' net2 <- bngeneplot(results = exampleEaRes,
#'         exp = exampleGeneExp, pathNum = 1, R = 10, returnNet=TRUE)
#' res <- compareBNs(list(net1$av, net2$av))
#' @export
compareBNs <- function(listOfNets){
    if (length(listOfNets)<2){
        stop("please provide multiple networks from bnlearn")
    }
    fms <- c()
    cmb <- combn(seq_len(length(listOfNets)), m=2)
    for (i in seq_len(ncol(cmb))){
        compareRes <- bnlearn::compare(listOfNets[[cmb[,i][1]]],
            listOfNets[[cmb[,i][2]]])
        prec <- compareRes$tp/(compareRes$fp+compareRes$tp)
        recall <- compareRes$tp/(compareRes$fn+compareRes$tp)
        fm <- 2 * prec * recall / (recall + prec)
        fms <- c(fms, fm)
    }
    return(fms)
}

#' bnpathtest
#'
#' Testing various R for bayesian network between pathways
#'
#' @param results the enrichment analysis result
#' @param exp gene expression matrix
#' @param expRow the type of the identifier of rows of expression matrix
#' @param expSample candidate rows to be included in the inference
#'                  default to all
#' @param algo structure learning method used in boot.strength()
#'             default to "hc"
#' @param algorithm.args parameters to pass to
#'                       bnlearn structure learnng function
#' @param cl cluster object from parallel::makeCluster()
#' @param qvalueCutOff the cutoff value for qvalue
#' @param adjpCutOff the cutoff value for adjusted pvalues
#' @param nCategory the number of pathways to be included
#' @param Rrange the sequence of R values to be tested
#' @param scoreType return the specified scores
#' @param orgDb perform clusterProfiler::setReadable
#'              based on this organism database
#' @param bypassConverting bypass symbol converting
#' @return list of graphs and scores
#' @examples
#' data("exampleEaRes");data("exampleGeneExp")
#' res <- bnpathtest(results = exampleEaRes, exp = exampleGeneExp,
#'        algo="hc", Rrange=seq(10, 30, 10), expRow = "ENSEMBL",
#'        scoreType="bge")
#' @export
#'
bnpathtest <- function (results, exp, expSample=NULL, algo="hc",
                        algorithm.args=NULL, expRow="ENSEMBL", cl=NULL,
                        orgDb=org.Hs.eg.db, bypassConverting=FALSE,
                        qvalueCutOff=0.05, adjpCutOff=0.05, nCategory=15,
                        Rrange=seq(2,40,2), scoreType="aic-g")
{
    if (!bypassConverting) {
        if (!is.null(orgDb)){
            results <- clusterProfiler::setReadable(results, OrgDb=orgDb)
        }
    }
    if (is.null(expSample)) {expSample <- colnames(exp)}
    # if (results@keytype == "kegg"){
    #     resultsGeneType <- "ENTREZID"
    # } else {
    #     resultsGeneType <- results@keytype
    # }

    res <- results@result
    if (!is.null(qvalueCutOff)) {
        res <- subset(res, res$qvalue < qvalueCutOff) }
    if (!is.null(adjpCutOff)) {
        res <- subset(res, res$p.adjust < adjpCutOff) }
    if (nCategory) {
        res <- res[seq_len(nCategory),]
        res <- res[!is.na(res$ID),]
    }


    pcs <- c()
    pwayNames <- c()
    for (i in seq_len(length(rownames(res)))) {
        genesInPathway <- strsplit(res[i, ]$geneID, "/")[[1]]
        if (!bypassConverting) {
            if (!is.null(orgDb)){
                genesInPathway <- clusterProfiler::bitr(genesInPathway,
                                                        fromType="SYMBOL",
                                                        toType=expRow,
                                                        OrgDb=orgDb)[expRow][,1]
            }
        }
        pathwayMatrix <- exp[ intersect(rownames(exp), genesInPathway),
                            expSample ]
        if (dim(pathwayMatrix)[1]==0) {
            message("no gene in the pathway present in expression data")
        } else {
            pathwayMatrixPca <- prcomp(t(pathwayMatrix), scale. = FALSE)$x[,1]
            avExp <- apply(pathwayMatrix, 2, mean)
            corFlag <- cor(pathwayMatrixPca, avExp)
            if (corFlag < 0){pathwayMatrixPca <- pathwayMatrixPca*-1}
            # pathwayMatrixSum <- apply(pathwayMatrix, 2, sum)
            pwayNames <- c(pwayNames, res[i,]$Description)
            pcs <- cbind(pcs, pathwayMatrixPca)
        }
    }
    colnames(pcs) <- pwayNames
    pcs <- data.frame(pcs, check.names=FALSE)

    strList <- list()
    graphList <- list()
    scoreList <- list()

    # strengthBf <- bf.strength(hc(pcs), pcs)
    # averageBf <- averaged.network(strengthBf)
    # averageBf <- cextend(averageBf, strict=FALSE)

    for (r in Rrange){
        # cat(paste("performing R:", r, "\n"))
        strength <- boot.strength(pcs, algorithm=algo,
            R=r, cluster=cl, algorithm.args=NULL)
        strList[[paste0("R",r)]] <- strength
        av <- averaged.network(strength)
        # Infer the edge direction by default
        # av <- chooseEdgeDir(av, pcs, scoreType)
        av <- cextend(av, strict=FALSE)

        if (is.dag(bnlearn::as.igraph(av))){
            graphList[[paste0("R",r)]] <- av
            scoreList[[paste0("R",r)]] <- score(av, pcs, scoreType)
        } else {
            graphList[[paste0("R",r)]] <- av
            scoreList[[paste0("R",r)]] <- NA
        }
    }

    return(list("df"=pcs, "str"=strList, "graph"=graphList, "score"=scoreList))
}


#' bngenetest
#'
#' Testing various R for bayesian network between genes
#'
#' @param results the enrichment analysis result
#' @param exp gene expression matrix
#' @param expRow the type of the identifier of rows of expression matrix
#' @param expSample candidate rows to be included in the inference
#'                  default to all
#' @param algo structure learning method used in boot.strength()
#'             default to "hc"
#' @param algorithm.args parameters to pass to
#'                       bnlearn structure learnng function
#' @param cl cluster object from parallel::makeCluster()
#' @param Rrange the sequence of R values to be tested
#' @param scoreType return the specified scores
#' @param convertSymbol whether the label of resulting network is
#'                      converted to symbol, default to TRUE
#' @param pathNum the pathway number (the number of row of the original result,
#'                                    ordered by p-value)
#' @param orgDb perform clusterProfiler::setReadable
#'              based on this organism database
#' @param bypassConverting bypass symbol converting
#'
#' @return list of graphs and scores
#' @examples
#' data("exampleEaRes");data("exampleGeneExp")
#' res <- bngenetest(results = exampleEaRes, exp = exampleGeneExp,
#' algo="hc", Rrange=seq(10, 30, 10), pathNum=1, scoreType="bge")
#' @export
#'
bngenetest <- function (results, exp, expSample=NULL, algo="hc",
                        Rrange=seq(2,40,2), cl=NULL, algorithm.args=NULL,
                        pathNum=NULL, convertSymbol=TRUE, expRow="ENSEMBL",
                        scoreType="aic-g", orgDb=org.Hs.eg.db,
                        bypassConverting=FALSE)
{
    if (bypassConverting) {convertSymbol <- FALSE}
    # if (results@keytype == "kegg"){
    #     resultsGeneType <- "ENTREZID"
    # } else {
    #     resultsGeneType <- results@keytype
    # }
    if (!bypassConverting) {
        if (!is.null(orgDb)){
            results <- setReadable(results, OrgDb=orgDb)
        }
    }
    if (is.null(expSample)) {expSample <- colnames(exp)}
    res <- results@result

    genesInPathway <- unlist(strsplit(res[pathNum, ]$geneID, "/"))
    if (!bypassConverting) {
        if (!is.null(orgDb)){
            genesInPathway <- clusterProfiler::bitr(genesInPathway,
                                                    fromType="SYMBOL",
                                                    toType=expRow,
                                                    OrgDb=orgDb)[expRow][,1]
        }
    }

    pcs <- exp[ intersect(rownames(exp), genesInPathway), expSample ]

    if (convertSymbol) {
            matchTable <- clusterProfiler::bitr(rownames(pcs), fromType=expRow,
                                toType="SYMBOL", OrgDb=orgDb)
            if (sum(duplicated(matchTable[,1])) >= 1) {
                message("Removing expRow that matches the multiple symbols")
                matchTable <- matchTable[
                    !matchTable[,1] %in% matchTable[,1][
                        duplicated(matchTable[,1])],]
            }
            rnSym <- matchTable["SYMBOL"][,1]
            rnExp <- matchTable[expRow][,1]
            pcs <- pcs[rnExp, ]
            rownames(pcs) <- rnSym
    }

    pcs <- data.frame(t(pcs))

    strList <- list()
    graphList <- list()
    scoreList <- list()

    # strengthBf <- bf.strength(hc(pcs), pcs)
    # averageBf <- averaged.network(strengthBf)
    # averageBf <- cextend(averageBf, strict=FALSE)

    for (r in Rrange){
        # cat(paste("performing R:", r, "\n"))
        strength <- boot.strength(pcs, algorithm=algo,
            algorithm.args=algorithm.args, R=r, cluster=cl)
        strList[[paste0("R",r)]] <- strength
        av <- averaged.network(strength)
        # Infer the edge direction by default
        # av <- chooseEdgeDir(av, pcs, scoreType)
        av <- cextend(av, strict=FALSE)

        if (is.dag(bnlearn::as.igraph(av))){
            graphList[[paste0("R",r)]] <- av
            scoreList[[paste0("R",r)]] <- score(av, pcs, scoreType)
        } else {
            graphList[[paste0("R",r)]] <- av
            scoreList[[paste0("R",r)]] <- NA
        }
    }

    return(list("df"=pcs, "str"=strList, "graph"=graphList, "score"=scoreList))
}
