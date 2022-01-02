

#' queryCpDistLs
#'
#' produce a plot of bnlearn::cpdist by performing bnlearn::cpdist on specified node, evidence and level.
#'
#' @param fitted bn.fit object
#' @param candidate name of node
#' @param evidences the evidences
#'
#' @return list of dataframe containing raw values
#' @examples queryCpDistLs(fitted, candidate="Mitotic Spindle Checkpoint", evidences=c("TP53<0.5","TP53>0.5","TP53>0.8"), n=5000)
#' @export
#'
queryCpDistLs <- function (fitted, candidate, evidences, discPalette="Set2",...) {
    evRenamed <- c()
    for (e in evidences){
      if (grepl("<=",e,fixed = TRUE)) {tmp <- unlist(strsplit(e,"<=")); newe <- paste0("`",tmp[1],"`<= ",tmp[2]); evRenamed <- c(evRenamed, newe)}
      else if (grepl(">=",e,fixed = TRUE)) {tmp <- unlist(strsplit(e,">=")); newe <- paste0("`",tmp[1],"`>= ",tmp[2]); evRenamed <- c(evRenamed, newe)}
      else if (grepl("<",e,fixed = TRUE)) {tmp <- unlist(strsplit(e,"<")); newe <- paste0("`",tmp[1],"`< ",tmp[2]); evRenamed <- c(evRenamed, newe)}
      else if (grepl(">",e,fixed = TRUE)) {tmp <- unlist(strsplit(e,">")); newe <- paste0("`",tmp[1],"`> ",tmp[2]); evRenamed <- c(evRenamed, newe)}
      else if (grepl("==",e,fixed = TRUE)) {tmp <- unlist(strsplit(e,"==")); newe <- paste0("`",tmp[1],"`== ",tmp[2]); evRenamed <- c(evRenamed, newe)}
    }

    ## If event is discrete
    if ("prob" %in% names(fitted[[candidate]])) {
        conc <- c()
        for (e in evRenamed){
            tmp <- eval(parse(text=paste0("cpdist(fitted, nodes='",candidate,"', evidence=",e,", method='ls', ...)")))
            conc <- rbind(conc, data.frame(tmp, rep(e, nrow(tmp))))
        }
        colnames(conc) <- c(candidate, "level")
        returnList <- list()
        returnList[["df"]] <- conc
        groupedDf <- conc %>% group_by(level, !!!syms(candidate)) %>%
                      summarise(n = n()) %>%
                      mutate(freq = n / sum(n))
        groupedPlot <- groupedDf %>%  ggplot(aes(level,freq,fill=eval(parse(text=candidate))))+
                  geom_bar(stat="identity",position='dodge')+
                  theme_bw()+scale_fill_brewer(palette=discPalette)

        returnList[["grouped"]] <- groupedDf
        returnList[["plot"]] <- groupedPlot
        return(returnList)
    }

    ## Otherwise
    conc <- c()
    for (e in evRenamed){
        tmp <- eval(parse(text=paste0("cpdist(fitted, nodes='",candidate,"', evidence=",e,", method='ls', ...)")))
        conc <- rbind(conc, data.frame(tmp, rep(e, nrow(tmp))))
    }
    colnames(conc) <- c(candidate, "level")
    
    returnList <- list()
    returnList[["df"]] <- conc

    if (length(candidate)==1){
        disMean <- conc %>% group_by(level) %>% summarise(mean=mean(get(candidate)), n=n())
        disMean$label <- paste0(disMean$level, " (mean=", round(disMean$mean,3), ")")
        lb <- c()
        for (nm in unique(conc$level)) {
            lab <- as.character(disMean %>% filter(level==nm) %>% select(label))
            num <- as.numeric(disMean %>% filter(level==nm) %>% select(n))
            lb <- c(lb, rep(lab, num))
        }
        conc$label <- lb
        p <- conc %>% ggplot(aes(y=level, x=get(candidate), color=label, fill=label)) +
            ggdist::stat_dotsinterval()+
            theme_bw() + xlab(candidate)
        returnList[["plot"]] <- p
    }
    return(returnList)
}


#' queryCpDistLw
#'
#' produce a plot of bnlearn::cpdist by performing bnlearn::cpdist on specified node, evidence and level.
#'
#'
#' @param fitted bn.fit object
#' @param candidate name of node
#' @param evidence evidence variable name
#' @param level level to be listed
#' @param alpha whether to reflect the weights by alpha (TRUE) or color (FALSE)
#' @param point geom_point the weighted mean
#' @importFrom dplyr summarise
#' @importFrom ggdist geom_dots
#' @importFrom ggdist scale_fill_ramp_continuous
#' @importFrom bnlearn cpdist
#' @return list of dataframe containing raw values
#' @examples queryCpDistLw(fitted, candidate="COL25A1", evidence = "Treatment", level=rownames(fitted$Treatment$prob))
#' @export
#'
queryCpDistLw <- function (fitted, candidate, evidence, level, point=FALSE, pointSize=5, alpha=TRUE, ...) {

    if ("prob" %in% names(fitted[[candidate]])) {stop("the cpdist(lw) for the categorical event node is currently not supported.")}

    conc <- c()
    if(length(evidence)!=1){stop("evidence must be one node")}
    for (l in level){
        ll <- list()
        ll[[evidence]] <- l
        tmp <- cpdist(fitted, nodes=c(candidate), evidence=ll, method="lw", ...)
        conc <- rbind(conc, data.frame(tmp, attributes(tmp)$weights, rep(l, nrow(tmp))))
    }
    
    colnames(conc) <- c(candidate, "weights", "level")
    conc$level <- factor(conc$level)
    wm <- conc %>% group_by(level) %>% dplyr::summarise_at(candidate, ~ weighted.mean(., weights))
    returnList <- list()
    returnList[["df"]] <- conc
    returnList[["wm"]] <- wm
    if (length(candidate)==1){
        ## Points are colored | set alpha by their weights
        if (alpha){
            p <- conc %>% ggplot()+
                  ggdist::geom_dots(aes(y=level, x=get(candidate), fill=level, alpha=weights), color=NA)+
                  theme_bw()+xlab(candidate)+ylab(evidence)
            if (point) {p <- p + geom_point(wm, mapping=aes(x=get(candidate), y=level, fill=level), shape=21, size=pointSize) }
        } else {
            p <- conc %>% ggplot()+
                  ggdist::geom_dots(aes(y=level, x=get(candidate), fill=weights), color=NA)+
                  scale_fill_ramp_continuous()+
                  theme_bw()+xlab(candidate)+ylab(evidence)
            if (point) {p <- p + geom_point(wm, mapping=aes(x=get(candidate), y=level), color="tomato", size=pointSize) }
        }
        p 
        returnList[["plot"]] <- p
    }
    return(returnList)
}

#' depKStest
#'
#' Perform KS-test recursively
#' @param type "cell_line" or "lineage"
#' @return pvalues and adjusted pvalues
#' @examples depKStest(pway, type="cell_line", cellLineName="253J_URINARY_TRACT", dep=depmap::depmap_crispr(), depMeta=depmap::depmap_metadata())
#' @export
#'
depKStest <- function(results, type, cellLineName=NULL, lineageName=NULL, adjMethod="bonferroni", dep=NULL, depMeta=NULL){

    # if (results@keytype == "ENTREZID" | results@keytype == "kegg"){
    #     q <- "entrez_id"
    # } else if (results@keytype == "SYMBOL"){
    #     q <- "gene_name"
    # } else {
    #     stop("provide entrezid, symbol or kegg")
    # }
    q <- "gene_name"

    if(is.null(dep)){stop("please provide depmap data")}
    if (type == "lineage"){
        if (is.null(depMeta)){stop("please provide metadata")}
        if (is.null(lineageName)){stop("please provide lineage name")}
        specificDep <- depMeta %>%
            dplyr::select(depmap_id, lineage) %>%
            dplyr::full_join(dep, by = "depmap_id") %>%
            filter(lineage == lineageName)
        y <- specificDep$dependency
    } else if (type == "cell_line"){
        if (is.null(cellLineName)){stop("please provide cell line name")}
        specificDep <- dep %>% filter(cell_line == cellLineName)
        y <- specificDep$dependency
    }

    ksp <- c()
    for (i in seq_len(length(rownames(results@result)))){
        genesInPathway <- strsplit(results@result[i, "geneID"], "/")[[1]]
        x <- (specificDep %>% filter((!!sym(q)) %in% genesInPathway) %>% select(dependency))$dependency
        if (length(x) == 0){
            ksp <- c(ksp, NA)
        } else {
            ksp <- c(ksp, ks.test(x, y)$p.value)
        }
    }
    pway@result$depP <- ksp
    pway@result$depAdjP <- p.adjust(ksp, adjMethod)
    return(pway)
}




#' inferMS
#'
#' multiscale bootstrap-based inference of Bayesian network
#'
#' @return object of class bn.strength
#'
#'
inferMS <- function(data, algo, algorithm.args, R, cl=NULL, r=seq(0.5, 1.5, 0.1))
{
    #stop("Currently cannot be used until the release of bnlearn 4.7")
    #threshold <- utils::getFromNamespace("threshold", "bnlearn")
    nList <- list()
    for (s in r){
        sb <- boot.strength(data=data, cluster=cl, algorithm=algo, R=R,
                            m=as.integer(nrow(data)*s), algorithm.args=algorithm.args)
        nList[[paste0("s",s)]] <- sb
    }
    md <- nList %>% purrr::reduce(dplyr::left_join, by = c("from","to"))

    stcol <- c("from", "to", colnames(md)[grepl("strength",colnames(md))])
    dicol <- c("from", "to", colnames(md)[grepl("direc",colnames(md))])

    st <- md[,stcol]
    ms <- c()
    for (i in seq_len(nrow(st))){
        tmpbp <- as.numeric(st[i,][3:length(st[i,])])
        ft <- pvclust::msfit(tmpbp, r, R)
        ms <- c(ms, as.numeric(1-pnorm(ft$coef[1]-ft$coef[2])))
    }
    st$MS <- ms

    di <- md[,dicol]
    dis <- c()
    for (i in seq_len(nrow(di))){
        tmpbp <- as.numeric(di[i,][3:length(di[i,])])
        ft <- pvclust::msfit(tmpbp, r, R)
        dis <- c(dis, as.numeric(1-pnorm(ft$coef[1]-ft$coef[2])))
    }
    di$DIR <- dis
    nn <- unique(c(di$from, di$to))
    res <- data.frame(di$from, di$to, st$MS, di$DIR)
    colnames(res) <- c("from","to","strength","direction")

    res = structure(res, method = "bootstrap", threshold = 0,
                    nodes = nn,
                    class = c("bn.strength", class(res)))

    ## Using bnlearn version 4.7.20210803 or above for function inclusion.threshold().
    if (packageVersion("bnlearn")=="4.7"){
        attributes(res)$threshold <- bnlearn::inclusion.threshold(res)
    }
    return(res)
}

#' bnpathtest
#'
#' Testing various R for bayesian network between pathways
#'
#' @param results the enrichment analysis result from clusterProfiler or ReactomePA
#' @param exp gene expression matrix
#' @param expRow the type of the identifier of rows of expression matrix
#' @param expSample candidate rows to be included in the inference, default to all
#' @param algo structure learning method used in boot.strength(), default to "hc"
#' @param cl cluster object from parallel::makeCluster()
#' @param qvalueCutoff the cutoff value for qvalue
#' @param adjpCutOff the cutoff value for adjusted pvalues
#' @param nCategory the number of pathways to be included
#' @param Rrange the sequence of R values to be tested
#' @param scoreType return the specified scores
#'
#' @return list of graphs and scores
#' @examples bnpathtest(results = pway, algo="hc", exp = vsted, Rrange=seq(10, 200, 5), expSample = rownames(subset(meta, Condition=="T")), expRow = "ENSEMBL", scoreType="bge")
#' @export
#'
bnpathtest <- function (results, exp, expSample=NULL, algo="hc", algorithm.args=NULL, expRow="ENSEMBL", cl=NULL,
                        qvalueCutOff=0.05, adjpCutOff=0.05, nCategory=15, Rrange=seq(2,40,2), scoreType="aic-g")
{
    if (is.null(expSample)) {expSample=colnames(exp)}
    # if (results@keytype == "kegg"){
    #     resultsGeneType <- "ENTREZID"
    # } else {
    #     resultsGeneType <- results@keytype
    # }

    res <- results@result
    if (!is.null(qvalueCutOff)) { res <- subset(res, qvalue < qvalueCutOff) }
    if (!is.null(adjpCutOff)) { res <- subset(res, p.adjust < adjpCutOff) }
    if (nCategory) {
        res <- res[seq_len(nCategory),]
        res <- res[!is.na(res$ID),]
    }


    pcs <- c()
    pwayNames <- c()
    for (i in seq_len(length(rownames(res)))) {
        genesInPathway <- strsplit(res[i, ]$geneID, "/")[[1]]
        genesInPathway <- clusterProfiler::bitr(genesInPathway,
                                                fromType="SYMBOL",
                                                toType=expRow,
                                                OrgDb=org.Hs.eg.db)[expRow][,1]

        pathwayMatrix <- exp[ intersect(rownames(exp), genesInPathway), expSample ]
        pathwayMatrixPca <- prcomp(t(pathwayMatrix))
        pwayNames <- c(pwayNames, res[i,]$Description)
        pcs <- cbind(pcs, pathwayMatrixPca$x[,"PC1"])
    }
    colnames(pcs) <- pwayNames
    pcs <- data.frame(pcs, check.names=FALSE)

    strList <- list()
    graphList <- list()
    scoreList <- list()

    # strengthBf <- bf.strength(hc(pcs), pcs)
    # averageBf <- averaged.network(strengthBf)
    # averageBf <- cextend(averageBf, strict=F)

    for (r in Rrange){
        strength <- boot.strength(pcs, algorithm=algo, R=r, cluster=cl, algorithm.args=NULL)
        strList[[paste0("R",r)]] <- strength
        av <- averaged.network(strength)
        # Infer the edge direction by default
        av <- chooseEdgeDir(av, pcs, scoreType)
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
#' @param results the enrichment analysis result from clusterProfiler or ReactomePA
#' @param exp gene expression matrix
#' @param expRow the type of the identifier of rows of expression matrix
#' @param expSample candidate rows to be included in the inference, default to all
#' @param algo structure learning method used in boot.strength(), default to "hc"
#' @param cl cluster object from parallel::makeCluster()
#' @param Rrange the sequence of R values to be tested
#' @param scoreType return the specified scores
#' @param convertSymbol whether the label of resulting network is converted to symbol, default to TRUE
#' @param pathNum the pathway number (the number of row of the original result, ordered by p-value)
#'
#' @return list of graphs and scores
#' @examples bngenetest(results = pway, algo="hc", exp = vsted, Rrange=seq(10, 200, 5), pathNum=15, expSample = rownames(subset(meta, Condition=="T")), expRow = "ENSEMBL", scoreType="bge")
#' @export
#'
bngenetest <- function (results, exp, expSample=NULL, algo="hc", Rrange=seq(2,40,2), cl=NULL, algorithm.args=NULL,
                        pathNum=NULL, convertSymbol=TRUE, expRow="ENSEMBL",scoreType="aic-g")
{
    # if (results@keytype == "kegg"){
    #     resultsGeneType <- "ENTREZID"
    # } else {
    #     resultsGeneType <- results@keytype
    # }

    res <- results@result

    genesInPathway <- unlist(strsplit(res[pathNum, ]$geneID, "/"))
    genesInPathway <- clusterProfiler::bitr(genesInPathway,
                                            fromType="SYMBOL",
                                            toType=expRow,
                                            OrgDb=org.Hs.eg.db)[expRow][,1]
    pcs <- exp[ intersect(rownames(exp), genesInPathway), expSample ]

    if (convertSymbol) {
          matchTable <- clusterProfiler::bitr(rownames(pcs), fromType=expRow, toType="SYMBOL", OrgDb=org.Hs.eg.db)
          if (sum(duplicated(matchTable[,1])) >= 1) {
            message("Removing expRow that matches the multiple symbols")
            matchTable <- matchTable[!matchTable[,1] %in% matchTable[,1][duplicated(matchTable[,1])],]
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
        strength <- boot.strength(pcs, algorithm=algo, algorithm.args=algorithm.args, R=r, cluster=cl)
        strList[[paste0("R",r)]] <- strength
        av <- averaged.network(strength)
        # Infer the edge direction by default
        av <- chooseEdgeDir(av, pcs, scoreType)
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


#' chooseEdgeDir
#' 
#' Recursively perform choose.direction() from bnlearn on undirected arcs.
#'
#' @noRd
chooseEdgeDir <- function(av, pcs, scoreType, debug=TRUE) {
    if (dim(undirected.arcs(av))[1]!=0){
        message("found undirected arc(s)")
        undir <- undirected.arcs(av)
        for (n in seq_len(dim(undir)[1])) {
            av <- choose.direction(av, c(undir[,1][n], undir[,2][n]),
                               pcs, criterion = scoreType, debug = debug)
        }
    }
    return(av)
}

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
    incPathG <- set.vertex.attribute(incPathG, "name", value=res[names(V(incPathG)), "Description"])
    undirG = graph_from_edgelist(as_edgelist(g), directed = FALSE)
    ovlG = as_edgelist(igraph::intersection(incPathG, undirG))
    ovlEdge = rbind(cbind(ovlG[,1], ovlG[,2]), cbind(ovlG[,2], ovlG[,1]))
    ovlEMat = c()
    for (i in seq_len(dim(ovlEdge)[1])){
        ovlEMat = c(ovlEMat, paste(ovlEdge[i,], collapse="_"))
    }
    edgeLabel = c()
    eList = as_edgelist(g)
    for (i in seq_len(dim(eList)[1])){
        if (paste(eList[i,], collapse="_") %in% ovlEMat){
            edgeLabel = c(edgeLabel, "solid")
        } else {
            edgeLabel = c(edgeLabel, "dotted")
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
#' @return subset of enrichment results
#' @examples obtainPath(results = pway, geneSymbol="E2F1")
#' @export

obtainPath <- function(res, geneSymbol) {
    # if (res@keytype=="kegg"){tot <- "ENTREZID"} else {tot <- res@keytype}
    # candidateGene <- clusterProfiler::bitr(geneSymbol, fromType = "SYMBOL", toType = tot, org.Hs.eg.db)[tot]
    res@result <- res@result[unlist(purrr::map(purrr::map(res@result$geneID,
                                                                  function (x) unlist(strsplit(x, "/"))),
                                                       function(x) geneSymbol %in% x)),]
    return(res)
}

#' bootReasonOne
#' 
#' bootstrapping probabilistic reasoning and return the mean, and the difference in mean from control level.
#' The sampling is performed by likelihood weighting. One evidence can be specified with multiple evidence levels.
#' If specified control level, the function additionally returns the mean differences.
#'
#' @param df data frame to perform structure learning
#' @param R bootstrap replicates number
#' @param node node names of the event
#' @param evidence node name of the evidence
#' @param level evidence levels
#' @param algo structure learning algorithm, default to hc
#' @param algorithm.args passed to algorithm function
#' @param cont control level
#' @importFrom bnlearn bn.fit
#' @importFrom dplyr summarize_at
#' @importFrom dplyr inner_join
#' @importFrom tidyr starts_with
#' @importFrom utils packageVersion
#' @return the mean, and the difference in mean from control level
#' @examples bootReasonOne(df, R=20, node=c("CENPH","NUP107"), evidence=c("TP53"), level=c(0, 0.25, 0.5, 0.75, 1), cont=0)
#' @export
bootReasonOne <- function (df, R, node, evidence, level, algo="hc", n=NULL, algorithm.args=NULL, cont=NULL) {
    bReason <- list()
    bDif <- list()
    
    for (i in seq_len(R)) {
        sampled <- sample(nrow(df), nrow(df), replace = TRUE)
        struc <- eval(parse(text=paste0(algo,"(df[sampled,]", algorithm.args,")")))
        fitted <- bn.fit(struc, df)
        cpdRes <- queryCpDistLw(fitted, candidate=node, evidence=evidence, level=level, n=n)
        cpdRes <- cpdRes$df
        tmpWm <- cpdRes %>% group_by(level) %>% summarize_at(node, ~ weighted.mean(., weights))
        bReason[[i]] <- tmpWm
        
        if (!is.null(cont)){
            difTbl <- c()
            cntTib <- tmpWm %>% filter(level==cont)
            for (lvl in level[!level %in% cont]) {
                lvlTib <- tmpWm %>% filter(level==lvl)
                difTib <- lvlTib[,node] - cntTib[,node]
                difTbl <- rbind(difTbl, difTib)
            }
            difTbl$level <- level[!level %in% cont]
            bDif[[i]] <- difTbl
        }
    }
    wmMean <- bReason %>%
        purrr::reduce(inner_join, by = "level")
    
    if (!is.null(cont)){
        wmDif <- bDif %>%
            purrr::reduce(inner_join, by = "level")
    }
    
    res <- c()
    difres <- c()
    
    for (n in node) {
        if (!is.null(cont)){
            difMean <- wmDif %>%
                select(starts_with(n)) %>%
                select(-1) %>% 
                summarize(Mean= rowMeans(., na.rm=TRUE), stdev=rowSds(as.matrix(.), na.rm=TRUE))
            difMean$level <- level[!level %in% cont]
            difMean$node <- n
            difres <- rbind(difres, difMean)
        }
        
      currentLevel <- wmMean$level
      levelMean <- wmMean %>%
        select(starts_with(n)) %>%
        select(-1) %>% 
        summarize(Mean= rowMeans(., na.rm=TRUE), stdev=rowSds(as.matrix(.), na.rm=TRUE))
      levelMean$level <- currentLevel
      levelMean$node <- n
      res <- rbind(res, levelMean)
    }
    resList <- list()
    resList[["distMean"]] <- res
    if (!is.null(cont)){
        resList[["control"]] <- paste(evidence,"=",cont)
        resList[["difMean"]] <- difres
    }
    return(resList)
}




#' bootReasonMultiple
#' 
#' bootstrapping probabilistic reasoning and return the mean of bootstrapped values.
#' Multiple events, evidences and levels can be specified. The likelihood weighting was used for sampling.
#' 
#' @param df data frame to perform structure learning
#' @param R bootstrap replicates number
#' @param nodes node names of the event
#' @param evidences node names of the evidence
#' @param level evidence levels
#' @param algo structure learning algorithm, default to hc
#' @param algorithm.args passed to algorithm function
#' @importFrom bnlearn bn.fit
#' @importFrom dplyr summarize_at
#' @importFrom dplyr inner_join
#' @importFrom tidyr starts_with
#' @return the mean of the bootstrapped values for each level
#' @examples bootReasonMultiple(df, R, nodes=c("A"), evidences=c("B","C"), levels=c(2,3))
#' @export
bootReasonMultiple <- function (df, R, nodes, evidences, levels, algo="hc", n=NULL, algorithm.args=NULL) {
    
    bReason <- list()
    
    if (length(evidences)!=length(levels)){stop("please provide same length for evidence and level")}
    
    evList <- list()
    for (n in seq_len(length(evidences))){
        evList[[evidences[n]]] <- levels[n]
    }
    

    for (i in seq_len(R)) {
        
        sampled <- sample(nrow(df), nrow(df), replace = TRUE)
        struc <- eval(parse(text=paste0(algo,"(df[sampled,]", algorithm.args,")")))
        fitted <- bn.fit(struc, df)
        
        levelDf <- data.frame()
        for (e in seq_len(length(evidences))) {
            newList <- list()
            newList[[evidences[e]]] <- levels[e]
            cpdRes <- cpdist(fitted, nodes=nodes, evidence=newList, method="lw", n=n)
            tmpWm <- cpdRes %>% summarise_all(function(x) weighted.mean(x, attributes(cpdRes)$weights))
            tmpWm$level <- paste0(evidences[e],"=",levels[e])
            levelDf <- rbind(levelDf, tmpWm)
        }
        bReason[[i]] <- levelDf
    }

    wmMean <- bReason %>%
        purrr::reduce(rbind) %>% group_by(level) %>% summarise_all(list(mean, sd))
    
    return(wmMean)
}


#' bootReasonOneDiscrete
#' 
#' Perform bootstrap reasoning for event of one discrete node if network is not DAG.
#' Logic sampling is used. Return the mean of frequency.
#'
#' @param df expression data
#' @param R bootstrap number
#' @param node candidate node (one only)
#' @param evidences expressions of evidence e.g. c("`A` < 3", "`B` >= 2")
#' @return mean and sd of frequency of event level
#' @importFrom dplyr summarise_all mutate ungroup
#' @importFrom bnlearn as.bn
#' @export
#' 
bootReasonOneDiscrete <- function (df, R, node, evidences, algo="hc", ref=NULL, n=NULL, algorithm.args=NULL, cont=NULL) {
    bReason <- list()
    bDif <- list()
    
    for (i in seq_len(R)) {
        sampled <- sample(nrow(df), nrow(df), replace = TRUE)
        struc <- eval(parse(text=paste0(algo,"(df[sampled,]", algorithm.args,")")))

        if (!is.null(ref)){
            struc <- as.bn(igraph::intersection(bnlearn::as.igraph(ref),bnlearn::as.igraph(struc)))
            if (!igraph::is.dag(as.igraph(struc))){
                message("intersection not DAG")
                break
            } else {
                message("intersection DAG")
            }
        }

        fitted <- bn.fit(struc, df)
        cpdRes <- queryCpDistLsWeb(fitted, candidate=node, evidences=evidences, n=n)
        cpdRes <- cpdRes$df
        cpdRes <- cpdRes %>% group_by(level, !!!syms(node)) %>%
                  summarise(n = n()) %>%
                  mutate(freq = n / sum(n))

        bReason[[i]] <- cpdRes
    }

    wmMean <- bReason %>%
            purrr::reduce(inner_join, by = c("level", node))
    wmMean <- cbind(
        wmMean %>% ungroup() %>% dplyr::select(1,2),
        wmMean %>% ungroup() %>% dplyr::select(starts_with("freq")) %>% summarize(freqMean= rowMeans(., na.rm=TRUE), freqStdev=rowSds(as.matrix(.), na.rm=TRUE)),
        wmMean %>% ungroup() %>% dplyr::select(starts_with("n")) %>% summarize(nMean= rowMeans(., na.rm=TRUE), nStdev=rowSds(as.matrix(.), na.rm=TRUE))
    )
    return(wmMean)
}

#' compareBNs
#' 
#' Take the list of networks and returns the F-measures
#' 
#' @param nets data frame to perform structure learning
#' @return F-measures of each combination of network
#' @examples compareBNs(nets)
#' @export
compareBNs <- function(listOfNets){
  if (length(listOfNets)<2){
    stop("please provide multiple networks from bnlearn")
  }
  fms <- c()
  cmb <- combn(seq_len(length(listOfNets)), m=2)
  for (i in seq_len(ncol(cmb))){
    compareRes <- bnlearn::compare(listOfNets[[cmb[,i][1]]], listOfNets[[cmb[,i][2]]])
    prec <- compareRes$tp/(compareRes$fp+compareRes$tp)
    recall <- compareRes$tp/(compareRes$fn+compareRes$tp)
    fm <- 2 * prec * recall / (recall + prec)
    fms <- c(fms, fm)
  }
  return(fms)
}