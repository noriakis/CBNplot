
#' bootReasonOneDiscrete
#' 
#' for shiny, perform bootstrap reasoning for event of one discrete node if network is not DAG.
#' Logic sampling is used.
#'
#' @param df expression data
#' @param R bootstrap number
#' @param node candidate node (one only)
#' @param evidences expressions of evidence e.g. c("`A` < 3", "`B` >= 2")
#' @return mean and sd of frequency of event level
#' @importFrom dplyr summarise_all mutate ungroup
#' @importFrom bnlearn as.bn
#' @noRd
bootReasonOneDiscrete <- function (df, R, node, evidences, algo="hc", ref=NULL, n=NULL, algorithm.args=NULL, cont=NULL) {
    bReason <- list()
    bDif <- list()
    
    for (i in seq_len(R)) {
        sampled <- sample(nrow(df), nrow(df), replace = TRUE)
        struc <- eval(parse(text=paste0(algo,"(df[sampled,]", algorithm.args,")")))

        if (!is.null(ref)){
            struc <- as.bn(igraph::intersection(bnlearn::as.igraph(ref),bnlearn::as.igraph(struc)))
            if (!igraph::is.dag(as.igraph(struc))){
                print("intersection not DAG")
                break
            } else {
                print("intersection DAG")
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


#' bootReasonOneWeb
#' 
#' for shiny, if event is numeric and the network is not DAG,
#' produce event distribution for one node from list of multiple evidences.
#'
#' @param df expression data
#' @param R bootstrap number
#' @param node candidate node (one only)
#' @param evidences list of evidences e.g. list(A=2, B=3)
#' @return mean and sd of weighted mean
#' 
#' @noRd
bootReasonOneWeb <- function (df, R, node, evidences, algo="hc", ref=NULL, n=NULL, algorithm.args=NULL, cont=NULL) {
    bReason <- list()
    bDif <- list()
    
    for (i in seq_len(R)) {
        sampled <- sample(nrow(df), nrow(df), replace = TRUE)
        struc <- eval(parse(text=paste0(algo,"(df[sampled,]", algorithm.args,")")))
        if (!is.null(ref)){
            struc <- as.bn(igraph::intersection(bnlearn::as.igraph(ref),bnlearn::as.igraph(struc)))
            if (!igraph::is.dag(as.igraph(struc))){
                print("intersection not DAG")
                break
            } else {
                print("intersection DAG")
            }
        }
        fitted <- bn.fit(struc, df)
        cpdRes <- queryCpDistLwWeb(fitted, candidate=node, evidenceList=evidences, n=n)
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


#' queryCpDistLsWeb
#'
#' for shiny, perform cpdist and produce plot if event is factor and network is DAG
#'
#' @param fitted bn.fit object
#' @param candidate One event node
#' @param evidences expressions of evidence e.g. c("`A` < 3", "`B` >= 2")
#' @param discPalette palette to be used in scale_fill_brewer
#'
#' @noRd
queryCpDistLsWeb <- function (fitted, candidate, evidences, discPalette="Set2",...) {
    conc <- c()
    for (e in evidences){
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


#' queryCpDistLwWeb
#'
#' for shiny, perform cpdist and produce plot if event is numeric and network is DAG
#'
#' @param fitted bn.fit object
#' @param candidate ONE event node
#' @param evidenceList list of evidences e.g. list(A=2, B=3)
#' @param point whether to show point of weighted mean on plot
#' @param pointSize point size
#' @param alpha The weight is shown by alpha (TRUE) or color (FALSE)
#'
#' @noRd
queryCpDistLwWeb <- function (fitted, candidate, evidenceList, point=FALSE, pointSize=5, alpha=TRUE, ...) {
    conc <- c() # Store the data per level
    for (l in seq_len(length(evidenceList))){
        ## Make label according to list
        labeller <- paste(names(evidenceList[[l]]), "=", evidenceList[[l]], collapse="_")
        tmp <- cpdist(fitted, nodes=c(candidate), evidence=evidenceList[[l]], method="lw", debug=FALSE, ...)
        conc <- rbind(conc, data.frame(tmp, attributes(tmp)$weights, rep(labeller, nrow(tmp))))
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
                  theme_bw()+xlab(candidate)
            if (point) {p <- p + geom_point(wm, mapping=aes(x=weightedMean, y=level, fill=level), shape=21, size=pointSize) }
        } else {
            p <- conc %>% ggplot()+
                  ggdist::geom_dots(aes(y=level, x=get(candidate), fill=weights), color=NA)+
                  scale_fill_ramp_continuous()+
                  theme_bw()+xlab(candidate)
            if (point) {p <- p + geom_point(wm, mapping=aes(x=weightedMean, y=level), color="tomato", size=pointSize) }
        }
        p 
        returnList[["plot"]] <- p
    }
    return(returnList)
}


#' bootReasonMultipleWeb
#' 
#' for shiny, reflect difference in color of plot when the network is not dag
#'
#' @param df expression data
#' @param R bootstrap number
#' @param nodes candidate node (accept multiple)
#' @param evidences list of evidences e.g. list(A=2, B=3)
#' @return mean and sd of weighted mean
#' @noRd
bootReasonMultipleWeb <- function (df, R, nodes, evidences, algo="hc", ref=NULL, n=NULL, algorithm.args=NULL) {
    mm <- list()
    bReason <- list()
    for (i in seq_len(R)) {
        sampled <- sample(nrow(df), nrow(df), replace = TRUE)
        struc <- eval(parse(text=paste0(algo,"(df[sampled,]", algorithm.args,")")))

        if (!is.null(ref)){
            struc <- as.bn(igraph::intersection(bnlearn::as.igraph(ref),bnlearn::as.igraph(struc)))
            if (!igraph::is.dag(as.igraph(struc))){
                print("intersection not DAG")
                break
            } else {
                print("intersection DAG")
                mm[[i]] <- struc
            }
        }

        fitted <- bn.fit(struc, df)
        
        labeller <- paste(names(evidences), "=", evidences, collapse="_")
        cpdRes <- cpdist(fitted, nodes=nodes, evidence=evidences, method="lw", n=n)
        tmpWm <- cpdRes %>% summarise_all(function(x) weighted.mean(x, attributes(cpdRes)$weights))
        tmpWm$level <- labeller
        bReason[[i]] <- tmpWm
    }

    wmMean <- bReason %>%
        purrr::reduce(rbind) %>% group_by(level) %>% summarise_all(list(mean, sd))
    
    return(list(wmMean,mm))
}
