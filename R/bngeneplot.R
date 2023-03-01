#' bngeneplot
#'
#' Plot gene relationship within the specified pathway
#' 
#'
#' @param results the enrichment analysis result
#' @param exp gene expression matrix
#' @param expRow the type of the identifier of rows of expression matrix
#' @param expSample candidate samples to be included in the inference
#'                  default to all
#' @param algo structure learning method used in boot.strength()
#'             default to "hc"
#' @param algorithm.args parameters to pass to bnlearn
#'                       structure learnng function
#' @param R the number of bootstrap
#' @param pathNum the pathway number
#'                (the number of row of the original result,
#'                 ordered by p-value)
#' @param convertSymbol whether the label of resulting network is 
#'                      converted to symbol, default to TRUE
#' @param bypassConverting bypass the symbol converting
#'                         If you use custom annotation databases that 
#'                         does not have SYMBOL listed in keys.
#'                         ID of rownames and those listed in EA result
#'                         must be same.
#' @param interactive whether to use bnviewer (default to FALSE)
#' @param cexCategory scaling factor of size of nodes
#' @param delZeroDegree delete zero degree nodes
#' @param disc discretize the expressoin data
#' @param tr Specify data.frame if one needs to discretize
#'           as the same parametersas the other dataset
#' @param remainCont Specify characters when perform discretization,
#'                   if some columns are to be remain continuous
#' @param cl cluster object from parallel::makeCluster()
#' @param showDir show the confidence of direction of edges
#' @param showDepHist whether to show depmap histogram
#' @param chooseDir if undirected edges are present,
#'                  choose direction of edges (default: FALSE)
#' @param scoreType score type to use on choosing direction
#' @param labelSize the size of label of the nodes
#' @param layout ggraph layout, default to "nicely"
#' @param clusterAlpha if specified multiple pathways,
#'                    the parameter is passed to geom_mark_hull()
#' @param strType "normal" or "ms" for multiscale implementation
#' @param sp query to graphite::pathways(), default to "hsapiens"
#' @param compareRef whether compare to the reference network
#' @param compareRefType "intersection" or "difference"
#' @param pathDb query to graphite::pathways(), default to "reactome"
#' @param dep the tibble storing dependency score from library depmap
#' @param sizeDep whether to reflect DepMap score to the node size
#' @param cellLineName the cell line name to be included
#' @param strengthPlot append the barplot depicting edges with high strength
#' @param nStrength specify how many edges are included in the strength plot
#' @param showLineage show the dependency score across the lineage
#' @param depMeta depmap::depmap_metadata(), needed for showLineage
#' @param strThresh the threshold for strength
#' @param hub visualize the genes with top-n hub scores
#' @param returnNet whether to return the network
#' @param otherVar other variables to be included in the inference
#' @param otherVarName the names of other variables
#' @param onlyDf return only data.frame used for inference
#' @param orgDb perform clusterProfiler::setReadable
#'              based on this organism database
#' @param shadowText whether to use shadow text for the better readability
#'                    default: TRUE
#' @param bgColor color for text background when shadowText is TRUE
#' @param textColor color for text when shadowText is TRUE
#' @param seed A random seed to make the analysis reproducible, default is 1.
#' @param useSiGN default to FALSE.
#' For using SiGN-BN in the function in Windows 10/11,
#' 1. Download the SiGN-BN HC+BS binary in WSL (https://sign.hgc.jp/signbn/download.html)
#' 2. Set PATH to executable (sign.1.8.3)
#' @return ggplot2 object
#'
#' @examples
#' data("exampleEaRes");data("exampleGeneExp")
#' res <- bngeneplot(results = exampleEaRes, exp = exampleGeneExp, pathNum = 1,
#'                   R = 10, convertSymbol = TRUE, expRow = "ENSEMBL")
#'
#' @import oaqc
#' @importFrom dplyr group_by summarize arrange n
#' @importFrom utils write.table
#' @importFrom graphite pathways convertIdentifiers pathwayGraph
#' @importFrom clusterProfiler setReadable
#' @importFrom rlang .data
#' @importFrom stats cor p.adjust prcomp weights
#' @export
#'

bngeneplot <- function (results, exp, expSample=NULL, algo="hc", R=20,
                        returnNet=FALSE, algorithm.args=NULL,
                        bypassConverting=FALSE,
                        pathNum=NULL, convertSymbol=TRUE, expRow="ENSEMBL",
                        interactive=FALSE, cexCategory=1, cl=NULL,
                        showDir=FALSE, chooseDir=FALSE, scoreType="bic-g",
                        labelSize=4, layout="nicely", clusterAlpha=0.2,
                        strType="normal", delZeroDegree=TRUE,
                        otherVar=NULL, otherVarName=NULL, onlyDf=FALSE,
                        disc=FALSE, tr=NULL, remainCont=NULL,
                        sp="hsapiens", compareRef=FALSE,
                        compareRefType="intersection", pathDb="reactome",
                        dep=NULL, depMeta=NULL, sizeDep=FALSE, showDepHist=TRUE,
                        cellLineName="5637_URINARY_TRACT",
                        showLineage=FALSE, orgDb=org.Hs.eg.db, shadowText=TRUE,
                        bgColor="white", textColor="black",
                        strengthPlot=FALSE, nStrength=10, strThresh=NULL,
                        hub=NULL, seed = 1, useSiGN=FALSE) {
    
    if (is.null(expSample)) {expSample <- colnames(exp)}
    if (compareRef & length(pathNum) > 1){
        stop("compareRef can be used with one pathNum or pathName.")}
    if (interactive & compareRef){
        stop("compareRef must be set to FALSE when use bnviewer.")}
    if (!is.numeric(pathNum)){
        stop("Please specify number(s) for pathNum.")}
    if (sizeDep & !convertSymbol){
        stop("sizeDep must be used with convertSymbol set to TRUE.")}
    if (sizeDep & length(pathNum) > 1){
        stop("sizeDep can be used with one pathNum.")}
    if (compareRef & is.null(pathDb)){
        stop("please specify which database to use as reference.")}
    if (showLineage & strengthPlot){
        stop("please specify one of showLineage or strengthPlot.")}

    # if (results@keytype == "kegg"){
    #     resultsGeneType <- "ENTREZID"
    # } else {
    #     resultsGeneType <- results@keytype
    # }

    if (bypassConverting){convertSymbol <- FALSE}
    if (!bypassConverting){
        if (!is.null(orgDb)){
            results <- setReadable(results, OrgDb=orgDb)
        }
    }
    ## The newer version of reactome.db
    results@result$Description <- gsub("Homo sapiens\r: ",
                                    "",
                                    results@result$Description)
    tmpCol <- colnames(results@result)
    tmpCol[tmpCol=="core_enrichment"] <- "geneID"
    tmpCol[tmpCol=="qvalues"] <- "qvalue"
    tmpCol[tmpCol=="setSize"] <- "Count"
    colnames(results@result) <- tmpCol

    if (showLineage) {
        if (is.null(dep)){dep <- depmap::depmap_crispr()}
        if (is.null(depMeta)){depMeta <- depmap::depmap_metadata()}
    }

    if (sizeDep) {
        if (is.null(cellLineName)){stop("Please specify cell line name.")}
        if (is.null(dep)){dep <- depmap::depmap_crispr()}
        filteredDep <- dep %>% filter(.data$cell_line==cellLineName)
        depHist <- ggplot(filteredDep, aes_(x=~dependency)) +
            geom_histogram(aes_(fill=~..count..), col="black") +
            scale_fill_gradient("Count", low = "blue", high = "red") +
            theme_minimal(base_family = "Arial Narrow") +
            ggtitle(cellLineName)+
            theme(plot.title = element_text(hjust=0.5, face="bold"),
                axis.text = element_text(size=10),
                axis.title = element_text(size=12))
    }
    if (sizeDep){
        ## when reflecting DepMap scores to node size
        scaleSizeLow <- 1
        scaleSizeHigh <- 10
    } else {
        scaleSizeLow <- 3
        scaleSizeHigh <- 8
    }

    res <- results@result

    genesInPathway <- unique(unlist(strsplit(res[pathNum, ]$geneID, "/")))
    if (!bypassConverting) {
        if (!is.null(orgDb)) {
            genesInPathway <- clusterProfiler::bitr(genesInPathway,
                                                    fromType="SYMBOL",
                                                    toType=expRow,
                                                    OrgDb= orgDb )[expRow][,1]
        }
    }

    pcs <- exp[ intersect(rownames(exp), genesInPathway), expSample ]

    if (!bypassConverting) {
        if (expRow!="SYMBOL"){
            if (convertSymbol && !is.null(orgDb)) {
                # rn <- clusterProfiler::bitr(rownames(pcs),
                #                             fromType=expRow,
                #                             toType="SYMBOL",
                #                             OrgDb=org.Hs.eg.db)["SYMBOL"][,1]
                ## Change expression matrix rownames to symbol
                ## If one "expRow" hit to multiple symbols,
                ## delete the ID from the subsequent analysis, showing warning.
                matchTable <- clusterProfiler::bitr(rownames(pcs), fromType=expRow,
                                                    toType="SYMBOL", OrgDb=orgDb)
                if (sum(duplicated(matchTable[,1])) >= 1) {
                    message("Removing IDs that matches the multiple symbols")
                    matchTable <- matchTable[
                        !matchTable[,1] %in% matchTable[,1][
                            duplicated(matchTable[,1])
                        ],
                    ]
                }
                rnSym <- matchTable["SYMBOL"][,1]
                rnExp <- matchTable[expRow][,1]
                pcs <- pcs[rnExp, ]
                rownames(pcs) <- rnSym
            }
        }
    }

    pcs <- data.frame(t(pcs))

    if (sizeDep){
        pcs <- pcs[,intersect(filteredDep$gene_name, colnames(pcs)),]
    }


    if (length(pathNum) >= 2) {
        clus <- c()

        for (num in pathNum){
            tmp <- data.frame(strsplit(res[num, ]$geneID, "/")[[1]])
            tmp$Pathway <- res[num, ]$Description
            clus <- rbind(clus, tmp)
        }

        colnames(clus) <- c("geneID", "Pathway")

        if (convertSymbol) {
            # m <- clusterProfiler::bitr(clus$geneID,
            #                            fromType=resultsGeneType,
            #                            toType="SYMBOL",
            #                            OrgDb=org.Hs.eg.db)
            # colnames(m) <- c("geneID","SYMBOL")
            mcl <- clus#merge(m, clus)
            colnames(mcl) <- c("SYMBOL","Pathway")
            cnt <- mcl %>% group_by(.data$SYMBOL) %>%
                arrange(.data$Pathway, .by_group = TRUE) %>%
                summarize(n=n(), Pathway=paste0(.data$Pathway, collapse = " + "))
            ovl <- cnt[cnt$n > 1,]
            mclSub <- subset(mcl, mcl$SYMBOL %in% cnt[cnt$n==1,]$SYMBOL)
            cls <- rbind(mclSub[,c("SYMBOL","Pathway")],
                        ovl[,c("SYMBOL","Pathway")])
            rownames(cls) <- cls$SYMBOL

        } else {
            if (!is.null(orgDb)){
                if (bypassConverting) {
                    m <- cbind(clus$geneID, clus$geneID)
                } else {
                    m <- clusterProfiler::bitr(clus$geneID,
                                                fromType="SYMBOL",
                                                toType=expRow,
                                                OrgDb=orgDb)
                }
                colnames(m) <- c("geneID", expRow)
                mcl <- merge(m, clus)
                cnt <- mcl %>% group_by_at(expRow) %>%
                    arrange(.data$Pathway, .by_group = TRUE) %>%
                    summarize(n=n(), Pathway=paste0(.data$Pathway, collapse = " + "))
                ovl <- cnt[cnt$n > 1,]
                mclSub <- subset(mcl,
                    mcl[expRow][,1] %in% as.character(
                        as.matrix(cnt[cnt$n==1,][,expRow])))
                cls <- data.frame(rbind(mclSub[,c(expRow,"Pathway")],
                    ovl[,c(expRow,"Pathway")]))
                rownames(cls) <- cls[, 1]
            }
        }
    }
    
    geneNames <- colnames(pcs)

    ## Insert other vars
    if (!is.null(otherVar)) {
        pcs <- cbind(pcs, otherVar)
        if (!is.null(otherVarName)){
            colnames(pcs) <- c(geneNames, otherVarName)
        }
    }

    if (disc){
        pcsRaw <- pcs
        pcs <- discDF(pcs, tr=tr, remainCont=remainCont)
    } else {
        pcsRaw <- pcs ## Hold pcsRaw anyway
    }

    if (onlyDf){
        return(pcs)
    }

    # print(dim(pcs))
    if (dim(pcs)[2]<=1){
        message("the number of gene is zero or one");return("error")}

    ## Bootstrap-based inference
    if (!useSiGN){
        if (strType == "normal"){
            strength <- withr::with_seed(seed = seed,
                boot.strength(pcs, algorithm=algo,
                    algorithm.args=algorithm.args, R=R, cluster=cl))
        } else if (strType == "ms"){
            strength <- withr::with_seed(seed = seed,
                inferMS(pcs, algo=algo, algorithm.args=algorithm.args, R=R, cl=cl))
        }
    } else {
        prefix <- gsub("\\.","",format(Sys.time(), "%Y%m%d%H%M%OS3"))
        tmpPath <- paste0(prefix,"tmpmat.txt")
        write.table(t(pcs), tmpPath, quote=FALSE,
            row.names=TRUE, col.names=FALSE, sep="\t")
        system(paste0('bash -c "signbn.1.8.3 --total-mem 1000 -N ',R,' -o ',
                      prefix,'_net.txt ',
                      tmpPath, '"'))
        unlink(tmpPath)
        net <- loadSign(paste0(prefix,'_net.txt'))
        strength <- net$str
    }

    ## Barplot of edge strength
    if (strengthPlot){
        strengthTop <- strength[order(strength$strength+strength$direction,
            decreasing = TRUE),][seq_len(nStrength),]
        strengthTop$label <- paste(strengthTop$from, "->", strengthTop$to)
        stp <- strengthTop %>%
            tidyr::pivot_longer(cols=c(.data$strength, .data$direction)) %>%
            ggplot(aes_(x=~label, y=~value, fill=~name))+
            geom_bar(position="dodge",stat="identity")+
            coord_flip(ylim=c(min(strengthTop$strength,
                strengthTop$direction)-0.05,1.0))+xlab("edges")+
            theme_bw()+scale_fill_manual(values = c("tomato","dodgerblue")) +
            scale_x_discrete(labels = function(x) stringr::str_wrap(x,
                                                                    width = 25))
    }

    ## Average by specified threshold
    if (!useSiGN) {
        if (!is.null(strThresh)){
            av <- averaged.network(strength, threshold=strThresh)
        } else {
            av <- averaged.network(strength)
        }
    } else {
        av <- net$av
    }

    # if (chooseDir){
    #     av <- chooseEdgeDir(av, pcs, scoreType)
    # }

    av <- cextend(av, strict=FALSE)

    if (interactive) {
        bnviewer::strength.viewer(bayesianNetwork = av,
            bayesianNetwork.boot.strength = strength)
    } else {

        g <- bnlearn::as.igraph(av)
        e <- as_edgelist(g, names = TRUE)
        if (dim(e)[1]==0){message("no edge present in graph");return("error")}
        eName <-paste0(e[,1], "_", e[,2])
        colnames(e) <- c("from","to")
        eDf <- merge(e, strength)
        rownames(eDf) <- paste0(eDf$from, "_", eDf$to)
        eDf <- eDf[eName, ]
        g <- set.edge.attribute(g, "color", index=E(g), eDf$strength)
        if (showDir){
            g <- set.edge.attribute(g, "label", index=E(g),
                round(eDf$direction,2))
        } else {
            g <- set.edge.attribute(g, "label", index=E(g), NA)
        }
        
        ## Hub genes
        hScore <- hub.score(g, scale = TRUE, weights = E(g)$color)$vector

        if (!is.null(hub)){
            defHub <- hScore[order(hScore, decreasing=TRUE)][seq_len(hub)]
            nodeShape <- names(V(g)) %in% names(defHub)
            nodeShape <- ifelse(nodeShape, 19, 21)
            V(g)$shape <- nodeShape
        } else {
            V(g)$shape <- rep(19, length(V(g)))
        }

        E(g)$width <- E(g)$color
        edgeWName <- "strength"

        if (sizeDep){
            sizeLab <- "-dependency"
            filteredDep <- filteredDep %>%
                filter(.data$gene_name %in% names(V(g))) %>%
                arrange(match(.data$gene_name, names(V(g))))
            depHistSub <- ggplot(filteredDep, aes_(x=~dependency)) +
                geom_histogram(binwidth=0.5,
                    aes_(fill=~..count..), col="black") +
                scale_fill_gradient("Count", low = "blue", high = "red")+
                theme_minimal(base_family = "Arial Narrow") +
                theme(axis.text=element_text(size=10),
                    axis.title=element_text(size=12))
            # Subset to those dependency scores are available
            if (!is.null(otherVar)){
                subV <- names(V(g)) %in% c(filteredDep$gene_name) |
                names(V(g)) %in% tail(colnames(pcs), n=dim(otherVar)[2])
            } else {
                subV <- names(V(g)) %in% c(filteredDep$gene_name)
            }
            depSubG <- V(g)[
                subV
            ]
            g <- igraph::induced_subgraph(g, depSubG)
            tmpSize <- vapply(names(V(g)),
            function(x) ifelse(
                x %in% filteredDep$gene_name,
                -1 * as.numeric(subset(filteredDep,
                    filteredDep$gene_name==x)$dependency),
                    NA), FUN.VALUE=1) #-1 * filteredDep$dependency
            # Size of metadata is set to mean of score
            tmpSize[is.na(tmpSize)] <- mean(tmpSize, na.rm=TRUE)
            V(g)$size <- tmpSize
            #meanExp <- apply(pcs[, names(V(g))], 2, mean)
            meanExpCol <- vapply(
                names(V(g)),
                function(x) ifelse(x %in% geneNames, mean(pcsRaw[, x]), NA),
                FUN.VALUE=1)
        } else {
            sizeLab <- "expression"
            #meanExp <- apply(pcs[, names(V(g))], 2, mean)
            meanExpCol <- vapply(names(V(g)),
                function(x) ifelse(x %in% geneNames,
                    mean(pcsRaw[, x]), NA), FUN.VALUE=1)
            meanExpSize <- meanExpCol
            meanExpSize[is.na(meanExpSize)] <- 1
            V(g)$size <- meanExpSize
        }

        ## Node color
        V(g)$color <- meanExpCol

        ## Cluster for multiple pathways
        if (length(pathNum) > 1) {
            V(g)$Pathway <- vapply(names(V(g)),
                function(x) ifelse(
                    x %in% geneNames, cls[x, ]$Pathway, "other variables"),
                FUN.VALUE="character")
        }

        ## Plot
        if (length(pathNum) == 1) {
            if (delZeroDegree){
                delG <- delete.vertices(g, igraph::degree(g)==0)
            } else {
                delG <- g
            }
            p <- ggraph(delG, layout=layout) + 
                geom_edge_diagonal(edge_alpha=1,
                                    position="identity",
                                    aes_(edge_colour=~color,
                                        width=~width, label=~label),
                                    label_size=3*(labelSize/4),
                                    label_colour=NA,
                                    angle_calc = "along",
                                    label_dodge=unit(3,'mm'),
                                    arrow=arrow(length=unit(4, 'mm')),
                                    end_cap=circle(5, 'mm'))+
                geom_node_point(aes_(color=~color, size=~size, shape=~shape),
                    show.legend=TRUE)+
                scale_color_continuous(low="blue", high="red",
                    name="expression") +
                scale_size(range=c(scaleSizeLow, scaleSizeHigh) * cexCategory,
                    name=sizeLab)+
                scale_edge_width(range=c(1, 3), guide="none")+
                scale_edge_color_continuous(low="dodgerblue",
                    high="tomato", name="strength")+
                guides(edge_color = guide_edge_colorbar(title.vjust = 3))+
                # geom_node_text(aes_(label=~name),
                # check_overlap=TRUE, repel=TRUE, size = labelSize) +
                scale_shape_identity()+
                theme_graph() + ggtitle(res[pathNum, "Description"])

            if (shadowText){
                p <- p + geom_node_text(
                            aes_(label=~stringr::str_wrap(name, width = 25)
                        ),
                    check_overlap=TRUE, repel=TRUE, size = labelSize,
                    color = textColor,
                    bg.color = bgColor, segment.color="black",
                    bg.r = .15)
            } else {
                p <- p + geom_node_text(
                            aes_(label=~stringr::str_wrap(name, width = 25)),
                    check_overlap=TRUE, repel=TRUE, size = labelSize)
            }

            if (sizeDep & !compareRef){
                if (showDepHist){
                    layoutSizedep <- "
                            ACC
                            BCC
                            "
                    p <- depHist + depHistSub + p +
                        plot_layout(design=layoutSizedep)
                }
            }

            if (compareRef){
                pathName <- res[pathNum, "Description"]
                graphiteP <- pathways(sp, pathDb)[[pathName]]
                if (!convertSymbol) {
                    graphiteP <- convertIdentifiers(graphiteP, expRow)
                } else {
                    graphiteP <- convertIdentifiers(graphiteP, "symbol")
                }

                refG <- pathwayGraph(graphiteP)
                refG <- igraph.from.graphNEL(refG, name=TRUE)
                refG <- set.vertex.attribute(refG, "name",
                    value=vapply(strsplit(names(V(refG)), ":"),
                        "[", 2, FUN.VALUE=character(1)))

                refV <- names(V(refG))
                intNodeLen <- length(intersect(refV, names(V(g))))

                refVSubset <- V(refG)[intersect( refV, names(V(g)))]
                avVSubset <- V(g)[intersect( refV, names(V(g)))]

                refSubG <- igraph::induced_subgraph(refG, refVSubset)
                avSubG <- igraph::induced_subgraph(g, avVSubset)

                refSubGUnd <- as.undirected(refSubG)
                avSubGUnd <- as.undirected(avSubG)

                difG <- difference(g, refSubG)
                difG <- delete.vertices(difG, igraph::degree(difG)==0)

                refELen <- length(E(refSubG))
                avELen <- length(E(avSubG))

                intG <- intersection(refSubG, avSubG, keep.all.vertices = FALSE)
                intG <- delete.vertices(intG, igraph::degree(intG)==0)


                if (compareRefType == "intersection"){
                    refPlot <- intG
                    t <- "Overlapping"
                } else if (compareRefType == "difference"){
                    refPlot <- difG
                    t <- "Different"
                }

                ovlELen <- length(E(refPlot))
                intP <- ggraph(refPlot, layout=layout) + 
                    geom_edge_diagonal(edge_alpha=1,
                                    position="identity",
                                    aes_(edge_colour=~color,
                                        width=~width, label=~label),
                                    label_size=3*(labelSize/4),
                                    label_colour=NA,
                                    angle_calc = "along",
                                    label_dodge=unit(3,'mm'),
                                    arrow=arrow(length=unit(4, 'mm')),
                                    end_cap=circle(5, 'mm'))+
                    geom_node_point(aes_(color=~color, size=~size,
                        shape=~shape), show.legend=TRUE)+
                    scale_color_continuous(low="blue", high="red",
                        name = "expression") +
                    scale_size(range=c(scaleSizeLow, 
                                        scaleSizeHigh) * cexCategory,
                        name=sizeLab)+
                    scale_edge_width(range=c(1, 3), guide="none")+
                    scale_edge_color_continuous(low="dodgerblue",
                        high="tomato", name="strength")+
                    guides(edge_color = guide_edge_colorbar(title.vjust = 3))+
                    # geom_node_text(aes_(label=~name), check_overlap=TRUE,
                    # repel=TRUE, size = labelSize) +
                    theme_graph() +
                    scale_shape_identity()+
                    ggtitle(paste(paste("Overlapping nodes =", intNodeLen),
                        paste("Edge number in reference network =", refELen),
                        paste("Edge number in inferred network =", avELen),
                        paste(paste0(t, " edges ="), ovlELen),
                        sep="\n"))
                if (shadowText){
                    intP <- intP + geom_node_text(
                                aes_(label=~stringr::str_wrap(name, width = 25)
                            ),
                        check_overlap=TRUE, repel=TRUE, size = labelSize,
                        color = "white",
                        bg.color = "black", segment.color="black",
                        bg.r = .15)
                } else {
                    intP <- intP + geom_node_text(
                        aes_(label=~stringr::str_wrap(name, width = 25)),
                        check_overlap=TRUE, repel=TRUE, size = labelSize)
                }

                p2 <- p + theme(legend.position="none")
                if (sizeDep & showDepHist){
                    layoutDep <- "
                            ACCDD
                            BCCDD
                            "
                    p <- depHist + depHistSub + p2 + intP +
                        plot_layout(design=layoutDep)

                } else {
                    p <- p2 + intP
                }
            }

            if (showLineage){
                lineageP <- depMeta %>%
                    dplyr::select(.data$depmap_id, .data$lineage) %>%
                    dplyr::full_join(dep, by = "depmap_id") %>%
                    dplyr::filter(.data$gene_name %in% names(V(g))) %>%
                    ggplot(aes_(x=~lineage, y=~dependency, fill=~lineage)) +
                    geom_boxplot(outlier.alpha = 0.1) +
                    theme_bw()+
                    theme(axis.text.x = element_text(angle = 45, hjust=1),
                        axis.text = element_text(size=12),
                        axis.title = element_text(size=14),
                        legend.position = "none")
                p <- p / lineageP

            }
        } else if (length(pathNum) > 1) {
            if (delZeroDegree){
                delG <- delete.vertices(g, igraph::degree(g)==0)
            } else {
                delG <- g
            }
            xy <- graphlayouts::layout_as_backbone(
                igraph::as.undirected(delG))$xy
            p <- ggraph(delG, layout="manual", x=xy[,1], y=xy[,2]) + 
                    geom_edge_diagonal(edge_alpha=1,
                                position="identity",
                                aes_(edge_colour=~color, width=~width,
                                    label=~label),
                                label_size=3*(labelSize/4),
                                label_colour=NA,
                                angle_calc = "along",
                                label_dodge=unit(3,'mm'),
                                arrow=arrow(length=unit(4, 'mm')),
                                end_cap=circle(5, 'mm'))+
                geom_node_point(aes_(color=~color, size=~size, shape=~shape),
                    show.legend=TRUE)+
                ggforce::geom_mark_hull(
                    aes_(xy[,1], xy[,2], group = ~Pathway,
                        fill = ~Pathway, label= ~Pathway),
                    concavity = 8, expand = unit(2, "mm"),
                    alpha = clusterAlpha, label.fill="transparent",
                    show.legend = FALSE, label.fontsize=12) +
                scale_fill_discrete(guide="none")+
                scale_color_continuous(low="blue", high="red",
                    name = "expression") +
                scale_size(range=c(3, 8) * cexCategory, name=sizeLab)+
                scale_edge_width(range=c(1, 3), guide="none")+
                scale_edge_color_continuous(low="dodgerblue", high="tomato",
                    name="strength")+
                guides(edge_color = guide_edge_colorbar(title.vjust = 3))+
                scale_shape_identity()+
                geom_node_text(aes_(label=~name), check_overlap=TRUE,
                    repel=TRUE, size = labelSize) +
                theme_graph()
        }

        if (strengthPlot){
            p <- p / stp + plot_layout(nrow=2, ncol=1, heights=c(0.8, 0.2))
        }

        if (returnNet){
            returnList <- list()
            returnList[["plot"]] <- p
            returnList[["str"]] <- strength
            returnList[["av"]] <- av
            returnList[["df"]] <- pcs
            return(returnList)
        } else {
            return(p)
        }
    }
}
