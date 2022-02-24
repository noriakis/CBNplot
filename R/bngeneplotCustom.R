#' bngeneplotCustom
#'
#' Plot gene relationship within the specified pathway using customized theme
#'
#' @param results the enrichment analysis result from clusterProfiler or ReactomePA
#' @param exp gene expression matrix
#' @param expRow the type of the identifier of rows of expression matrix
#' @param expSample candidate rows to be included in the inference, default to all
#' @param algo structure learning method used in boot.strength(), default to "hc"
#' @param algorithm.args parameters to pass to bnlearn structure learnng function
#' @param otherVar other variables to be included in the inference
#' @param otherVarName the names of other variables
#' @param R the number of bootstrap
#' @param onlyDf return only data.frame used for inference
#' @param pathNum the pathway number (the number of row of the original result, ordered by p-value)
#' @param convertSymbol whether the label of resulting network is converted to symbol, default to TRUE
#' @param cexCategory scaling factor of size of nodes
#' @param cl cluster object from parallel::makeCluster()
#' @param showDir show the confidence of direction of edges
#' @param chooseDir if undirected edges are present, choose direction of edges
#' @param labelSize the size of label of the nodes
#' @param layout ggraph layout, default to "nicely"
#' @param strType "normal" or "ms" for multiscale implementation
#' @param dep the tibble storing dependency score from library depmap
#' @param sizeDep whether to reflect DepMap score to the node size
#' @param cellLineName the cell line name to be included
#' @param strengthPlot append the barplot depicting edges with high strength
#' @param nStrength specify how many edges are included in the strength plot
#' @param strThresh the threshold for strength
#' @param hub visualize the genes with top-n hub scores
#' @param fontFamily font family name to be used for plotting
#' @param glowEdgeNum edges with top-n confidence of direction are highlighted
#' @param nodePal vector of coloring of nodes (low, high)
#' @param edgePal vector of coloring of edges (low, high)
#' @param textCol color of texts in network plot
#' @param titleCol color of title in network plot
#' @param backCol color of background in network plot
#' @param barTextCol text color in barplot
#' @param barPal bar color
#' @param barBackCol background color in barplot
#' @param barLegendKeyCol legend key color in barplot
#' @param barAxisCol axis color in barplot
#' @param barPanelGridCol panel grid color in barplot
#' @param returnNet whether to return the network
#' @param titleSize the size of title
#' @param scoreType score type to use on inference
#' @param orgDb perform clusterProfiler::setReadable based on this organism database
#' @param interactive whether to use bnviewer (default to FALSE)
#' @param disc discretize the expressoin data
#' @param tr Specify data.frame if one needs to discretize as the same parameters as the other dataset
#' @param remainCont Specify characters when perform discretization, if some columns are to be remain continuous

#' 
#' @examples
#' data("exampleEaRes");data("exampleGeneExp")
#' res <- bngeneplotCustom(results=exampleEaRes, exp=exampleGeneExp, chooseDir=TRUE, pathNum=1, glowEdgeNum=NULL, hub=3, R=40, fontFamily="sans")
#' @return ggplot2 object
#'
#' @export
#'

bngeneplotCustom <- function (results, exp, expSample=NULL, algo="hc", R=20,
                             pathNum=NULL, convertSymbol=TRUE, expRow="ENSEMBL", interactive=FALSE,
                             cexCategory=1, cl=NULL, showDir=FALSE, chooseDir=FALSE, algorithm.args=NULL,
                             labelSize=4, layout="nicely", strType="normal", returnNet=FALSE,
                             otherVar=NULL, otherVarName=NULL, onlyDf=FALSE, disc=FALSE, tr=NULL, remainCont=NULL,
                             dep=NULL, sizeDep=FALSE, orgDb=org.Hs.eg.db, cellLineName="5637_URINARY_TRACT", fontFamily="sans",
                             strengthPlot=FALSE, nStrength=10, strThresh=NULL, hub=NULL, glowEdgeNum=NULL,
                             nodePal=c("blue","red"), edgePal=c("blue","red"), textCol="black", titleCol="black", backCol="white",
                             barTextCol="black", barPal=c("red","blue"), barBackCol="white", scoreType="bic-g",
                             barLegendKeyCol="white", barAxisCol="black", barPanelGridCol="black", titleSize=24) {
    # if (is.null(hub) || is.null(glowEdgeNum) || hub <= 0 || glowEdgeNum <= 0) {stop("please specify number >= 1 for hub, glowEdgeNum")}
    if (length(nodePal)!=2 || length(edgePal)!=2 || length(barPal)!=2){stop("Please pick two colors for nodePal, edgePal and barPal.")}
    if (!is.numeric(pathNum)){stop("Please specify number(s) for pathNum.")}
    if (is.null(expSample)) {expSample=colnames(exp)}
    
    ## Deprecated
    # if (results@keytype == "kegg"){
    #     resultsGeneType <- "ENTREZID"
    # } else {
    #     resultsGeneType <- results@keytype
    # }
    results <- setReadable(results, OrgDb=orgDb)
    tmpCol <- colnames(results@result)
    tmpCol[tmpCol=="core_enrichment"] <- "geneID"
    tmpCol[tmpCol=="qvalues"] <- "qvalue"
    tmpCol[tmpCol=="setSize"] <- "Count"
    colnames(results@result) <- tmpCol

    # if (showLineage) {
    #     if (is.null(dep)){dep = depmap::depmap_crispr()}
    #     if (is.null(depMeta)){depMeta = depmap::depmap_metadata()}
    # }

    if (sizeDep) {
        if (is.null(cellLineName)){stop("Please specify cell line name.")}
        if (is.null(dep)){dep = depmap::depmap_crispr()}
        filteredDep <- dep %>% filter(cell_line==cellLineName)
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

    if (sizeDep){
        pcs <- pcs[,intersect(filteredDep$gene_name, colnames(pcs)),]
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

    if (strType == "normal"){
      strength <- boot.strength(pcs, algorithm=algo, algorithm.args=algorithm.args, R=R, cluster=cl)
    } else if (strType == "ms"){
      strength <- inferMS(pcs, algo=algo, algorithm.args=algorithm.args, R=R, cl=cl)
    }

    if (strengthPlot){
        strengthTop <- strength[order(strength$strength+strength$direction, decreasing = TRUE),][seq_len(nStrength),]
        strengthTop$label <- paste(strengthTop$from, "to", strengthTop$to)
        stp <- strengthTop %>% tidyr::pivot_longer(cols=c(strength, direction)) %>%
            ggplot(aes(x=label, y=value, fill=name))+
            geom_bar(position="dodge",stat="identity",alpha=0.7)+
            xlab("edges")+
            coord_flip(y=c(min(strengthTop$strength, strengthTop$direction)-0.05,1.0))+
            scale_fill_manual(values=c(barPal[1], barPal[2]), na.value="transparent")+
            scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 25))+
            theme(
                text=element_text(size=16,  family=fontFamily, color=barTextCol),
                plot.background = element_rect(fill = barBackCol, colour = NA),
                legend.background = ggplot2::element_blank(),
                legend.key = element_rect(fill = barLegendKeyCol),
                axis.text=element_text(color=barAxisCol),
                axis.ticks = ggplot2::element_blank(),
                axis.line = ggplot2::element_blank(),
                axis.title.y = element_text(vjust=-0.5),
                panel.grid.minor = ggplot2::element_blank(),
                panel.background = element_blank(),
                panel.grid = ggplot2::element_line(linetype = 3, color = barPanelGridCol, size = 0.2))
    }

    if (!is.null(strThresh)){
        av <- averaged.network(strength, threshold=strThresh)
    } else {
        av <- averaged.network(strength)
    }

    if (chooseDir){
        av <- chooseEdgeDir(av, pcs, scoreType)
    }

    av <- cextend(av, strict=FALSE)

    g <- bnlearn::as.igraph(av)
    e <- as_edgelist(g, names = TRUE)
    eName <-paste0(e[,1], "_", e[,2])
    colnames(e) <- c("from","to")
    eDf <- merge(e, strength)
    rownames(eDf) <- paste0(eDf$from, "_", eDf$to)
    eDf <- eDf[eName, ]
    g <- set.edge.attribute(g, "color", index=E(g), eDf$strength)
    if (showDir){
        g <- set.edge.attribute(g, "label", index=E(g), round(eDf$direction,2))
    } else {
        g <- set.edge.attribute(g, "label", index=E(g), NA)
    }

    ## Define hub genes
    hScore <- hub.score(g, scale = TRUE, weights = E(g)$color)$vector

    ## Make color for glowing edges
    if (!is.null(glowEdgeNum)){
        highEdge <- E(g)[order(E(g)$label, decreasing=TRUE)][seq_len(glowEdgeNum)]
        glowEdge <- E(g) %in% highEdge
        E(g)$glowEdge <- E(g)$color
        E(g)$glowEdge[!glowEdge] <- NA
        E(g)$alphaEdge <- rep(1, length(E(g)))
        E(g)$alphaEdge[!glowEdge] <- NA
    } else {
        E(g)$glowEdge <- NA
        E(g)$alphaEdge <- NA
    }

    E(g)$width <- E(g)$color
    edgeWName <- "strength"

    if (sizeDep){
        sizeLab <- "-dependency"
        filteredDep <- filteredDep %>% filter(gene_name %in% names(V(g))) %>%
            arrange(match(gene_name, names(V(g))))
        # Subset to those dependency scores are available
        depSubG <- V(g)[names(V(g)) %in% c(filteredDep$gene_name) || names(V(g)) %in% tail(colnames(pcs), n=dim(otherVar)[2])]
        g <- igraph::subgraph(g, depSubG)
        V(g)$size <- vapply(names(V(g)), function(x) ifelse(x %in% filteredDep$gene_name, -1 * as.numeric(subset(filteredDep, gene_name==x)$dependency), NA), FUN.VALUE=1)#-1 * filteredDep$dependency
        #meanExp <- apply(pcs[, names(V(g))], 2, mean)
        meanExpCol <- vapply(names(V(g)), function(x) ifelse(x %in% geneNames, mean(pcs[, x]), mean(apply(pcs[,geneNames], 2, mean))), FUN.VALUE=1) # colorize the other nodes later
    } else {
        sizeLab <- "expression"
        #meanExp <- apply(pcs[, names(V(g))], 2, mean)
        meanExpCol <- vapply(names(V(g)), function(x) ifelse(x %in% geneNames, mean(pcs[, x]), mean(apply(pcs[,geneNames], 2, mean))), FUN.VALUE=1) # colorize the other nodes later
        meanExpSize <- meanExpCol
        meanExpSize[is.na(meanExpSize)] <- 2
        V(g)$size <- meanExpSize
    }

    V(g)$color <- meanExpCol

    if (!is.null(hub)){
        ## Make color for glowing nodes
        defHub <- hScore[order(hScore, decreasing=TRUE)][seq_len(hub)]
        nodeShape <- names(V(g)) %in% names(defHub)

        V(g)$cyberColor <- V(g)$color
        V(g)$cyberColor[!nodeShape] <- NA

        V(g)$cyberAlpha <- rep(1, length(V(g)))
        V(g)$cyberAlpha[!nodeShape] <- NA
    } else {
        V(g)$cyberColor <- NA
        V(g)$cyberAlpha <- NA
    }

    V(g)$shape <- rep(19, length(V(g)))

    ## Cluster for multiple pathways
    if (length(pathNum) > 1) {
        V(g)$Pathway <- vapply(names(V(g)), function(x) ifelse(x %in% geneNames, cls[x, ]$Pathway, "other variables"), FUN.VALUE="character")
    }

    ## Plot
    if (length(pathNum) == 1) {
        delG <- delete.vertices(g, igraph::degree(g)==0)
        p <- ggraph(delG, layout=layout) + geom_edge_diagonal(edge_alpha=0.3,
                                                              position="identity",
                                                              aes_(edge_colour=~I(color), width=~I(width), label=~I(label)),
                                                              label_size=3,
                                                              label_colour=NA,
                                                              family=fontFamily,
                                                              angle_calc = "along",
                                                              label_dodge=unit(3,'mm'),
                                                              arrow=arrow(length=unit(4, 'mm')),
                                                              end_cap=circle(5, 'mm'))+
            geom_node_point(aes_(color=~color, size=~size, shape=~shape),
                            show.legend=TRUE, alpha=0.4)+
            scale_color_continuous(low=nodePal[1], high=nodePal[2], name="expression",
                                   na.value="transparent") +
            scale_size(range=c(scaleSizeLow, scaleSizeHigh) * cexCategory, name=sizeLab)+
            scale_edge_width(range=c(1, 3), guide="none")+
            scale_edge_color_continuous(low=edgePal[1], high=edgePal[2], name="strength",
                                        na.value="transparent")+
            guides(edge_color = guide_edge_colorbar(title.vjust = 3))+
            geom_node_text(aes_(label=~name, color=~color), check_overlap=TRUE, nudge_y=0.2, repel=TRUE, fontface="bold",
                           family=fontFamily, size = labelSize) +
            scale_shape_identity()+
            theme_graph() +
            ggtitle(res[pathNum, "Description"])+
            theme(
                text=element_text(size=16,  family=fontFamily, color=textCol),
                plot.title = element_text(size=titleSize, family=fontFamily, color = titleCol),
                panel.background = element_blank(),
                plot.background = element_rect(fill = backCol, colour = NA))

        ## Glowing phase
        layers <- 10
        size <- 1
        edgeSize <- 0.01
        glow_size <- 1.5
        glow_edge_size <- 0.5

        if (!is.null(hub)){
            for (i in seq_len(layers+1)) {
                p <- p + geom_node_point(
                    aes_(
                        color=~cyberColor,
                        alpha=~cyberAlpha
                    ),
                    fill=NA,
                    size=size+(glow_size*i))
            }
        }

        if (!is.null(glowEdgeNum)){
            for (i in seq_len(layers+1)) {
                p <- p + geom_edge_diagonal(
                    position="identity",
                    aes_(edge_colour=~I(glowEdge),
                         edge_alpha=~I(alphaEdge)),
                    width=edgeSize+(glow_edge_size*i),
                    arrow=arrow(length=unit(i*0.1, 'mm')),
                    end_cap=circle(5, 'mm')
                )
            }
        }

        p <- p + scale_alpha(range = c(0.01, 0.1), guide="none")+
            scale_edge_alpha(range=c(0.05, 0.1),guide="none")

    } else if (length(pathNum) > 1) {
        stop("The plotting across multiple pathways is currently not supported in custom.")
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
