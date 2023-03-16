#' bnpathplotCustom
#'
#' Plot pathway relationship using customized theme
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
#' @param R the number of bootstrap
#' @param color color of node, default to adjusted p-value
#' @param cexCategory scaling factor of size of nodes
#' @param cexLine scaling factor of width of edges
#' @param cl cluster object from parallel::makeCluster()
#' @param showDir show the confidence of direction of edges
#' @param chooseDir if undirected edges are present, choose direction of edges
#' @param labelSize the size of label of the nodes
#' @param layout ggraph layout, default to "nicely"
#' @param qvalueCutOff the cutoff value for qvalue
#' @param adjpCutOff the cutoff value for adjusted pvalues
#' @param nCategory the number of pathways to be included
#' @param strType "normal" or "ms" for multiscale implementation
#' @param strThresh threshold for strength, automatically determined if NULL
#' @param compareRef whether compare to the reference network between pathway
#' @param hub change the shape of node according to hub scores (default NULL)
#' @param dep the tibble storing dependency score from library depmap
#' @param sizeDep whether to reflect DepMap score to the node size
#' @param cellLineName the cell line name to be included
#' @param strengthPlot append the barplot depicting edges with high strength
#' @param nStrength specify how many edges are included in the strength plot
#' @param fontFamily font family name to be used for plotting
#' @param glowEdgeNum edges with top-n confidence of direction are highlighted
#' @param nodePal vector of coloring of nodes (low, high)
#' @param edgePal vector of coloring of edges (low, high)
#' @param textCol color of texts in network plot
#' @param backCol color of background in network plot
#' @param barTextCol text color in barplot
#' @param barPal bar color
#' @param edgeLink use geom_edge_link() instead of geom_edge_diagonal()
#' @param barBackCol background color in barplot
#' @param barLegendKeyCol legend key color in barplot
#' @param barAxisCol axis color in barplot
#' @param barPanelGridCol panel grid color in barplot
#' @param returnNet whether to return the network
#' @param scoreType score type to use on inference
#' @param disc discretize the expressoin data
#' @param tr Specify data.frame if one needs to discretize
#'           as the same parameters as the other dataset
#' @param remainCont Specify characters when perform discretization,
#'                   if some columns are to be remain continuous
#' @param otherVar other variables to be included in the inference
#' @param otherVarName the names of other variables
#' @param onlyDf return only data.frame used for inference
#' @param orgDb perform clusterProfiler::setReadable
#'              based on this organism database
#' @param bypassConverting bypass the symbol converting
#'                         ID of rownames and those listed in EA result
#'                         must be same
#' @param seed A random seed to make the analysis reproducible, default is 1.
#' @examples 
#' data("exampleEaRes");data("exampleGeneExp")
#' res <- bnpathplotCustom(results=exampleEaRes, exp=exampleGeneExp,
#'                         fontFamily="sans", glowEdgeNum=3, hub=3)
#' @return ggplot2 object
#'
#' @export
#'
bnpathplotCustom <- function (results, exp, expSample=NULL, algo="hc",
                            R=20, expRow="ENSEMBL", color="p.adjust",
                            cexCategory=1, cl=NULL, showDir=FALSE,
                            chooseDir=FALSE, labelSize=4, layout="nicely",
                            strType="normal", compareRef=FALSE, disc=FALSE,
                            tr=NULL, remainCont=NULL, qvalueCutOff=0.05,
                            adjpCutOff=0.05, nCategory=15, cexLine=1,
                            returnNet=FALSE, dep=NULL, sizeDep=FALSE,
                            cellLineName="5637_URINARY_TRACT",
                            fontFamily="sans", otherVar=NULL, otherVarName=NULL,
                            onlyDf=FALSE, algorithm.args=NULL,
                            strengthPlot=FALSE, nStrength=10, edgeLink=FALSE,
                            strThresh=NULL, hub=NULL, glowEdgeNum=NULL,
                            nodePal=c("blue","red"), edgePal=c("blue","red"),
                            textCol="black", backCol="white",
                            barTextCol="black", barPal=c("red","blue"),
                            barBackCol="white", scoreType="bic-g",
                            barLegendKeyCol="white", orgDb=org.Hs.eg.db,
                            barAxisCol="black", barPanelGridCol="black",
                            seed = 1, bypassConverting=FALSE) {
    # if (is.null(hub) ||
    # is.null(glowEdgeNum) ||
    # hub <= 0 || glowEdgeNum <= 0) {
    #     stop("please specify number >= 1 for hub, glowEdgeNum")}
    if (length(nodePal)!=2 ||
        length(edgePal)!=2 ||
        length(barPal)!=2){
        stop("Please pick two colors for nodePal, edgePal and barPal.")}
    if (compareRef){stop("compareRef is currently not supported.")}
    if (is.null(expSample)) {expSample <- colnames(exp)}
    if (length(results)>1) {
        stop("For multiple databases, please use bnpathplot().")}
    # if (results@keytype == "kegg"){
    #     resultsGeneType <- "ENTREZID"
    # } else {
    #     resultsGeneType <- results@keytype
    # }
    attributes(results)$result$Description <- gsub("Homo sapiens\r: ",
                                    "",
                                    attributes(results)$result$Description)
    if (attributes(results)$class[1]=="enrichResult"){
        typeOfOntology <- attributes(results)$ontology
    } else if (attributes(results)$class[1]=="gseaResult"){
        typeOfOntology <- attributes(results)$setType
    }
    if (!bypassConverting) {
        if (!is.null(orgDb)){
            results <- setReadable(results, OrgDb = orgDb)
        }
    }

    tmpCol <- colnames(attributes(results)$result)
    ## Make comparable
    tmpCol[tmpCol=="core_enrichment"] <- "geneID"
    tmpCol[tmpCol=="qvalues"] <- "qvalue"
    tmpCol[tmpCol=="setSize"] <- "Count"

    colnames(attributes(results)$result) <- tmpCol
    if (!"enrichmentScore" %in% colnames(attributes(results)$result)){
        attributes(results)$result["enrichmentScore"] <- 0
    }

    if (sizeDep) {
        if (is.null(dep)){dep <- depmap::depmap_crispr()}
        filteredDep <- dep %>% filter(.data$cell_line==cellLineName)
    }

    res <- attributes(results)$result

    pcs <- c()
    pwayNames <- c()

    ## Filtering according to qval, pval and nCategory
    if (!is.null(qvalueCutOff)) {
        res <- subset(res, res$qvalue < qvalueCutOff) }
    if (!is.null(adjpCutOff)) {
        res <- subset(res, res$p.adjust < adjpCutOff) }
    if (nCategory) {
        res <- res[seq_len(nCategory),]
        res <- res[!is.na(res$ID),]
    }

    pathDep <- c()
    for (i in seq_len(length(rownames(res)))) {
        genesInPathway <- strsplit(res[i, ]$geneID, "/")[[1]]

        if (sizeDep){
            pathDep <- c(pathDep,
                -1 * mean((filteredDep %>%
                    filter(.data$gene_name %in% genesInPathway))$dependency))
        }
        if (!bypassConverting) {
            if (!is.null(orgDb)){
                genesInPathway <- clusterProfiler::bitr(genesInPathway,
                                                        fromType="SYMBOL",
                                                        toType=expRow,
                                                        OrgDb=orgDb)[expRow][,1]
            }
        }
        pathwayMatrix <- exp[ intersect(rownames(exp),
                                        genesInPathway), expSample ]
        if (dim(pathwayMatrix)[1]==0) {
            message("no gene in the pathway present in expression data")
        } else {
            pathwayMatrixPca <- prcomp(t(pathwayMatrix), scale. = FALSE)$x[,1]
            avExp <- apply(pathwayMatrix, 2, mean)
            corFlag <- cor(pathwayMatrixPca, avExp)
            if (corFlag < 0){pathwayMatrixPca <- pathwayMatrixPca*-1}
            pwayNames <- c(pwayNames, res[i,]$Description)
            pcs <- cbind(pcs, pathwayMatrixPca)
            # pathwayMatrixSum <- apply(pathwayMatrix, 2, sum)
            # pwayNames <- c(pwayNames, res[i,]$Description)
            # pcs <- cbind(pcs, pathwayMatrixSum)
        }
    }

    if (sizeDep) {names(pathDep) <- res$Description}
    colnames(pcs) <- pwayNames

    pcs <- data.frame(pcs, check.names=FALSE)
    if (dim(pcs)[1]==0){return("error")}
    
    if (!is.null(otherVar)) {
        pcs <- cbind(pcs, otherVar)
        if (!is.null(otherVarName)){
            colnames(pcs) <- c(pwayNames, otherVarName)
        }
    }

    if (disc){
        pcs <- discDF(pcs, tr=tr, remainCont=remainCont)
    }

    if (onlyDf){
        return(pcs)
    }

    if (strType == "normal"){
        strength <- withr::with_seed(seed = seed, boot.strength(pcs,
            algorithm=algo, algorithm.args=algorithm.args, R=R, cluster=cl))
    } else if (strType == "ms"){
        strength <- withr::with_seed(seed = seed,
            inferMS(pcs, algo=algo, algorithm.args=algorithm.args, R=R, cl=cl))
    }

    if (strengthPlot){
        strengthTop <- strength[order(strength$strength+strength$direction,
            decreasing = TRUE),][seq_len(nStrength),]
        strengthTop$label <- paste(strengthTop$from, "to", strengthTop$to)
        stp <- strengthTop %>%
            tidyr::pivot_longer(cols=c(.data$strength, .data$direction)) %>%
            ggplot(aes(x=.data$label, y=.data$value, fill=.data$name))+
            geom_bar(position="dodge",stat="identity",alpha=0.7)+
            xlab("edges")+
            coord_flip(ylim=c(min(strengthTop$strength,
                strengthTop$direction)-0.05,1.0))+
            scale_fill_manual(values=c(barPal[1], barPal[2]),
                na.value="transparent")+
            scale_x_discrete(labels = function(x) stringr::str_wrap(x,
                                                                width = 25))+
            theme(
                text=element_text(size=16,  family=fontFamily,
                    color=barTextCol),
                plot.background = element_rect(fill = barBackCol, colour = NA),
                legend.background = ggplot2::element_blank(),
                legend.key = element_rect(fill = barLegendKeyCol),
                axis.text=element_text(color=barAxisCol),
                axis.ticks = ggplot2::element_blank(),
                axis.line = ggplot2::element_blank(),
                axis.title.y = element_text(vjust=-0.5),
                panel.grid.minor = ggplot2::element_blank(),
                panel.background = element_blank(),
                panel.grid = ggplot2::element_line(linetype = 3,
                    color = barPanelGridCol, size = 0.2))
    }

    if (!is.null(strThresh)){
        av <- averaged.network(strength, threshold=strThresh)
    } else {
        av <- averaged.network(strength)
    }

    # if (chooseDir){
    #     av <- chooseEdgeDir(av, pcs, scoreType)
    # }

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

    hScore <- hub.score(g, scale = TRUE, weights = E(g)$color)$vector


    ## Make color for glowing edges
    if (!is.null(glowEdgeNum)){
        highEdge <- E(g)[order(E(g)$label,
                                decreasing=TRUE)][seq_len(glowEdgeNum)]
        glowEdge <- E(g) %in% highEdge
        E(g)$glowEdge <- E(g)$color
        E(g)$glowEdge[!glowEdge] <- NA
        E(g)$alphaEdge <- rep(1, length(E(g)))
        E(g)$alphaEdge[!glowEdge] <- NA
    } else {
        E(g)$glowEdge <- NA
        E(g)$alphaEdge <- NA
    }


    ## Edge width
    if (is.null(otherVar)) {
        if (length(attributes(results)$termsim) != 0) {
            w <- attributes(results)$termsim[ names(V(g)), names(V(g)) ]
            wd <- reshape2::melt(w)
            wd <- wd[wd[,1] != wd[,2],]
            wd <- wd[!is.na(wd[,3]),]

            rownames(wd) <- paste0(wd[,1], "_", wd[,2])
            wd <- wd[ eName, ]
            rawWidth <- wd[,3]
            if (showDir){
                rawWidth[is.na(rawWidth)] <- 0
            }
            E(g)$width <- sqrt(rawWidth * 5) * cexLine
            edgeWName <- "similarity"
        } else {
            E(g)$width <- E(g)$color
            edgeWName <- "strength"
        }
    } else {
        E(g)$width <- E(g)$color
        edgeWName <- "strength"
    }## IFELSE OTHERVAR


    if (sizeDep) {
        sizeLab <- "-mean dependency"
        scaleSizePathLow <- 1
        scaleSizePathHigh <- 10
        V(g)$size <- vapply(names(V(g)),
            function(x) ifelse(x %in% names(pathDep),
                as.numeric(pathDep[x]), 1), FUN.VALUE=1)
    } else {
        sizeLab <- "gene count"
        scaleSizePathLow <- 3
        scaleSizePathHigh <- 8
        V(g)$size <- vapply(names(V(g)),
            function(x) ifelse(x %in% res$Description,
                as.numeric(subset(res, res$Description==x)["Count"]), 3),
            FUN.VALUE=1)
    }


    V(g)$color <- vapply(names(V(g)),
        function(x) ifelse(x %in% res$Description,
            as.numeric(subset(res, res$Description==x)[color]),
            as.numeric(apply(res[color], 2, mean))), FUN.VALUE=1)


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


    ## Plot
        delG <- delete.vertices(g, igraph::degree(g)==0)
        if (edgeLink) {
            p <- ggraph(delG, layout=layout) + 
                    geom_edge_link(edge_alpha=0.3,
                        position="identity",
                        aes_(edge_colour=~color, width=~width,
                            label=~label),
                        label_size=3,
                        label_colour=NA,
                        family=fontFamily,
                        angle_calc = "along",
                        label_dodge=unit(3,'mm'),
                        arrow=arrow(length=unit(4, 'mm')),
                        end_cap=circle(5, 'mm'))
            } else {
            p <- ggraph(delG, layout=layout) + 
                    geom_edge_diagonal(edge_alpha=0.3,
                        position="identity",
                        aes_(edge_colour=~color, width=~width,
                            label=~label),
                        label_size=3,
                        label_colour=NA,
                        family=fontFamily,
                        angle_calc = "along",
                        label_dodge=unit(3,'mm'),
                        arrow=arrow(length=unit(4, 'mm')),
                        end_cap=circle(5, 'mm'))                
            }
            p <- p + geom_node_point(aes_(color=~color, size=~size, shape=~shape),
                            show.legend=TRUE, alpha=0.4)+
            scale_color_continuous(low=nodePal[1], high=nodePal[2], name=color,
                                    na.value="transparent") +
            scale_size(range=c(scaleSizePathLow,
                scaleSizePathHigh) * cexCategory, name=sizeLab)+
            scale_edge_width(range=c(1, 3), guide="none")+
            scale_edge_color_continuous(low=edgePal[1], high=edgePal[2],
                                        name="strength",
                                        na.value="transparent")+
            guides(edge_color = guide_edge_colorbar(title.vjust = 3))+
            geom_node_text(aes_(label=~stringr::str_wrap(name, width = 25),
                        color=~color), check_overlap=TRUE,
                        nudge_y=0.2, repel=TRUE, fontface="bold",
                        family=fontFamily, size = labelSize) +
            scale_shape_identity()+
            theme_graph() +
            theme(
                text=element_text(size=16,  family=fontFamily, color=textCol),
                # plot.title = element_text(size=24,
                # family=fontFamily, color = titleCol),
                panel.background = element_blank(),
                plot.background = element_rect(fill = backCol, colour = NA))

    ## Glowing phase
    layers <- 10
    size <- 8
    edgeSize <- 0.01
    glow_size <- 1.5
    glow_edge_size <- 1.1

    if (!is.null(hub)){
        for (i in seq_len(layers+1)){
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
        for (i in seq_len(layers+1)){
            if (edgeLink) {
                p <- p + geom_edge_link(
                        position="identity",
                        aes_(edge_colour=~glowEdge,
                            edge_alpha=~alphaEdge),
                        width=edgeSize+(glow_edge_size*i),
                        arrow=arrow(length=unit(i*0.1, 'mm')),
                        end_cap=circle(5, 'mm')
                    )
            } else {
                p <- p + geom_edge_diagonal(
                        position="identity",
                        aes_(edge_colour=~glowEdge,
                            edge_alpha=~alphaEdge),
                        width=edgeSize+(glow_edge_size*i),
                        arrow=arrow(length=unit(i*0.1, 'mm')),
                        end_cap=circle(5, 'mm')
                    )                
            }
        }
    }

    p <- p+scale_alpha(range = c(0.01, 0.1), guide="none")+
        scale_edge_alpha(range=c(0.01, 0.1),guide="none")

    if (strengthPlot){
        p <- p / stp + plot_layout(nrow=2, ncol=1, heights=c(0.6, 0.4))
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
