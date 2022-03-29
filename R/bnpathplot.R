#' bnpathplot
#'
#' Plot pathway relationship
#'
#' @param results the enrichment analysis result from clusterProfiler or ReactomePA
#' @param exp gene expression matrix
#' @param expRow the type of the identifier of rows of expression matrix
#' @param expSample candidate rows to be included in the inference, default to all
#' @param algo structure learning method used in boot.strength(), default to "hc"
#' @param algorithm.args parameters to pass to bnlearn structure learnng function
#' @param R the number of bootstrap
#' @param interactive whether to use bnviewer (default to FALSE)
#' @param color color of node, default to adjusted p-value
#' @param cexCategory scaling factor of size of nodes
#' @param delZeroDegree delete zero degree nodes
#' @param disc discretize the expressoin data
#' @param tr Specify data.frame if one needs to discretize as the same parameters as the other dataset
#' @param remainCont Specify characters when perform discretization, if some columns are to be remain continuous
#' @param cexLine scaling factor of width of edges
#' @param cl cluster object from parallel::makeCluster()
#' @param showDir show the confidence of direction of edges
#' @param chooseDir if undirected edges are present, choose direction of edges
#' @param scoreType score type to use on choosing edge direction
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
#' @param returnNet whether to return the network
#' @param otherVar other variables to be included in the inference
#' @param otherVarName the names of other variables
#' @param onlyDf return only data.frame used for inference
#' @param edgeLink whether to set edge to geom_edge_link(), FALSE to use geom_edge_diagonal()
#' @param databasePal palette to be used in scale_color_brewer, when the multiple results are to be shown
#' @param orgDb perform clusterProfiler::setReadable based on this organism database
#' @param shadowText whether to use shadow text for the better readability (default: TRUE)
#' @param bgColor color for text background when shadowText is TRUE
#' @param textColor color for text when shadowText is TRUE
#' @param seed A random seed to make the analysis reproducible, default is 1.
#' @return ggplot2 object
#'
#' @import ggraph ggplot2 patchwork igraph org.Hs.eg.db enrichplot ExperimentHub Rmpfr
#' @importFrom bnlearn score boot.strength averaged.network cextend undirected.arcs choose.direction
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by_at filter select
#' @importFrom utils read.csv
#' @importFrom reshape2 melt
#' @importFrom stringr str_starts
#' @examples
#' data("exampleEaRes");data("exampleGeneExp")
#' res <- bnpathplot(results = exampleEaRes, exp = exampleGeneExp, R = 10, expRow = "ENSEMBL")
#'
#' @export
#'
bnpathplot <- function (results, exp, expSample=NULL, algo="hc", algorithm.args=NULL,
                      expRow="ENSEMBL", cl=NULL, returnNet=FALSE, otherVar=NULL, otherVarName=NULL,
                      qvalueCutOff=0.05, adjpCutOff=0.05, nCategory=15, R=20, interactive=FALSE,
                      color="p.adjust", cexCategory=1, cexLine=0.5, chooseDir=FALSE, showDir=FALSE, delZeroDegree=TRUE,
                      labelSize=4, layout="nicely", onlyDf=FALSE, disc=FALSE, tr=NULL, remainCont=NULL,
                      shadowText=TRUE, bgColor="white", textColor="black",
                      compareRef=FALSE, strThresh=NULL, strType="normal", hub=NULL, scoreType="bic-g", databasePal="Set2",
                      dep=NULL, sizeDep=FALSE, orgDb=org.Hs.eg.db, edgeLink=TRUE, cellLineName="5637_URINARY_TRACT", strengthPlot=FALSE, nStrength=10, seed = 1)
{
    # Set type of ontology and rename
    if (length(results)>1){
        typeOfOntologies <- list()
        for (i in seq_len(length(results))){
            if (attributes(results[[i]])$class[1]=="enrichResult"){
                typeOfOntologies[[i]]=results[[i]]@ontology
            } else if (attributes(results[[i]])$class[1]=="gseaResult"){
                typeOfOntologies[[i]]=results[[i]]@setType
            }
            if (!is.null(orgDb)){
                results[[i]] <- setReadable(results[[i]], OrgDb = orgDb)
            }
        }
    } else {
        if (attributes(results)$class[1]=="enrichResult"){
            typeOfOntology=results@ontology
        } else if (attributes(results)$class[1]=="gseaResult"){
            typeOfOntology=results@setType
        }
        if (!is.null(orgDb)){
            results <- setReadable(results, OrgDb = orgDb)
        }
    }

    if (length(results)>1){
        for (i in seq_len(length(results))){
            tmpCol <- colnames(results[[i]]@result)
            ## Make comparable
            tmpCol[tmpCol=="core_enrichment"] <- "geneID"
            tmpCol[tmpCol=="qvalues"] <- "qvalue"
            tmpCol[tmpCol=="setSize"] <- "Count"
            colnames(results[[i]]@result) <- tmpCol
            if (!"enrichmentScore" %in% colnames(results[[i]]@result)){
                results[[i]]@result["enrichmentScore"] <- 0
            }
        }
    } else {
        tmpCol <- colnames(results@result)
        ## Make comparable
        tmpCol[tmpCol=="core_enrichment"] <- "geneID"
        tmpCol[tmpCol=="qvalues"] <- "qvalue"
        tmpCol[tmpCol=="setSize"] <- "Count"

        colnames(results@result) <- tmpCol
        if (!"enrichmentScore" %in% colnames(results@result)){
            results@result["enrichmentScore"] <- 0
        }
    }

    if (nCategory==0){stop("category must be > 0")}
    if (is.null(expSample)) {expSample=colnames(exp)}
    if (interactive & compareRef){stop("compareRef must be set to FALSE when use bnviewer.")}
    if (length(results) > 1 & compareRef) {stop("cannot specify compareRef when multiple databases")}
    if (length(results)==1){
        if (compareRef & typeOfOntology!="Reactome"){stop("compareRef for the pathways is currently available for reactome only.")}
    }

    ## Make dict of results
    if (length(results)>1){
        nc <- list()
        for (i in seq_len(length(results))) {
            for (p in results[[i]]@result$Description){
                nc[[p]] <- typeOfOntologies[[i]]
            }
        }
        #if (length(unique(vapply(results, function(x) gsub("kegg", "ENTREZID", x@keytype), FUN.VALUE="character")))!=1) {stop("if specify multiple databases, keytype must be same")}
    }

    if (sizeDep) {
        ## Filter to specified cell line
        if (is.null(dep)){dep = depmap::depmap_crispr()}
        filteredDep <- dep %>% filter(cell_line==cellLineName)
        depHist <- ggplot(filteredDep, aes(x=dependency)) +
            geom_histogram(aes(fill=..count..), col="black") +
            scale_fill_gradient("Count", low = "blue", high = "red") +
            theme_minimal(base_family = "Arial Narrow") +
            ggtitle(cellLineName)+
            theme(plot.title = element_text(hjust=0.5, face="bold"),
                  axis.text = element_text(size=10),
                  axis.title = element_text(size=12))
    }

    if (length(nCategory)!=length(results)){
        nCategory <- rep(nCategory, length(results))
    }

    if (length(results)>1){
        res <- c()
        for (n in seq_len(length(results))) {
            tmpres <- results[[n]]@result
            ## cutoff for q-value and adjusted p-value is same for all data
            if (!is.null(qvalueCutOff)) { tmpres <- subset(tmpres, qvalue < qvalueCutOff) }
            if (!is.null(adjpCutOff)) { tmpres <- subset(tmpres, p.adjust < adjpCutOff) }
            if (nCategory[n]) {
                # nCategory is defined per data
                tmpres <- tmpres[seq_len(nCategory[n]),]
                tmpres <- tmpres[!is.na(tmpres$ID),]
            }
            if (!is.null(colnames(res))){
                tmpres <- tmpres[,intersect(colnames(tmpres),colnames(res))]
                res <- res[,intersect(colnames(tmpres),colnames(res))]
            }
            res <- rbind(res, tmpres)
        }
    } else {
            res <- results@result
            if (!is.null(qvalueCutOff)) { res <- subset(res, qvalue < qvalueCutOff) }
            if (!is.null(adjpCutOff)) { res <- subset(res, p.adjust < adjpCutOff) }
            if (nCategory) {
                # nCategory is defined per data
                res <- res[seq_len(nCategory),]
                res <- res[!is.na(res$ID),]
            }
    }

    ## Some pathway databases have duplicate names
    if (sum(duplicated(res$Description))>0){
        dnames <- res$Description[duplicated(res$Description)]
        dupind <- which(res$Description %in% dnames)
        res[dupind, "Description"] <- paste0(res[dupind, "Description"],"_",res[dupind, "ID"])
        for (desc in res$Description){
            if (!desc %in% names(nc)){
                nc[[desc]] <- "Duplicated"
            }
        }
    }

    pcs <- c()
    pwayNames <- c()
    pathDep = c()
    for (i in seq_len(length(rownames(res)))) {
        genesInPathway <- strsplit(res[i, ]$geneID, "/")[[1]]

        if (sizeDep) {
            pathDep = c(pathDep, -1 * mean((filteredDep %>% filter(gene_name %in% genesInPathway))$dependency))
        }
        if (!is.null(orgDb)){
            genesInPathway <- suppressMessages(clusterProfiler::bitr(genesInPathway,
                                                    fromType="SYMBOL",
                                                    toType=expRow,
                                                    OrgDb=orgDb)[expRow][,1])
        }
        pathwayMatrix <- exp[ intersect(rownames(exp), genesInPathway), expSample ]
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
    
    if (sizeDep) {names(pathDep) = res$Description}
    colnames(pcs) <- pwayNames
    pcs <- data.frame(pcs, check.names=FALSE)
    if (dim(pcs)[1]==0){return("error")}


    if (!is.null(otherVar)) {
        pcs <- cbind(pcs, otherVar)
        if (!is.null(otherVarName)){
          checkColor <- otherVarName
          colnames(pcs) <- c(pwayNames, otherVarName)
        } else {
          checkColor <- colnames(otherVar)
        }
    } else {
        checkColor <- NULL
    }

    if (disc){
        pcs <- discDF(pcs, tr=tr, remainCont=remainCont)
    }


    if (onlyDf){
        return(pcs)
    }

    if (strType == "normal"){
      strength <- withr::with_seed(seed = seed ,boot.strength(pcs, algorithm=algo, algorithm.args=algorithm.args, R=R, cluster=cl))
    } else if (strType == "ms"){
      strength <- withr::with_seed(seed = seed, inferMS(pcs, algo=algo, algorithm.args=algorithm.args, R=R, cl=cl))
    }

    
    if (strengthPlot){
        strengthTop <- strength[order(strength$strength+strength$direction, decreasing = TRUE),][seq_len(nStrength),]
        strengthTop$label <- paste(strengthTop$from, "->", strengthTop$to)
        stp <- strengthTop %>% tidyr::pivot_longer(cols=c(strength, direction)) %>%
            ggplot(aes(x=label, y=value, fill=name))+geom_bar(position="dodge",stat="identity")+
            coord_flip(y=c(min(strengthTop$strength, strengthTop$direction)-0.05,1.0))+xlab("edges")+
            theme_bw()+scale_fill_manual(values = c("tomato","dodgerblue")) +
            scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 25))
    }

    if (dim(strength)[1]==0) {return("error")}
    if (!is.null(strThresh)){
        av <- averaged.network(strength, threshold=strThresh)
    } else {
        av <- averaged.network(strength)
    }

    if (chooseDir){
      av <- chooseEdgeDir(av, pcs, scoreType)
    }

    av <- cextend(av, strict=FALSE)

    if (interactive) {
        strength.viewer(bayesianNetwork = av, bayesianNetwork.boot.strength = strength)
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
          g <- set.edge.attribute(g, "label", index=E(g), round(eDf$direction,2))
        } else {
          g <- set.edge.attribute(g, "label", index=E(g), NA)
        }
        hScore <- hub.score(g, scale = TRUE, weights = E(g)$color)$vector

        if (!is.null(hub)){
          defHub <- hScore[order(hScore, decreasing=TRUE)][seq_len(hub)]
          nodeShape <- names(V(g)) %in% names(defHub)
          nodeShape <- ifelse(nodeShape, 19, 21)
          V(g)$shape <- nodeShape
        } else {
          V(g)$shape <- rep(19, length(V(g)))
        }

        ## Edge width
        if (is.null(otherVar)) {
            if (length(results)==1) {
                if (length(results@termsim) != 0) {
                    w <- results@termsim[ names(V(g)), names(V(g)) ]
                    wd <- melt(w)
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
            }## IFELSE RESULTS1
        } else {
            E(g)$width <- E(g)$color
            edgeWName <- "strength"
        }## IFELSE OTHERVAR

        ## Node size
        if (sizeDep) {
            sizeLab = "-mean dependency"
            scaleSizePathLow = 1
            scaleSizePathHigh = 10
            ## Other vars are set to 1.
            # V(g)$size <- ifelse(names(V(g)) %in% names(pathDep), pathDep, 1)
            tmpSize <- vapply(names(V(g)), function(x) ifelse(x %in% names(pathDep), as.numeric(pathDep[x]), NA), FUN.VALUE=1)
            tmpSize[is.na(tmpSize)] <- mean(tmpSize, na.rm=TRUE)
            V(g)$size <- tmpSize
        } else {
            sizeLab = "gene count"
            scaleSizePathLow = 3
            scaleSizePathHigh = 8
            tmpSize <- vapply(names(V(g)), function(x) ifelse(x %in% res$Description, as.numeric(subset(res, Description==x)["Count"]), NA), FUN.VALUE=1)
            tmpSize[is.na(tmpSize)] <- mean(tmpSize, na.rm=TRUE)
            V(g)$size <- tmpSize
            #V(g)$size <- ifelse(names(V(g)) %in% res$Description, res$Count, 3)
        }
        ## Node color
        #V(g)$color <- subset(res, res$Description %in% names(V(g)))[color][,1]
        if (length(results)>1){
            V(g)$color <- vapply(names(V(g)), function(x) ifelse(x %in% checkColor, "Metadata", nc[[x]]), FUN.VALUE="character")
        } else {
            V(g)$color <- vapply(names(V(g)), function(x) ifelse(x %in% res$Description, as.numeric(subset(res, Description==x)[color]), NA), FUN.VALUE=1)
        }

        ## Insert reference mapping
        if (length(results)==1) {
            if (typeOfOntology == "Reactome" & compareRef) {
                edgeLabel <- returnReactomeIntersection(res, g)
                E(g)$refovl = edgeLabel
            } else {
                E(g)$refovl = rep("solid", length(E(g)))
            }
        } else {E(g)$refovl = rep("solid", length(E(g)))}

        ## Plot
        if (delZeroDegree){
                delG <- delete.vertices(g, igraph::degree(g)==0)
            } else {
                delG <- g
        }
        if (edgeLink){
            p <- ggraph(delG, layout=layout) +
                geom_edge_link(edge_alpha=1,
                               position="identity",
                               angle_calc="along",
                               aes_(width=~I(width),
                                    edge_colour=~I(color),
                                    edge_linetype=~I(refovl),
                                    label=~I(label)),
                               label_dodge = unit(3, 'mm'),
                               label_colour = NA,
                               label_size = 3*(labelSize/4),
                               arrow = arrow(length=unit(4, 'mm')),
                               end_cap=circle(5, 'mm'))
        } else {
            ## If not link, diagonal is used
            p <- ggraph(delG, layout=layout) +
                 geom_edge_diagonal(edge_alpha=1,
                               position="identity",
                               angle_calc="along",
                               aes_(width=~I(width),
                                    edge_colour=~I(color),
                                    edge_linetype=~I(refovl),
                                    label=~I(label)),
                               label_dodge = unit(3, 'mm'),
                               label_colour = NA,
                               label_size = 3*(labelSize/4),
                               arrow = arrow(length=unit(4, 'mm')),
                               end_cap=circle(5, 'mm'))
        }
        p <- p + geom_node_point(aes_(color=~color, size=~size, shape=~shape))+
                scale_edge_color_continuous(low="dodgerblue", high="tomato", name="strength")+
                guides(edge_color = guide_edge_colorbar(title.vjust = 3))+
                scale_size(range=c(scaleSizePathLow, scaleSizePathHigh) * cexCategory, name=sizeLab)+
                scale_edge_width_continuous(range=c(1, 5) * cexLine, name=edgeWName) +
                scale_shape_identity()+
                theme_graph()

        ## Use shadowtext or not
        if (shadowText){
            p <- p + geom_node_text(aes_(label=~stringr::str_wrap(name, width = 25)),
                check_overlap=TRUE, repel=TRUE, size = labelSize,
                color = textColor,
                bg.color = bgColor, segment.color="black",
                bg.r = .15)
        } else {
            p <- p + geom_node_text(aes_(label=~stringr::str_wrap(name, width = 25)),
                check_overlap=TRUE, repel=TRUE, size = labelSize)
        }

        if (length(results)>1){
            p <- p + scale_colour_brewer(palette=databasePal, name="Database")
        } else {
            p <- p + scale_color_continuous(low="red", high="blue", name=color)
        }

        if (compareRef){
            p <- p + scale_edge_linetype_manual(values=c("dotted","solid"), labels=c("No", "Yes"), name="In reference")
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

