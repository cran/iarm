#' Item Characteristic Curves
#'
#'Plots Item Characteristic Curves for dichotomous and polytomous items. The plot can display observed scores as
#'total scores (method="score") or as average scores within adjacent class intervals (method="cut").
#'Class intervals can be useful when the sample size is not large enough to contain an adequate
#'number of respondents with the same total score for each possible total score. The function includes
#'the option to plot observed scores according to values of an exogenous variable to evaluate differential item functioning
#'(dif="yes").
#'
#' @param data An object of class "data.frame" containing the items (include all items present in the model).
#' The variables need to be numeric.
#' @param itemnumber A numeric vector indicating the columns of the data (the items) which ICCs are going to be plotted.
#' Maximum of four items per plot.
#' @param pallete An object of class "character". Choose a pre-made color pallete from package RColorBrewer.
#' Only available for dif="no".
#' @param xticks A numeric scalar. Specify x-axis tick values.
#' @param yticks A numeric scalar. Specify y-axis tick values.
#' @param thetain A numeric scalar. Specify minimum theta values for person parameters.
#' @param thetaend A numeric scalar. Specify maximum theta values for person parameters.
#' @param method The method for displaying observed scores. Choose "score" to plot total scores.
#' Choose "cut" to plot class intervals.
#' @param grid Chooses whether the background grid should be displayed. Options are "yes" or "no".
#' @param cinumber A numeric scalar. The number of adjacent class intervals in which participants will be divided.
#' Notice that the number of class intervals cannot be higher than the number of total scores.
#' @param itemdescrip A character vector indicating the description of the plotted items. Maximum of four descriptions
#' (one description per item plotted).
#' @param axis.rumm Configures whether the plot should display the entire trait range or solely the trait range close
#' to the observed scores (similar to private software RUMM2030). Options are "yes" or "no".
#' @param dif Configures whether the observed scores will be plotted according to values of an exogenous variable
#' to evaluate differential item function. Options are "yes" or "no".
#' @param difvar Chooses the variable which will be used to evaluate differential item functioning. Only
#' necessary when dif="yes".
#' @param diflabels A character vector indicating the labels to values of the variable chosen to evaluate differential item functioning.
#' Only necessary when dif="yes".
#' @param difstats Displays the partial gamma coefficient to indicate the magnitude of differential item
#' functioning. Options are "yes" or "no". Only necessary when dif="yes".
#' @param title A character vector. The title of the plot.
#' @param icclabel Displays the labels of Expected Item Score and Observed Item Score. Options are "yes"
#' or "no".
#' @param xaxistitle A character vector. The x-axis title.
#' @param yaxistitle A character vector. The y-axis title.
#' @importFrom psychotools pcmodel raschmodel personpar
#' @importFrom Hmisc cut2
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom stats aggregate predict
#' @importFrom utils capture.output
#' @export
#' @author Pedro Henrique Ribeiro Santiago \email{pedro.ribeirosantiago@@adelaide.edu.au}, Marianne Mueller
#' @examples
#' \dontrun{# Creates a plot for Item 1 using total scores
#' ICCplot(desc2[,5:13], itemnumber=1, method="score", itemdescrip="Item 1")
#'
#' # Creates a plot for Item 1 using 8 class intervals
#' ICCplot(desc2[,5:13], itemnumber=1, method="cut", cinumber=8, itemdescrip="Item 1")
#'
#' # Creates a plot for Item 1 using 8 class intervals without RUMM style axis
#' ICCplot(desc2[,5:13], itemnumber=1, method="cut", cinumber=8, itemdescrip="Item 1", axis.rumm="no")
#'
#' # Creates a plot for Item 3 using 8 class intervals and evaluating DIF according to gender
#' ICCplot(desc2[,5:13], itemnumber=3, method="cut", cinumber=8, itemdescrip="Item 3",
#' dif="yes", difvar=desc2$gender, diflabels=c("Men", "Women"))
#'
#' # Creates a plot with three items using 5 class intervals and evaluating DIF according to gender
#' ICCplot(desc2[,5:13], itemnumber=1:3, method="cut", cinumber=5,
#' itemdescrip=c("Item 1","Item 2","Item 3"), dif="yes"
#' difvar=desc2$gender, diflabels=c("Men", "Women"))
#'}
ICCplot <- function(data, itemnumber, pallete='Paired', xticks=1.0, yticks=0.5,
                    thetain=-6.000, thetaend=6.000, method="score", grid="yes", cinumber=6, itemdescrip="",
                    axis.rumm="yes", dif="no", difvar=NA, diflabels=c("Group1", "Group 2", "Group 3", "Group 4", "Group5"),
                    difstats="yes", title="Item Characteristic Curve", icclabel="yes",
                    xaxistitle="Theta", yaxistitle="Item Score") {

  pltC <- function() {

    if (grid=="no") {
      background <- element_blank()
      gridy <- element_blank()
      gridx <- element_blank()
      panelgrid <- element_blank()
    } else {
      background <- element_rect(fill = "white", colour="black")
      gridy <- element_blank()
      gridx <- element_blank()
      panelgrid <- element_line(colour="grey87", size=0.25)
    }

    if (icclabel=="yes") {
      icclabels <- "bottom"
    } else {icclabels <- "none"}
    xpos <- ypos <- annotateText <- hjustvar <- vjustvar <- NULL
    annotations <- data.frame(
      xpos = c(-Inf,-Inf,Inf,Inf),
      ypos =  c(-Inf, Inf,-Inf,Inf),
      annotateText = itmn,
      hjustvar = c(-0.5) ,
      vjustvar = c(3.0))

    if (any(lapply(data,class)=="character")==TRUE) {
      stop(' "Input variables must be numeric" ')
    } else
      if (any(lapply(data,class)=="factor")==TRUE) {
        stop(' "Input variables must be numeric" ')
      }
    else {}

    seqres <- NA
    for (i in 1:ncol(data)) {
      seqres[i] <- all(abs(diff(sort(unique(data[,i])))) == 1)
      seqres
    }
    if (any(seqres=="FALSE")==TRUE) {
      stop(' "You need to provide the number of responses for all items categories.
              If there were zero responses to one category, please include this information in the data" ')
    } else {}

    maxr <- max(data, na.rm=TRUE)
    minr <- min(data, na.rm=TRUE)
    if(minr>0) {data2 <- (data - minr)} #Checks whether the minimum value is 0
    if(minr==0) {data2 <- data}
    else {}
    adat <- as.data.frame(data2)
    adat <- adat[complete.cases(adat),]

    cati <- NA
    for (i in 1:ncol(adat)) {
      cati[i] <- max(adat[,i])
      mxsc <- sum(cati)
    } #Calculates the number of categories

    ttsc <- rowSums(adat, na.rm=TRUE)
    if (ncol(adat)>2) {
      rtsc <- rowSums(adat[,-itmc], na.rm=TRUE)
    } else {
      rtsc <- adat[,-itmc]
    }
    adat <- cbind(adat, ttsc, rtsc)
    adat <- adat[order(ttsc),]
    itms <- ncol(adat)-2
    pop <- seq(from=thetain, to=thetaend, by=0.01)
    if ((maxr-minr)>1) {
      pcmo <- pcmodel(adat[,1:itms])
      prcu <- predict(pcmo, newdata=pop, type="cumprobability")
    } else {
      ramo <- raschmodel(adat[,1:itms])
      prcu <- predict(ramo, newdata=pop, type="cumprobability")
    }
    type <- 1
    cure <- cbind(prcu, pop, type)
    cure <- as.data.frame(cure)
    cat <- table(adat[,itmc])
    ncat <- length(cat)
    pron=c()
    for (i in 1:(ncol(cure)-2)) {
      pron[i]=(sum(cure[,i])==nrow(cure))
    }
    rcat <- which(pron == TRUE)
    nrem <- length(rcat)
    ini <- rcat[itmc] + 1
    fin <- rcat[itmc] + (ncat-1)
    mxcl <- ncat*itmc
    mncl <- (ncat*(itmc-1))+2
    if (ncat>2) {
      cend <- cbind(rowSums(cure[,ini:fin]), cure$pop, cure$type, rowSums(cure[,1:(length(cure)-2)])-nrem)
    } else {
      cend <- cbind(cure[,mxcl], cure$pop, cure$type, rowSums(cure[,1:(length(cure)-2)])-nrem)
    }
    cend <- as.data.frame(cend)
    names(cend) <- c("score", "theta", "type", "ttsc")
    if ((maxr-minr)>1) {
      ppar <- person_estimates(pcmo)
      psav <- ppar
      ppar <- ppar[2:(nrow(ppar)-1),3]
    } else {
      ppar <- personpar(ramo)
    }
    mxpp <- length(ppar)
    ppar <- ppar[1:mxpp]

    if (method=="cut"&cinumber >= max(adat$ttsc)) {
      stop('Number of class intervals need to be smaller than the total scores')
    } else {}

    if (method=="cut"){
      pcrv <- adat
      pmsc <- which(pcrv$ttsc==0)
      if (length(pmsc)==0) {
        pcrv <- pcrv
      } else {
        pcrv <- pcrv[-which(pcrv$ttsc==0),]
      }
      pxsc <- which(pcrv$ttsc==mxsc)
      if (length(pxsc)==0) {
        pcrv <- pcrv
      } else {
        pcrv <- pcrv[-which(pcrv$ttsc==mxsc),]
      }
      pcrv$class <- cut2(pcrv$ttsc, g=cinumber, oneval=TRUE, levels.mean = TRUE)
    } else {}

    if (method=="cut")
    {if (cinumber > length(unique(pcrv$class)))
    { stop('There are not enough subjects in each total score to produce this number of class intervals')}
      else{}
    }
    else {}

    if(method=="cut" & dif!="yes") {
      obs <- aggregate(pcrv[,itmc], by=list(pcrv$class), FUN=mean)
      obs[,1] <- as.numeric(levels(obs[,1]))
      type <- 2
      x <- rep(NA, cinumber)
      for (i in 1:cinumber) {
        x[i] = cend$theta[which(abs(cend$ttsc-obs[i,1])==min(abs(cend$ttsc-obs[i,1])))]
      }
      x <- x[1:cinumber]
      obs <- cbind(obs[,2], x, type)
      obs <- as.data.frame(obs)
      names(obs) <- c("score", "theta", "type")
      cend <- cend[,1:3]
      cend <- as.data.frame(cend)
      names(cend) <- c("score", "theta", "type")
      cend <- rbind(cend, obs)
      cend[,1:3] <- as.numeric(unlist(cend[,1:3]))

      if(axis.rumm=="yes") {
        classbreak=c(
          max(min(subset(cend$theta, cend$type==1)),
              min((subset(cend$theta, cend$type!=1))-
                    1.0)),
          min(max((subset(cend$theta, cend$type==1))),
              max((subset(cend$theta, cend$type!=1)))+
                1.0))
      } else {
        classbreak=c(NA,NA)
      }

      theta <- NULL
      score <- NULL
      myICCplot <- ggplot(cend, aes(x=theta, y=score, col=as.factor(type))) + #The graph itself
        ggtitle(title) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5), legend.position=icclabels,
              panel.grid = panelgrid,
              panel.grid.major.y = gridy,
              panel.grid.major.x = gridx,
              panel.border = element_rect(colour="black", size=0.25, fill=NA),
              panel.background = background) +
        scale_x_continuous(breaks = round(seq(min(cend$theta), max(cend$theta), by = xticks),1),
                           limits=classbreak) +
        scale_y_continuous(breaks = round(seq(min(0), max(cend$score)+0.5, by = yticks),1)) +
        scale_color_brewer(palette = pallete, name="", labels=c("Expected Item Score", "Average Observed Item Score")) +
        labs(y = yaxistitle, x = xaxistitle) +
        geom_point(na.rm=TRUE) +
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="black")
    }

    else if (method!="cut"&dif!="yes") {

      obs <- aggregate(adat[,itmc], by=list(adat$ttsc), FUN=mean)
      zero <- which(obs[,1]==0)
      if (length(zero)==0) { obs <- obs}
      else {obs <- obs[-which(obs[,1]==0),] }
      ext <- which(obs==mxsc)
      if (length(ext)==0) { obs <- obs}
      else {obs <- obs[-which(obs==mxsc),] }
      rest <- aggregate(adat[,itmc], by=list(adat$rtsc), FUN=mean)
      type <- 2
      tlsc <- as.vector(unlist(obs[1]))
      names(obs) <- c("score", "obsmean")

      if ((maxr-minr)>1) {
        pseq <- seq(from=psav[2,1], to=psav[(nrow(psav)-1),1], by=1)
      }
      else {
        pseq <- as.numeric(unlist(regmatches(capture.output(ppar)
                                             [c(seq(1, length(capture.output(ppar)), 2))], gregexpr("[[:digit:]]+",
                                                                                                    capture.output(ppar)[c(seq(1, length(capture.output(ppar)), 2))]))))
      }
      ppnw <- cbind(ppar, pseq)
      ppnw <- as.data.frame(ppnw)
      names(ppnw) <- c("theta", "score")
      ptmp <- ppnw[ppnw$score %in% tlsc,]
      ptmp <- as.data.frame(ptmp)
      names(ptmp) <- c("theta", "score")
      obs <- merge(obs, ptmp, by="score")
      obs <- cbind(obs[,2:3], type)
      obs <- as.data.frame(obs)
      names(obs) <- c("score", "theta", "type")
      cend <- cend[,1:3]
      cend <- rbind(cend, obs)

      if(axis.rumm=="yes") {
        classbreak=c(
          max(min(subset(cend$theta, cend$type==1)),
              min((subset(cend$theta, cend$type!=1))-
                    1.0)),
          min(max((subset(cend$theta, cend$type==1))),
              max((subset(cend$theta, cend$type!=1)))+
                1.0))
      }
      else {classbreak=c(NA,NA)}

      myICCplot <- ggplot(cend, aes(x=theta, y=score, col=as.factor(type))) +
        ggtitle(title) + #Choose title
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5), legend.position=icclabels,
              panel.grid = panelgrid,
              panel.grid.major.y = gridy,
              panel.grid.major.x = gridx,
              panel.border = element_rect(colour="black", size=0.25, fill=NA),
              panel.background = background) +
        scale_x_continuous(breaks = round(seq(min(cend$theta), max(cend$theta), by = xticks),1),
                           limits=classbreak) +
        scale_y_continuous(breaks = round(seq(min(0), max(cend$score)+0.5, by = yticks),1)) +
        scale_color_brewer(palette = pallete, name="", labels=c("Expected Item Score", "Average Observed Item Score")) +
        labs(y = yaxistitle, x = xaxistitle) +
        geom_point(na.rm=TRUE) +
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="black")
    }

    else if (method=="cut"&dif=="yes") {

      difd <- cbind(data2, difvar)
      difd <- difd[complete.cases(difd),]
      difd[,1:ncol(difd)] <- as.numeric(unlist(difd[,1:ncol(difd)]))
      difs <- min(difd$difvar)
      if (difs>1) {difd$difvar <- ((difd$difvar)-(difs-1))}
      else if (difs==1) {}
      else if (difs<1) {difd$difvar <- ((difd$difvar)+1)}
      else {}
      ldif <- list()
      ttsc <- list()
      obs <- list()
      zero <- list()
      ext <- list()
      class <- list()
      x <- list()
      for (i in 1:length(unique(difd$difvar))){
        x[[i]]<- vector()
      }

      for (i in 1:length(unique(difd$difvar))) {
        ldif[[i]] <- subset(difd, difd$difvar==i)
        ldif[[i]] <- ldif[[i]][complete.cases(ldif[[i]]),]
        ttsc[[i]] <- rowSums(ldif[[i]][,1:(ncol(ldif[[i]])-1)], na.rm=TRUE)
        ldif[[i]] <- cbind(ldif[[i]][,1:(ncol(ldif[[i]])-1)], ttsc[[i]])
        ldif[[i]] <- ldif[[i]][order(ldif[[i]][,ncol(ldif[[i]])]),]

        zero[[i]] <- which(ldif[[i]]$`ttsc[[i]]`==0)
        if (length(zero[[i]])==0) { ldif[[i]] <- ldif[[i]]}
        else {ldif[[i]] <- ldif[[i]][-which(ldif[[i]]$`ttsc[[i]]`==0),] }
        ext[[i]] <- which(ldif[[i]]==mxsc)
        if (length(ext[[i]])==0) { ldif[[i]] <- ldif[[i]]}
        else {ldif[[i]] <- ldif[[i]][-which(ldif[[i]]$`ttsc[[i]]`==mxsc),] }
        ldif[[i]]$class <- cut2(ldif[[i]][,ncol(ldif[[i]])], g=cinumber, oneval=TRUE, levels.mean = TRUE)
        obs[[i]] <- aggregate(ldif[[i]][,itmc], by=list(ldif[[i]][,ncol(ldif[[i]])]), FUN=mean)
        type[i] <- i+1
        obs[[i]][,1] <- as.numeric(levels(obs[[i]][,1]))
      }

      for (i in 1:length(unique(difd$difvar))) {
        for (j in 1:nrow(obs[[i]])) {
          x[[i]][j] = cend$theta[which(abs(cend$ttsc-obs[[i]][j,1])==min(abs(cend$ttsc-obs[[i]][j,1])))]
        }}

      for (i in 1:length(unique(difd$difvar))) {
        obs[[i]] <- cbind(obs[[i]][,2], x[[i]], type[i])
        obs[[i]] <- as.data.frame(obs[[i]])
        names(obs[[i]]) <- c("score", "theta", "type")}

      cend <- cend[,1:3]
      cend <- as.data.frame(cend)
      names(cend) <- c("score", "theta", "type")
      big_data = do.call(rbind, obs)
      cend <- rbind(cend, big_data)

      allg <- partgam_DIF(data2, difvar)
      pgmm <- allg[itmc,3]
      pgmm <- format(round(pgmm, digits=2), nsmall=2)

      if(difstats=="yes") {

        annotateDIF <- data.frame(
          xpos = c(-Inf,-Inf,Inf,Inf),
          ypos =  c(-Inf, Inf,-Inf,Inf),
          annotateText = paste("gamma == ", pgmm),
          hjustvar = c(-0.35) ,
          vjustvar = c(4.0)) }

      else {annotateDIF <- data.frame(
        xpos = c(-Inf,-Inf,Inf,Inf),
        ypos =  c(-Inf, Inf,-Inf,Inf),
        annotateText = c(""),
        hjustvar = c(-0.35) ,
        vjustvar = c(4.0)) }

      if(axis.rumm=="yes") {
        classbreak=c(
          max(min(subset(cend$theta, cend$type==1)),
              min((subset(cend$theta, cend$type!=1))-
                    1.0)),
          min(max((subset(cend$theta, cend$type==1))),
              max((subset(cend$theta, cend$type!=1)))+
                1.0))
      }
      else {classbreak=c(NA,NA)}

      myICCplot <- ggplot(cend, aes(x=theta, y=score, col=as.factor(type))) +
        ggtitle(title) + #Choose title
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5), legend.position=icclabels,
              panel.grid = panelgrid,
              panel.grid.major.y = gridy,
              panel.grid.major.x = gridx,
              panel.border = element_rect(colour="black", size=0.25, fill=NA),
              panel.background = background) +
        scale_x_continuous(breaks = round(seq(min(cend$theta), max(cend$theta), by = xticks),1),
                           limits=classbreak) +
        scale_y_continuous(breaks = round(seq(min(0), max(cend$score)+0.5, by = yticks),1)) +
        scale_color_manual(values=c("#F0F0F0","#1F78B4","#B2DF8A","#33A02C","#FB9A99",
                                    "#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99",
                                    "#B15928"), name="", labels=c("Expected Item Score", diflabels)) +
        labs(y = yaxistitle, x = xaxistitle) +
        geom_point(na.rm=TRUE) +
        geom_line(na.rm=TRUE) +
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="black") +
        geom_text(data=annotateDIF,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="black", parse=TRUE)
    }

    else if (method!="cut"&dif=="yes") {

      difd <- cbind(data2, difvar)
      difd <- difd[complete.cases(difd),]
      difd[,1:ncol(difd)] <- as.numeric(unlist(difd[,1:ncol(difd)]))
      difs <- min(difd$difvar)
      if (difs>1) {difd$difvar <- ((difd$difvar)-(difs-1))}
      else if (difs==1) {}
      else if (difs<1) {difd$difvar <- ((difd$difvar)+1)}
      else {}

      ldif <- list()
      ttsc <- list()
      obs <- list()
      zero <- list()
      ext <- list()
      tlsc <- list()
      ptmp <- list()

      for (i in 1:length(unique(difd$difvar))) {
        ldif[[i]] <- subset(difd, difd$difvar==i)
        ldif[[i]] <- ldif[[i]][complete.cases(ldif[[i]]),]
        ttsc[[i]] <- rowSums(ldif[[i]][,1:(ncol(ldif[[i]])-1)], na.rm=TRUE)
        ldif[[i]] <- cbind(ldif[[i]][,1:(ncol(ldif[[i]])-1)], ttsc[[i]])
        ldif[[i]] <- ldif[[i]][order(ldif[[i]][,ncol(ldif[[i]])]),]
        obs[[i]] <- aggregate(ldif[[i]][,itmc], by=list(ldif[[i]][,ncol(ldif[[i]])]), FUN=mean)
        zero[[i]] <- which(obs[[i]][,1]==0)
        if (length(zero[[i]])==0) { obs[[i]] <- obs[[i]]}
        else {obs[[i]] <- obs[[i]][-which(obs[[i]][,1]==0),] }
        ext[[i]] <- which(obs[[i]]==mxsc)
        if (length(ext[[i]])==0) { obs[[i]] <- obs[[i]]}
        else {obs[[i]] <- obs[[i]][-which(obs[[i]]==mxsc),] }
        type[i] <- i+1
        names(obs[[i]]) <- c("score", "obsmean")
        tlsc[[i]] <- as.vector(unlist(obs[[i]][1]))
        ptmp[[i]] <- cbind(ppar[1:length(tlsc[[i]])], tlsc[[i]])
        ptmp[[i]] <- as.data.frame(ptmp[[i]])
        names(ptmp[[i]]) <- c("theta", "score")
        obs[[i]] <- merge(obs[[i]], ptmp[[i]], by="score")
        obs[[i]] <- cbind(obs[[i]][,2:3], type[i])
        obs[[i]] <- as.data.frame(obs[[i]])
        names(obs[[i]]) <- c("score", "theta", "type")
      }

      cend <- cend[,1:3]
      cend <- as.data.frame(cend)
      names(cend) <- c("score", "theta", "type")
      big_data = do.call(rbind, obs)
      cend <- rbind(cend, big_data)

      allg <- partgam_DIF(data2, difvar)
      pgmm <- allg[itmc,3]
      pgmm <- format(round(pgmm, digits=2), nsmall=2)

      if(difstats=="yes") {

        annotateDIF <- data.frame(
          xpos = c(-Inf,-Inf,Inf,Inf),
          ypos =  c(-Inf, Inf,-Inf,Inf),
          annotateText = paste("gamma == ", pgmm),
          hjustvar = c(-0.35) ,
          vjustvar = c(4.0)) }

      else {annotateDIF <- data.frame(
        xpos = c(-Inf,-Inf,Inf,Inf),
        ypos =  c(-Inf, Inf,-Inf,Inf),
        annotateText = c(""),
        hjustvar = c(-0.35) ,
        vjustvar = c(4.0)) }

      if(axis.rumm=="yes") {
        classbreak=c(
          max(min(subset(cend$theta, cend$type==1)),
              min((subset(cend$theta, cend$type!=1))-
                    1.0)),
          min(max((subset(cend$theta, cend$type==1))),
              max((subset(cend$theta, cend$type!=1)))+
                1.0))
      }
      else {classbreak=c(NA,NA)}

      myICCplot <- ggplot(cend, aes(x=theta, y=score, col=as.factor(type))) +
        ggtitle(title) + #Choose title
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5), legend.position=icclabels,
              panel.grid = panelgrid,
              panel.grid.major.y = gridy,
              panel.grid.major.x = gridx,
              panel.border = element_rect(colour="black", size=0.25, fill=NA),
              panel.background = background) +
        scale_x_continuous(breaks = round(seq(min(cend$theta), max(cend$theta), by = xticks),1),
                           limits= classbreak) +
        scale_y_continuous(breaks = round(seq(min(0), max(cend$score)+0.5, by = yticks),1)) +
        scale_color_manual(values=c("#F0F0F0","#1F78B4","#B2DF8A","#33A02C","#FB9A99",
                                    "#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99",
                                    "#B15928"), name="", labels=c("Expected Item Score", diflabels)) +
        labs(y = yaxistitle, x = xaxistitle) +
        geom_point(na.rm=TRUE) +
        geom_line(na.rm=TRUE) +
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="black") +
        geom_text(data=annotateDIF,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="black", parse=TRUE)
    }

  }

  if (length(itemnumber)>4)

  { stop(' "The function plots only a maximum of 4 items simultaneously" ')}

  else if (length(itemnumber)==1) {
    itmc=itemnumber
    itmn=itemdescrip[1]
    plot <- pltC()
    plot
  }

  else {
    plst <- list()
    for (i in 1:length(itemnumber)) {
      itmc=itemnumber[i]
      itmn=itemdescrip[i]
      plst[[i]] <- pltC()
    }

    do.call(grid.arrange, args=(c(plst, nrow=2, ncol=2)))
    paste("Please press Zoom on the Plots window to see the plot")
  }
}

