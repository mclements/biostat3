## Forked from Matt Cooper's pastebin
## http://pastebin.com/FjkWnCWm
## Described in this article:
## https://mcfromnz.wordpress.com/2011/11/06/kaplan-meier-survival-plot-with-at-risk-table/
## Author: Matt Cooper, 2012-06-21
## Edited: Andreas Karlsson, 2015-02-17
##
## Do not run
## data(colon)
## fit <- survfit(Surv(time,status)~rx, data=colon)
## ggkm(fit, timeby=500, ystratalabs=c("Obs","Lev","Lev+5FU"), table=F, pval=F)
## colon$adhere <- factor(colon$adhere,labels =c("Yes","No"))
## fit <- survfit(Surv(time,status)~rx + adhere, data=colon)
## ggkm(fit, timeby=500, ystratalabs=c("Obs","Lev","Lev+5FU"), subs="No", main="Survival curve for those that don't adhere")
## ggkm(fit, timeby=500, ystratalabs=c("Obs","Lev","Lev+5FU"), subs="Yes", main="Survival curve for those that do adhere")



ggkm <- function(sfit,
		 table = TRUE,
		 returns = FALSE,
		 xlabs = "Time",
		 ylabs = "Survival Probability",
		 xlims = c(0,max(sfit$time)),
		 ylims = c(0,1),
		 ystratalabs = NULL,
		 ystrataname = NULL,
		 timeby = 100,
		 main = "Kaplan-Meier Plot",
		 pval = TRUE,
		 subs = NULL,
		 ...) {

#############
### libraries
#############

    require(ggplot2)
    require(survival)
    require(dplyr)
    require(gridExtra)

#################################
### sorting the use of subsetting
#################################

    times <- seq(0, max(sfit$time), by = timeby)

    if(is.null(subs)){
        subs1 <- 1:length(levels(summary(sfit)$strata))
        subs2 <- 1:length(summary(sfit,censored=T)$strata)
        subs3 <- 1:length(summary(sfit,times = times,extend = TRUE)$strata)
    } else{
        for(i in 1:length(subs)){
            if(i==1){
                ssvar <- paste("(?=.*\\b=",subs[i],sep="")
            }
            if(i==length(subs)){
                ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],"\\b)",sep="")
            }
            if(!i %in% c(1, length(subs))){
                ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],sep="")
            }
            if(i==1 & i==length(subs)){
                ssvar <- paste("(?=.*\\b=",subs[i],"\\b)",sep="")
            }
        }
        subs1 <- which(regexpr(ssvar,levels(summary(sfit)$strata), perl=T)!=-1)
        subs2 <- which(regexpr(ssvar,summary(sfit,censored=T)$strata, perl=T)!=-1)
        subs3 <- which(regexpr(ssvar,summary(sfit,times = times,extend = TRUE)$strata, perl=T)!=-1)
    }

    if( !is.null(subs) ) pval <- FALSE

##################################
### data manipulation pre-plotting
##################################

    if(is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*","",names(sfit$strata))) #[subs1]
    if(is.null(ystrataname)) ystrataname <- "Strata"
    m <- max(nchar(ystratalabs))
    times <- seq(0, max(sfit$time), by = timeby)

    .df <- data.frame(                      # data to be used in the survival plot
                      time = sfit$time[subs2],
                      n.risk = sfit$n.risk[subs2],
                      n.event = sfit$n.event[subs2],
                      surv = sfit$surv[subs2],
                      strata = factor(summary(sfit, censored = T)$strata[subs2]),
                      upper = sfit$upper[subs2],
                      lower = sfit$lower[subs2]
                      )

    levels(.df$strata) <- ystratalabs       # final changes to data for survival plot
    zeros <- data.frame(time = 0, surv = 1,
                        strata = factor(ystratalabs, levels=levels(.df$strata)),
                        upper = 1, lower = 1)
    .df <- rbind_list(zeros, .df)
    d <- length(levels(.df$strata))

###################################
### specifying plot parameteres etc
###################################

    p <- ggplot( .df, aes(time, surv)) +
        geom_step(aes(linetype = strata), size = 0.7) +
            theme_bw() +
            theme(axis.title.x = element_text(vjust = 0.5)) +
                scale_x_continuous(xlabs, breaks = times, limits = xlims) +
                    scale_y_continuous(ylabs, limits = ylims) +
                        theme(panel.grid.minor = element_blank()) +
                            theme(legend.position = c(ifelse(m < 10, .28, .35),ifelse(d < 4, .25, .35))) +    # MOVE LEGEND HERE [first is x dim, second is y dim]
                                theme(legend.key = element_rect(colour = NA)) +
                                    labs(linetype = ystrataname) +
                                        theme(plot.margin = unit(c(0, 1, .5,ifelse(m < 10, 1.5, 2.5)),"lines")) +
                                            ggtitle(main)

    ## Create a blank plot for place-holding
    ## .df <- data.frame()
    blank.pic <- ggplot(.df, aes(time, surv)) +
        geom_blank() + theme_bw() +
            theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
                 axis.title.x = element_blank(),axis.title.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid.major = element_blank(),panel.border = element_blank())

#####################
### p-value placement #
#####################a

    if(pval) {
        sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
        pval <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE)
        pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
        p <- p + annotate("text",x = 0.6 * max(sfit$time),y = 0.1,label = pvaltxt)
    }

###################################################
### Create table graphic to include at-risk numbers
###################################################

    if(table) {
        risk.data <- data.frame(
            strata = factor(summary(sfit,times = times,extend = TRUE)$strata[subs3]),
            time = summary(sfit,times = times,extend = TRUE)$time[subs3],
            n.risk = summary(sfit,times = times,extend = TRUE)$n.risk[subs3]
            )
        risk.data$strata <- factor(risk.data$strata, levels=rev(levels(risk.data$strata)))

        data.table <- ggplot(risk.data,aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
                                        #, color = strata)) +
            geom_text(size = 3.5) + theme_bw() +
                scale_y_discrete(breaks = as.character(levels(risk.data$strata)),
                                 labels = rev(ystratalabs)) +
                                        # scale_y_discrete(#format1ter = abbreviate,
                                        # breaks = 1:3,
                                        # labels = ystratalabs) +
                                     scale_x_continuous("Numbers at risk", limits = xlims) +
                                         theme(axis.title.x = element_text(size = 10, vjust = 1),
                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                              panel.border = element_blank(),axis.text.x = element_blank(),
                                              axis.ticks = element_blank(),axis.text.y = element_text(face = "bold",hjust = 1))

        data.table <- data.table +
            theme(legend.position = "none") + xlab(NULL) + ylab(NULL)

        data.table <- data.table +
            theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 2.5, 3.5) - 0.28 * m), "lines")) # ADJUST POSITION OF TABLE FOR AT RISK

#######################
### Plotting the graphs
#######################

        ## p <- ggplotGrob(p)
        ## p <- addGrob(p, textGrob(x = unit(.8, "npc"), y = unit(.25, "npc"), label = pvaltxt,
        ## gp = gpar(fontsize = 12)))
        grid.arrange(p, blank.pic, data.table, clip = FALSE, nrow = 3,
                     ncol = 1, heights = unit(c(2, .1, .25),c("null", "null", "null")))

        if(returns) {
            a <- arrangeGrob(p, blank.pic, data.table, clip = FALSE, nrow = 3,
                             ncol = 1, heights = unit(c(2, .1, .25), c("null", "null", "null")))
            return(a)
        }
    } else {
        ## p <- ggplotGrob(p)
        ## p <- addGrob(p, textGrob(x = unit(0.5, "npc"), y = unit(0.23, "npc"),
        ## label = pvaltxt, gp = gpar(fontsize = 12)))

        if(returns) return(p)
    }
}
