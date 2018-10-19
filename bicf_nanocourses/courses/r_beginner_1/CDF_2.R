tbl=read.csv(file="AGO-ANKRD52_RNAseq.csv",header=T)
cdf_AGO2_nontargets=ecdf(tbl$AGOKO_notup_ANKRD52KO_logFC)
cdf_AGO2_targets=ecdf(tbl$AGOKO_up_ANKRD52KO_logFC)
pdf("AGO2_targets.pdf",width=6,height=6)
plot(cdf_AGO2_nontargets,
     verticals=TRUE,
     do.points=FALSE,
     lwd=2,
     xlim=c(-2,2),
     ylim=c(0,1),
     main="AGO2 targets vs. non-targets",
     xlab="ANKRD52 KO/WT (log2 FC)",ylab="Cumulative Fraction",
     xaxs="i",yaxs="i",
     col.01line = "WHITE",
     frame.plot=FALSE)
lines(cdf_AGO2_targets,verticals=TRUE,do.points=FALSE,lwd=2,col="red")
ks.test=ks.test(tbl$AGOKO_up_ANKRD52KO_logFC,tbl$AGOKO_notup_ANKRD52KO_logFC,
                alternative="less")
plot.text <- paste("K-S p-value",round(ks.test$p.value,digits=4),sep="=")
text(-1.35,0.8,label=plot.text)
legend("topleft",c("AGO2 targets","AGO2 non-targets"),
       lty=c(1,1),lwd=c(2,2),col=c("red","black"),bty="n")
dev.off()