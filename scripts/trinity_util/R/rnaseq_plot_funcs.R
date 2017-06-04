plot_MA = function(logCounts, logFoldChange, FDR, logFC_cutoff, FDR_cutoff, xlab="logCounts", ylab="logFC", title="MA plot", pch=20) {

    plot(logCounts, logFoldChange, col=ifelse((FDR<FDR_cutoff & abs(logFoldChange)>=logFC_cutoff), "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);;

}


plot_Volcano = function(logFoldChange, FDR, logFC_cutoff, FDR_cutoff, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20) {

   plot(logFoldChange, -1*log10(FDR), col=ifelse((FDR<FDR_cutoff & abs(logFoldChange)>=logFC_cutoff), "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);

}


plot_MA_and_Volcano = function(logCounts, logFoldChange, FDR, logFC_cutoff, FDR_cutoff, xlab="logCounts", ylab="logFC", title="MA plot") {

    def.par = par(no.readonly = TRUE) # save default, for resetting...

    gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
    layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 

    plot_MA(logCounts, logFoldChange, FDR, logFC_cutoff, FDR_cutoff);
    plot_Volcano(logFoldChange, FDR, logFC_cutoff, FDR_cutoff);

    # draw again, but use a smaller dot for data points
    plot_MA(logCounts, logFoldChange, FDR, logFC_cutoff, FDR_cutoff, pch='.');
    plot_Volcano(logFoldChange, FDR, logFC_cutoff, FDR_cutoff, pch='.');
    

    par(def.par)   
        
    
}
