library(ProFuse)
library(ProFit)
library(ProFound)
library(LaplacesDemon)

library(magicaxis)
library(RColorBrewer)

library(glue)

profitImageScale = function(z, zlim, col = heat.colors(12),
                             breaks, axis.pos=1, axis.padj=0, add.axis=TRUE, axis.tck = 0, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos, padj=axis.padj, tck=axis.tck)}
}


multiImageMakePlots = function(datalist, imagestack, segim, parm, plottext,
                               whichcomponents=list(sersic="all",moffat="all",ferrer="all",pointsource="all"),
                               cmap = rev(colorRampPalette(brewer.pal(9,'RdYlBu'))(100)),
                               errcmap = rev(c("#B00000",colorRampPalette(brewer.pal(9,'RdYlBu'))(100)[2:99],"#0000B0")),
                               errischisq=FALSE, dofs=0){
  
  magZP = datalist[[1]]$magzero
  contlw = 3
  maxsigma = 5

  ndofs = ifelse(missing(dofs), 0, length(dofs)) # degrees of freedom
  if(ndofs>0) stopifnot(length(dofs) <= 2)
  
  parmar = c(0.5,0.5,0.25,0)
  par(mar=parmar, oma=c(1,1,1,1), mgp=c(0,0,0))
  
  # Create the plotting layout
  nrows = 5
  ncols = datalist$Nim + 1
  plotArr = cbind(c(1,2,3,4+(ncols-1)*nrows,0) ,array(seq(4, 4+(ncols-1)*nrows-1), dim=c(nrows,ncols-1)))
  layout(plotArr, widths=c(0.19,0.19,0.19,0.19,0.19,0.05),heights=rep(1/5,5))
  
  ## Segmentation map plot
  if (missing(segim)){
    magimage(image=imagestack)
  } else {  
    profoundSegimPlot(image=imagestack, segim=segim)
  }
  
  for (ii in seq(1,datalist$Nim)){ # loop over all exposures
    image = datalist[[ii]]$image
    region = datalist[[ii]]$region
    sigma = datalist[[ii]]$sigma
    
    
    fitModellist = profitRemakeModellist(parm = parm, Data = datalist[[ii]], offset = datalist[[ii]]$offset)$modellist
    modelimage = profitMakeModel(modellist = fitModellist,
                               magzero = magZP, psf = datalist[[ii]]$psf, dim = dim(image),
                               psfdim = dim(datalist[[ii]]$psf), whichcomponents = whichcomponents)$z
    
    data = list(x=1:dim(image)[1],y=1:dim(image)[2],z=image)
    residual = image - modelimage
    
    medimg = median(abs(image[region]))/2
    maximg = max(abs(image[region]))
    
    zlims = c(0,1)
    stretch="asinh"
    stretchscale = 1/medimg
    
    if (ii==1){
      ## Model image plot
      magimage(modelimage, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap)
      
      # get the stacked calc region
      imagedim = dim(image)
      seg = datalist[[ii]]$segim==datalist[[ii]]$segim[ceiling(imagedim[1]/2),ceiling(imagedim[2]/2)]
      
      segcon = magimage(1-seg, add=T, col=NA)
      contour(segcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
      legend('topleft', legend='Model', cex=1.25)
      
      ## Colorbar
      par(mar = parmar + c(0,12,0,0), mgp=c(0,0,0))
      breaks = seq(-maxsigma, maxsigma, length.out=length(cmap)+1)
      profitImageScale(zlim=c(-maxsigma,maxsigma), col=errcmap, breaks = breaks, axis.pos=2, axis.padj=1)
      
    }
    
    par(mar=parmar, mgp=c(0,0,0))
    
    ## Data image plot
    magimage(data$z, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap)
    
    tempcon=magimage(1-region, add=T, col=NA)
    contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    legend('topleft', legend='Data', cex=1.25)
    
    ## Data - Model image plot
    magimage(residual, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap)
    contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    legend('topleft', legend='Data-Model', cex=1.25)
    
    
    ## Chi-squared plots
    par(mar=parmar, mgp=c(0,0,0))
    
    errsign = (-1)^(image > modelimage)
    if(!errischisq){
      error = residual/sigma
    } else {
      error = errsign*sqrt(abs(sigma))
    }
    
    errmap = error
    error = error[region]
    maxerr = max(abs(error))
    stretcherr = 1/median(abs(error))
    errmap[!region & (errmap>maxerr)] = maxerr
    minerr = -maxerr
    
    errmap[!region & (errmap<minerr)] = minerr
    errmap[errmap > maxsigma] = maxsigma
    errmap[errmap < -maxsigma] = -maxsigma
    
    magimage(errmap, magmap=FALSE, zlim=c(-maxsigma,maxsigma), col=errcmap)
    contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    legend('topleft',legend=bquote(chi*" = (Data-Model)"/sigma), cex=1.25)
    
    dx = 0.1
    xlims = c(-4,4)
    x = seq(xlims[1],xlims[2],dx)
    y = hist(error,breaks=c(-(maxerr+dx),x,maxerr+dx),plot=FALSE)$count[2:length(x)]/sum(region)/dx
    
    ylims = c(min(y[y>0]), 0.5)
    y[y<=0] = ylims[1]/10
    
    vardata = var(error)
    tdof=2*vardata/(vardata-1)
    tdof=interval(tdof,0,Inf)
    
    magplot(x[1:(length(x)-1)]+dx,y, xlim=xlims, ylim=ylims, xlab="",ylab="", xaxs="i", type="s",log="y")
    lines(x, dnorm(x), col="blue", xaxs="i")
    lines(x, dt(x,tdof), col="red", xaxs="i")
    
    labs = c(expression(chi),bquote(norm(1)),bquote(Student-T(.(signif(tdof,5)))))
    cols = c("black","blue","red")
    ltys = c(1,1,1)
    legend("bottom",legend=labs,col=cols,lty=ltys, cex=1.25)    
    
    error = log10(error^2)
    xr = range(error)
    xlims = c(-3,min(2,max(xr)))
    x = seq(xlims[1],xlims[2],dx)
    dxbin = 10^x[2:length(x)]-10^x[1:(length(x)-1)]
    y = hist(error,breaks=c(xr[1]-dx,x,xr[2]+dx), plot=FALSE)$count[2:length(x)]
    y = y/sum(y)/dxbin
    
    ylims = c(min(y[y>0]),10)
    y[y<=0] = ylims[1]/10
    
    magplot(x[1:(length(x)-1)]+dx, y, xlim=xlims,ylim=ylims, xlab="",ylab="", xaxs="i", type="s",log="y")
    
    lines(x, dchisq(10^x,1), col="blue", xaxs="i")
    labs = c(bquote(chi^2),expression(chi^2*(1)))
    cols = c("black","blue")
    
    if(ndofs > 0){
      dofcols = c("red", "darkgreen")
      for(i in 1:length(dofs)){
        dofstr = sprintf("%.3e",dofs[i])
        lines(x, dchisq(xp,dofs[i]), col=dofcols[i], xaxs="i")
        labs = c(labs,bquote(chi^2 (.(dofstr))))
        ltys = c(ltys,1)
        cols = c(cols, dofcols[i])
      }
    }
    
    abline(v=0,lty=2,col='red')
    legend("bottomleft",legend=labs,col=cols,lty=ltys, cex=1.25)
    
  }
  
  ### Plot text
  if (!missing(plottext)){
    par(mar = c(0,0,0,0), mgp = c(3,1,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.0, y = 0.9, plottext, adj=c(0,1),
         cex = 1.7, col = "black")
  }
  
}


plotDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Plots/Pilot'
resultsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/ProFuse_Outputs/Pilot'
sampleFilename = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Catalogues/Subcats/DEVILS_HST_pilot_sample.csv'

dfSample = read.csv(sampleFilename)
dfSample['UID_fmt'] = format(dfSample['UID'], scientific=FALSE)

unfitCount = 0
for (ii in seq(nrow(dfSample))){
  sourceName = dfSample[['UID']][ii]
  resultFilename = glue('{resultsDir}/{sourceName}/{sourceName}_single_output.rds')
  if (file.exists(resultFilename)){
    if (!file.exists(glue('{plotDir}/{sourceName}/{sourceName}_single_multifit.png'))){
      if (!dir.exists(glue('{plotDir}/{sourceName}'))){
        dir.create(glue('{plotDir}/{sourceName}'))
      }
      result = readRDS(resultFilename)
      
      print(glue('{plotDir}/{sourceName}/{sourceName}_single_multifit.png'))
      png(filename=glue('{plotDir}/{sourceName}/{sourceName}_single_multifit.png'), width=1500, height=1500, pointsize=20)
      multiImageMakePlots(datalist=result$dataList, imagestack=result$imgStack, segim=result$seg$segim, parm=result$highFit$parm,
                          plottext=glue("UID:\n{result$sourceName}\n\n # Exposures: {result$dataList$Nim}\n log(Mstar): {format(round(log10(dfSample[['StellarMass']][ii]), 2), nsmall = 2)}\n z: {format(round(dfSample[['zBest']][ii], 3), nsmall = 2)}"))
      dev.off()
      
      rm(result)
    }
    
  } else {
    unfitCount = unfitCount + 1
    print(glue('{sourceName} .rds file does not exist: {resultFilename}'))
  }
}

print(unfitCount)
