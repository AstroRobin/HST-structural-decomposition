### Plot the fitting results of multiple exposures in a single figure
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


multiImageMakePlots = function(datalist, imagestack, segim, parm, plottext, magzps, offsets,
                               whichcomponents=list(sersic="all",moffat="all",ferrer="all",pointsource="all"),
                               cmap = rev(colorRampPalette(brewer.pal(9,'RdYlBu'))(100)),
                               errcmap = rev(c("#B00000",colorRampPalette(brewer.pal(9,'RdYlBu'))(100)[2:99],"#0000B0")),
                               errischisq=FALSE, dofs=0){
  
  contlw = 3
  maxsigma = 5
  
  nexps = sum(!grepl('mon.names|parm.names|N|Nim', names(datalist)))
  ndofs = ifelse(missing(dofs), 0, length(dofs)) # degrees of freedom
  if(ndofs>0) stopifnot(length(dofs) <= 2)
  
  parmar = c(1, 1, 0.5, 0.5)
  parmgp = c(3, 0.2, 0)
  par(mar=parmar, oma=c(2,1,1,1))
  
  # Create the plotting layout
  ncols = nexps + 1
  if (ncols > 3){
    nrows = 3
    plotArr = rbind( array(c(seq(1, ncols*2), 0, ncols*2+1), dim=c(2, ncols+1)), c(seq(ncols*2+2, ncols*2 + ncols + 2)) )
  } else if (ncols == 3){
    nrows = 4
    plotArr = rbind( array(c(seq(1, ncols*2), 0, ncols*2+1), dim=c(2, ncols+1)), c(seq(ncols*2+2, ncols*2+ncols), 0, 0), c(seq(ncols*2+ncols+1, ncols*2+ncols+2), 0, 0) )
  } else {
    nrows = 4
    plotArr = rbind( array(c(seq(1, ncols*2), 0, ncols*2+1), dim=c(2, ncols+1)), c(seq(ncols*2+2, ncols*2+ncols+1), 0), c(seq(ncols*2+ncols+2, ncols*2+ncols+3), 0) )
  }
  
  layout(plotArr, widths=c(rep(1/(ncols+0.2), ncols), 0.2/(ncols+0.2)), heights=rep(1/nrows, nrows))
  
  ## Segmentation map plot
  if (missing(segim)){
    magimage(image=imagestack, mgp=parmgp)
  } else {
    profound = profoundProFound(image=imagestack, segim=segim)
    targetIdx = which(profound$segstats$segID == profound$segim[dim(imagestack)[1]%/%2,dim(imagestack)[2]%/%2])
    cutBox = rep(min(profound$segstats$R90[[targetIdx]]*10, dim(imagestack)[1]),2)
    imageCut = magcutout(image=imagestack, loc=c(profound$segstats$xcen[[targetIdx]],profound$segstats$ycen[[targetIdx]]), box=cutBox)$image
    segimCut = magcutout(image=segim, loc=c(profound$segstats$xcen[[targetIdx]],profound$segstats$ycen[[targetIdx]]), box=cutBox)$image
    profoundSegimPlot(image=imageCut, segim=segimCut, mgp=parmgp)
  }
  
  if (missing(offsets)){
    offsets = list()
    for (ii in seq(nexps)){
      offsets[[ii]] = datalist[[ii]]$offset
    }
  }
  
  toterror = c() # the total error accumulated across all pixels in all exposures
  
  for (ii in seq(1, nexps)){ # loop over all exposures
    image = datalist[[ii]]$image
    region = datalist[[ii]]$region
    sigma = datalist[[ii]]$sigma
    magZP = magzps[[ii]]
    
    fitModellist = profitRemakeModellist(parm = parm, Data = datalist[[ii]], offset = offsets[[ii]])$modellist
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
      magimage(modelimage, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap, mgp=parmgp)
      
      targetseg = datalist[[ii]]$segim # get the segmentation map of just the target galaxy
      targetseg[targetseg != targetseg[dim(image)[1]%/%2,dim(image)[2]%/%2]] = 0
      
      segcon = magimage(1-targetseg, add=T, col=NA)
      contour(segcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
      legend('topleft', legend='Model', cex=1.25) # change region to segim
    
      
      ## Bulge-to-total ratio plot
      if (datalist[[ii]]$Nmod >= 2){ # Are there two components?
        if (names(datalist[[ii]]$modellist)[[1]] == "pointsource"){ # psf + sersic
          bulgecomp = list('pointsource'=1)
          diskcomp = list('sersic'=1)
        } else if (names(datalist[[ii]]$modellist)[[1]] == "sersic") { # sersic + sersic
          bulgecomp = list('sersic'=1)
          diskcomp = list('sersic'=2)
        }
        
        bulgemodel = profitMakeModel(modellist = fitModellist, 
                                     magzero = magZP, psf = NULL, dim = dim(image), whichcomponents = bulgecomp)$z
        diskmodel = profitMakeModel(modellist = fitModellist, 
                                    magzero = magZP, psf = NULL, dim = dim(image), whichcomponents = diskcomp)$z
        
        b2timage = bulgemodel/(bulgemodel + diskmodel)
        
      }
    
    }
    
    ## Data image plot
    magimage(data$z, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap, mgp=parmgp)
    
    tempcon=magimage(1-region, add=T, col=NA)
    contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    legend('topleft', legend='Data', cex=1.25)
    
    ## Data - Model image plot
    # magimage(residual, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap)
    # contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    # legend('topleft', legend='Data-Model', cex=1.25)
    
    
    ## Chi-squared plots
    errsign = (-1)^(image > modelimage)
    if(!errischisq){
      error = residual/sigma
    } else {
      error = errsign*sqrt(abs(sigma))
    }
    
    errmap = error
    error = error[region]
    toterror = c(toterror, error)
    
    # modify error map scale and limits for plotting
    maxerr = max(abs(error))
    stretcherr = 1/median(abs(error))
    errmap[!region & (errmap>maxerr)] = maxerr
    minerr = -maxerr
    
    errmap[!region & (errmap<minerr)] = minerr
    errmap[errmap > maxsigma] = maxsigma
    errmap[errmap < -maxsigma] = -maxsigma
    
    magimage(errmap, magmap=FALSE, zlim=c(-maxsigma,maxsigma), col=errcmap, mgp=parmgp)
    contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    legend('topleft',legend=bquote(chi*"=(Data-Model)"/sigma), cex=1.25)

  }
  
  ## Colorbar
  #par(mar = parmar + c(0,max(2,22-2*ncols),0,0), mgp=c(0,0,0))
  breaks = seq(-maxsigma, maxsigma, length.out=length(cmap)+1)
  profitImageScale(zlim=c(-maxsigma,maxsigma), col=errcmap, breaks = breaks, axis.pos=4, axis.padj = -1, mgp=c(0,0,0))

  
  ### Bulge-to-total plot
  if (datalist[[ii]]$Nmod >= 2){ # Are there two components?
    b2tcmap = rev(colorRampPalette(brewer.pal(9,'BrBG'))(100))
    magimage(b2timage, zlim=c(0,1), col=b2tcmap, stretch='lin', lo=0, hi=1, mgp=parmgp)
    legend('topleft',legend=bquote("Bulge-to-Total ratio"), cex=1.25)
  } else {
    plot.new()
  }
  
  ### chi plot
  dx = 0.1
  xlims = c(-4,4)
  x = seq(xlims[1],xlims[2],dx)
  maxerr = max(abs(toterror))
  y = hist(toterror, breaks=c(-(maxerr+dx), x, maxerr+dx), plot=FALSE)$count[2:length(x)]/length(toterror)/dx
  
  ylims = c(min(y[y>0]), 0.5)
  y[y<=0] = ylims[1]/10
  
  vardata = var(error)
  tdof=2*vardata/(vardata-1)
  tdof=LaplacesDemon::interval(tdof,0,Inf)
  
  magplot(x[1:(length(x)-1)]+dx, y, xlim=xlims, ylim=ylims, xlab=expression(chi), ylab="", xaxs="i", type="s", log="y", mgp=parmgp)
  lines(x, dnorm(x), col="blue", xaxs="i")
  lines(x, dt(x,tdof), col="red", xaxs="i")
  
  labs = c(expression(chi),bquote(norm(1)),bquote(Student-T(.(signif(tdof,5)))))
  cols = c("black","blue","red")
  ltys = c(1,1,1)
  legend("bottom",legend=labs,col=cols,lty=ltys, cex=1.25)
  
  
  ### log10 (error^2)
  chisq = sum(toterror^2)/length(toterror)
  
  logerrorsq = log10(toterror^2)
  xr = range(logerrorsq)
  xlims = c(-3, min(2, max(xr)))
  x = seq(xlims[1], xlims[2], dx)
  dxbin = 10^x[2:length(x)]-10^x[1:(length(x)-1)]
  y = hist(logerrorsq, breaks=c(xr[1]-dx,x,xr[2]+dx), plot=FALSE)$count[2:length(x)]
  y = y/sum(y)/dxbin
  ylims = c(min(y[y>0]),10)
  y[y<=0] = ylims[1]/10
  
  magplot(x[1:(length(x)-1)]+dx, y, xlim=xlims, ylim=ylims, xlab=expression(log10(chi^2)), ylab="", xaxs="i", type="s", log="y", mgp=parmgp)
  xp=10^x
  lines(x, dchisq(xp,1), col="blue", xaxs="i")
  labs = c(bquote(chi^2),expression(chi^2*(1)))
  cols = c("black","blue")
  ltys = c(1,1,1)

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
  legend("top",legend=bquote(chi[nu]^2*.(sprintf("=%.2f",chisq))), lty=NULL, cex=1.25)
  
  ### Plot text
  if (!missing(plottext)){
    par(mar = c(0,0,0,0), mgp = c(3,1,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.0, y = 0.9, plottext, adj=c(0,1),
         cex = 1.5, col = "black")
    
  }
  
}



ellipsePlot = function(Data, modellist, model1='sersic', model2='sersic', loc1=1, loc2=2, pixscale=1, FWHM=1, SBlim=26, df=100, raw=FALSE){
  
  if(missing(Data)){stop('Data object of class profit.data must be provided!')}
  if(class(Data)!="profit.data"){stop("Data must be of class profit.data, as output by profitSetupData!")}
  if(missing(modellist)){modellist=Data$modellist}
  
  region = Data$region
  if(is.null(region)) {
    region = numeric(length(Data$image)) + 1
  }
  
  if (!is.null(model2)){
    comp1 = profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(model1 = loc1), dim = Data$imagedim, psf=Data$psf)
    comp2 = profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(model2 = loc2), dim = Data$imagedim, psf=Data$psf)
    total = profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(model1 = 'all'), dim = Data$imagedim, psf=Data$psf)
    
    imageellipse=profitEllipse(Data$image*region, xcen=modellist[[model2]]$xcen[loc2], ycen=modellist[[model2]]$ycen[loc2], ang=modellist[[model2]]$ang[loc2], axrat=modellist[[model2]]$axrat[loc2], box=modellist[[model2]]$box[loc2])
    sigmaellipse=profitEllipse(Data$sigma*region, xcen=modellist[[model2]]$xcen[loc2], ycen=modellist[[model2]]$ycen[loc2], ang=modellist[[model2]]$ang[loc2], axrat=modellist[[model2]]$axrat[loc2], box=modellist[[model2]]$box[loc2])
    comp1ellipse=profitEllipse(comp1$z, xcen=modellist[[model1]]$xcen[loc1], ycen=modellist[[model1]]$ycen[loc1], ang=modellist[[model1]]$ang[loc1], axrat=modellist[[model1]]$axrat[loc1], box=modellist[[model1]]$box[loc1])
    comp2ellipse=profitEllipse(comp2$z, xcen=modellist[[model2]]$xcen[loc2], ycen=modellist[[model2]]$ycen[loc2], ang=modellist[[model2]]$ang[loc2], axrat=modellist[[model2]]$axrat[loc2], box=modellist[[model2]]$box[loc2])
    totalellipse=profitEllipse(total$z, xcen=modellist$sersic$xcen[loc2], ycen=modellist$sersic$ycen[loc2], ang=modellist$sersic$ang[loc2], axrat=modellist$sersic$axrat[loc2], box=modellist$sersic$box[loc2])
    psfellipse=cbind(seq(0,10*FWHM,len=1e3), dnorm(seq(0,10*FWHM,len=1e3),sd=FWHM/(2*sqrt(2*log(2)))))
    
    imageellipse[,1]=imageellipse[,1]*pixscale
    sigmaellipse[,1]=sigmaellipse[,1]*pixscale
    comp1ellipse[,1]=comp1ellipse[,1]*pixscale
    comp2ellipse[,1]=comp2ellipse[,1]*pixscale
    totalellipse[,1]=totalellipse[,1]*pixscale
    
    sigmaellipse=cbind(sigmaellipse,imageellipse[,2]-sigmaellipse[,2])
    sigmaellipse=cbind(sigmaellipse,imageellipse[,2]+sigmaellipse[,2])
    imageellipse[,2]=-2.5*suppressWarnings(log10(imageellipse[,2]))+5*log10(pixscale)+Data$magzero
    sigmaellipse[,2:4]=-2.5*suppressWarnings(log10(sigmaellipse[,2:4]))+5*log10(pixscale)+Data$magzero
    comp1ellipse[,2]=-2.5*suppressWarnings(log10(comp1ellipse[,2]))+5*log10(pixscale)+Data$magzero
    comp2ellipse[,2]=-2.5*suppressWarnings(log10(comp2ellipse[,2]))+5*log10(pixscale)+Data$magzero
    totalellipse[,2]=-2.5*suppressWarnings(log10(totalellipse[,2]))+5*log10(pixscale)+Data$magzero
    psfellipse[,2]=-2.5*suppressWarnings(log10(psfellipse[,2]))+5*log10(pixscale)+Data$magzero
    
    imageellipse[is.na(imageellipse)]=SBlim
    imageellipse[is.infinite(imageellipse)]=SBlim
    sigmaellipse[is.na(sigmaellipse)]=SBlim
    sigmaellipse[is.infinite(sigmaellipse)]=SBlim
    comp1ellipse[is.na(comp1ellipse)]=SBlim+5
    comp1ellipse[is.infinite(comp1ellipse)]=SBlim+5
    comp2ellipse[is.na(comp2ellipse)]=SBlim
    comp2ellipse[is.infinite(comp2ellipse)]=SBlim
    totalellipse[is.na(totalellipse)]=SBlim
    totalellipse[is.infinite(totalellipse)]=SBlim
    psfellipse[is.na(psfellipse)]=SBlim+5
    psfellipse[is.infinite(psfellipse)]=SBlim+5
    
    if(raw){
      predict.image=list(x=imageellipse[,1],y=imageellipse[,2])
      predict.sigma=list(x=sigmaellipse[,1],y=sigmaellipse[,2])
      predict.sigma.lo=list(x=sigmaellipse[,1],y=sigmaellipse[,3])
      predict.sigma.hi=list(x=sigmaellipse[,1],y=sigmaellipse[,4])
      predict.comp1=list(x=comp1ellipse[,1],y=comp1ellipse[,2])
      predict.comp2=list(x=comp2ellipse[,1],y=comp2ellipse[,2])
      predict.total=list(x=totalellipse[,1],y=totalellipse[,2])
      predict.psf=list(x=psfellipse[,1],y=psfellipse[,2])
    }else{
      
      smooth.image=smooth.spline(imageellipse,df=df)
      smooth.sigma.mid=smooth.spline(sigmaellipse[,c(1,2)],df=df)
      smooth.sigma.lo=smooth.spline(sigmaellipse[,c(1,3)],df=df)
      smooth.sigma.hi=smooth.spline(sigmaellipse[,c(1,4)],df=df)
      smooth.comp1=smooth.spline(comp1ellipse,df=df)
      smooth.comp2=smooth.spline(comp2ellipse,df=df)
      smooth.total=smooth.spline(totalellipse,df=df)
      
      predict.image=predict(smooth.image,imageellipse[,1])
      predict.sigma.mid=predict(smooth.sigma.mid,imageellipse[,1])
      predict.sigma.lo=predict(smooth.sigma.lo,imageellipse[,1])
      predict.sigma.hi=predict(smooth.sigma.hi,imageellipse[,1])
      predict.comp1=predict(smooth.comp1,imageellipse[,1])
      predict.comp2=predict(smooth.comp2,imageellipse[,1])
      predict.total=predict(smooth.total,imageellipse[,1])
    }
    
    psfellipse[,2]=psfellipse[,2]+min(predict.comp1$y)-min(psfellipse[1,2])
    
    smooth.total=smooth.spline(totalellipse,df=df)
    refpredict=predict(smooth.total,imageellipse[,1])
    
    xhi=refpredict$x[min(which(refpredict$y>SBlim))]
    yhi=min(predict.image$y,predict.total$y)
    
    sigma.polygon=rbind(cbind(predict.sigma.lo$x,predict.sigma.lo$y),
                        cbind(rev(predict.sigma.hi$x),rev(predict.sigma.hi$y))
    )
    
    layout(rbind(1,2), heights=c(0.7,0.4))
    par(oma=c(3.1,3.6,1.1,1.1))
    par(mar=c(0,0,0,0))
    
    if(pixscale==1){
      xlab='Project Major Axis / pix'
      ylab1=expression(mu*' (mag/pix'*''^2*')')
      ylab2=expression(Delta*mu*' (mag/pix'*''^2*')')
    }else{
      xlab='Project Major Axis / asec'
      ylab1=expression(mu*' (mag/asec'*''^2*')')
      ylab2=expression(Delta*mu*' (mag/asec'*''^2*')')
    }
    
    magplot(0,0,type='n',ylim=c(SBlim+1,yhi),xlim=c(0,xhi),col='black',xlab='', ylab=ylab1, grid=T,labels = c(F,T))
    polygon(sigma.polygon[,1],sigma.polygon[,2],col = hsv(v=0,alpha = 0.1), border=NA)
    lines(predict.image,col='black')
    lines(predict.comp1,col='red')
    lines(predict.comp2,col='blue')
    lines(predict.total,col='darkgreen')
    lines(psfellipse,col='purple')
    abline(h=SBlim,lty=3)
    abline(v=xhi,lty=3)
    abline(v=FWHM/2,lty=3)
    legend('topright',legend=c('Data','ProFit Total','ProFit Bulge','ProFit Disk','PSF'),lty=1,col=c('black','darkgreen','red','blue','purple'),bg='white')
  
  } else {
    total = profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(model1 = loc1), dim = Data$imagedim, psf=Data$psf)
  }
  
  magplot(predict.image$x, predict.image$y-predict.total$y, type='l', xlim=c(0,xhi), ylim=c(-0.5,0.5), col='darkgreen', xlab=xlab, ylab=ylab2, grid=T, labels=c(T,T))
  polygon(sigma.polygon[,1], sigma.polygon[,2]-c(predict.total$y,rev(predict.total$y)), col = hsv(v=0,alpha = 0.1), border=NA)
  abline(h=0)
  abline(v=FWHM/2,lty=3)
}