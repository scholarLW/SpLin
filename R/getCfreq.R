getCfreq = function(x, binGroup, nsg, gmane)
{
  x = as.numeric(x)
  xmean = rep(0, length(nsg))
  for(ikk in 1:length(nsg))
  {
    ix = binGroup$X2 == nsg[ikk]
    nasee = binGroup$X1[ix]
    index = which(gmane %in% nasee)
    if(length(index)>0)
    {
      allcellnum = x[index]
      xmean[ikk] = log2(1+mean(allcellnum) * length(allcellnum>0)/length(allcellnum))
    }
  }
  return(xmean)
}
