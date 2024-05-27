library(speaq)

peakList <- detectSpecPeaks(Re(spec[c(1:16),]),
                            nDivRange = c(128),                
                            scales = seq(1, 16, 2),
                            baselineThresh = 5E6,
                            SNR.Th = -1,
                            verbose=T
);



resFindRef<- findRef(peakList)
refInd <- resFindRef$refInd

maxShift = 50;

Y <- dohCluster(Re(spec[c(1:16),]),
                peakList = peakList,
                refInd = refInd,
                maxShift  = maxShift,
                acceptLostPeak = TRUE, verbose=T)
