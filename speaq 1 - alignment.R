library(speaq)

peakList <- detectSpecPeaks(Re(spec),
                            nDivRange = c(128),                
                            scales = seq(1, 16, 2),
                            baselineThresh = 5E6,
                            SNR.Th = -1,
                            verbose=FALSE
);



resFindRef<- findRef(peakList)
refInd <- resFindRef$refInd

maxShift = 50;

Y <- dohCluster(Re(spec),
                peakList = peakList,
                refInd = refInd,
                maxShift  = maxShift,
                acceptLostPeak = TRUE, verbose=T);
