library(ggplot2)

spectra <- fullres[-1]
teste <- as.matrix(as.numeric(colnames(spectra)))
spec=as.matrix(spectra)

wine.filled <- speaq::PeakFilling(Y.grouped = wine.grouped, 
                                  Y.spec = as.matrix(spectra),  
                                  max.index.shift = 50,
                                  nCPU = 2) # nCPU set to 1 for the vignette build


ROI.ppm <- 1.330
roiWidth.ppm <- 0.025

speaq::ROIplot(Y.spec = as.matrix(spectra), 
               X.ppm = t, 
               ungrouped.peaks = wine.peaks,
               grouped.peaks = wine.grouped, 
               ROI.ppm = ROI.ppm,
               roiWidth.ppm = roiWidth.ppm)


wine.peaks <- speaq::getWaveletPeaks(Y.spec=as.matrix(spectra), 
                                      X.ppm=teste, 
                                      baselineThresh = 3E6,
                                      SNR.Th = -1, 
                                      nCPU = 2, 
                                 include_nearbyPeaks = TRUE)
