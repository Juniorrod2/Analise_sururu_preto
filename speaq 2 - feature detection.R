spectra <- Re(spec)
ppm_data <- as.numeric(colnames(spectra))

peaks <- speaq::getWaveletPeaks(Y.spec=spectra, 
                                     X.ppm=ppm_data, 
                                     baselineThresh = 5E6,
                                     SNR.Th = -1, 
                                     nCPU = 4, 
                                     include_nearbyPeaks = F) 


wine.grouped <- speaq::PeakGrouper(Y.peaks =peaks,  
                                   min.samp.grp = 5, 
                                   grouping.window.width = 200)
ROI.ppm <- 1.330
roiWidth.ppm <- 0.025

speaq::ROIplot(Y.spec = spectra, 
               X.ppm = ppm_data, 
               ungrouped.peaks = peaks,
               grouped.peaks = wine.grouped, 
               ROI.ppm = ROI.ppm,
               roiWidth.ppm = roiWidth.ppm)



filled <- speaq::PeakFilling(Y.grouped = wine.grouped, 
                             Y.spec = spectra,  
                             max.index.shift = 50,
                             nCPU = 4)

features <- BuildFeatureMatrix(filled)

plot_interactive_Spectra(features,limit_n_points = F)
