function (spectra, driver_peak, mode = "cor", spectrum_resolution = 0.5,ref_spectrum = 1) 
{
  
  driver_peak <- 1.040949
  mode <- "cor"
  ref_spectrum <- 2
  spectra <- spec
  Spectrum <- NAPRMN::pepsMatrixToDF(spectra)
  ppm_range <- as.numeric(colnames(Spectrum[-1]))
  
  driver <- which.min(abs(ppm_range - driver_peak))
  driver_signal <- Spectrum[driver]
  
  cor <- cor(driver_signal[1], Spectrum[-1])
  cov <- cov(driver_signal[1], Spectrum[-1])
  ref <- Spectrum[ref_spectrum, ]
  ref <- tidyr::gather(ref, -Sample, value = "int", key = "ppm")
  cor <- tidyr::gather(as.data.frame(cor), value = "corr", 
    key = "ppm")
  cov <- tidyr::gather(as.data.frame(cov), value = "cov", 
    key = "ppm")
  ref$adjInt <- ref$int * ((cor$corr^2) * (cor$corr/abs(cor$corr)) + 
    0.01)
  ref$corr <- round((cor$corr^2) * (cor$corr/abs(cor$corr)), 
    digits = 2)
  ref$cov <- cov$cov
  plot_cov <- ggplot2::ggplot(ref, ggplot2::aes(as.numeric(ppm), 
    cov, group = Sample)) + ggplot2::geom_path() + ggplot2::geom_point(ggplot2::aes(color = corr), 
    size = 0.05) + ggplot2::scale_color_gradient2(low = "blue", 
    high = "red") + ggplot2::scale_x_reverse()
  if (mode == "cov") {
    Stocsy_plot <- ggplot2::ggplot(ref, ggplot2::aes(as.numeric(ppm), 
      cov, group = Sample, color = corr, xend = dplyr::lead(as.numeric(ppm)), 
      yend = dplyr::lead(cov))) + ggplot2::geom_segment() + 
      ggplot2::scale_color_gradient2(low = "green", mid = "blue", 
        high = "red", midpoint = 0, limits = c(-1, 1)) + 
      ggplot2::scale_x_reverse(breaks = seq(10, 0, -1)) + 
      ggplot2::labs(x = "Chemical Shift")
  }
  else {
    Stocsy_plot <- ggplot2::ggplot(ref, ggplot2::aes(as.numeric(ppm), 
      adjInt, group = Sample, color = corr, xend = dplyr::lead(as.numeric(ppm)), 
      yend = dplyr::lead(adjInt))) + ggplot2::geom_segment() + 
      ggplot2::scale_color_gradient2(low = "green", mid = "blue", 
        high = "red", midpoint = 0, limits = c(-1, 1)) + 
      ggplot2::scale_x_reverse(breaks = seq(10, 0, -1)) + 
      ggplot2::labs(x = "Chemical Shift")
  }
  plotly::ggplotly(Stocsy_plot)
}
