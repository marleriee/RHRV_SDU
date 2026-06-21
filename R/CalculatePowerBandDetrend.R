############################## CalculatePowerBandDetrend ######################
#' Spectral Power per Band with Configurable Detrending
#' @description
#' A modified version of \code{\link[RHRV]{CalculatePowerBand}} that adds an
#' explicit \emph{detrend} argument controlling what is removed from each
#' window before the FFT is taken.
#' @details
#' \code{CalculatePowerBand} always subtracts only the window mean (DC
#' offset) before applying the Hamming taper; no linear or higher-order
#' trend is removed. Slow drift within a window has nowhere to go except
#' into the lowest frequency bins, which can inflate ULF/VLF estimates.
#' This function adds the option to additionally fit and remove a linear
#' trend from each window before tapering:
#' \itemize{
#'   \item{\code{"none"}: window is used as-is, no preprocessing.}
#'   \item{\code{"mean"}: DC removal only, identical to
#'     \code{CalculatePowerBand}'s existing behaviour.}
#'   \item{\code{"linear"}: mean removal plus a fitted straight-line trend
#'     removed via \code{lm(window ~ t)} residuals.}
#' }
#' Windowing, the Hamming taper, zero-padding, and band-power summation are
#' otherwise identical to \code{CalculatePowerBand}, so the returned
#' \emph{FreqAnalysis} structure is compatible with
#' \code{\link[RHRV]{SplitPowerBandByEpisodes}} and the rest of the
#' frequency-analysis pipeline without further changes.
#'
#' Default band limits are the 1996 Task Force values (Camm et al. 1996)
#' rather than \code{CalculatePowerBand}'s package defaults, since these are
#' typically the convention being tested against. Override with the usual
#' \code{ULFmin}/\code{ULFmax}/etc. arguments if needed.
#' @param HRVData Data structure that stores the beats register and
#' information related to it.
#' @param indexFreqAnalysis An integer referencing the data structure that
#' will contain the frequency analysis.
#' @param size Size of window, in seconds.
#' @param shift Displacement of window for successive calculations, in
#' seconds.
#' @param sizesp Size of window for the FFT, in samples (with zero padding).
#' If \code{NULL}, set to the next power of two of the window size.
#' @param ULFmin,ULFmax Frequency limits (Hz) for the ULF band.
#' @param VLFmin,VLFmax Frequency limits (Hz) for the VLF band.
#' @param LFmin,LFmax Frequency limits (Hz) for the LF band.
#' @param HFmin,HFmax Frequency limits (Hz) for the HF band.
#' @param detrend One of \code{"linear"} (default), \code{"mean"}, or
#' \code{"none"}. See Details.
#' @param taper Logical. If \code{TRUE} (default), apply a Hamming taper as
#' in \code{CalculatePowerBand}. If \code{FALSE}, no taper is applied.
#' @param verbose Deprecated; see \code{\link[RHRV]{SetVerbose}}.
#' @return The \emph{HRVData} structure with the \emph{FreqAnalysis} entry at
#' \code{indexFreqAnalysis} populated with \code{ULF}, \code{VLF}, \code{LF},
#' \code{HF}, \code{HRV}, \code{LFHF}, and \code{Time}, plus the
#' \code{detrend} setting used.
#' @examples
#' \dontrun{
#' data(HRVData)
#' HRVData <- BuildNIHR(HRVData)
#' HRVData <- FilterNIHR(HRVData)
#' HRVData <- InterpolateNIHR(HRVData)
#' HRVData <- CreateFreqAnalysis(HRVData)
#' HRVData <- CalculatePowerBandDetrend(
#'   HRVData, size = 300, shift = 60, detrend = "linear"
#' )
#' }
#' @seealso \code{\link[RHRV]{CalculatePowerBand}},
#' \code{\link[RHRV]{SplitPowerBandByEpisodes}}
#' @export

CalculatePowerBandDetrend <- function(HRVData,
                                       indexFreqAnalysis = length(HRVData$FreqAnalysis),
                                       size, shift, sizesp = NULL,
                                       ULFmin = 0,      ULFmax = 0.0033,
                                       VLFmin = 0.0033, VLFmax = 0.04,
                                       LFmin  = 0.04,   LFmax  = 0.15,
                                       HFmin  = 0.15,   HFmax  = 0.4,
                                       detrend = c("linear", "mean", "none"),
                                       taper = TRUE,
                                       verbose = NULL) {

  detrend <- match.arg(detrend)

  HRVData <- HandleVerboseArgument(HRVData, verbose)
  VerboseMessage(HRVData$Verbose,
                 paste0("Calculating power per band (detrend = ", detrend, ")"))
  CheckAnalysisIndex(indexFreqAnalysis, length(HRVData$FreqAnalysis), "frequency")

  if (max(ULFmin, ULFmax, VLFmin, VLFmax, LFmin, LFmax, HFmin, HFmax) >
      (HRVData$Freq_HR / 2)) {
    stop("Some frequency in the band's limits is bigger than the Nyquist frequency (",
         HRVData$Freq_HR / 2, " Hz).")
  }

  signal    <- 1000.0 / (HRVData$HR / 60.0)   # ms, same as CalculatePowerBand
  signalLen <- length(signal)
  shiftsamples <- shift * HRVData$Freq_HR
  sizesamples  <- floor(size * HRVData$Freq_HR)

  if (is.null(sizesp)) {
    sizesp <- 2^ceiling(log2(sizesamples))
  }
  lenZeroPadding <- if (sizesp <= sizesamples) 0 else (sizesp - sizesamples)

  hamming <- if (taper) {
    0.54 - 0.46 * cos(2 * pi * (0:(sizesamples - 1)) / (sizesamples - 1))
  } else {
    rep(1, sizesamples)
  }
  hammingfactor <- if (taper) 1.586 else 1

  power <- function(spec_arg, freq_arg, fmin, fmax) {
    band <- spec_arg[freq_arg >= fmin & freq_arg < fmax]
    hammingfactor * sum(band) / (2 * length(spec_arg)^2)
  }

  # Number of windows -- identical logic to CalculatePowerBand
  nw <- 1
  begnw <- 1
  repeat {
    begnw <- begnw + shiftsamples
    if ((begnw + sizesamples - 1) >= signalLen) break
    nw <- nw + 1
  }
  VerboseMessage(HRVData$Verbose, paste("Windowing signal...", nw, "windows"))

  freqs <- seq(from = 0, to = HRVData$Freq_HR / 2,
               length.out = (sizesamples + lenZeroPadding) / 2)

  ULF <- VLF <- LF <- HF <- HRV <- numeric(nw)
  t_index <- seq_len(sizesamples)   # reused for the linear fit each window

  for (i in 1:nw) {
    beg <- 1 + (shiftsamples * (i - 1))
    window <- signal[beg:(beg + sizesamples - 1)]
    window[is.na(window)] <- 0

    # ---- the only substantive difference from CalculatePowerBand ----
    window <- switch(detrend,
      "none"   = window,
      "mean"   = window - mean(window),
      "linear" = residuals(lm(window ~ t_index))
    )
    # -------------------------------------------------------------------

    window <- window * hamming
    window <- c(window, rep(0, len = lenZeroPadding))

    spec_tmp <- Mod(fft(window))^2
    spec <- spec_tmp[1:(length(spec_tmp) / 2)]

    HRV[i] <- power(spec, freqs, 0.0, 0.5 * HRVData$Freq_HR)
    ULF[i] <- power(spec, freqs, ULFmin, ULFmax)
    VLF[i] <- power(spec, freqs, VLFmin, VLFmax)
    LF[i]  <- power(spec, freqs, LFmin,  LFmax)
    HF[i]  <- power(spec, freqs, HFmin,  HFmax)
  }

  realShiftTime <- shiftsamples / HRVData$Freq_HR
  realSizeTime  <- sizesamples / HRVData$Freq_HR

  HRVData$FreqAnalysis[[indexFreqAnalysis]]$HRV  <- HRV
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF  <- ULF
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF  <- VLF
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF   <- LF
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF   <- HF
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFHF <- LF / HF
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$Time <- seq(
    realSizeTime / 2, by = realShiftTime, length.out = length(HF)
  )

  HRVData$FreqAnalysis[[indexFreqAnalysis]]$size    <- size
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$shift   <- shift
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$sizesp  <- sizesp
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$detrend <- detrend

  HRVData$FreqAnalysis[[indexFreqAnalysis]]$type   <- "fourier"
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULFmin <- ULFmin
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULFmax <- ULFmax
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLFmin <- VLFmin
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLFmax <- VLFmax
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFmin  <- LFmin
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFmax  <- LFmax
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$HFmin  <- HFmin
  HRVData$FreqAnalysis[[indexFreqAnalysis]]$HFmax  <- HFmax

  VerboseMessage(HRVData$Verbose, "Power per band calculated")

  return(HRVData)
}


############################## compare_detrend_methods ##########################
#' Compare Detrending Methods for Windowed Power-Band Estimates
#' @description
#' Runs \code{\link{CalculatePowerBandDetrend}} on the same \emph{HRVData}
#' object once for each of \code{"none"}, \code{"mean"}, and \code{"linear"}
#' detrending, and returns the mean band power per method in a single table.
#' Intended as a quick diagnostic for whether trend leakage is contributing
#' to inflated ULF/VLF estimates in short windows.
#' @param HRVData Data structure that stores the beats register and
#' information related to it.
#' @param size Size of window, in seconds.
#' @param shift Displacement of window for successive calculations, in
#' seconds.
#' @param ULFmin,ULFmax Frequency limits (Hz) for the ULF band.
#' @param VLFmin,VLFmax Frequency limits (Hz) for the VLF band.
#' @param LFmin,LFmax Frequency limits (Hz) for the LF band.
#' @param HFmin,HFmax Frequency limits (Hz) for the HF band.
#' @return A \code{data.frame} with one row per detrend method, containing
#' \code{detrend}, \code{mean_ULF}, \code{mean_VLF}, \code{mean_LF},
#' \code{mean_HF}, and \code{n_windows}.
#' @examples
#' \dontrun{
#' cmp <- compare_detrend_methods(HRVData, size = 300, shift = 60)
#' print(cmp)
#' }
#' @seealso \code{\link{CalculatePowerBandDetrend}}
#' @export
compare_detrend_methods <- function(HRVData, size, shift,
                                     ULFmin = 0,      ULFmax = 0.0033,
                                     VLFmin = 0.0033, VLFmax = 0.04,
                                     LFmin  = 0.04,   LFmax  = 0.15,
                                     HFmin  = 0.15,   HFmax  = 0.4) {

  methods <- c("none", "mean", "linear")
  out <- data.frame()

  for (m in methods) {
    hd <- CreateFreqAnalysis(HRVData)
    idx <- length(hd$FreqAnalysis)
    hd <- CalculatePowerBandDetrend(
      hd, indexFreqAnalysis = idx, size = size, shift = shift,
      ULFmin = ULFmin, ULFmax = ULFmax,
      VLFmin = VLFmin, VLFmax = VLFmax,
      LFmin  = LFmin,  LFmax  = LFmax,
      HFmin  = HFmin,  HFmax  = HFmax,
      detrend = m
    )

    fa <- hd$FreqAnalysis[[idx]]
    out <- rbind(out, data.frame(
      detrend   = m,
      mean_ULF  = mean(fa$ULF, na.rm = TRUE),
      mean_VLF  = mean(fa$VLF, na.rm = TRUE),
      mean_LF   = mean(fa$LF,  na.rm = TRUE),
      mean_HF   = mean(fa$HF,  na.rm = TRUE),
      n_windows = length(fa$ULF)
    ))
  }

  out
}
