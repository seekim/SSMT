# SSMT
MATLAB codes for the paper entitled "State-space multitaepr time-frequency analysis"

This repository contains implementations of the algorithms published in PNAS paper entitled "State-space multitaepr time-frequency analysis".

Please cite the following paper if these codes are helpful for your research:

S.-E. Kim, M. Behr, D. Ba, and E. N. Brown, "State-Space Multitaper Time-Frequency Analysis," Proceedings of the National Academy of Sciences, vol. 115, no.1, pp. E5-E14, Jan 2018. (http://www.pnas.org/content/early/2017/12/15/1702877115)

Codes:
  1. main.m: Main code    
  2. EM_parameters.m: Compute noise & state variance using EM algorithm 
  3. periodogram.m: Compute periodogram
  4. multitaper.m: Compute multitaper spectrogram
  5. SS_ST.m: Compute SS periodogram
  6. SS_MT.m: Compute SS mutitaper spectrogram

Instructions: Download all the codes in a directory and run main.m. The code will show four spectrograms for EEG data recorded undere general anesthesia. 
