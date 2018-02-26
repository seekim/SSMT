%   This code implements the state-space multitaper spectrogram  
%   described in Kim et al., 2018 PNAS. 
%
%   Usage:
%   main.m: Main code    
%   EM_parameters.m: Compute noise & state variance using EM algorithm 
%   periodogram.m: Compute periodogram
%   multitaper.m: Compute multitaper spectrogram
%   SS_ST.m: Compute SS periodogram
%   SS_MT.m: Compute SS mutitaper spectrogram
%
%   From the paper:
%  "State-space multitpaer time-freqeuncy analysis"
%   Kim, S-E, Behr, MK, Ba, D & Brown, EN
%   PNAS, 2018
%
%   Copyright 2018 The General Hospital Coporation, authored by Seong-Eun Kim, Ph.D.
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%   Last modified 1/11/2018
%
%************************************************************************** 

function [spect, results] = SS_MT(yy, fs, TW, K, stateNoise, observationNoise, ...
                                       initialState, initialVariance)

[nw, N] = size(yy);

[tapers,concentrations]=dpss(nw,TW,K);

y_ex = repmat(yy,1,1,K);
y_ex = permute(y_ex,[1 3 2]);

for i = 1:N
    mtY(:,:,i) = tapers.*y_ex(:,:,i);
end
mtFrequencyY = fft(mtY,nw,1);

mtStateEstimate = zeros(nw,K,N);
mtStatePrediction = zeros(nw,K,N);
mtVarianceEstimate = zeros(nw,K,N);
mtVariancePrediction = zeros(nw,K,N);
mtKalmanGain = zeros(nw,K,N);

mtStateEstimate(:,:,1) = initialState;
mtVarianceEstimate(:,:,1) = initialVariance;
mtSystemNoise = stateNoise;
mtObservationNoise = observationNoise;

for i = 2:N
    mtStatePrediction(:,:,i) = mtStateEstimate(:,:,i-1);
    mtVariancePrediction(:,:,i) = mtVarianceEstimate(:,:,i-1) + mtSystemNoise;
    mtKalmanGain(:,:,i) = mtVariancePrediction(:,:,i) ./ ...
        (mtObservationNoise + mtVariancePrediction(:,:,i));
    mtStateEstimate(:,:,i) = mtStatePrediction(:,:,i) + ...
        mtKalmanGain(:,:,i) .* (mtFrequencyY(:,:,i) - mtStatePrediction(:,:,i));
    mtVarianceEstimate(:,:,i) = (1-mtKalmanGain(:,:,i)) .* mtVariancePrediction(:,:,i);
end

mtSpect = abs(mtStateEstimate).^2;
spect = squeeze(mean(mtSpect,2))/fs;
mtSpect = mtSpect/fs;

results = struct(...
    'spect', spect, ...
    'stateEstimate', mtStateEstimate, ...
    'statePrediction', mtStatePrediction, ...
    'varianceEstimate', mtVarianceEstimate, ...
    'variancePrediction', mtVariancePrediction, ...
    'kalmanGain', mtKalmanGain,...
    'systemNoise', mtSystemNoise, ...
    'observationNoise', mtObservationNoise,...
    'mtSpect',mtSpect);
