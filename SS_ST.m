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

function [spect, results] = SS_ST(yy, fs, systemNoise, observationNoise, ...
                                     initialState, initialVariance)

[nw, N] = size(yy);

rtaper = rectwin(nw)/sqrt(nw);

kyy = zeros(nw,N);
for i = 1 : N
    kyy(:,i) = rtaper.*yy(:,i);
end
frequencyY = fft(kyy,nw,1);

% Initialization
statePrediction = zeros(nw,N); 
stateEstimate = zeros(nw,N);
variancePrediction = zeros(nw,N); 
varianceEstimate = zeros(nw,N);
kalmanGain = zeros(nw,N);

stateEstimate(:,1) = initialState;
varianceEstimate(:,1) = initialVariance;

% Kalman filtering
for i = 2:N 
    statePrediction(:,i) = stateEstimate(:,i-1);
    variancePrediction(:,i) = varianceEstimate(:,i-1) + systemNoise;
    kalmanGain(:,i) = variancePrediction(:,i) ./ ...
        (observationNoise + variancePrediction(:,i));
    stateEstimate(:,i) = statePrediction(:,i) + ...
        kalmanGain(:,i) .* (frequencyY(:,i) - statePrediction(:,i));
    varianceEstimate(:,i) = (1-kalmanGain(:,i)) .* variancePrediction(:,i);
end

spect = abs(stateEstimate).^2/fs;

results = struct(...
    'spect', spect, ...
    'stateEstimate', stateEstimate, ...
    'statePrediction', statePrediction, ...
    'varianceEstimate', varianceEstimate, ...
    'variancePrediction', variancePrediction, ...
    'kalmanGain',kalmanGain,...
    'systemNoise', systemNoise, ...
    'observationNoise', observationNoise);
