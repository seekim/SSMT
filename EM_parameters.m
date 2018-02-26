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

function [sn, on, is, iv, logLikelihoods] = ...
        EM_parameters(alpha, beta, frequencyY, observationNoise, ...
                      convergence_criterion, obsnoiseCutoff, max_iteration)
    % we will estimate observation noise variance using frequencies
    % in (0, obsnoiseCutoff] (where obsnoiseCutoff is an index, not
    % Hz). if obsnoiseCutoff is not provided, instead use all
    % frequencies (including 0)

%% Handle single-taper case

% In the multitaper case, our data frequencyY comes with dimensions
% (frequency, taper, time). In the single-taper case, it comes with
% dimensions (frequency, time). Therefore, if we see two
% dimensions, resize it as (frequency, 1, time).

fromSingleTaper = false;
if (ndims(frequencyY) == 2)
    fromSingleTaper = true;
    [nw, N] = size(frequencyY);
    frequencyY = reshape(frequencyY, nw, 1, N);
end

%% Set up variables

assert (ndims(frequencyY) == 3)
[nw, K, N] = size(frequencyY);

logLikelihoods = [];


%% Hyperparameters

sn = ones(nw,K);
on = ones(nw,K)*observationNoise;
is = frequencyY(:,:,1);
iv = is.*conj(is);

if fromSingleTaper
    sn = reshape(sn, nw, 1);
    on = reshape(on, nw, 1);
    is = reshape(is, nw, 1);
    iv = reshape(iv, nw, 1);
end

d = inf;
nIterations = 0;
while mean(d(:)) > convergence_criterion && nIterations < max_iteration
    [nsn, non, nis, niv, ll] = em_step(sn, on, is, iv);
    % don't limit the frequency spectrum here; we still care
    % about high frequencies converging
    for freq_i = 1:nw
        for taper_i = 1:K
            d(freq_i,taper_i) = norm((sn(freq_i, taper_i, :)) - ...
                                     (nsn(freq_i, taper_i, :))) / ...
                norm((sn(freq_i, taper_i, :)) + ...
                     (nsn(freq_i, taper_i, :)));
        end
    end
    sn = nsn;
    on = non;
    is = nis;
    iv = niv;
    nIterations = nIterations + 1;
    logLikelihoods(end+1) = ll;
end

if d <= convergence_criterion
    fprintf('EM converged in %d iterations\n', nIterations);
else
    fprintf('EM terminated after %d iterations\n', nIterations);
end


function [nextSystemNoise, nextObsNoise, ...
          nextInitialState, nextInitialVariance, ll] = ...
        em_step(systemNoise, observationNoise, ...
                initialState, initialVariance)

    %% Expectation step

    % forward pass

    statePrediction = zeros(nw,K,N);
    stateEstimate = zeros(nw,K,N);
    variancePrediction = zeros(nw,K,N);
    varianceEstimate = zeros(nw,K,N);
    kalmanGain = zeros(nw,K,1);

    stateEstimate(:,:,1) = initialState;
    varianceEstimate(:,:,1) = initialVariance;

    for i = 2:N
        statePrediction(:,:,i) = stateEstimate(:,:,i-1);
        variancePrediction(:,:,i) = varianceEstimate(:,:,i-1) + systemNoise;
        kalmanGain = variancePrediction(:,:,i) ./ ...
            (observationNoise + variancePrediction(:,:,i));
        stateEstimate(:,:,i) = statePrediction(:,:,i) + ...
            kalmanGain .* (frequencyY(:,:,i) - statePrediction(:,:,i));
        varianceEstimate(:,:,i) = (1-kalmanGain) .* variancePrediction(:,:,i);
    end

    % backward pass (remember that the state transition matrix is the
    % identity)

    smoothState = zeros(nw,K,N);
    smoothVariance = zeros(nw,K,N);
    % lagVariance(j,k) is sigma^2_{k,k-1|K}(omega_j)
    % (which is the conjugate of sigma^2_{k-1,k|K}(omega_j) )
    lagVariance = zeros(nw,K,N);


    smoothState(:,:,N) = stateEstimate(:,:,N);
    smoothVariance(:,:,N) = varianceEstimate(:,:,N);

    for i = (N-1):-1:1
        smoothGain = ...
            varianceEstimate(:,:,i) ./ variancePrediction(:, :, i+1);
        smoothState(:,:,i) = stateEstimate(:,:,i) + ...
            smoothGain .* (smoothState(:,:,i+1) - statePrediction(:,:,i+1));
        smoothVariance(:,:,i) = varianceEstimate(:,:,i) + ...
            smoothGain.^2 .* (smoothVariance(:,:,i+1) - ...
                           variancePrediction(:,:,i+1));
        % smoothGain is real, so no need to take the conjugate
        lagVariance(:,:,i+1) = smoothGain .* smoothVariance(:,:,i+1);
    end

    % Run the FIS for one more step to get our initial estimates.
    % Under our model, our estimates at time 0 are our predictions at time
    % 1.
    smoothGain = ...
        varianceEstimate(:,:,1) ./ variancePrediction(:, :, 2);
    smoothInitialState = statePrediction(:,:,1) + ...
        smoothGain .* (smoothState(:,:,1) - statePrediction(:,:,1));
    smoothInitialVariance = variancePrediction(:,:,1) + ...
            smoothGain.^2 .* (smoothVariance(:,:,1) - ...
                           variancePrediction(:,:,1));
    lagVariance(:,:,1) = smoothGain .* smoothVariance(:,:,1);

    %% Maximization step

    systemNoiseStateTerm = 2 * (abs(smoothInitialState).^2 ...
                                + sum(abs(smoothState(:,:,1:end-1)).^2, 3)) ...
        + abs(smoothState(:,:,end)).^2;
    systemNoiseVarianceTerm = 2 * (smoothInitialVariance + ...
                                   sum(smoothVariance(:,:,1:end-1), 3)) ...
        + smoothVariance(:,:,end);
    sampleCovariance = smoothInitialState .* conj(smoothState(:,:,1)) + ...
        sum(smoothState(:,:,1:end-1) .* ...
            conj(smoothState(:,:,2:end)), 3);
    lagCovSums = sum(lagVariance,3);

    nextSystemNoise = (systemNoiseStateTerm + systemNoiseVarianceTerm ...
                       - 2*real(sampleCovariance) - 2*real(lagCovSums) + beta) ...
        ./ (alpha + N);

    if exist('obsnoiseCutoff','var')
        whitenoiseIndices = 1:obsnoiseCutoff;
        wnw = length(whitenoiseIndices);
    else
        whitenoiseIndices = 1:nw;
    end

    stateNoiseSum = sum(abs(smoothState(whitenoiseIndices,:,:)).^2, 3);
    stateNoiseVarianceSum = sum(smoothVariance(whitenoiseIndices,:,:), 3);
    obsSum = sum(abs(frequencyY(whitenoiseIndices,:,:)).^2,3);
    obsStateSum = 2*sum(real(conj(frequencyY(whitenoiseIndices,:,:)).*smoothState(whitenoiseIndices,:,:)),3);
    
    nextObsNoise = (beta+sum(obsSum+stateNoiseSum+stateNoiseVarianceSum-obsStateSum))/(N*wnw+alpha-1);
    nextObsNoise = repmat(nextObsNoise,nw,1);
    
    nextInitialState = initialState; %smoothInitialState;
    nextInitialVariance = initialVariance; %smoothVariance(:,:,1);
    
    temp_obs = repmat(observationNoise,1,1,N);
    lls = - (abs(frequencyY - statePrediction).^2) ./ (variancePrediction + temp_obs) ...
          - log(variancePrediction + temp_obs) - log(pi);
    ll = sum(lls(~isnan(lls(:))));

end


end