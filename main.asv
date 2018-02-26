% Test the state-space multitaper spectral estimation algorithm on
% SED10.mat (EEG data on a human subject under general anesthesia)
%
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

close all; clear all; clc;

disp('Loading data...');
subjName = 'SED10'; 
subj = [subjName,'.mat']; 
data = load(subj);

data=data.data;
fs = 250; %Sampling rate
channel=1; % Frontal Channel

% Get meaningful data from EEG data
y = data(channel,600*fs+1:1880*fs)';

Nt = length(y);
clear data;

win = 2; % length of window (second)
sf = 1/win; % one step of frequency
nw = win*fs; % the number of elements in a window
N = floor(Nt/nw); % the number of window

% Matrix form of data according to the size of window
yy = reshape(y(1:nw*N),nw,N); 

rtaper = rectwin(nw);
rtaper = rtaper/sqrt(nw);

kyy = zeros(nw,N);
for i = 1 : N
    kyy(:,i) = rtaper.*yy(:,i);
end
% Fourier transfrom of data
frequencyY = fft(kyy,nw,1);

%% EM Algorithm
% First we can limit the frequency range to 0 to 30 Hz and we can adjust
% the max level depending on the EEG data for greater denoising. 
OBSNOISE_CUTOFF = 25*win; % 25 Hz

% Initially we can set the alpha and beta as 1 
alpha = 1;
beta = 1;

% Initial guess for the observation noise
observationNoise = 100;
GUESS_WINDOW_LENGTH = 150; % EM estimatin: 5 min = 300 sec = 300/win = 150

% Estimation of parameters using the EM algorithm for non-tapered data
[sn, on, is, iv, lls] = EM_parameters(alpha, beta, ...
    frequencyY(:,1:GUESS_WINDOW_LENGTH), ...
    observationNoise, 1e-5, OBSNOISE_CUTOFF, 500);
                                  
% Multitapering
TW = 2; % Time-bandwidth production
K = 3; % The Number of tapers
[tapers,concentrations]=dpss(nw,TW,K); % Get the optimal tapers

y_ex = repmat(yy,1,1,K);
y_ex = permute(y_ex,[1 3 2]);
mtY = y_ex;
for i = 1:N
    mtY(:,:,i) = tapers.*y_ex(:,:,i);
end
mtFrequencyY = fft(mtY,nw,1);

% Estimation of parameters using the EM algorithm for multtitapered data
[mtSn, mtOn, mtIs, mtIv, mtLls] = EM_parameters(alpha, beta, ...
    mtFrequencyY(:,:,1:GUESS_WINDOW_LENGTH), ...
    observationNoise, 1e-5, OBSNOISE_CUTOFF, 500);

%% Spectral estimation (periodogram, multitaper, SS-P, SS-MT)

% periodogram
spect1 = periodogram(yy, fs);
% multitaper spectrogram
spect2 = multitaper(yy, fs, TW, K);
% state-space periodogram
spect3 = SS_ST(yy, fs, sn, on, is, iv);
% statate-space multitaper spectrogram
spect4 = SS_MT(yy, fs, TW, K, mtSn, mtOn, mtIs, mtIv);

%% Plot for spectrogram comparisons

fig = figure('color','w','units','normalized','position',[0 0 0.7 0.9]); clf;

fmax = 30;
cmin = -15;
cmax = 10;

colormap jet

ax(1) = subplot(411);
imagesc((1:N)*win/60, (0:fmax*win)*sf, pow2db(spect1(1:fmax*win,:)));
axis xy;
set(gca,'clim',[cmin cmax])
set(ax(1),'xticklabel',[]);
ylabel('Frequency (Hz)');
colorbar
title('Periodogram')
drawnow

ax(2) = subplot(412);
imagesc((1:N)*win/60, (0:fmax*win)*sf, pow2db(spect2(1:fmax*win,:)));
axis xy;
set(gca,'clim',[cmin cmax])
set(ax(2),'xticklabel',[]);
ylabel('Frequency (Hz)');
colorbar
title('Multitaper Spectrogram')
drawnow

ax(3) = subplot(413);
imagesc( (1:N)*win/60, (0:fmax*win)*sf, pow2db(spect3(1:fmax*win,:)));
axis xy;
set(gca,'clim',[cmin cmax])
set(ax(3),'xticklabel',[]);
ylabel('Frequency (Hz)');
colorbar
title('State-Space Periodogram')
drawnow

ax(4) = subplot(414);
imagesc( (1:N)*win/60, (0:fmax*win)*sf, pow2db(spect4(1:fmax*win,:)));
axis xy;
set(gca,'clim',[cmin cmax])
ylabel('Frequency (Hz)');
xlabel('Time (min)');
colorbar
title('State-Space Multitaper Spectrogram')
drawnow


