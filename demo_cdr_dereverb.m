%DEMO_CDR_DEREVERB
%
% Demonstration of CDR-based noise and reverberation suppression.
%
% To use this with your own recordings:
% 1. Change wave filename
% 2. Adapt microphone spacing (cfg.d)
% 2. Adapt cfg.TDOA, or use the DOA-independent estimator (estimate_cdr_nodoa)
%
% Reference:
% Andreas Schwarz, Walter Kellermann, "Coherent-to-Diffuse Power Ratio
% Estimation for Dereverberation", IEEE/ACM Trans. on Audio, Speech and
% Lang. Proc., 2015 (under review); preprint available: arXiv:1502.03784
% PDF: http://arxiv.org/pdf/1502.03784
%
% Andreas Schwarz (schwarz@lnt.de)
% Multimedia Communications and Signal Processing
% Friedrich-Alexander-Universitaet Erlangen-Nuernberg (FAU)
% Cauerstr. 7, 91058 Erlangen, Germany

addpath(genpath('lib'));

%% filterbank initialization
cfg.K = 512; % FFT size
cfg.N = 128; % frame shift
cfg.Lp = 1024; % prototype filter length
%p=IterLSDesign(cfg.Lp,cfg.K,cfg.N);
load('lib/filterbank/prototype_K512_N128_Lp1024.mat');

%% algorithm and scenario configuration
cfg.fs = 16000;      % sampling rate [Hz]
cfg.c = 342;         % speed of sound [m/s]
cfg.d_mic = 0.08;   % mic spacing [m]

% all estimators except estimate_cdr_nodoa require the TDOA of the signal; make sure
% to adapt this when loading another wave file
cfg.TDOA = 2.15e-04; % ground truth for wav/roomC-2m-75deg.wav

cfg.nr.lambda = 0.68; % smoothing factor for PSD estimation
cfg.nr.mu = 1.3;     % noise overestimation factor
cfg.nr.floor = 0.1;  % minimum gain
%cfg.nr.alpha = 1; cfg.nr.beta = 1; % power subtraction
cfg.nr.alpha = 2; cfg.nr.beta = 0.5; % magnitude subtraction
%cfg.nr.alpha = 2; cfg.nr.beta = 1; % Wiener filter

%cfg.estimator = @estimate_cdr_unbiased;           % unbiased estimator (CDRprop1)
cfg.estimator = @estimate_cdr_robust_unbiased;    % unbiased, "robust" estimator (CDRprop2)
%cfg.estimator = @estimate_cdr_nodoa;              % DOA-independent estimator (CDRprop3)
%cfg.estimator = @estimate_cdr_nodiffuse;          % noise coherence-independent estimator (CDRprop4; does not work for TDOA -> 0!)

%% preparation
[x,fs_in] = audioread('wav/roomC-2m-75deg.wav');
x = resample(x,cfg.fs,fs_in);

%% Signal processing
% The algorithm itself is real-time capable, i.e., no processing of the entire
% utterance is necessary. Here however, for efficiency of the MATLAB implementation,
% the entire signal is processed at once.

fprintf('Performing signal enhancement... ');tic;

% analysis filterbank
X=DFTAnaRealEntireSignal(x,cfg.K,cfg.N,p);

% estimate PSD and coherence
Pxx = estimate_psd(X,cfg.nr.lambda);
Cxx = estimate_cpsd(X(:,:,1),X(:,:,2),cfg.nr.lambda)./sqrt(Pxx(:,:,1).*Pxx(:,:,2));

frequency = linspace(0,cfg.fs/2,cfg.K/2+1)'; % frequency axis

% define coherence models
Css = exp(1j * 2 * pi * frequency * cfg.TDOA);              % target signal coherence; not required for estimate_cdr_nodoa
Cnn = sinc(2 * frequency * cfg.d_mic/cfg.c); % diffuse noise coherence; not required for estimate_cdr_nodiffuse

% apply CDR estimator (=SNR)
SNR = cfg.estimator(Cxx, Cnn, Css);
SNR = max(real(SNR),0);

weights = spectral_subtraction(SNR,cfg.nr.alpha,cfg.nr.beta,cfg.nr.mu);
weights = max(weights,cfg.nr.floor);
weights = min(weights,1);

% postfilter input is computed from averaged PSDs of both microphones
Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));

% apply postfilter
Processed = weights .* Postfilter_input;

% synthesis filterbank
y = DFTSynRealEntireSignal(Processed,cfg.K,cfg.N,p);
fprintf('done (%.2fs).\n', toc);

%% output
% write output file
audiowrite('wav/out.wav',y,cfg.fs);

%% visualization
figure(1)
subplot(211)
imagesc(10*log10(SNR))
set(gca,'YDir','normal')
caxis([-15 15])
colorbar
title('Estimated CDR (=SNR) [dB]')
xlabel('frame index')
ylabel('subband index')
subplot(212)
imagesc(weights)
set(gca,'YDir','normal')
caxis([0 1])
colorbar
title('Filter gain')
xlabel('frame index')
ylabel('subband index')
