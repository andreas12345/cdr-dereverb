%ESTIMATE_DELAY_GCC_PHAT Delay estimation using GCC-PHAT.
%
% Estimates the delay between two signals using the GCC-PHAT method.
% Returns the delay
% in samples.
%
% Andreas Schwarz (andreas.schwarz@fau.de), Dec. 2014

function shift = estimate_delay_gcc_phat(x1,x2)
factor =  20; % oversampling (padding) factor to increase time resolution
nfft=1024;
window=hanning(nfft);
n_overlap = nfft/2;
X1 = specgram(x1,nfft,1,window,n_overlap);
X2 = specgram(x2,nfft,1,window,n_overlap);

delta = 10; % regularization constant
norm_CPSD = mean(X1.*conj(X2),2)./mean((abs(X1.*conj(X2))+delta),2);
% padding to increase time resolution of the cross-correlation function
c = fftshift(ifft(norm_CPSD,factor*nfft,'symmetric'));
%plot(c);

[~,shift] = max(c);
shift = (shift-(length(c)/2+1))/factor;
