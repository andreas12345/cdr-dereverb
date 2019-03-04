%ESTIMATE_CDR_ROBUST_UNBIASED
% Unbiased estimation of the Coherent-to-Diffuse Ratio (CDR) from the complex
% coherence of a mixed (noisy) signal, using knowledge of both signal and noise
% coherence. This is a variation of estimate_cdr_unbiased which shows better
% performance in practice. Equivalent to CDRprop2 in [1].
%
% CDR = estimate_cdr_nodiffuse(X, N, S)
%       X: complex coherence of mixed (noisy) signal
%       N: coherence of noise component (real-valued)
%       S: coherence of signal component (magnitude one)
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
function CDR = estimate_cdr_robust_unbiased(Cxx,Cnn,Css)
Css = bsxfun(@times, ones(size(Cxx)), Css);
Cnn = bsxfun(@times, ones(size(Cxx)), Cnn);

% limit the magnitude of Cxx to prevent numerical problems
magnitude_threshold = 1-1e-10;
critical = abs(Cxx)>magnitude_threshold;
Cxx(critical) = magnitude_threshold .* Cxx(critical) ./ abs(Cxx(critical));

CDR = 1./(-abs(Cnn-exp(1j*angle(Css)))./(Cnn.*cos(angle(Css))-1)).*abs((exp(-1j*angle(Css)).*Cnn - (exp(-1i*angle(Css)).*Cxx))./(real(exp(-1i*angle(Css)).*Cxx) - 1));

% Ensure we don't get any negative or complex results due to numerical effects
CDR = max(real(CDR),0);
end
