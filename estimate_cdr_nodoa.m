%ESTIMATE_CDR_NODOA 
% Blind (DOA-independent), unbiased estimation of the Coherent-to-Diffuse Ratio (CDR)
% from the complex coherence of a mixed (noisy) signal. Equivalent to CDRprop3 in [1].
%
% CDR = estimate_cdr_nodoa(Cxx, Cnn)
%       Cxx: complex coherence of mixed (noisy) signal
%       Cnn: coherence of noise component (real-valued)
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
function CDR = estimate_cdr_nodoa(Cxx,Cnn,~)
Cnn = bsxfun(@times, ones(size(Cxx)), Cnn); % extend to dimension of Cxx

% limit the magnitude of Cxx to prevent numerical problems
magnitude_threshold = 1-1e-10;
critical = abs(Cxx)>magnitude_threshold;
Cxx(critical) = magnitude_threshold .* Cxx(critical) ./ abs(Cxx(critical));

CDR =  (-(abs(Cxx).^2 + Cnn.^2.*real(Cxx).^2 - Cnn.^2.*abs(Cxx).^2 - 2.*Cnn.*real(Cxx) + Cnn.^2).^(1/2) - abs(Cxx).^2 + Cnn.*real(Cxx))./(abs(Cxx).^2-1);

% Ensure we don't get any negative or complex results due to numerical effects
CDR = max(real(CDR),0);
end
