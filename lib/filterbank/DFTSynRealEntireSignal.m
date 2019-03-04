function x = DFTSynRealEntireSignal(X_in,K,N,p)

Lp = length(p);
x = zeros(1, N*size(X_in,2),'like',X_in);

if all(X_in(:)==0)
    return
end

tdl = zeros(1,length(p),'like',X_in);

% restore conjugate symmetric spectrum and apply iFFT
Y = real(ifft([X_in; conj(X_in(end-1:-1:2,:))],[],1));

tmp = bsxfun(@times,repmat(Y,[ceil(Lp/K) 1]),p.');

Nzeros = zeros(1,N,'like',X_in);

for i = 1:size(X_in,2)
    k = N*(i - 1) + 1;
   
    tdl = [Nzeros tdl(1:Lp-N)] + tmp(:,i).';
    x(k:k+N-1) = (K*N)*tdl(Lp:-1:Lp-N+1);
end
