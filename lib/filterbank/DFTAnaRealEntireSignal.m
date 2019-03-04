function X_out=DFTAnaRealEntireSignal(x_in,K,N,p)
Lp = length(p);
Lx = size(x_in,1);
n_ch = size(x_in, 2);
n_blocks = ceil(Lx/N);

X_out = zeros([K/2+1, n_blocks, n_ch],'like',x_in);
for ch_ix=1:n_ch
    x_tmp = x_in(:,ch_ix);
    x_tmp = x_tmp(:).';
    x_buffer = buffer([zeros(1,N-1) x_tmp],Lp,Lp-N);
    x_buffer = x_buffer(Lp:-1:1,1:n_blocks);
    x_tmp = [];
    U = reshape(bsxfun(@times, x_buffer, p.'),K,ceil(Lp/K),n_blocks);
    x_buffer = [];
    V = squeeze(fft(sum(U,2),[],1));
    U = [];
    X_out(:,:,ch_ix) = V(1:K/2+1,:);
end
