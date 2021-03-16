function output = sameconv(A, B,T)
%  G = sameconv(A, B);
%   
%  Causally filters A with B, giving a column vector with same height as
%  A.  (B not flipped as in standard convolution).
%
%  Convolution performed efficiently in (zero-padded) Fourier domain.
B = [0;B];
[am, an] = size(A);
[bm, bn] = size(B);
nn = T+bm-1;
output = zeros(am,an);
for i = 1:ceil(am/T)
    AA = A( (i*T-T+1):(i*T) );
    G = ifft(sum(fft(AA,nn).*fft(B,nn),2));
    output( (i*T-T+1):(i*T) ) = G(1:T,:);
end
