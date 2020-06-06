function G = sameconv(A, B)
%  G = sameconv(A, B);
%   
%  Causally filters A with B, giving a column vector with same height as
%  A.  (B not flipped as in standard convolution).
%
%  Convolution performed efficiently in (zero-padded) Fourier domain.


[am, an] = size(A);
[bm, bn] = size(B);
nn = am+bm-1;
G = zeros(nn,an);
for j = 1:an
    As = A(:,j);
    G(:,j) = ifft(sum(fft(As,nn).*fft(B,nn),2));
end

G = G(1:am,:);
G = reshape(G,[],1);