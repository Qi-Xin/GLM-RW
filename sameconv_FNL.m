function G = sameconv_FNL(A, B, n)
%  G = sameconv(A, B);
%   
%  Causally filters A with B, giving a column vector with same height as
%  A.  (B not flipped as in standard convolution).
%
%  Convolution performed efficiently in (zero-padded) Fourier domain.
B = [0;B];

[am, an] = size(A);
[bm, bn] = size(B);
nn = am+bm-1;
G = zeros(nn,1);

sp = find(A==1);
for i = 1:(length(sp)-n)
    temp = min([sp(i+n),sp(i)+bm-1]);
    G(sp(i+n-1)+1:temp) = B(2+(sp(i+n-1)-sp(i)):(temp-sp(i+n-1)+1)+(sp(i+n-1)-sp(i)));
end

G = G(1:am,:);
