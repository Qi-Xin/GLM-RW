function G = sameconv_OnlyLastOne(A, B)
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
    sp = find(As==1);
    for i = 1:length(sp)
        if i == length(sp)
            temp = sp(i)+bm-1;
        else
            temp = min([sp(i+1),sp(i)+bm-1]);
        end
        G(sp(i)+1:temp,j) = B(2:(temp-sp(i)+1));
    end
end

G = G(1:am,:);
G = reshape(G,[],1);