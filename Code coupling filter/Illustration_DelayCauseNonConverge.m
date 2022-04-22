% This code illustrates the "don't converge" problem

clearvars;

len = 5;
X = [1;1;1;0;1];
Y = [0;1;1;0;1];
C1 = [1];

max_iter = 1e2;
x_rcd = zeros(len,max_iter);
y_rcd = zeros(len,max_iter);
x = X;
y = Y;
para_rcd = zeros(4,2*max_iter);
tempy = [0;0];
tempx = [0;0];
x = zeros(len,1);
y = zeros(len,1);

for iter = 1:max_iter
    tempy = glmfit([sameconv(x,C1)],Y,'poisson');
    y = sameconv(x,C1)*tempy(2)+tempy(1);
    y_rcd(:,iter) = y;
    
    para_rcd(:,2*iter-1) = [tempy;tempx];
    
    tempx = glmfit([sameconv(y,C1)],X,'poisson');
    x = sameconv(y,C1)*tempx(2)+tempx(1);
    x_rcd(:,iter) = x;
    
    para_rcd(:,2*iter) = [tempy;tempx];

end

figure;
hold on
plot((1:(2*max_iter))/2, para_rcd', '-');
title('Parameters');
legend('\beta_2','k_{12}','\beta_1','k_{21}');
xlabel('iteration (i)');


function G = sameconv(A, B)
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

    G = ifft(sum(fft(A,nn).*fft(B,nn),2));
    G = G(1:am,:);
end