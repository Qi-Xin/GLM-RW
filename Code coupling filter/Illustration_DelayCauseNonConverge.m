% This code illustrates the "don't converge" problem

clearvars;

len = 5;
Y1 = [1;1;1;0;1];
Y2 = [0;1;1;0;1];
C1 = [1];   % Coupling filter is very simple, just consider one time bin back

max_iter = 1e2;
log_fr_1_rcd = zeros(len,max_iter);   % "rcd" means "record"
log_fr_2_rcd = zeros(len,max_iter);
para_rcd = zeros(4,2*max_iter);
tempy = [0;0];
tempx = [0;0];

% Initialization
log_fr_1 = zeros(len,1);
log_fr_2 = zeros(len,1);

for iter = 1:max_iter
    tempy = glmfit([sameconv(log_fr_1,C1)],Y2,'poisson');
    log_fr_2 = sameconv(log_fr_1,C1)*tempy(2)+tempy(1);
    log_fr_2_rcd(:,iter) = log_fr_2;
    
    para_rcd(:,2*iter-1) = [tempy;tempx];
    
    tempx = glmfit([sameconv(log_fr_2,C1)],Y1,'poisson');
    log_fr_1 = sameconv(log_fr_2,C1)*tempx(2)+tempx(1);
    log_fr_1_rcd(:,iter) = log_fr_1;
    
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