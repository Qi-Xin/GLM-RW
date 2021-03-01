% This code demonstrate normal Gibbs regression between two neurons doesn't
% have the problem. The problem is caused by the forced causality - delay. 
clearvars;

for i = 1:1e5
    X = [2;3;1;0] + random('normal',0,1,4,1);
    Y = [1;3;4;2] + random('normal',0,1,4,1);
    C1 = [1];
    C2 = [0;1];

    max_iter = 1e1;
    x_rcd = zeros(4,max_iter);
    y_rcd = zeros(4,max_iter);
    x = X;
    y = Y;
    for iter = 1:max_iter
        temp = regress(Y, [[1;1;1;1],x]);
        y = [[1;1;1;1],x]*temp;
        y_rcd(:,iter) = y;
        temp = regress(X, [[1;1;1;1],y]);
        x = [[1;1;1;1],y]*temp;
        x_rcd(:,iter) = x;
    end

    if iter>2 & max(x_rcd(:,iter)-x_rcd(:,iter-1))>1e-3
        break;
    end
    
end

%%
figure;
subplot(2,1,1)
hold on
plot(x_rcd');
subplot(2,1,2)
hold on
plot(y_rcd');