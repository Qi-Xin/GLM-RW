% This code illustrates the "don't converge" problem: forced causality -
% delay is responsible for it. 

clearvars;
addpath(genpath('D:/code'))
X = [0.7;5;0.6];
Y = [2;1;4.5];
C1 = [1];

max_iter = 1e2;
x_rcd = zeros(3,max_iter);
y_rcd = zeros(3,max_iter);
x = X;
y = Y;
for iter = 1:max_iter
    temp = regress(Y, [[1;1;1],sameconv(x,C1)]);
    y = [[1;1;1],sameconv(x,C1)]*temp;
    y_rcd(:,iter) = y;
    temp = regress(X, [[1;1;1],sameconv(y,C1)]);
    x = [[1;1;1],sameconv(y,C1)]*temp;
    x_rcd(:,iter) = x;
end
%% 
figure;
hold on
mycolor = linspace(0,1,7);
for iter = 1:7
    temp = x_rcd(:,iter);
    hx{iter} = plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'*-','linewidth',2,'color',[mycolor(iter) 0 1-mycolor(iter)]);
    mylgd{iter} = ['x',num2str(iter)];
end
iter = iter+1;
temp = X;
hx{iter} = plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'k*-','linewidth',3);
mylgd{iter} = 'X';
legend(mylgd);
view(3)
grid on
%%
figure;
hold on
mycolor = linspace(0,1,7);
for iter = 1:7
    temp = y_rcd(:,iter);
    hy{iter} =plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'*-','linewidth',2,'color',[mycolor(iter) 0 1-mycolor(iter)]);
    mylgd{iter} = ['y',num2str(iter)];
end
iter = iter+1;
temp = Y;
hx{iter} = plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'k*-','linewidth',3);
mylgd{iter} = 'Y';
legend(mylgd);
view(3)
grid on
%%
figure;
subplot(2,1,1)
hold on
plot(x_rcd');
title('x^{(i)}');
legend('x_1^{(i)}','x_2^{(i)}','x_3^{(i)}');
xlabel('iteration (i)');
subplot(2,1,2)
hold on
plot(y_rcd');
title('y^{(i)}');
legend('y_1^{(i)}','y_2^{(i)}','y_3^{(i)}');
xlabel('iteration (i)');
%%

% Show getting y from x and constant
iter = 90;
figure
hold on

% old x
temp = x_rcd(:,iter-1);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'bo-','linewidth',1.5);

% convolution with old x
temp = sameconv(x_rcd(:,iter-1),C1);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'c*-','linewidth',1.5);

% constant
temp = [1;1;1];
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'g*-','linewidth',1.5);

% target Y
temp = Y;
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'kx-','linewidth',1.5);

% new y
temp = y_rcd(:,iter);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'ro-','linewidth',1.5);

v1 = [1;1;1];
v2 = sameconv(x_rcd(:,iter-1),C1);
fill3( [0,v1(1),v1(1)+v2(1),v2(1)] , [0,v1(2),v1(2)+v2(2),v2(2)] ,[0,v1(3),v1(3)+v2(3),v2(3)],'r','FaceAlpha',0.3, ...
    'EdgeAlpha',0);
grid on

legend('x','convolution of x and C','bias','target Y','fitted y','plane of bias and convolution');
title('Getting y of the 11 th iteration')
view(3);
xlim([0 3]);
ylim([0 6]);
zlim([0 5]);
%%

% Show getting x from y and constant
iter = 90;
figure
hold on

% old y
temp = y_rcd(:,iter);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'bo-','linewidth',1.5);

% convolution with old y
temp = sameconv(y_rcd(:,iter),C1);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'c*-','linewidth',1.5);

% constant
temp = [1;1;1];
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'g*-','linewidth',1.5);

% target X
temp = X;
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'kx-','linewidth',1.5);

% new x
temp = x_rcd(:,iter);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'ro-','linewidth',1.5);

v1 = [1;1;1];
v2 = sameconv(y_rcd(:,iter),C1);
fill3( [0,v1(1),v1(1)+v2(1),v2(1)] , [0,v1(2),v1(2)+v2(2),v2(2)] ,[0,v1(3),v1(3)+v2(3),v2(3)],'r','FaceAlpha',0.3, ...
    'EdgeAlpha',0);
grid on

legend('y','convolution of y and C','bias','target X','fitted x','plane of bias and convolution');
title('Getting x of the 11 th iteration')
view(3);
xlim([0 3]);
ylim([0 6]);
zlim([0 5]);
%%

% Show getting y from x and constant
iter = 91;
figure
hold on

% old x
temp = x_rcd(:,iter-1);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'bo-','linewidth',1.5);

% convolution with old x
temp = sameconv(x_rcd(:,iter-1),C1);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'c*-','linewidth',1.5);

% constant
temp = [1;1;1];
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'g*-','linewidth',1.5);

% target Y
temp = Y;
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'kx-','linewidth',1.5);

% new y
temp = y_rcd(:,iter);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'ro-','linewidth',1.5);

v1 = [1;1;1];
v2 = sameconv(x_rcd(:,iter-1),C1);
fill3( [0,v1(1),v1(1)+v2(1),v2(1)] , [0,v1(2),v1(2)+v2(2),v2(2)] ,[0,v1(3),v1(3)+v2(3),v2(3)],'r','FaceAlpha',0.3, ...
    'EdgeAlpha',0);
grid on

legend('x','convolution of x and C','bias','target Y','fitted y','plane of bias and convolution');
title('Getting y of the 12 th iteration')
view(3);
xlim([0 3]);
ylim([0 6]);
zlim([0 5]);
%%

% Show getting y from x and constant
iter = 91;
figure
hold on

% old y
temp = y_rcd(:,iter);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'bo-','linewidth',1.5);

% convolution with old y
temp = sameconv(y_rcd(:,iter),C1);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'c*-','linewidth',1.5);

% constant
temp = [1;1;1];
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'g*-','linewidth',1.5);

% target X
temp = X;
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'kx-','linewidth',1.5);

% new x
temp = x_rcd(:,iter);
plot3([0,temp(1)],[0,temp(2)],[0,temp(3)],'ro-','linewidth',1.5);

v1 = [1;1;1];
v2 = sameconv(y_rcd(:,iter),C1);
fill3( [0,v1(1),v1(1)+v2(1),v2(1)] , [0,v1(2),v1(2)+v2(2),v2(2)] ,[0,v1(3),v1(3)+v2(3),v2(3)],'r','FaceAlpha',0.3, ...
    'EdgeAlpha',0);
grid on

legend('y','convolution of y and C','bias','target X','fitted x','plane of bias and convolution');
title('Getting x of the 12 th iteration')
view(3);
xlim([0 3]);
ylim([0 6]);
zlim([0 5]);

