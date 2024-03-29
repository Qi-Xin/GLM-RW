function [I_per] = get_signal(signalType,maxSig,T,jitter)
% 1 for no signal, 2 for square wave, 3 for gamma, 4 for gaussian

I_per = 0*ones(1,T);
if signalType == 2
    I_per(ceil(0.1*T+1):ceil(0.4*T)) = ones(1,length(I_per(ceil(0.1*T+1):ceil(0.4*T))));
    I_per(ceil(0.6*T+1):ceil(0.8*T)) = ones(1,length(I_per(ceil(0.6*T+1):ceil(0.8*T))));
    I_per = I_per/max(I_per)*maxSig;
end
if signalType == 3
    pd = makedist('InverseGaussian','mu',0.4*T *(1+jitter*(rand(1,1)-0.5)),'lambda',1.2*T *(1+jitter*(rand(1,1)-0.5)));
    I_per = pdf(pd,[1:2*T]);
    I_per = I_per/max(I_per)*maxSig *(1+jitter*(rand(1,1)-0.5));
end
if signalType == 4
    pd = makedist('Normal','mu',0.5*T *(1+0*jitter*(rand(1,1)-0.5)),'sigma',0.1*T *(1+0*jitter*(rand(1,1)-0.5)));
    I_per = pdf(pd,[1:T]);
    I_per = I_per/max(I_per)*maxSig*(1+jitter*(rand(1,1)-0.5));
end

end

