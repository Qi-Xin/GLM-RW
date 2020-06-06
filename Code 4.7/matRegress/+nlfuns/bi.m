function [f,df,ddf] = bi(x)
%  [f,df,ddf] = expfun(x)
%
%  replacement for 'exp' that returns 3 arguments (value, 1st & 2nd deriv)

f = exp(x)./(1+exp(x));
df = exp(x)./(exp(x)+1).^2;
ddf = exp(x).*(exp(x)-1)./(1+exp(x)).^3;
