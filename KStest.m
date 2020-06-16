function [pvalue, rate, h_output, k_output] = KStest(y, h, I, k, bias, PlotFlag)
% Plot KS test

[Ntrials,T] = size(y);

h_length = length(h);
k_length = length(k);

z = [];
rate = zeros(Ntrials,T);
h_output = zeros(Ntrials,T);
k_output = zeros(Ntrials,T);

for i = 1:Ntrials
    rate(i,1)=exp(bias);
    h_output(i,:) = sameconv(y(i,:)',h');
    k_output(i,:) = sameconvSti(I(i,:)',flipud(k'));
    rate(i,:) = exp( h_output(i,:) + k_output(i,:) + bias );
    sp = find(y(i,:) == 1);
    z = [z,sum(rate(i,1:sp(1)))];
    for j = 2:length(sp)
        z = [z,sum(rate( i , (sp(j-1)+1) : sp(j) ))];
    end
	
end

if PlotFlag==1
    figure
    [eCDF, zvals] = ecdf(z);
    mCDF = 1-exp(-zvals); 
    plot(mCDF,eCDF) 
    hold on 
    plot([0 1], [0 1]+1.36/sqrt(length(z)),'k') 
    plot([0 1], [0 1]-1.36/sqrt(length(z)),'k') 
    hold off 
    xlabel('Model CDF','FontSize',13) 
    ylabel('Empirical CDF','FontSize',13)
    title('Goodness-of-fit test','FontSize',16);
    axis([0 1 0 1]);
end

test_cdf = makedist('exp', 'mu', 1);
[h, pvalue] = kstest(z, 'CDF', test_cdf);

figure;
histogram(z,0:0.1:10,'Normalization','pdf');
xlim([0 10]);


end

