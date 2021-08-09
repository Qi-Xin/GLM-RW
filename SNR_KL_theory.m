clearvars;

t_tot = 1e3;
maxSig = 7e-2;
minSig = 2e-2;

sum_range = 10;
EPEfr_v = zeros(1,25);
EPEconst_v = zeros(1,25);
EPEfr_kl = zeros(1,25);
EPEconst_kl = zeros(1,25);

j = 0;
width_range = 10:10:250;
for width = width_range
    j = j + 1;
    pd = makedist('Normal','mu',t_tot/2,'sigma',width);
    normalcurve = pdf(pd,[1:1e3]);
    fr = minSig*ones(1,t_tot) + (maxSig-minSig)*normalcurve/max(normalcurve);
    meanfr = mean(fr)*ones(1,t_tot);

    for i = 0:sum_range
        EPEfr_v(j) = EPEfr_v(j) + sum( fr.^i.*exp(-fr)./factorial(i) .* (i-fr).^2 );
        EPEconst_v(j) = EPEconst_v(j) + sum( fr.^i.*exp(-fr)./factorial(i) .* (i-meanfr).^2 );
        EPEfr_kl(j) = EPEfr_kl(j) + sum( fr.^i.*exp(-fr)./factorial(i) .* log(fr.^i.*exp(-fr)./factorial(i)) );
        EPEconst_kl(j) = EPEconst_kl(j) + sum( fr.^i.*exp(-fr)./factorial(i) .* log(meanfr.^i.*exp(-meanfr)./factorial(i)) );
    end
    
end

SNR_v = (EPEconst_v - EPEfr_v)./EPEfr_v;
SNR_kl = (EPEconst_kl - EPEfr_kl)./EPEfr_kl;

%%
figure
yyaxis left
plot(width_range,SNR_kl);

yyaxis right
plot(width_range,SNR_v);

legend('KL based SNR','Variance based SNR');
xlabel('Bump width (ms)');
%%
figure
plot(1:1e3,fr)
