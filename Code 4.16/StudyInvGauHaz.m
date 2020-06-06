clearvars

ddt = 1e-2;
max1 = 0.5;
pd = makedist('InverseGaussian','mu',3,'lambda',0.1);
pdfIG = pdf(pd,0:ddt:max1);
cdfIG = cdf(pd,0:ddt:max1);
svvIG = 1 - cdfIG;
hazard = pdfIG./svvIG;
figure
hold on
plot(0:ddt:max1, hazard);
plot(0:ddt:max1, pdfIG);
plot(0:ddt:max1, svvIG);
legend('Hazard','pdf','svv');