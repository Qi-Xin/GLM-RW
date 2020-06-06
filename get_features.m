function [peakVal,peakLoc,dipVal,dipLoc] = get_features(f_rcd,fLU_rcd,sti_or_his,tau_E_range,tau_M_range)

[le,lm] = size(f_rcd);
if sti_or_his == 1
    for i = 1:le
        for j = 1:lm
            f_rcd{i,j} = flipud(f_rcd{i,j});
            fLU_rcd{i,j} = flipud(fLU_rcd{i,j});
        end
    end
end

% peak
peakVal = struct;
peakVal.m = NaN*zeros(le,lm);
peakVal.l = NaN*zeros(le,lm);
peakVal.u = NaN*zeros(le,lm);
peakLoc = struct;
peakLoc.m = NaN*zeros(le,lm);
peakLoc.l = NaN*zeros(le,lm);
peakLoc.u = NaN*zeros(le,lm);
for i = 1:le
    for j = 1:lm
        f = f_rcd{i,j};
        fLU = fLU_rcd{i,j};
        fL = fLU(:,1);
        fU = fLU(:,2);
        peakVal.m(i,j) = max(f);
        peakLoc.m(i,j) = find(f == peakVal.m(i,j));
        peakVal.l(i,j) = fL(peakLoc.m(i,j));
        peakVal.u(i,j) = fU(peakLoc.m(i,j));
        [~, peakLoc.l(i,j)] = min(abs(fU(1:peakLoc.m(i,j))-peakVal.m(i,j)));
        [~, peakLoc.u(i,j)] = min(abs(fU(peakLoc.m(i,j):end)-peakVal.m(i,j)));
        peakLoc.u(i,j) = peakLoc.u(i,j) + peakLoc.m(i,j)-1;
    end
end

% dip
dipVal = struct;
dipVal.m = NaN*zeros(le,lm);
dipVal.l = NaN*zeros(le,lm);
dipVal.u = NaN*zeros(le,lm);
dipLoc = struct;
dipLoc.m = NaN*zeros(le,lm);
dipLoc.l = NaN*zeros(le,lm);
dipLoc.u = NaN*zeros(le,lm);
for i = 1:le
    for j = 1:lm
        f = f_rcd{i,j};
        fLU = fLU_rcd{i,j};
        fL = fLU(:,1);
        fU = fLU(:,2);
        if i == 1
            f = f(26:50);
            fL = fL(26:50);
            fU = fU(26:50);
        end
        dipVal.m(i,j) = min(f);
        dipLoc.m(i,j) = find(f == dipVal.m(i,j));
        dipVal.l(i,j) = fL(dipLoc.m(i,j));
        dipVal.u(i,j) = fU(dipLoc.m(i,j));
        [~, dipLoc.l(i,j)] = min(abs(fL(1:dipLoc.m(i,j))-dipVal.m(i,j)));
        [~, dipLoc.u(i,j)] = min(abs(fL(dipLoc.m(i,j):end)-dipVal.m(i,j)));
        dipLoc.u(i,j) = dipLoc.u(i,j) + dipLoc.m(i,j)-1;
        if i == 1
            dipLoc.m(i,j) = dipLoc.m(i,j)+26;
            dipLoc.u(i,j) = dipLoc.u(i,j)+26;
            dipLoc.l(i,j) = dipLoc.l(i,j)+26;
        end
    end
end

figure
subplot(2,2,1);hold on;
target = peakVal;
xx = log(tau_E_range);
for i = 1:lm
    curve = target.m(:,i);
    curvel = target.m(:,i)-target.l(:,i);
    curveu = target.u(:,i)-target.m(:,i);
    errorbar(xx+0.05*(i-2),curve,curvel,curveu);
end
set(gca,'xtick',xx,'xticklabel',tau_E_range);
legend('\tau_M=10','\tau_M=25','\tau_M=50');
title('Peak Value');
xlabel('\tau_{input}');

subplot(2,2,2);hold on;
target = dipVal;
xx = log(tau_E_range);
for i = 1:lm
    curve = target.m(:,i);
    curvel = target.m(:,i)-target.l(:,i);
    curveu = target.u(:,i)-target.m(:,i);
    errorbar(xx+0.05*(i-2),curve,curvel,curveu);
end
set(gca,'xtick',xx,'xticklabel',tau_E_range);
%legend('\tau_M=10','\tau_M=25','\tau_M=50');
title('Dip Value');
xlabel('\tau_{input}');

subplot(2,2,3);hold on;
target = peakLoc;
xx = log(tau_E_range);
for i = 1:lm
    curve = target.m(:,i);
    curvel = target.m(:,i)-target.l(:,i);
    curveu = target.u(:,i)-target.m(:,i);
    errorbar(xx+0.05*(i-2),curve,curvel,curveu);
end
set(gca,'xtick',xx,'xticklabel',tau_E_range);
%legend('\tau_M=10','\tau_M=25','\tau_M=50');
title('Peak Location');
xlabel('\tau_{input}');
ylim([0 length(f)]);

subplot(2,2,4);hold on;
target = dipLoc;
xx = log(tau_E_range);
for i = 1:lm
    curve = target.m(:,i);
    curvel = target.m(:,i)-target.l(:,i);
    curveu = target.u(:,i)-target.m(:,i);
    errorbar(xx+0.05*(i-2),curve,curvel,curveu);
end
set(gca,'xtick',xx,'xticklabel',tau_E_range);
%legend('\tau_M=10','\tau_M=25','\tau_M=50');
title('Dip Location');
xlabel('\tau_{input}');
ylim([0 length(f)]);

end

