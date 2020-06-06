function [] = plotraster(spikeMat, tVec, tt)
hold all;
for trialCount = 1:size(spikeMat,1)
    spikePos = tVec(find(spikeMat(trialCount, :)));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount)/1e4 spikePos(spikeCount)/1e4], ...
            [trialCount-0.4-0.5 trialCount+0.4-0.5], 'k');
    end
end
ylim([0 size(spikeMat, 1)+1]);
axis([0 tVec(length(tVec))/1e4 0 size(spikeMat,1)]);
%xlabel('t/s','FontSize',13);
%ylabel('Trial','FontSize',13);
%title(tt,'FontSize',16);