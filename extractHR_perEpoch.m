
%% cut up time-frequency result

startTime   = 18; % time (s) when people start interoceptively counting their heart rate
timing      = [25 35 45];
for tIdx = 1:length(timing)
    tStart = find(tData>startTime+sum(timing(1:tIdx-1)),1,'first');
    tEnd = find(tData>startTime+sum(timing(1:tIdx)),1,'first');
    if length(tEnd)==0 % no end time
        tEnd = tData(end);
    end
    disp(['median HR across time for epoch ' num2str(tIdx) '; '  num2str(startTime+sum(timing(1:tIdx-1))) ' - ' num2str(startTime+sum(timing(1:tIdx))) 's; Duration ' num2str(timing(tIdx)) 's: ' num2str(median(60*fvec2(ceil(maxIdx_filt(tStart:tEnd)))))])
end
tStart = find(tData>startTime+sum(timing(1:1-1)),1,'first');
tEnd = find(tData>startTime+sum(timing(1:3)),1,'first');
disp(['median HR across selectively the three epochs; ' num2str(startTime) '-' num2str(startTime+sum(timing)) 's: ' num2str(median(60*fvec2(ceil(maxIdx_filt(tStart:tEnd)))))])
