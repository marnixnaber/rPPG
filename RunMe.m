
%% READ ME

% This script performs a remote photoplethysmography (rPPG) analysis on
% videos with algorithms as described in the following publication:
% 
% van der Kooij & Naber (2018). Standardized procedures for the testing and 
% reporting of remote heart rate imaging. Behavior Research Methods
%
% Below you can vary the parameters for the signal processing steps (e.g.
% frequency filtering).
% 
% In the "extractFaceFromVideo.m" file you will find more parameters that
% can be adjusted (e.g., sensitivity to detect faces, number of points to
% track the face, and method to detect skin pixels)


% --------LICENSE & ACKNOWLEDGMENT-----------
% 
% Copyright © 2017 Marnix Naber, The Netherlands
% 
% This program is distributed under the terms of the GNU General Public
% License (see gpl.txt)
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% When you have used this script for scientific purposes, please acknowledge 
% and cite the following reference: 
% 
% van der Kooij & Naber (2018). Standardized procedures for the testing and 
% reporting of remote heart rate imaging. Behavior Research Methods
%
%
%
% -----------------CONTACT-------------------
% 
% For questions about, remarks on, or problems with the code, 
% please contact: marnixnaber@gmail.com
%
% This script was tested in Matlab version 2014b. In case this script
% does not run because of an error reporting a missing function, then
% please check your matlab version and installed toolboxes. To run this
% script succesfully, the image processing toolbox and computer vision
% system toolbox should be installed.

%% set parameters

videoFileName                               = 'rPPG_video.mp4';

signalProcessing = struct();
signalProcessing.samplingRate               = 60;       % [frames per second] sampling rate: temporal resolution of pixel value signal will increase with interpolation to X Hz
signalProcessing.interpMethod               = 'pchip';  % rPPG signal is always interpolated to a frequency of the sampling rate
signalProcessing.FOI                        = [45 165]; % range: [min max] frequency of interest of heart rate in Beats per minute (BPM)
signalProcessing.highPassPixelFilter.active = 1;        % [0 1]; 1 = apply low pass filter to pixel values ... to remove artifacts by movement or illumination
signalProcessing.highPassPixelFilter.params = [6 (signalProcessing.FOI(1)/60)/(signalProcessing.samplingRate/2)];     % [int 0.01-0.10] butterworth parameters --> [6 0.04] is ideal for frame rate of 30

signalProcessing.lowPassPSDFilter.active    = 1;        % [0 1]; 1 = apply low pass filter to power density spectrum to remove spurious peaks due to noise
signalProcessing.lowPassPSDFilter.params    = [8 0.2];  % [int 0.01-0.10] butterworth parameters [Xth_order cutoff_freq]
signalProcessing.lowPassHRTimeFilter.active = 1;        % [0 1]; 1 = apply low pass filter to heart rate over time to remove spurious changes in HR due to noise
signalProcessing.lowPassHRTimeFilter.params = [8 0.03]; % [int 0.01-0.10] butterworth parameters [Xth_order cutoff_freq]
            
[signalProcessing.highPassPixelFilter.NthOrder,signalProcessing.highPassPixelFilter.cutOffFreq] = butter(signalProcessing.highPassPixelFilter.params(1),signalProcessing.highPassPixelFilter.params(2));
[signalProcessing.lowPassPSDFilter.NthOrder,signalProcessing.lowPassPSDFilter.cutOffFreq]       = butter(signalProcessing.lowPassPSDFilter.params(1),signalProcessing.lowPassPSDFilter.params(2));
[signalProcessing.lowPassHRTimeFilter.NthOrder,signalProcessing.lowPassHRTimeFilter.cutOffFreq] = butter(signalProcessing.lowPassHRTimeFilter.params(1),signalProcessing.lowPassHRTimeFilter.params(2));

signalProcessing.ica.nComps                  = 3;        %
signalProcessing.ica.nIte                    = 2000;     %
signalProcessing.ica.stab                    = 'on';     %
signalProcessing.ica.verbose                 = 'off';    %

%%

[pixelValPerFrame,faceMap,vidInfo,faceDetection,faceTracking] = extractFaceFromVideo(videoFileName,'all', 1);


%% cutoff beginning and end if missing values, and fill missing values, and resample 

startIdx        = find(isfinite(pixelValPerFrame(:,1)),1,'first');
endIdx          = find(isfinite(pixelValPerFrame(:,1)),1,'last');

xdata           = vidInfo.tStamp(startIdx:endIdx)-vidInfo.tStamp(startIdx);
resampledXdata  = linspace(xdata(1),xdata(end),ceil(signalProcessing.samplingRate*(xdata(end)-xdata(1))));
resampledYdata = NaN(length(resampledXdata),3);

for c = 1:3    
    ydata = pixelValPerFrame(startIdx:endIdx,c);
    vect = isfinite(ydata);                
    resampledYdata(:,c) = interp1(xdata(vect),ydata(vect),resampledXdata,signalProcessing.interpMethod);
end

%% Filter out low frequency changes

pixelVal_filt = resampledYdata;
for c = 1:3
    if signalProcessing.highPassPixelFilter.active
        pixFilter = filtfilt(signalProcessing.highPassPixelFilter.NthOrder,signalProcessing.highPassPixelFilter.cutOffFreq,resampledYdata(:,c));
        pixelVal_filt(:,c) = resampledYdata(:,c)-pixFilter;
    else
        pixelVal_filt(:,c) = resampledYdata(:,c);
    end
end

%% ICA on different color channels

finalSignal = struct();
finalSignal.comp                = fastica(pixelVal_filt','numOfIC',signalProcessing.ica.nComps,'maxNumIterations',signalProcessing.ica.nIte,'stabilization',signalProcessing.ica.stab,'verbose',signalProcessing.ica.verbose);

%% Fourier transform

signalProcessing.fft.L                  = signalProcessing.samplingRate*length(resampledXdata);
signalProcessing.fft.NFFT               = 2^nextpow2(signalProcessing.fft.L); % Next power of 2 from length of y
signalProcessing.fft.freq               = signalProcessing.samplingRate/signalProcessing.fft.NFFT*(0:signalProcessing.fft.NFFT-1);
signalProcessing.fft.freqInterestRange  = signalProcessing.FOI/60;
signalProcessing.fft.fRange2            = find(signalProcessing.fft.freq>signalProcessing.fft.freqInterestRange(1) & signalProcessing.fft.freq<signalProcessing.fft.freqInterestRange(2));
signalProcessing.fft.HRRange            = 60*signalProcessing.fft.freq(signalProcessing.fft.fRange2);

finalSignal.powerVal     = [];
for co = 1:signalProcessing.ica.nComps

    Y                 = fft(finalSignal.comp(co,:),signalProcessing.fft.NFFT); % calculate frequency spectrum
    finalSignal.powerVal(:,co)    = Y.*conj(Y)/signalProcessing.fft.NFFT;

end

%% Low-pass filter fourier spectrum

finalSignal.coherenceVal = [];
for co = 1:signalProcessing.ica.nComps
    if signalProcessing.lowPassPSDFilter.active
        finalSignal.powerVal(:,co) = filtfilt(signalProcessing.lowPassPSDFilter.NthOrder,signalProcessing.lowPassPSDFilter.cutOffFreq,finalSignal.powerVal(:,co));
    end

    finalSignal.coherenceVal(:,co) = finalSignal.powerVal(:,co)./sqrt(sum(finalSignal.powerVal(signalProcessing.fft.fRange2,co).^2));
%     finalSignal.coherenceVal(:,co) = finalSignal.powerVal(:,co)./sqrt(sum(finalSignal.powerVal(:,co).^2));

end

%% Extract HR

finalSignal.maxPower        = [];
finalSignal.maxCoherence    = [];
tempMaxPowerIdx             = [];
tempMaxCoherenceIdx         = [];
for co = 1:signalProcessing.ica.nComps
    [finalSignal.maxPower(co),tempMaxPowerIdx(co)] = max(finalSignal.powerVal(signalProcessing.fft.fRange2,co));
    [finalSignal.maxCoherence(co),tempMaxCoherenceIdx(co)] = max(finalSignal.coherenceVal(signalProcessing.fft.fRange2,co));
end
[~,tempMaxPowerCompIdx] = max(finalSignal.maxPower);
[~,tempMaxCoherenceCompIdx] = max(finalSignal.maxCoherence);

finalSignal.HRfreqIdx_powerBased    = tempMaxPowerIdx(tempMaxPowerCompIdx);
finalSignal.HR_powerBased           = signalProcessing.fft.freq(signalProcessing.fft.fRange2(tempMaxPowerIdx(tempMaxPowerCompIdx)))*60;
finalSignal.HR_powerBased_PerComp   = signalProcessing.fft.freq(signalProcessing.fft.fRange2(tempMaxPowerIdx))*60;
finalSignal.HRBestCompIdx_powerBased    = tempMaxPowerCompIdx;

finalSignal.HRfreqIdx_coherenceBased    = tempMaxCoherenceIdx(tempMaxCoherenceCompIdx);
finalSignal.HR_coherenceBased           = signalProcessing.fft.freq(signalProcessing.fft.fRange2(tempMaxCoherenceIdx(tempMaxCoherenceCompIdx)))*60;
finalSignal.HR_PerComp_coherenceBased   = signalProcessing.fft.freq(signalProcessing.fft.fRange2(tempMaxCoherenceIdx))*60;
finalSignal.HRBestCompIdx_coherenceBased    = tempMaxCoherenceCompIdx;

%% plot frequency spectrum

figure();
subplot(2,2,1);
plot(signalProcessing.fft.HRRange,finalSignal.powerVal(signalProcessing.fft.fRange2,:));
hold on
line([repmat(finalSignal.HR_powerBased,1,2)],[0 max(finalSignal.powerVal(signalProcessing.fft.fRange2,finalSignal.HRBestCompIdx_powerBased))],'Color','k','LineStyle',':','LineWidth',2)
xlabel('HR')
ylabel('Power')

subplot(2,2,2)
plot(signalProcessing.fft.HRRange,finalSignal.coherenceVal(signalProcessing.fft.fRange2,:));
hold on
line([repmat(finalSignal.HR_coherenceBased,1,2)],[0 max(finalSignal.coherenceVal(signalProcessing.fft.fRange2,finalSignal.HRBestCompIdx_coherenceBased))],'Color','k','LineStyle',':','LineWidth',2)
xlabel('HR')
ylabel('Coherence')
legend('Component 1','Component 2','Component 3')

subplot(2,2,4);
plot(resampledXdata,finalSignal.comp(finalSignal.HRBestCompIdx_coherenceBased,:),'k')
title('Best component - coherence based')
ylabel('Pixel val')
xlabel('Time [s]')




%% Time frequency analysis

computerCoherence = 1; % 0 = raw data, 1 = coherence transformed

winSize     = 10;   % X seconds
tempRes     = 240;  % temporal resolution, number of time points - THIS NEEDS TO BE CHANGED TO NUMBER OF POINTS PER SECOND
freqRes     = 120;  % frequency resolution, number of frequencies

winSize = winSize*ceil(signalProcessing.samplingRate);
tBins   = floor(linspace(1,length(finalSignal.comp)-winSize,tempRes));
fvec    = linspace(signalProcessing.FOI(1)/60,signalProcessing.FOI(2)/60,freqRes);

timeFreqData = NaN(length(fvec),length(tBins));
tData = NaN(1,length(tBins));
countT = 0;
for tt = tBins
    countT = countT+1;
    t = tt:tt+winSize;
    [timeFreqData(:,countT),fvec2] = plomb(finalSignal.comp(finalSignal.HRBestCompIdx_powerBased,t),resampledXdata(t),fvec);
%     [timeFreqData(:,countT),fvec2] = plomb(finalSignal.comp(finalSignal.HRBestCompIdx_coherenceBased,t),resampledXdata(t),fvec);
    tData(countT) = mean(resampledXdata(t));
    
    if computerCoherence
        timeFreqData(:,countT) = timeFreqData(:,countT)./sqrt(sum(timeFreqData(:,countT)).^2);
    end
    
end

% signalProcessing.lowPassHRTimeFilter.params = [8 0.03];  % [int 0.01-0.10] butterworth parameters [Xth_order cutoff_freq]
% [signalProcessing.lowPassHRTimeFilter.NthOrder,signalProcessing.lowPassHRTimeFilter.cutOffFreq] = butter(signalProcessing.lowPassHRTimeFilter.params(1),signalProcessing.lowPassHRTimeFilter.params(2));
        
figure();
imagesc(timeFreqData)
hold on
[~,maxIdx] = max(timeFreqData);
line(1:length(maxIdx),maxIdx,'Color','g','LineWidth',3,'LineStyle',':')

if signalProcessing.lowPassHRTimeFilter.active
    maxIdx_filt = filtfilt(signalProcessing.lowPassHRTimeFilter.NthOrder,signalProcessing.lowPassHRTimeFilter.cutOffFreq,maxIdx);
end
line(1:length(maxIdx_filt),maxIdx_filt,'Color','b','LineWidth',3)

set(gca,'xtick',round(linspace(1,length(tBins),5)))
set(gca,'xticklabel',round(tData(round(linspace(1,length(tBins),5)))))

set(gca,'ytick',round(linspace(1,length(fvec),10)))
set(gca,'yticklabel',round(60*fvec(round(linspace(1,length(fvec),10)))))

xlabel('Time (s)')
ylabel('Heart rate (bpm)')
colormap('hot')

