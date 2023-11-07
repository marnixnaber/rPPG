
%% READ ME

% This script performs a remote photoplethysmography (rPPG) analysis on
% videos with algorithms as described in the following publication:
% 
% van der Kooij, K.M., Naber, M. An open-source remote heart rate imaging
% method with practical apparatus and algorithms.
% Behav Res 51, 2106–2119 (2019).
% https://doi.org/10.3758/s13428-019-01256-8

% Below you can vary the parameters for the signal processing steps (e.g.
% frequency filtering).
% 
% In the "extractFaceFromVideo.m" file you will find more parameters that
% can be adjusted (e.g., sensitivity to detect faces, number of points to
% track the face, and method to detect skin pixels)
%
% IMPORTANT: to get an accurate measurement you need to ensure that the 
% following conditions are met:
% - Enough luminance (face a window with daylight when recording a video)
% - Position as close to the camera as possible (<40cm from the camera)
% - High camera frame rate (>20 frames per second). A high resolution is less but also important
% - No or high quality video encoding (lossless, no compression)

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
% and cite the reference above.
%
%
% -----------------CONTACT-------------------
% 
% For questions about, remarks on, or problems with the code, 
% please contact: marnixnaber@gmail.com
%
% This script was tested in Matlab version 2014b and 2019b. In case this script
% does not run because of an error reporting a missing function, then
% please check your matlab version and installed toolboxes. To run this
% script succesfully, the image processing toolbox and computer vision
% system toolbox should be installed. This script further uses a fast
% independent component analysis designed by 
% Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen 

%% set parameters

videoFileName                               = 'rPPG_video.mp4';

signalProcessing = struct();
signalProcessing.HrDetectionMethod          = 'fastica';

% 'fastica' (see Van der Kooij & Naber, 2019); 
% 'ica' (see Poh, M. Z., McDuff, D. J., & Picard, R. W., 2010); 
% 'pos' (see Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G., 2017; state-of-the-art rppg)


signalProcessing.samplingRate               = 60;       % [frames per second] sampling rate: temporal resolution of pixel value signal will increase with interpolation to X Hz
signalProcessing.interpMethod               = 'pchip';  % rPPG signal is always interpolated to a frequency of the sampling rate
signalProcessing.FOI                        = [45 165]; % range: [min max] frequency of interest of heart rate in Beats per minute (BPM)

signalProcessing.highPassPixelFilter.active = 1;        % [0 1]; 1 = apply low pass filter to pixel values ... to remove artifacts by movement or illumination
signalProcessing.highPassPixelFilter.params = [6 (signalProcessing.FOI(1)/60)/(signalProcessing.samplingRate/2)];     % [int 0.01-0.10] butterworth parameters --> [6 0.04] is ideal for frame rate of 30

signalProcessing.lowPassPSDFilter.active    = 1;        % [0 1]; 1 = apply low pass filter to power density spectrum to remove spurious peaks due to noise
signalProcessing.lowPassPSDFilter.params    = [5 0.001];  % [int 0.01-0.10] butterworth parameters [Xth_order cutoff_freq]

signalProcessing.lowPassHRTimeFilter.active = 1;        % [0 1]; 1 = apply low pass filter to heart rate over time to remove spurious changes in HR due to noise
signalProcessing.lowPassHRTimeFilter.params = [8 0.02]; % [int 0.01-0.10] butterworth parameters [Xth_order cutoff_freq]
            
[signalProcessing.highPassPixelFilter.NthOrder,signalProcessing.highPassPixelFilter.cutOffFreq] = butter(signalProcessing.highPassPixelFilter.params(1),signalProcessing.highPassPixelFilter.params(2));
[signalProcessing.lowPassPSDFilter.NthOrder,signalProcessing.lowPassPSDFilter.cutOffFreq]       = butter(signalProcessing.lowPassPSDFilter.params(1),signalProcessing.lowPassPSDFilter.params(2));
[signalProcessing.lowPassHRTimeFilter.NthOrder,signalProcessing.lowPassHRTimeFilter.cutOffFreq] = butter(signalProcessing.lowPassHRTimeFilter.params(1),signalProcessing.lowPassHRTimeFilter.params(2));

signalProcessing.ica.nComps                  = 3;        %
signalProcessing.ica.nIte                    = 2000;     %
signalProcessing.ica.stab                    = 'on';     %
signalProcessing.ica.verbose                 = 'off';    %

signalProcessing.tfa.method                 = 'plomb';  % time-frequency analysis method: 'stft' or 'plomb'
signalProcessing.tfa.computeCoherence       = 1; % 0 = raw data, 1 = coherence transformed
signalProcessing.tfa.winSize                = 10;   % X seconds
signalProcessing.tfa.tempRes                = 240;  % temporal resolution, number of time points - THIS NUMBER NEEDS TO BE CONVERTED FROM FRAMES TO FRAMES/SECOND IN A NEXT RELEASE
signalProcessing.tfa.freqRes                = 120;  % frequency resolution, number of frequencies


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

%% SKIN COLOR TO PULSE SIGNAL 

finalSignal = struct();
if strcmp(signalProcessing.HrDetectionMethod,'fastica') % fastICA on different color channels
    
    finalSignal.comp                = fastica(pixelVal_filt','numOfIC',signalProcessing.ica.nComps,'maxNumIterations',signalProcessing.ica.nIte,'stabilization',signalProcessing.ica.stab,'verbose',signalProcessing.ica.verbose);

elseif strcmp(signalProcessing.HrDetectionMethod,'ica') 
    [W,finalSignal.comp] = ica(pixelVal_filt',signalProcessing.ica.nComps); % JADE ICA - J. F. Cardoso 1997, G. D. Clifford, MIT, 2004

elseif strcmp(signalProcessing.HrDetectionMethod,'pos') % The Plane Orthogonal to Skin-Tone (POS) Method 
    
    nFrames                     = size(pixelVal_filt,1);
    pixelVal_filt_colorRotated  = zeros(1,nFrames); 
    nSecWindow                  = 1.6;
    nWindowFrames = ceil(nSecWindow*signalProcessing.samplingRate); 
    for startFrameIdx = 1:nFrames-nWindowFrames
        endFrameIdx = startFrameIdx + nWindowFrames;
        
        tempPixelVal_filt_detrended     = (pixelVal_filt(startFrameIdx:endFrameIdx,:)./mean(pixelVal_filt(startFrameIdx:endFrameIdx,:)))'; % divide by mean to center around zero. This will bias to smaller numbers though .. maybe minus mean() is better?
        tempColorRotated                = [0, 1, -1; -2, 1, 1] * tempPixelVal_filt_detrended; % matrix multiplication to rotate color space to strengthen blood induced pulses ... ends up with 2 instead of 3 dimensions
        tempColorRotated_balanced       = tempColorRotated(1,:) + ((std(tempColorRotated(1,:)) / std(tempColorRotated(2,:))) * tempColorRotated(2,:)); % equalize variance in both dimensions and add them up
        pixelVal_filt_colorRotated(startFrameIdx:endFrameIdx) = pixelVal_filt_colorRotated(startFrameIdx:endFrameIdx) + (tempColorRotated_balanced - mean(tempColorRotated_balanced)); % detrend resulting signal and add up to signal

    end

    finalSignal.comp = [pixelVal_filt_colorRotated; pixelVal_filt_colorRotated; pixelVal_filt_colorRotated];
end

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

subplot(2,2,3)
plot(signalProcessing.fft.HRRange,sum(finalSignal.powerVal(signalProcessing.fft.fRange2,:)')./max(sum(finalSignal.powerVal(signalProcessing.fft.fRange2,:)')),'k');
hold on
plot(signalProcessing.fft.HRRange,sum(finalSignal.coherenceVal(signalProcessing.fft.fRange2,:)')/max(sum(finalSignal.coherenceVal(signalProcessing.fft.fRange2,:)')),'b');
line([repmat(finalSignal.HR_coherenceBased,1,2)],[0 1],'Color','k','LineStyle',':','LineWidth',2)
xlabel('HR')
ylabel('Sum across components')
legend('Power','Coherence')

subplot(2,2,4);
plot(resampledXdata,finalSignal.comp(finalSignal.HRBestCompIdx_coherenceBased,:),'k')
title('Best component - coherence based')
ylabel('Pixel val')
xlabel('Time [s]')

disp(['HR in power spectrum: ' num2str(finalSignal.HR_powerBased)])
disp(['HR in coherence spectrum: ' num2str(finalSignal.HR_coherenceBased)])



%% Time frequency analysis

signalProcessing.tfa.winSizeFr = signalProcessing.tfa.winSize*ceil(signalProcessing.samplingRate); % transform from seconds to frames

if signalProcessing.tfa.winSizeFr >= length(finalSignal.comp)
    disp(['ERROR: variable signalProcessing.tfa.winSize is set to ' num2str(signalProcessing.tfa.winSize) ' seconds,'])
    disp(['which is as equal to or larger than the duration of the video'])
    disp('RECOMMENDATION: reduce the sliding window size to half the duration or even smaller (in seconds)')

end

win     = blackman(signalProcessing.tfa.winSizeFr, 'periodic')';
tBins   = floor(linspace(1,length(finalSignal.comp)-signalProcessing.tfa.winSizeFr,signalProcessing.tfa.tempRes));
fvec    = linspace(signalProcessing.FOI(1)/60,signalProcessing.FOI(2)/60,signalProcessing.tfa.freqRes);
nfft = signalProcessing.tfa.freqRes*(signalProcessing.samplingRate/2);
signalProcessing.tfa.freq       = signalProcessing.samplingRate/nfft*(0:nfft-1);
signalProcessing.tfa.fRange2    = find(signalProcessing.tfa.freq>=signalProcessing.fft.freqInterestRange(1) & signalProcessing.tfa.freq<signalProcessing.fft.freqInterestRange(2));


timeFreqData = NaN(length(fvec),length(tBins));
tData = NaN(1,length(tBins));
countT = 0;
for tt = tBins
    countT = countT+1;
    t = tt:tt+signalProcessing.tfa.winSizeFr-1;

    if strcmp(signalProcessing.tfa.method,'plomb')    
        [timeFreqData(:,countT),fvec2] = plomb(finalSignal.comp(finalSignal.HRBestCompIdx_powerBased,t),resampledXdata(t),fvec);
    elseif strcmp(signalProcessing.tfa.method,'stft')    
        
        Y = fft(finalSignal.comp(finalSignal.HRBestCompIdx_powerBased,t).*win,nfft);
        Y = Y(1:ceil(end/2));
        Y = Y(signalProcessing.tfa.fRange2);
        
        fvec2   = signalProcessing.tfa.freq(signalProcessing.tfa.fRange2);
        
        timeFreqData(:,countT) = Y.*conj(Y)/signalProcessing.tfa.freqRes;
        
    end
    
    tData(countT) = mean(resampledXdata(t));

    if signalProcessing.tfa.computeCoherence
        timeFreqData(:,countT) = timeFreqData(:,countT)./sqrt(sum(timeFreqData(:,countT)).^2);
    end

end

% signalProcessing.lowPassHRTimeFilter.params = [8 0.03];  % [int 0.01-0.10] butterworth parameters [Xth_order cutoff_freq]
% [signalProcessing.lowPassHRTimeFilter.NthOrder,signalProcessing.lowPassHRTimeFilter.cutOffFreq] = butter(signalProcessing.lowPassHRTimeFilter.params(1),signalProcessing.lowPassHRTimeFilter.params(2));
       %%
       
figure();
imagesc(timeFreqData)
hold on
[~,maxIdx] = max(timeFreqData);
line(1:length(maxIdx),maxIdx,'Color','g','LineWidth',3,'LineStyle',':')

if signalProcessing.lowPassHRTimeFilter.active
    maxIdx_filt = filtfilt(signalProcessing.lowPassHRTimeFilter.NthOrder,signalProcessing.lowPassHRTimeFilter.cutOffFreq,maxIdx);
end
line(1:length(maxIdx_filt),maxIdx_filt,'Color','b','LineStyle','--','LineWidth',3)

set(gca,'xtick',round(linspace(1,length(tBins),5)))
set(gca,'xticklabel',round(tData(round(linspace(1,length(tBins),5)))),'FontSize',12)

set(gca,'ytick',round(linspace(1,length(fvec),10)))
set(gca,'yticklabel',round(60*fvec(round(linspace(1,length(fvec),10)))),'FontSize',12)

xlabel('Time (s)','FontSize',16)
ylabel('Heart rate (bpm)','FontSize',16)
colormap('hot')

disp(['median HR across time  - max fit: ' num2str(median(60*fvec2(ceil(maxIdx))))])
disp(['median HR across time  - max fit smoothed: ' num2str(median(60*fvec2(ceil(maxIdx_filt))))])


tempData = zeros(prod(size(timeFreqData)),3);
for k = 1:size(tempData,1)
    tempData(k,:) = [ceil(k/size(timeFreqData,1)) mod(k-1,size(timeFreqData,1))+1 timeFreqData(k)];
end

selectData = tempData(:,3) > prctile(tempData(:,3),95);
tempData = tempData(selectData,:);

weights = sqrt(tempData(:,3));

%     tempData(:,1) = tData(tempData(:,1));
%     tempData(:,2) = fvec(tempData(:,2));
%     P = polyfitweighted(tempData(:,1),tempData(:,2),1,weights);
%     smoothx = linspace(min(tempData(:,1)),max(tempData(:,1)),100);
%     smoothy = polyval(P,smoothx);
%
tempData(:,1) = tData(tempData(:,1));
tempData(:,2) = fvec(tempData(:,2));
P = polyfitweighted(tempData(:,1),tempData(:,2),1,weights);
smoothx = linspace(min(tempData(:,1)),max(tempData(:,1)),size(timeFreqData,2));
smoothy = polyval(P,smoothx);

% figure()
% plot(smoothx,60*smoothy,'w','LineWidth',2)
plot([1:length(smoothy)],(smoothy-fvec(1))*(size(timeFreqData,1)/(fvec(end)-fvec(1))),'w','LineWidth',2)
legend('\color{green} Max','\color{blue} Smooth max','\color{white} SNR weighted','Location','SouthEast','Box','off','FontSize',16)
% legend boxoff
disp(['median HR across time - SNR weighted fit: ' num2str(median(60*smoothy))])

%%
figure();
plot(tData,60*fvec2(ceil(maxIdx)),'g:')
hold on
plot(tData,60*fvec2(ceil(maxIdx_filt)),'b--')
%set(gca,'yticklabel',[50:10:120],'FontSize',12)
plot(tData,60*smoothy,'r')
ylabel('Heart Rate (BPM)','FontSize',16)
xlabel('Time (s)','FontSize',16)
legend('Max','Smooth max','SNR weighted','FontSize',16)


%% 
