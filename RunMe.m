
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

mkdir('Plots')
mkdir('Raw data')
mkdir('results')
mkdir('ResultsExcel')
mkdir('videos')
    
%% get script folder
clearvars -except scriptRun
close all; clc;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename)); clear tmp
scriptFolder = pwd;
%% set parameters

chooseSettings
% videosFolder = uigetdir;
files = dir(fullfile(videosFolder, '*.mp4'));

finalVarNames = {'sbjName',  ...
    'HR_powerBased', 'HR_coherenceBased', 'maxIdx_filt','smoothy','s2nrat', ...
    'nFrames', 'vidDuration', 'frameRate', 'threshold', 'warningMess', ...
    'nFacesDetected', 'enoughPointsVisible', 'luminance', 'flag', 'nPeaksAboveThreshold', ...
    'maxFreqDiff'};
%% Start main loop
thisPlot = 0;
for thisVideo = 1:length(files)
    try
        videoFileName = [videosFolder filesep files(thisVideo).name];
        tic
        
        filename = {files(thisVideo).name};
        sbjName = erase(files(thisVideo).name, '.avi');
        sbjName = erase(files(thisVideo).name, '.mp4');
        sbjName = char(sbjName);
        
        %     condition = extractAfter(videoFileName, '-'); condition = extractBefore(condition, '-');
        
        cprintf(-[1,0,1], strcat('This is iteration: ', num2str(thisVideo), '\n'));
        disp(['Processing Subject: ' sbjName])
        
        %     disp(['This is iteration: ' num2str(thisVideo)])
        %     disp(['Processing Subject: ' videoFileName])
        
        % videoFileName                               = 'rPPG_video.mp4';
        
        signalProcessing = struct();
        signalProcessing.HrDetectionMethod          = 'pos';    % 'ica' (see Van der Kooij & Naber, 2019); 'pos' (see Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G., 2017)
        
        signalProcessing.samplingRate               = 60;       % [frames per second] sampling rate: temporal resolution of pixel value signal will increase with interpolation to X Hz
        signalProcessing.interpMethod               = 'pchip';  % rPPG signal is always interpolated to a frequency of the sampling rate
        signalProcessing.FOI                        = [45 165]; % range: [min max] frequency of interest of heart rate in Beats per minute (BPM)
        
        signalProcessing.highPassPixelFilter.active = 1;        % [0 1]; 1 = apply low pass filter to pixel values ... to remove artifacts by movement or illumination
        signalProcessing.highPassPixelFilter.params = [6 (signalProcessing.FOI(1)/60)/(signalProcessing.samplingRate/2)];     % [int 0.01-0.10] butterworth parameters --> [6 0.04] is ideal for frame rate of 30
        
        signalProcessing.lowPassPSDFilter.active    = 0;        % [0 1]; 1 = apply low pass filter to power density spectrum to remove spurious peaks due to noise
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
        signalProcessing.tfa.winSize                = winSize;   % X seconds
        signalProcessing.tfa.tempRes                = 240;  % temporal resolution, number of time points - THIS NUMBER NEEDS TO BE CONVERTED FROM FRAMES TO FRAMES/SECOND IN A NEXT RELEASE
        signalProcessing.tfa.freqRes                = 120;  % frequency resolution, number of frequencies
        
        
        %% extract RGB values of facial skin per frame
        
        if overwriteRGB | ~exist([scriptFolder filesep 'Raw Data' filesep expPlotsFolder filesep sbjName '.mat'], 'file')
            [pixelValPerFrame,faceMap,vidInfo,faceDetection,faceTracking] = v2_extractFaceFromVideo(videoFileName,'all', makeFaceVid, showVid, sbjName, expPlotsFolder);
            save([scriptFolder filesep 'Raw data' filesep expPlotsFolder filesep sbjName '.mat'],'pixelValPerFrame','vidInfo')
        else
            load([scriptFolder filesep 'Raw data' filesep expPlotsFolder filesep sbjName '.mat'],'pixelValPerFrame','vidInfo')
        end
        
        
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
        if strcmp(signalProcessing.HrDetectionMethod,'ica')
            for c = 1:3
                if signalProcessing.highPassPixelFilter.active
                    pixFilter = filtfilt(signalProcessing.highPassPixelFilter.NthOrder,signalProcessing.highPassPixelFilter.cutOffFreq,resampledYdata(:,c));
                    pixelVal_filt(:,c) = resampledYdata(:,c)-pixFilter;
                else
                    pixelVal_filt(:,c) = resampledYdata(:,c);
                end
            end
        end
        
        %% SKIN COLOR TO PULSE SIGNAL
        
        if strcmp(signalProcessing.HrDetectionMethod,'ica') % ICA on different color channels
            finalSignal = struct();
            finalSignal.comp                = fastica(pixelVal_filt','numOfIC',signalProcessing.ica.nComps,'maxNumIterations',signalProcessing.ica.nIte,'stabilization',signalProcessing.ica.stab,'verbose',signalProcessing.ica.verbose);
            
            %% POS on different color channels
            
        elseif strcmp(signalProcessing.HrDetectionMethod,'pos') % The Plane Orthogonal to Skin-Tone (POS) Method
            
            WinSec = 1.6; % Window in seconds
            
            %lines and comments correspond to pseudo code algorithm on reference page 7
            N = size(pixelVal_filt,1);%line 0 - RGB is of length N frames
            H = zeros(1,N);%line 1 - initialize to zeros of length of video sequence
            l = ceil(WinSec*signalProcessing.samplingRate);%line 1 - window length equivalent to reference:
            C = zeros(length(l),3);
            for n = 1:N-1%line 2 - loop from first to last frame in video sequence
                %line 3 - spatial averaging was performed when video was read
                m = n - l + 1;%line 4 condition
                if(m > 0)%line 4
                    Cn = ( pixelVal_filt(m:n,:) ./ mean(pixelVal_filt(m:n,:)) )';%line 5 - temporal normalization
                    S = [0, 1, -1; -2, 1, 1] * Cn;%line 6 - projection
                    h = S(1,:) + ((std(S(1,:)) / std(S(2,:))) * S(2,:));%line 7 - tuning
                    H(m:n) = H(m:n) + (h - mean(h));%line 8 - overlap-adding
                end%line 9 - end if
            end%line 10 - end for
            
            finalSignal.comp=[H; H; H];
            
            if signalProcessing.highPassPixelFilter.active
                
                for c = 1:3
                    pixFilter = filtfilt(signalProcessing.highPassPixelFilter.NthOrder,signalProcessing.highPassPixelFilter.cutOffFreq,finalSignal.comp(c,:)');
                    finalSignal.comp(c,:) = finalSignal.comp(c,:)-pixFilter';
                end
            end
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
        
        %% compute coherence values (takes into account noise in spectrum)
        
        finalSignal.coherenceVal = finalSignal.powerVal;
        
        for co = 1:signalProcessing.ica.nComps
            finalSignal.coherenceVal(:,co) = finalSignal.powerVal(:,co)./sqrt(sum(finalSignal.powerVal(signalProcessing.fft.fRange2,co).^2));
        end
        
        %% Low-pass filter fourier spectrum
        
        finalSignal.powerVal_filt = finalSignal.powerVal;
        finalSignal.coherenceVal_filt = finalSignal.coherenceVal;
        for co = 1:signalProcessing.ica.nComps
            if signalProcessing.lowPassPSDFilter.active
                finalSignal.powerVal_filt(:,co) = filtfilt(signalProcessing.lowPassPSDFilter.NthOrder,signalProcessing.lowPassPSDFilter.cutOffFreq,finalSignal.powerVal(:,co));
                finalSignal.coherenceVal_filt(:,co) = filtfilt(signalProcessing.lowPassPSDFilter.NthOrder,signalProcessing.lowPassPSDFilter.cutOffFreq,finalSignal.coherenceVal(:,co));
            end
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
        [~,tempMaxPowerCompIdx]         = max(finalSignal.maxPower);
        [~,tempMaxCoherenceCompIdx]     = max(finalSignal.maxCoherence);
        
        finalSignal.HRfreqIdx_powerBased        = tempMaxPowerIdx(tempMaxPowerCompIdx);
        finalSignal.HR_powerBased               = signalProcessing.fft.freq(signalProcessing.fft.fRange2(tempMaxPowerIdx(tempMaxPowerCompIdx)))*60;
        finalSignal.HR_powerBased_PerComp       = signalProcessing.fft.freq(signalProcessing.fft.fRange2(tempMaxPowerIdx))*60;
        finalSignal.HRBestCompIdx_powerBased    = tempMaxPowerCompIdx;
        
        finalSignal.HRfreqIdx_coherenceBased        = tempMaxCoherenceIdx(tempMaxCoherenceCompIdx);
        finalSignal.HR_coherenceBased               = signalProcessing.fft.freq(signalProcessing.fft.fRange2(tempMaxCoherenceIdx(tempMaxCoherenceCompIdx)))*60;
        finalSignal.HR_PerComp_coherenceBased       = signalProcessing.fft.freq(signalProcessing.fft.fRange2(tempMaxCoherenceIdx))*60;
        finalSignal.HRBestCompIdx_coherenceBased    = tempMaxCoherenceCompIdx;
        
        %% compute signal to noise ratio
        % see De Haan, G., & Jeanne, V. (2013). Robust pulse rate from chrominance-based rPPG. IEEE Transactions on Biomedical Engineering, 60(10), 2878-2886.
        
        binWidthSpec    = 0.5;
        NyquistFreq     = signalProcessing.samplingRate/2;
        NFFT            = (NyquistFreq*2*60)/binWidthSpec; % number of bins in power spectrum
        
        [powerSpec,freq] = periodogram(finalSignal.comp(tempMaxPowerCompIdx,:),hamming(length(resampledXdata)),NFFT,signalProcessing.samplingRate);
        
        selectRange = (freq >= signalProcessing.fft.freqInterestRange(1))&(freq <= signalProcessing.fft.freqInterestRange(2));
        % plot(freq(selectRange)*60,powerSpec(selectRange))
        
        hrFreq          = finalSignal.HR_powerBased/60;
        sumPowerHR      = sum(powerSpec((freq >= hrFreq-0.1)&(freq <= hrFreq+0.1)));
        sumPowerSpec    = sum(powerSpec((freq >= signalProcessing.fft.freqInterestRange(1))&(freq <= signalProcessing.fft.freqInterestRange(2))));
        s2nrat          = pow2db(sumPowerHR/(sumPowerSpec-sumPowerHR));
        
        %% plot frequency spectrum
        checkIfWantToPlot(plotYesNo)
        
        subplot(2,2,1);
        plot(signalProcessing.fft.HRRange,finalSignal.powerVal_filt(signalProcessing.fft.fRange2,:));
        hold on
        plot(signalProcessing.fft.HRRange,finalSignal.powerVal(signalProcessing.fft.fRange2,:),'LineStyle',':');
        line([repmat(finalSignal.HR_powerBased,1,2)],[0 max(finalSignal.powerVal(signalProcessing.fft.fRange2,finalSignal.HRBestCompIdx_powerBased))],'Color','k','LineStyle',':','LineWidth',2)
        xlabel('HR')
        ylabel('Power')
        
        subplot(2,2,2)
        ha = plot(signalProcessing.fft.HRRange,finalSignal.coherenceVal_filt(signalProcessing.fft.fRange2,:));
        hold on
        plot(signalProcessing.fft.HRRange,finalSignal.coherenceVal(signalProcessing.fft.fRange2,:),'LineStyle',':');
        line([repmat(finalSignal.HR_coherenceBased,1,2)],[0 max(finalSignal.coherenceVal(signalProcessing.fft.fRange2,finalSignal.HRBestCompIdx_coherenceBased))],'Color','k','LineStyle',':','LineWidth',2)
        xlabel('HR')
        ylabel('Coherence')
        legend(ha,{'Component 1','Component 2','Component 3'})
        
        subplot(2,2,3)
        plot(signalProcessing.fft.HRRange,sum(finalSignal.powerVal_filt(signalProcessing.fft.fRange2,:)')./max(sum(finalSignal.powerVal_filt(signalProcessing.fft.fRange2,:)')),'k');
        hold on
        plot(signalProcessing.fft.HRRange,sum(finalSignal.powerVal(signalProcessing.fft.fRange2,:)')./max(sum(finalSignal.powerVal(signalProcessing.fft.fRange2,:)')),'k:');
        plot(signalProcessing.fft.HRRange,sum(finalSignal.coherenceVal_filt(signalProcessing.fft.fRange2,:)')/max(sum(finalSignal.coherenceVal_filt(signalProcessing.fft.fRange2,:)')),'b');
        plot(signalProcessing.fft.HRRange,sum(finalSignal.coherenceVal(signalProcessing.fft.fRange2,:)')/max(sum(finalSignal.coherenceVal(signalProcessing.fft.fRange2,:)')),'b:');
        line([repmat(finalSignal.HR_coherenceBased,1,2)],[0 1],'Color','k','LineStyle',':','LineWidth',2)
        xlabel('HR')
        ylabel('Sum across components')
        legend('Power','Coherence')
        
        subplot(2,2,4);
        plot(resampledXdata,finalSignal.comp(finalSignal.HRBestCompIdx_coherenceBased,:),'k')
        title('Best component - coherence based')
        ylabel('Pixel val')
        xlabel('Time [s]')
        
        if ~(exist('plots', 'dir'))
            mkdir('plots');
        end
        cd 'plots';
        
        if ~(exist(expPlotsFolder, 'dir'))
            mkdir(expPlotsFolder);
        end
        cd(expPlotsFolder);
        
        if ~(exist(sbjName, 'dir'))
            mkdir(sbjName);
        end
        cd(sbjName)
        thisPlot = thisPlot+1;
        saveas(gcf,[sbjName '_' num2str(thisPlot) '_frequency_spectrum'],'jpg')
        cd(scriptFolder)
        
        disp(['HR in power spectrum: ' num2str(finalSignal.HR_powerBased)])
        disp(['HR in coherence spectrum: ' num2str(finalSignal.HR_coherenceBased)])
        
        %% check if the first and second peak have similar power and are distant from one another
        checkIfWantToPlot(plotYesNo)
        powerSpectrum = finalSignal.powerVal_filt(signalProcessing.fft.fRange2,1);
        
        plot(signalProcessing.fft.HRRange, powerSpectrum)
        
        [PkAmp, PkFrequency] = findpeaks(abs(powerSpectrum), signalProcessing.fft.HRRange); %find amplitude ant time of all peaks
        [~,idx] = sort(PkAmp,'descend'); % get the index of the peaks by amplitude
        
        maxPkAmp = max(PkAmp); % maximum amplitude peak
        PkAmp1Quarter = maxPkAmp/4; % divide the max amplitude peak by 4/4
       
        nPeaksAboveThreshold = sum(PkAmp > PkAmp1Quarter*3); % find number of peaks above threshold

        peaksAboveThreshold = PkAmp(PkAmp > PkAmp1Quarter*3); % select only the peaks above threshold
        freqAboveThreshold = PkFrequency(PkAmp > PkAmp1Quarter*3)';
        
        peaksFreq = [peaksAboveThreshold freqAboveThreshold];
        
        [~, rowMaxPeak] = max(peaksFreq(:,1));
        highestFreq = peaksFreq(rowMaxPeak,2);
        maxFreqDiff = max(abs(peaksFreq(:,2)-highestFreq));
        
        
%         [sortedPeaksAboveThreshold, sortedPeakOrder] = sort(peaksAboveThreshold,'descend'); % sort from highest to smallest
        [sortedFreqAboveThreshold, ~] = sort(freqAboveThreshold,'descend'); % time sorted based on amplitude order
        
        % Plot peaks above threshold
        hold on
        for thisPeak = 1:nPeaksAboveThreshold
            plot([PkFrequency(idx(thisPeak)) PkFrequency(idx(thisPeak))], ...
                [0 PkAmp(idx(thisPeak))],'Color','r','LineStyle', '--', 'LineWidth',2)
        end
        
        % Plot lines dividing max peak amplitude in 4/4
        plot([min(signalProcessing.fft.HRRange) max(signalProcessing.fft.HRRange)],...
            [PkAmp1Quarter PkAmp1Quarter],'Color','b','LineStyle', '--', 'LineWidth',1)
        
        plot([min(signalProcessing.fft.HRRange) max(signalProcessing.fft.HRRange)], ...
            [PkAmp1Quarter*2 PkAmp1Quarter*2],'Color','b','LineStyle', '--', 'LineWidth',1)
                
        plot([min(signalProcessing.fft.HRRange) max(signalProcessing.fft.HRRange)],...
            [PkAmp1Quarter*3 PkAmp1Quarter*3],'Color','b','LineStyle', '--', 'LineWidth',1)
        
        % Save plots
        cd(['plots' filesep expPlotsFolder filesep sbjName])
        thisPlot = thisPlot+1;
        saveas(gcf,[sbjName '_' num2str(thisPlot) '_frequency_spectrum_peakDiff'],'jpg')
        cd(scriptFolder)
        
        % Check if this power spectrum has to be flagged
       diffMostDistancePeak = abs(sortedFreqAboveThreshold(1) - sortedFreqAboveThreshold(end));

       if nPeaksAboveThreshold > 1
           flag = 1;
           if maxFreqDiff > 10
               flag = 2;
           end
       elseif nPeaksAboveThreshold <= 1
           flag = 0;
       end
       
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
        %     if strcmp(plotYesNo, 'Yes')
        checkIfWantToPlot(plotYesNo)
        
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
        %     end
        disp(['median HR across time: ' num2str(median(60*fvec2(ceil(maxIdx_filt))))])
        
        %         cd([scriptFolder filesep 'plots' filesep expPlotsFolder filesep sbjName])
        %         saveas(gcf,[sbjName '_2 spectrum_2'],'jpg')
        %
        %         cd(scriptFolder)
        
        %%
        %         checkIfWantToPlot(plotYesNo)
        tempData = zeros(prod(size(timeFreqData)),3);
        for k = 1:size(tempData,1)
            tempData(k,:) = [ceil(k/size(timeFreqData,1)) mod(k-1,size(timeFreqData,1))+1 timeFreqData(k)];
        end
        
        selectData = tempData(:,3) > prctile(tempData(:,3),95);
        tempData = tempData(selectData,:);
        
        weights = sqrt(tempData(:,3));
        
        %     tempData(:,1) = tData(tempData(:,1));
        %     tempData(:,2) = fvec(tempData(:,2));
        P = polyfitweighted(tempData(:,1),tempData(:,2),1,weights);
        smoothx = linspace(min(tempData(:,1)),max(tempData(:,1)),100);
        smoothy = polyval(P,smoothx);
        
        plot(smoothx,smoothy,'w','LineWidth',2)
        
        tempData(:,1) = tData(tempData(:,1));
        tempData(:,2) = fvec(tempData(:,2));
        P = polyfitweighted(tempData(:,1),tempData(:,2),1,weights);
        smoothx = linspace(min(tempData(:,1)),max(tempData(:,1)),100);
        smoothy = polyval(P,smoothx);
        legend('MaxIdx', 'MaxIdxFilt', 'Smoothy', 'Location', 'southeast')
        cd(['plots' filesep expPlotsFolder filesep sbjName])
        saveas(gcf,[sbjName '_smooth'],'jpg')
        cd(scriptFolder)
        
        disp(['median HR across time - alternative calculation: ' num2str(median(60*smoothy))])
        
        %     p = polyfit(tempData(:,1),tempData(:,2),1);
        
        %% Plot HR
        checkIfWantToPlot(plotYesNo)
        
        plot(tData,60*fvec2(ceil(maxIdx_filt)))
        hold on
        plot(smoothx,60*smoothy)
        ylabel('Heart Rate (BPM)')
        xlabel('Time (s)')
        
        cd(['plots' filesep expPlotsFolder filesep sbjName])
        saveas(gcf,[sbjName '_3 Hr_Time'],'jpg')
        cd(scriptFolder)
        
        %% GET LUMINOSITY
        pixelPerFrameNoNaN = pixelValPerFrame;
        pixelPerFrameNoNaN(any(isnan(pixelPerFrameNoNaN), 2), :) = [];
        
        normalizedPixelVal = pixelPerFrameNoNaN/255;
        meanLuminancePerChan = mean(normalizedPixelVal,1);
        meanLuminanceTot = mean(meanLuminancePerChan,2);
        
        
        %% Write results
        if ~exist('faceDetection', 'var')
            faceDetection.threshold = NaN;
            faceDetection.warningMess = NaN;
            faceDetection.nFacesDetected = NaN;
            faceTracking.enoughPointsVisible = NaN;
        end
        
        if thisVideo == 1 || ~exist('allSbjResults','var')
            % if you chose to overwrite the data, create a new
            % allSbjResults
            thisSbjResults = table(convertCharsToStrings(sbjName), ...
                finalSignal.HR_powerBased, finalSignal.HR_coherenceBased, median(60*fvec2(ceil(maxIdx_filt))),...
                median(60*smoothy), s2nrat, vidInfo.nFrames, vidInfo.vidDuration, vidInfo.frameRate, ...
                faceDetection.threshold, faceDetection.warningMess, faceDetection.nFacesDetected, ...
                faceTracking.enoughPointsVisible, meanLuminanceTot, flag, nPeaksAboveThreshold, maxFreqDiff);
            thisSbjResults.Properties.VariableNames = finalVarNames;
            
            allSbjResults = table(convertCharsToStrings(sbjName), ...
                finalSignal.HR_powerBased, finalSignal.HR_coherenceBased, median(60*fvec2(ceil(maxIdx_filt))), ...
                median(60*smoothy), s2nrat, vidInfo.nFrames, vidInfo.vidDuration, vidInfo.frameRate, ...
                faceDetection.threshold, faceDetection.warningMess, faceDetection.nFacesDetected, ...
                faceTracking.enoughPointsVisible, meanLuminanceTot, flag, nPeaksAboveThreshold, maxFreqDiff);
            allSbjResults.Properties.VariableNames = finalVarNames;
            
        else
            thisSbjResults = table(convertCharsToStrings(sbjName), ...
                finalSignal.HR_powerBased, finalSignal.HR_coherenceBased, median(60*fvec2(ceil(maxIdx_filt))),...
                median(60*smoothy), s2nrat, vidInfo.nFrames, vidInfo.vidDuration, vidInfo.frameRate, ...
                faceDetection.threshold, faceDetection.warningMess, faceDetection.nFacesDetected, ...
                faceTracking.enoughPointsVisible, meanLuminanceTot, flag, nPeaksAboveThreshold, maxFreqDiff);
            thisSbjResults.Properties.VariableNames = finalVarNames;
            
            allSbjResults = [allSbjResults; thisSbjResults];
        end
        
        
        % Save each sbj individually
        cd([scriptFolder filesep 'results'])
        if ~exist([sbjName '.mat'], 'file')
            %             thisSbjResults.Properties.VariableNames =finalVarNames;
            save([sbjName '.mat'], 'thisSbjResults');
        end
        
        
        % save "ECG" traces reconstructed from video
        if ~(exist('ECGs ExtractedFromVideo', 'dir'))
            mkdir('ECGs ExtractedFromVideo');
        end
        
        cd(scriptFolder)
        
        thisSbjElapsedTime = toc;
        disp(['This sbj took ' num2str(thisSbjElapsedTime) ' seconds'])
        
        %     clearvars -except files plotYesNo thisVideo allSbjResults videosFolder overwriteRGB
        
        close all; close all hidden
        
    catch ME
        cprintf('*red', strcat('Script crushed on sbj: ', sbjName, '\n')); % Show the name of the file that's not working
        cd(videosFolder) % go into the folder with the videos
        if ~(exist('Not Working', 'dir')) % make sure tha there's a folder called "not working"
            mkdir('Not Working');
        end
        % Move the file that's not working to the Not Working folder
        movefile(char(filename), char([videosFolder filesep 'Not Working']))
        cd([scriptFolder filesep 'results'])
        
        %         allSbjResults.Properties.VariableNames = finalVarNames;
        %         % Save excel file
        %         writetable(allSbjResults, ['allSbjResults_' date '.xlsx']);
        %
        %         % Save entire workspace
        %         save(['allSbjResults_'  date '.mat']);
        cd(scriptFolder)

        cprintf(-[0,0,1], strcat('There was an error. ', '\n'));
        %         Gluce_RunMe_MN_GF_11_Jun_2021
        %         rethrow(ME)
    end
end

%% Save data
disp('All done. Saving the data. Let''s publish on on Nature')
allSbjResults.Properties.VariableNames = finalVarNames;

timeNow = char(datetime('now','Format','HH:mm:ss')); % get minutes and seconds
timeNow = strrep(timeNow,':','_');

if ~(exist('ResultsExcel', 'dir'))
    mkdir('ResultsExcel');
end

cd([scriptFolder filesep 'ResultsExcel'])

writetable(allSbjResults, [resultsFilename '_' date '_' timeNow '.xlsx']); % save all results to xlsx
writetable(allSbjResults, [resultsFilename '_' date '_' timeNow]) % save all results to txt
save([resultsFilename '_' date '_' timeNow]); % save entire workspace
%% Additional info
scriptRunning = mfilename('fullpath');
additionalInfo = table(string(scriptRunning), string(date));
additionalInfo = [additionalInfo, scriptSettings];
additionalInfo.Properties.VariableNames{1} = 'scriptUsedForTheseData';
additionalInfo.Properties.VariableNames{2} = 'date';
% Write the info the the results xlsx file
writetable(additionalInfo,[resultsFilename '_' date '_' timeNow '.xlsx'], 'Sheet',2);





