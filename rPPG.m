function [icaMaxPower, icaBestHR, rgbMaxPower, rgbBestHR, power_ica, phase_ica, power_rgb, phase_rgb, meanPVperFrame, meanPVperFrame_filt, avePixelValue, aveMotion] = rPPG(fileName,TOI,FOI,lpFiltParams,lpFreqFiltParams,reselectROI,reselectSkinColor,hueSatMapResolution,moveDetReso,outputValuesAndPlots)


% [icaMaxPower, icaBestHR, rgbMaxPower, rgbBestHR, power_ica, phase_ica, power_rgb, phase_rgb, avePixelValue, aveMotion] = rPPG(fileName,TOI,FOI,lpFiltParams,lpFreqFiltParams,reselectROI,reselectSkinColor,hueSatMapResolution,moveDetReso,outputValuesAndPlots)
% 
% This script performs a remote photoplethysmography (rPPG) analysis on
% videos with algorithms as described in the following publication:
% 
% van der Kooij & Naber (2017). Standardized procedures for the testing and 
% reporting of remote heart rate imaging. Journal title, etc.
%
% 
% ---------------- INPUT -------------------
% 
% fileName: string. File name of the video including the extension (e.g., *.avi)
% 
% TOI: 'all', '', or [A B] in seconds. Time of interest --> 'all' and '' (empty char)
% will analyze all available frames. When entered two doubles [A B], this
% script will only analyze the frames that fall within the period from A 
% seconds till B seconds. Default = 'all';
%
% FOI: 'all', '', or [A B] in beats per minute. Frequency of interest --> 
% 'all' and '' (empty char) will analyze all available heart rate
% frequencies. When entered two doubles [A B], this script will only
% analyze the frequencies that fall within the range from A BPM till B BPM. 
% Default = [45 165];
%
% lpFiltParams: [A B]. Low pass filter on raw pixel value traces to remove
% artifacts caused by motion or changes in light. A is the Ath order of the 
% low pass butterworth filter and has to be an integer between 0 and 10. 
% B is the cut-off frequency and has to be a fraction. Default: [6 0.04]
%
% lpFreqFiltParams: [A B]. Low pass filter on power frequency spectrum 
% to take into account changes in heart rate and thus multiple, individual 
% peaks in the power frequency spectrum. A is the Ath order of the 
% low pass butterworth filter and has to be an integer between 0 and 10. 
% B is the cut-off frequency and has to be a fraction. Default: [3 0.2]
% 
% reselectROI: 0 or 1. A region of interest within the video is chosen 
% (frame cropping) to remove background and irrelevant image parts. The
% coordinates of the selected polygon as ROI is saved as "fileName_face.mat". 
% When running this script for the second time on the same video, it will
% automatically load the selected ROI, unless reselectROI = 1. Default: 0;
% 
% reselectSkinColor: 0 or 1. Only the pixels with skin color will be analyzed. 
% Skin color will be determined by selecting a hue (color) and saturation
% range. The selected range is saved as "fileName_skincolor.mat". When
% running this script for the second time on the same video, it will
% automatically load the selected range, unless reselectSkinColor = 1. Default: 0;
% 
% hueSatMapResolution: integer. Dimensions in pixels of the hue and saturation map 
% that is used to select skin color tone in the video. Default = 500;
% 
% moveDetReso: uneven integer. If higher than 0, the amount of 
% movement across each frame sequence is calculated. The number indicates
% the size of each analyzed image part. For example, in a frame of
% 1100x1100 pixels and with moveDetReso = 11, movement is
% calculated within 100x100 image regions of 11x11 pixels.
%
% outputValuesAndPlots: 0 or 1. If 1, show power and HR values. Default: 0;
%
% ---------------- OUTPUT -------------------
%
% icaMaxPower: power of the largest peak in the frequency power spectrum
% after an independent component analysis. This number represents the 
% pulsation strength at the detected heart rate per component
% 
% icaBestHR: detected pulse (heart) rate per component
%
% rgbMaxPower: power of the largest peak in the frequency power spectrum
% before an independent component analysis. This number represents the 
% pulsation strength at the detected heart rate per color channel [R G B]
%
% rgbBestHR: detected pulse (heart) rate per color channel [R G B]
%
% power_ica: power spectrum after ICA
%
% phase_ica: phase spectrum after ICA
%
% power_rgb: power spectrum before ICA 
%
% phase_rgb: phase spectrum before ICA 
%
% meanPVperFrame: original raw pixel values as a function of video time per
% color channel [R G B]
%
% meanPVperFrame_filt: low-pass frequency filtered pixel values
% (i.e., same as meanPVperFrame but filtered as specified with variable lpFiltParams)
%
% avePixelValue: average pixel value per color channel [R G B]
%
% aveMotion: average motion across all frames and segmented image regions
% per color channel [R G B]
%
%
% 
% 
% --------LICENSE & ACKNOWLEDGMENT-----------
% 
% Copyright © 2017 Marnix Naber, The Netherlands
% 
% This program is distributed under the terms of the GNU Genreal Public
% % License (see gpl.txt)
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
% van der Kooij & Naber (2017). Standardized procedures for the testing and 
% reporting of remote heart rate imaging. Journal title, etc.
%
%
%
%
% -----------------CONTACT-------------------
% 
% For questions, remarks, or problems with the code, 
% please contact: marnixnaber@gmail.com
%
% This script was tested in Matlab version 2014b. In case this script
% does not run because of an error reporting a missing function, then
% please check your matlab version and installed toolboxes. To run this
% script succesfully, the image processing toolbox should be installed.


%% input variables

switch nargin,
    case 0
        fileName                = 'rPPG_video.mp4';
        TOI                     = 'all';
        FOI                     = [45 165]; 
        lpFiltParams            = [6 0.04];
        lpFreqFiltParams        = [3 0.2];
        hueSatMapResolution     = 500;
        reselectROI             = 0; 
        reselectSkinColor       = 0;
        moveDetReso             = 0;
        outputValuesAndPlots    = 1;
    case 1
        TOI                     = 'all';
        FOI                     = [45 165]; 
        lpFiltParams            = [6 0.04];
        lpFreqFiltParams        = [3 0.2];
        hueSatMapResolution     = 500;
        reselectROI             = 0; 
        reselectSkinColor       = 0;
        moveDetReso             = 0;
        outputValuesAndPlots    = 0;
    case 2
        FOI                     = [45 165]; 
        lpFiltParams            = [6 0.04];
        lpFreqFiltParams        = [3 0.2];
        hueSatMapResolution     = 500;
        reselectROI             = 0; 
        reselectSkinColor       = 0;
        moveDetReso             = 0;
        outputValuesAndPlots    = 0;
    case 3
        lpFiltParams            = [6 0.04];
        lpFreqFiltParams        = [3 0.2];
        hueSatMapResolution     = 500;
        reselectROI             = 0; 
        reselectSkinColor       = 0;
        moveDetReso             = 0;
        outputValuesAndPlots    = 0;
    case 4
        lpFreqFiltParams        = [3 0.2];
        hueSatMapResolution     = 500;
        reselectROI             = 0; 
        reselectSkinColor       = 0;
        moveDetReso             = 0;
        outputValuesAndPlots    = 0;
    case 5
        hueSatMapResolution     = 500;
        reselectROI             = 0; 
        reselectSkinColor       = 0;
        moveDetReso             = 0;
        outputValuesAndPlots    = 0;
    case 6
        reselectROI             = 0; 
        reselectSkinColor       = 0;
        moveDetReso             = 0;
        outputValuesAndPlots    = 0;
    case 7
        reselectSkinColor       = 0;
        moveDetReso             = 0;
        outputValuesAndPlots    = 0;
    case 8
        moveDetReso             = 0;
        outputValuesAndPlots    = 0;
    case 9 
        outputValuesAndPlots    = 0;
end


%% start of rPPG analysis


disp('Loading movie: ');
disp(fileName);

vidFile                 = VideoReader(fileName);

nFrames                 = vidFile.NumberOfFrames;
vidDuration             = vidFile.Duration;
frameRate               = vidFile.NumberOfFrames/vidDuration;

disp(['Number of frames: ' num2str(nFrames)]);
disp(['Frame rate: ' num2str(frameRate)]);
disp(['Video duration: ' num2str(vidDuration)]);
disp(['-------------------------------------------------------']);
disp(['-------------------------------------------------------']);
disp(['-------------------------------------------------------']);

if ischar(TOI) & strcmp(lower(TOI),'all') | ischar(TOI) & length(TOI) == 0
    TOI = [0 vidDuration];
    frameNums_sel  = 1:nFrames;
elseif ischar(TOI)
   ERROR('Unknown argument for variable time of interest (TOI)') 
else
    frameNums_sel   = floor(TOI(1)*frameRate)+1:floor(TOI(2)*frameRate);
end

nFrames_sel     = length(frameNums_sel);

if frameNums_sel(end) > nFrames % requested to analyze more video time than available
    disp('time of interest (TOI) indicates a time range exceeds the available video time');
    disp([num2str(nFrames_sel-find(frameNums_sel==nFrames,1,'first')) ' more frames required.']);
    frameNums_sel   = floor(TOI(1)*frameRate)+1:nFrames;
    nFrames_sel     = length(frameNums_sel);
end
    
%% Select or automatically load regions of interest in video image
% Select face (remove head rest) ( or automatically load if alrea dy
% selected)
if ~exist([fileName '_face.mat']) | reselectROI
    temp.cdata = read(vidFile, frameNums_sel(1));
    ha = figure(1);
    imshow(temp.cdata)

    title('Select the region of interest by selecting multiple locations with the left mouse button');
    
    drawnow;
    
    % select face with mouse
    [selROIBW, xiFace, yiFace] = roipoly(uint8(temp.cdata));
    
    save([fileName '_face.mat'],'selROIBW');
else
    load([fileName '_face.mat']);
end

[x,y] = find(selROIBW);
selROIBW = selROIBW(min(x):max(x),min(y):max(y));

if ~exist([fileName '_skincolor.mat']) | reselectSkinColor
    
    temp.cdata = read(vidFile, frameNums_sel(1));
    targetIm = temp.cdata(min(x):max(x),min(y):max(y),:);
    
    for k = 1:3
        tempTargetIm = targetIm(:,:,k);
        tempTargetIm(~selROIBW) = NaN;
        targetIm(:,:,k)=tempTargetIm;
    end
    
    % Select HSV skin color
    hsvChannels = rgb2hsv(targetIm);
    huemap = hsvChannels(:,:,1);
    satmap = hsvChannels(:,:,2);
    intmap = hsvChannels(:,:,3);

    figure();
    subplot(2,2,1)
    imshow(targetIm);
    title(fileName);
    
    subplot(2,2,2)
    imagesc(huemap);
    colormap('jet')
    colorbar
    axis off
    
    subplot(2,2,3);
    hist(huemap(:),linspace(0,1,256));
    title('Hue Histogram')

    subplot(2,2,4);
    hist(satmap(:),linspace(0,1,256));
    title('Saturation Histogram')

    R = satmap(:);
    TH = huemap(:)*2*pi;

    % make polar plot of colors in image
    [X,Y] = pol2cart(TH,R);
    
    imSize = hueSatMapResolution;
    tempHueIm = zeros(imSize,imSize);
    tempSatIm = zeros(imSize,imSize);
    tempIntIm = zeros(imSize,imSize);
    
    Xn = round(X*(imSize/2-1)+imSize/2+1);
    Yn = round(Y*(imSize/2-1)+imSize/2+1);
    zeroCoor = imSize/2;
    
    % fill maps
    tempHueIm((Xn-1).*imSize+Yn) = huemap(:);
    tempSatIm((Xn-1).*imSize+Yn) = satmap(:);
    tempIntIm((Xn-1).*imSize+Yn) = intmap(:);
    
    colorspaceIm = reshape([tempHueIm tempSatIm tempIntIm],size(tempIntIm,1),size(tempIntIm,1),3);
    colorspaceIm = hsv2rgb(colorspaceIm);
    
    dontstop = 1;
    while dontstop
        
        figure(999);
        imshow(colorspaceIm);
        hold on
        plot(zeroCoor,zeroCoor,'go');
        title({'Select the wedge range of skin colors by clicking '; 'the left mouse button at top left and bottom right corner'});
        drawnow;

        pause(1);
        
        xi = [];
        yi = [];
        for i = 1:2
            [xi(i), yi(i)] = ginput(1);

            % draw selection point
            plot(xi(i),yi(i),'gx');

            % draw red line showing the edges of the wedges
            plot(linspace(zeroCoor,imSize,2),polyval(polyfit([zeroCoor xi(i)],[zeroCoor yi(i)],1),linspace(zeroCoor,imSize,2)),'r') 

            drawnow;
        end
        [TH,R] = cart2pol(xi-zeroCoor,yi-zeroCoor);

        for i = 1:2

            % draw green line
            plot(linspace(xi(1),xi(2),2),polyval(polyfit([zeroCoor xi(i)],[zeroCoor yi(i)],1),linspace(xi(1),xi(2),2)),'g') 
            [xic,yic] = pol2cart(linspace(TH(1),TH(2),10),zeros(1,10)+R(i));

            plot(xic+zeroCoor,yic+zeroCoor,'g') 
        end
        drawnow;
        hold off

        [~,R] = cart2pol((xi-zeroCoor-1)/(imSize/2-1),(yi-zeroCoor-1)/(imSize/2-1));

        % recompute range of skin color hue and saturation values
        TH(TH<=0) = TH(TH<=0)+2*pi;
        TH = TH./(2*pi);

        disp(['Range of selected hue values: ' num2str(TH)]);
        disp(['Range of selected saturation values: ' num2str(R)]);

        if R(1) >= R(2)
            error('Radius of second point should be larger than first point!');
        end
        
        selectSkin = huemap>TH(1) & satmap>R(1) & satmap<R(2) | huemap<TH(2) & satmap>R(1) & satmap<R(2);

        selectIm = targetIm;
        selectIm(reshape(repmat(~selectSkin,1,3),size(targetIm,1),size(targetIm,2),3)) = NaN;
        
        h = figure(998);
        imshow(selectIm)
        hold on
        text(10,10,['Selected hue range: ' num2str(TH)],'Color','g')
        text(10,40,['Selected saturation range: ' num2str(R)],'Color','g')
        title('Press "A" for accept, press "R" to reject and reselect skin colors');
        drawnow;
        hold off
        
        keyspressed = [];

        k=0;
        while ~k
            k = waitforbuttonpress;
            currkey = get(gcf,'currentcharacter');
            if strcmp(currkey,'a') | strcmp(currkey,'r')
                k = 1;
            else
                k = 0;
            end
        end
        
        if strcmp(currkey,'a')
            dontstop = 0;
        else
            close(h);
        end
    end
    
    save([fileName '_skincolor.mat'],'TH','R');
else
    load([fileName '_skincolor.mat']);
end


%%
if moveDetReso > 0
    hbm = vision.BlockMatcher( ...
        'ReferenceFrameSource', 'Input port', 'BlockSize', [moveDetReso moveDetReso]);
    hbm.OutputValue = ...
        'Horizontal and vertical components in complex form';
end

meanPVperFrame = NaN(nFrames_sel,3);
aveMotionPerFrame = NaN(nFrames_sel,3);
for f = 1:round(nFrames_sel) % loop through each image frame to calculate average pixel value in region of interest
    if sum(f == round(linspace(1,round(nFrames_sel),11)))
        disp(['Extracting pixel values in video: ' num2str(round((f/nFrames_sel)*100)) '%'])
    end
    
    temp.cdata      = read(vidFile, frameNums_sel(f));      % read frame
%     targetIm        = temp.cdata(rect(2):rect(2)+rect(4)-1,rect(1):rect(1)+rect(3)-1,:);

    targetIm = double(temp.cdata(min(x):max(x),min(y):max(y),:));
    for k = 1:3
        tempTargetIm = targetIm(:,:,k);
        tempTargetIm(~selROIBW) = NaN;
        targetIm(:,:,k)= tempTargetIm;
    end
    
    hsvChannels     = rgb2hsv(uint8(targetIm));
    huemap          = hsvChannels(:,:,1);
    satmap          = hsvChannels(:,:,2);
    
    selectSkin = huemap>TH(1) & satmap>R(1) & satmap<R(2) | huemap<TH(2) & satmap>R(1) & satmap<R(2);
    
    targetIm = double(targetIm);
    targetIm(reshape(repmat(~selectSkin,1,3),size(targetIm,1),size(targetIm,2),3)) = NaN;
    
    if f == 1 & outputValuesAndPlots
        h = figure();
        imagesc(uint8(targetIm));
        saveas(h,[fileName '_ExampleSkinColorSelection.jpg'])
    end
    
    % Compute motion for the two images
    if f > 1 & moveDetReso > 0
        for c = 1:3
            motion = step(hbm, targetIm(:,:,c), prevImage(:,:,c));  
            aveMotionPerFrame(f,c) = mean(sqrt(real(motion(isfinite(motion(:)))).^2 + imag(motion(isfinite(motion(:)))).^2));
        end
    end
    
    for c = 1:3 % loop through colors (RGB)
        selectColorChannel      = targetIm(:,:,c); %% select region of interest from image
        meanPVperFrame(f,c)      = mean(selectColorChannel(isfinite(selectColorChannel(:)))); % store mean value of all pixels in region of interest
    end
    
    prevImage = targetIm;
end

colorNames = {'Red channel','Green channel','Blue channel'};

%  [X Y] = meshgrid(1:motionResolution:size(targetIm, 2), 1:motionResolution:size(targetIm, 1));
% imshow(uint8(targetIm)); hold on;
% quiver(X(:), Y(:), real(motion(:)), imag(motion(:)), 0); hold off;
%    
%% average pixel value per color channel

avePixelValue = [];
for c = 1:3 % loop through colors (RGB)
    avePixelValue(c) = mean(meanPVperFrame(isfinite(meanPVperFrame(1,c)),c));
end

if outputValuesAndPlots
    figure();
    disp('Average pixel value per channel [R G B]: ')
    disp(avePixelValue);
    plot(avePixelValue);
    ylabel('Pixel Value')
    set(gca,'xtick',1:3)
    set(gca,'xticklabel',colorNames)
end

%% average motion per color channel

aveMotion = [];
for c = 1:3 % loop through colors (RGB)
    aveMotion(c) = mean(aveMotionPerFrame(isfinite(aveMotionPerFrame(1,c)),c));
end
if outputValuesAndPlots
    figure();
    disp('Average motion per channel [R G B]: ')
    disp(aveMotion);
    plot(aveMotion);
    ylabel('Motion')
    set(gca,'xtick',1:3)
    set(gca,'xticklabel',colorNames)
end

%% 

[bfilt,afilt]           = butter(lpFiltParams(1),lpFiltParams(2));

if outputValuesAndPlots
    figure();
    for c = 1:3
        subplot(2,2,c);
        plot(meanPVperFrame(:,c))
        hold on
        plot(filtfilt(bfilt,afilt,meanPVperFrame(:,c)),'r')
        title(colorNames{c})
        ylabel('Pixel values')
    end
end

meanPVperFrame_filt = NaN(nFrames_sel,3);

for c = 1:3
    meanPVperFrame_filt(:,c)   = meanPVperFrame(:,c)-filtfilt(bfilt,afilt,meanPVperFrame(:,c));
end

if outputValuesAndPlots
    figure();
    for c = 1:3
        subplot(2,2,c);
        plot(meanPVperFrame_filt(:,c))
        hold on
        title(colorNames{c})
        ylabel('Filtered pixel values')
    end
end

recordHz        = length(meanPVperFrame_filt(:,2))/diff(TOI);

%                 Fs      = 30;%allFrameRates(ppNum,runNum,bodyCondNum); % DIT MOET IK OPNIEUW LADEN
Fs      = recordHz;
L       = length(meanPVperFrame_filt);
NFFT    = 2^nextpow2(L); % Next power of 2 from length of y

freq    = Fs/NFFT*(0:NFFT-1);

freqInterestRange = FOI/60;
fRange2 = find(freq>freqInterestRange(1) & freq<freqInterestRange(2));

power_rgb = [];
phase_rgb = [];

for c = 1:3
    Y       = fft(meanPVperFrame_filt(:,c),NFFT); % calculate frequency spectrum

    power_rgb(:,c)  = Y.*conj(Y)/NFFT;
    phase_rgb(:,c)  = angle(Y(1:NFFT));
end

if outputValuesAndPlots
    figure();
    
    for c = 1:3
        subplot(2,2,c);
        plot(freq(fRange2)*60,power_rgb(fRange2,c))
    %     [pks,locs,w,p] = findpeaks(power_rgb(:,c));
        hold on

        title(colorNames{c})
        ylabel('Power')
    end
end

[rgbMaxPower,rgbBestHR] = max(power_rgb(fRange2,:));
rgbBestHR = freq(fRange2(rgbBestHR))*60;

if outputValuesAndPlots
    disp('Power per color channel [R G B]: ')
    disp(rgbMaxPower)
    disp('HR per color channel [R G B]: ')
    disp(rgbBestHR)
end

comp                = fastica(meanPVperFrame_filt','numOfIC',3,'maxNumIterations',2000,'stabilization','on','verbose','off');

componentNames = {'Component 1','Component 2','Component 3'};
if outputValuesAndPlots
    figure();
    for c = 1:3
        subplot(2,2,c);
        plot(comp(c,:))
        hold on
        title(componentNames{c})
        ylabel('Pixel value')

    end
end

power_ica           = [];
power_ica_lpFilt    = [];
phase_ica           = [];

[dfilt,cfilt]       = butter(lpFreqFiltParams(1),lpFreqFiltParams(2));

for c = 1:3
    Y       = fft(comp(c,:),NFFT); % calculate frequency spectrum
    power_ica(:,c)          = Y.*conj(Y)/NFFT;
    phase_ica(:,c)          = angle(Y(1:NFFT));
    
    power_ica_lpFilt(:,c)   = filtfilt(dfilt,cfilt,power_ica(:,c));
end

if outputValuesAndPlots
    h = figure();
    for c = 1:3
        subplot(2,2,c);
        plot(freq(fRange2)*60,power_ica(fRange2,c),'k')
        hold on
        plot(freq(fRange2)*60,power_ica_lpFilt(fRange2,c),'r')
    %     [pks,locs,w,p] = findpeaks(power_ica(:,c));
        hold on
        title(componentNames{c})
        ylabel('Power')
    end
    saveas(h,[fileName '_PowerSpectrum.jpg'])
end

[icaMaxPower,icaBestHR] = max(power_ica(fRange2,:));
icaBestHR = freq(fRange2(icaBestHR))*60;

if outputValuesAndPlots
disp('Power per ICA component: ')
disp(icaMaxPower)
disp('HR per ICA component: ')
disp(icaBestHR)
end



[icaMaxPower_lpFilt,icaBestHR_lpFilt] = max(power_ica_lpFilt(fRange2,:));
icaBestHR_lpFilt = freq(fRange2(icaBestHR_lpFilt))*60;

if outputValuesAndPlots
disp('Low pass filtered power per ICA component: ')
disp(icaMaxPower_lpFilt)
disp('HR after low pass filtering of power spectrum, per ICA component: ')
disp(icaBestHR_lpFilt)
end

%%

clear vidFile;
