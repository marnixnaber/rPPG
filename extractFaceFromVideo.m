function [pixelValPerFrame,faceMap,vidInfo,faceDetection,faceTracking] = extractFaceFromVideo(fileName,TOI, showSel)


%% set defaults
switch nargin,
    case 0
        fileName                = 'rPPG_video.mp4';
        TOI                     = 'all';
        showSel                 = 1;
    case 1
        TOI                     = 'all';
        showSel                 = 1;
    case 2
        showSel                 = 1;
end

%% parameters

faceDetection                   = struct();
faceDetection.active                = 1;        % 0 = no face detection, 1 = if face detection on
faceDetection.threshold             = 20;       % [a.u. cascade] threshold that determines how likely faces are detected (the lower, the more likely a face is detected but also the more likely a false alarm)
faceDetection.thresAdaptSpeed       = 0.9;      % [0 1] speed at which threshold is lowered if no face is detected
faceDetection.minFaceSizePrc        = 5;        % [0-100] face should cover a minimum percentage of area within the video frame
faceDetection.textSize              = 1/50;     % [0-1] proportion text size of warning messages; proportion of the y-scale of video images
faceDetection.noFaceNFramesReset    = 150;      % [integer] wait X of frames of no face detected until face detector threshold is reset 

faceTracking                    = struct;
faceTracking.active                 = 1;        % track the movement of the face rather than detecting it every X frames (X set with variable faceDetectFrameRate)
faceTracking.minPoints              = 10;       % minimal points that need to be detected in face to be able to track face

skinDetection                   = struct;
skinDetection.nFaceSectors          = 12;       % number of row/column sectors that can be selected in the face as skin (8x8 = 64)
skinDetection.method                = 1;        % [0 1 2 3]; 
% 0 = preset template segmentation of face to remove background pixels.
% Disadvantage: not every face has the same shape

% 1 = preset template segmentation of face WITHOUT eyes and WITHOUT background pixels.
% Disadvantage: not every face has the same shape
%
% 2 = automatic selection of skin pixels through k-means color clustering
% Advantage of detecting skin every frame is that skin is detected independent of head angle.
% Disadvantage of detecting skin every frame is that it takes time to calculate.
%
% 3 = manual selection of skin pixels through hue and saturation range selection;
% Disadvantage is that you have to manually set the hue and saturation range.
    skinDetection.nColorClusters    = 4;        % [#]; Number of color clusters to look for
    skinDetection.calcFrameRate     = 0;        % [# frames]; 0 = calculate template only for first frame, 1 = refresh template every frame, 2> every X frames --> the larger the number, the faster the computation of HR
    skinDetection.fracPixelsPresent = 0.8;      % [0-1]; fraction of colored pixels within a face grid sector that needs to be detected as skin color

%% DO NOT CHANGE THE FOLLOWING parameters

faceDetection.noFaceNFrames         = 0;        % number of frames that a face is not detected
faceDetection.nTimesRedetected      = 0;
faceDetection.warningMess           = 0;         % 1 = face too small, get closer; 2 = face partially outside image, move to center
faceDetection.nTimesDetected        = 0;
faceDetection.nFacesDetected        = 0;
faceDetection.startThreshold        = faceDetection.threshold;       % reference


faceTracking.Method                 = 1;        % 1 = MinEigenFeatures
faceTracking.numPts                 = 0;        % number of points detected
faceTracking.points                 = [];       % x and y coordinates of points
faceTracking.nTimesPointsRedetected = 0;

%% load point tracker

% Create the point tracker object.
faceTracking.pointTracker = vision.PointTracker('MaxBidirectionalError', 2);

enoughPointVisibleText = {'No','Yes'};

%% load movie
disp('Loading movie: ');
disp(fileName);

vidFile                 = VideoReader(fileName);

%% inspect video properties

vidInfo = struct();
vidInfo.nFrames     = vidFile.NumberOfFrames;
vidInfo.vidDuration = vidFile.Duration;
vidInfo.frameRate   = (vidFile.NumberOfFrames-1)/vidInfo.vidDuration;

disp(['Number of frames: ' num2str(vidInfo.nFrames)]);
disp(['Frame rate: ' num2str(vidInfo.frameRate)]);
disp(['Video duration: ' num2str(vidInfo.vidDuration)]);
disp(['-------------------------------------------------------']);
disp(['-------------------------------------------------------']);
disp(['-------------------------------------------------------']);

% break
%% check time of interest and select frame numbers

if ischar(TOI) & strcmp(lower(TOI),'all') | ischar(TOI) & length(TOI) == 0
    TOI = [0 vidInfo.vidDuration];
    vidInfo.frameNums_sel  = 1:vidInfo.nFrames;
elseif ischar(TOI)
   ERROR('Unknown argument for variable time of interest (TOI)') 
else
    vidInfo.frameNums_sel   = floor(TOI(1)*vidInfo.frameRate)+1:floor(TOI(2)*vidInfo.frameRate);
end

vidInfo.nFrames_sel     = length(vidInfo.frameNums_sel);

if vidInfo.frameNums_sel(end) > vidInfo.nFrames % requested to analyze more video time than available
    disp('time of interest (TOI) indicates a time range exceeds the available video time');
    disp([num2str(vidInfo.nFrames_sel-find(vidInfo.frameNums_sel==vidInfo.nFrames,1,'first')) ' more frames required.']);
    vidInfo.frameNums_sel   = floor(TOI(1)*vidInfo.frameRate)+1:vidInfo.nFrames;
    vidInfo.nFrames_sel     = length(vidInfo.frameNums_sel);
end

vidInfo.tStamp      = zeros(1,vidInfo.nFrames_sel);
faceMap = [];


%% loop through video frames and extract face region

tic
for f = 1:round(vidInfo.nFrames_sel) % loop through each image frame to calculate average pixel value in region of interest
    if sum(f == round(linspace(1,round(vidInfo.nFrames_sel),11)))
        disp(['Extracting pixel values in video: ' num2str(round((f/vidInfo.nFrames_sel)*100)) '% - ' num2str(toc) 's'])
    end
    
    temp.cdata          = read(vidFile, vidInfo.frameNums_sel(f));      % read frame
    vidInfo.tStamp(f)   = (vidInfo.frameNums_sel(f)-1)*(1/vidInfo.frameRate);
    Im                  = temp.cdata;
    ImGray              = rgb2gray(Im);
    
    % predefine variables
    if f == 1
        vidInfo.resolution = size(Im);
        if length(vidInfo.resolution) == 3 % color channels
            pixelValPerFrame = NaN(round(vidInfo.nFrames_sel),vidInfo.resolution(3));
        else
            pixelValPerFrame = NaN(1,round(vidInfo.nFrames_sel));
        end
        
        warningTextPosition = [vidInfo.resolution(2)*0.01 vidInfo.resolution(1)*0.8];
        warningTextFontSize = ceil(vidInfo.resolution(2)*faceDetection.textSize*2);
        if warningTextFontSize < 8
            warningTextFontSize = 8;
        end
    end
    
%% detect face

    faceDetection.nTimesRedetectedOld = faceDetection.nTimesRedetected;
    if faceDetection.active & faceDetection.nFacesDetected==0 | faceDetection.active & size(faceTracking.points,1) < faceTracking.minPoints % detect if not yet detected or if not enough features available for tracking
        
        faceDetection.nTimesRedetected = faceDetection.nTimesRedetected+1;
        
        %%%%
        faceDetection = detectFace(ImGray,faceDetection);
        %%%%
        
        if faceDetection.nFacesDetected == 0 % not detected, start counting how many frames no face
            faceDetection.noFaceNFrames = faceDetection.noFaceNFrames+1;
        else
            faceDetection.noFaceNFrames = 0;
        end
            
        % low threshold to make face detection more sensitive, but only if
        % face was never detected before. Otherwise wait 5s to reset detector
        if faceDetection.nFacesDetected == 0 & faceDetection.nTimesDetected == 0
            if faceDetection.noFaceNFrames < faceDetection.noFaceNFramesReset
                faceDetection.threshold = round(faceDetection.threshold*faceDetection.thresAdaptSpeed);
            else
                faceDetection.threshold = faceDetection.startThreshold;
                faceDetection.noFaceNFrames = 0;
            end
        elseif faceDetection.nFacesDetected == 0 & faceDetection.nTimesDetected > 0 % face is probably gone for a short moment
            if faceDetection.noFaceNFrames >= faceDetection.noFaceNFramesReset
                faceDetection.nTimesDetected = 0;
                faceDetection.noFaceNFrames = 0;
                faceDetection.threshold = faceDetection.startThreshold;
            end
        end
        
        if ~isempty(faceDetection.bbox)
            x = faceDetection.bbox(1, 1); 
            y = faceDetection.bbox(1, 2); 
            w = faceDetection.bbox(1, 3); 
            h = faceDetection.bbox(1, 4);
            faceDetection.bboxPolygon = [x, y, x+w, y, x+w, y+h, x, y+h];
        else
            faceDetection.bboxPolygon = [];
        end
        
        faceTracking.bbox           = faceDetection.bbox;
        faceTracking.bboxPolygon    = faceDetection.bboxPolygon;
        if ~isempty(faceTracking.bbox)
            faceTracking.bboxPoints     = bbox2points(faceTracking.bbox(1, :));
        end
    end
    
%% track face with unique features within bounding box around face

    if faceDetection.nFacesDetected > 0 & faceTracking.active
        
        if isempty(faceTracking.points)
            
            faceTracking.nTimesPointsRedetected = faceTracking.nTimesPointsRedetected+1;
            faceTracking.points = detectMinEigenFeatures(ImGray, 'ROI', faceTracking.bbox);

            faceTracking.points = faceTracking.points.Location;
            initialize(faceTracking.pointTracker, faceTracking.points, ImGray);
            
            faceTracking.oldPoints = faceTracking.points;
        end
        
        [faceTracking.points, faceTracking.isFound] = step(faceTracking.pointTracker, ImGray);
        faceTracking.visiblePoints   = faceTracking.points(faceTracking.isFound, :);
        faceTracking.oldInliers      = faceTracking.oldPoints(faceTracking.isFound, :);
        
        faceTracking.numPts = size(faceTracking.visiblePoints,1);
        
        if faceTracking.numPts >= 2 % need at least 2 points

            % Estimate the geometric transformation between the old points
            % and the new points and eliminate outliers
            [faceTracking.xform, faceTracking.oldInliers, faceTracking.visiblePoints] = estimateGeometricTransform(...
                faceTracking.oldInliers, faceTracking.visiblePoints, 'similarity', 'MaxDistance', 4);

            % Apply the transformation to the bounding box points
            faceTracking.bboxPoints = transformPointsForward(faceTracking.xform, faceTracking.bboxPoints);

            % Insert a bounding box around the object being tracked
            faceTracking.bboxPolygon = reshape(faceTracking.bboxPoints', 1, []);

            % Reset the points
            faceTracking.oldPoints = faceTracking.visiblePoints;
            setPoints(faceTracking.pointTracker, faceTracking.oldPoints);
        end
        
        faceTracking.validityTime(f) = sum(faceTracking.isFound);
        if sum(faceTracking.isFound) < 2 % not enough visible
            faceTracking.enoughPointsVisible = 0;
        else
            faceTracking.enoughPointsVisible = 1;
        end
        
        if faceTracking.numPts < faceTracking.minPoints
            faceTracking.points = [];
            faceTracking.oldPoints = [];
            release(faceTracking.pointTracker);
        end
        
        if showSel
            
            % Insert the bounding box around the object being tracked
            ImDisp = insertShape(Im, 'FilledPolygon', faceTracking.bboxPolygon,'Opacity',0.4);

            % Display tracked points 
            ImDisp = insertMarker(ImDisp, faceTracking.visiblePoints, '+','Color', 'white');
        end
    else
        faceTracking.enoughPointsVisible = 0;
        
        if showSel
            ImDisp = Im;
        end
    end
    
    if showSel & faceDetection.nFacesDetected > 0 & faceTracking.active
        textString = {
            ['# Faces detected: ' num2str(faceDetection.nFacesDetected)],...
            ['Threshold:' num2str(faceDetection.threshold)],...
            ['Frame #:' num2str(f)],...
            ['Frame Rate:' sprintf('%2.2f',vidInfo.frameRate) ' per second'],...
            ['Time:' sprintf('%4.2f',vidInfo.tStamp(f)) 's'],...
            ['Enough points visible:' enoughPointVisibleText{faceTracking.enoughPointsVisible+1}]};

        textPosition = zeros(size(textString,2),2);
        for i = 1:size(textString,2)
            textPosition(i,:) = [vidInfo.resolution(1)*faceDetection.textSize vidInfo.resolution(2)*faceDetection.textSize+2*vidInfo.resolution(2)*faceDetection.textSize*(i-1)];
        end
        textFontSize = ceil(vidInfo.resolution(2)*faceDetection.textSize);
        if textFontSize < 8
            textFontSize = 8;
        end
        ImDisp = insertText(ImDisp,textPosition,textString,'FontSize',textFontSize,'TextColor','r','BoxOpacity',0.4);
    end
    
%% check for size and location bounding box around face
    
    if faceDetection.nFacesDetected > 0 & faceTracking.active
        faceDetection.faceArea = max(faceDetection.bbox(:,3).*faceDetection.bbox(:,4));
        boxInsideImage  = faceTracking.bboxPolygon(1) > 0 & faceTracking.bboxPolygon(7) > 0 & faceTracking.bboxPolygon(2) > 0 & faceTracking.bboxPolygon(4) > 0 & faceTracking.bboxPolygon(3) <= size(ImGray,2) & faceTracking.bboxPolygon(5) <= size(ImGray,2) & faceTracking.bboxPolygon(6) <= size(ImGray,1) & faceTracking.bboxPolygon(8) <= size(ImGray,1);
        boxBigEnough    = faceDetection.faceArea >= prod(size(ImGray))*0.01*faceDetection.minFaceSizePrc;    
        if ~boxBigEnough
            faceDetection.warningMess = 1; % face too small, get closer to camera
            warningText = 'Face is too small. Move closer to the camera';
        elseif ~boxInsideImage
            faceDetection.warningMess = 2; % face outside frame, move to center view of camera
            warningText = 'Face is not visible. Center faces in the video';
        else
            faceDetection.warningMess = 0; % everything ok
        end
    else
        faceDetection.warningMess = 3; % no face detected
        warningText = 'Face is not detected.';
    end
    
    % show warning in image
    if faceDetection.warningMess > 0 & showSel
        ImDisp = insertText(ImDisp,warningTextPosition,{warningText},'FontSize',warningTextFontSize,'TextColor','r','BoxOpacity',0.2);
    end
    
%% select face based on template

    if faceDetection.nFacesDetected > 0 & boxInsideImage & boxBigEnough % detect skin every X frames, only if face is actually detected
        
        skinDetection.bboxPolygon = faceTracking.bboxPolygon;
        skinDetection.firstTimeDetected = faceDetection.nTimesRedetectedOld~=faceDetection.nTimesRedetected;
        skinDetection.nTimesDetected = faceDetection.nTimesDetected;
        skinDetection.f = f;
        
        if skinDetection.method == 0 % predefined template
            [faceIm,skinImB_3] = selectFaceTemplate(skinDetection, Im);
        elseif skinDetection.method == 1 % predefined template
            [faceIm,skinImB_3] = selectFaceTemplateNoEyes(skinDetection, Im);
        elseif skinDetection.method == 2 % template based on k-clustering of colors
            [skinDetection,faceIm,skinImB_3] = selectFaceTemplateColorCluster(skinDetection, Im);
        elseif skinDetection.method == 3 % template based on manual setting of hue and saturation ranges
            [skinDetection,faceIm,skinImB_3] = selectFaceTemplateManualHueSaturationSelection(skinDetection, Im);
        end
        
        faceIm(~skinImB_3) = 0;
  
%% attach cropped face image to original video image

        if showSel
            ImDisp(1:size(faceIm,1),end+1:end+size(faceIm,2),:) = faceIm;
        end
    
%% average pixel values per color channel
        
        if length(vidInfo.resolution) == 3 % color channels
            for c = 1:3
                faceImPerChannel = faceIm(:,:,c);
                pixelValPerFrame(f,c) = mean(faceImPerChannel(skinImB_3(:,:,c)));
            end
%             pixelValPerFrame(f,:)
        else
            pixelValPerFrame(f) = mean(faceIm(skinImB_3(:)));
        end 
        
    else
        if length(vidInfo.resolution) == 3 % color channels
            pixelValPerFrame(f,:) = [NaN NaN NaN];
        else
            pixelValPerFrame(f) = NaN;
        end
    end
    

%% show video frame 
    
    if showSel & faceDetection.nFacesDetected > 0 & faceTracking.active
        figure(1);
            
        imshow(ImDisp)
        drawnow;
        clear ImDisp
    end
    
end
toc

release(faceTracking.pointTracker);