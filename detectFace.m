function [faceDetection] = detectFace(Im,faceDetection)

FDetect = vision.CascadeObjectDetector('FrontalFaceLBP','MergeThreshold',faceDetection.threshold);
faceDetection.bbox    = step(FDetect,Im);

if size(faceDetection.bbox,1) > 1 % if more than one face detected, take largest
    faceDetection.nFacesDetected = size(faceDetection.bbox,1);
    [~,maxFFIdx] = max(faceDetection.bbox(:,3).*faceDetection.bbox(:,4));
    faceDetection.nTimesDetected = faceDetection.nTimesDetected+1;
    faceDetection.bbox = faceDetection.bbox(maxFFIdx,:);
elseif size(faceDetection.bbox,1) == 1 % just one face detected
    faceDetection.nFacesDetected = 1;
    maxFFIdx = 1;
    faceDetection.nTimesDetected = faceDetection.nTimesDetected+1;
    faceDetection.bbox = faceDetection.bbox(maxFFIdx,:);
elseif size(faceDetection.bbox,1) == 0 % no face detected
    faceDetection.nFacesDetected = 0;
end