function [skinDetection,faceIm,skinImB_3] = selectFaceTemplateColorCluster(skinDetection, Im) 

xyCoor = round([min(skinDetection.bboxPolygon(1:2:end)) min(skinDetection.bboxPolygon(2:2:end)) max(skinDetection.bboxPolygon(1:2:end)) max(skinDetection.bboxPolygon(2:2:end))]);
xyCoor(xyCoor<1) = 1;

if length(size(Im)) == 3
    faceIm = Im(xyCoor(2):xyCoor(4),xyCoor(1):xyCoor(3),:);
else
    faceIm = Im(xyCoor(2):xyCoor(4),xyCoor(1):xyCoor(3));
end

% reset face sectors
faceGridVal = zeros(skinDetection.nFaceSectors,skinDetection.nFaceSectors);
faceGridVal(:) = 1:skinDetection.nFaceSectors^2;

% Divide face region in face sectors and rotate depending on the angle of the face
boxSize         = round(sqrt((skinDetection.bboxPolygon(4)-skinDetection.bboxPolygon(2))^2+(skinDetection.bboxPolygon(3)-skinDetection.bboxPolygon(1))^2))-1;
faceGridVal     = imrotate(imresize(faceGridVal,[boxSize boxSize],'Method','Nearest'),atand((skinDetection.bboxPolygon(4)-skinDetection.bboxPolygon(2))/(skinDetection.bboxPolygon(1)-skinDetection.bboxPolygon(3))));
if size(faceIm,1) > size(faceGridVal,1) | size(faceIm,2) > size(faceGridVal,2)
    if length(size(Im)) == 3
        faceIm      = faceIm(1:size(faceGridVal,1),1:size(faceGridVal,2),:);
%         faceGridVal = repmat(faceGridVal(1:size(faceIm,1),1:size(faceIm,2)),1,1,3);
    else
        faceIm      = faceIm(1:size(faceGridVal,1),1:size(faceGridVal,2));
%         faceGridVal = faceGridVal(1:size(faceIm,1),1:size(faceIm,2));
    end
    faceGridVal = faceGridVal(1:size(faceIm,1),1:size(faceIm,2));
else
%     if length(size(Im)) == 3
%         faceGridVal = repmat(faceGridVal(1:size(faceIm,1),1:size(faceIm,2)),1,1,3);
%     else
%         faceGridVal = faceGridVal(1:size(faceIm,1),1:size(faceIm,2));
%     end
    faceGridVal = faceGridVal(1:size(faceIm,1),1:size(faceIm,2));
end
faceIm = double(faceIm);



% perform skin detection every X frames (or only first frame)
if skinDetection.calcFrameRate==0 & skinDetection.firstTimeDetected==1 | mod(skinDetection.f,skinDetection.calcFrameRate) & skinDetection.calcFrameRate > 1 | skinDetection.calcFrameRate==1

    % remove grid sectors that are probably not skin
    faceTemplate = zeros(size(faceIm,1),size(faceIm,2));
    faceTemplate = ~(faceGridVal(:)<11 | faceGridVal(:)> 53 | faceGridVal(:)== 14 | faceGridVal(:)== 15 | faceGridVal(:)== 16 | faceGridVal(:)== 49 | faceGridVal(:)== 50);
    faceTemplate = logical(faceTemplate);

    
    % select color cluster
    faceSize = size(faceIm,1);

    % convert to LAB space (luminance, saturation, hue)
    % faceIm_lab = rgb2lab(faceIm);
    faceIm_lab = rgb2lab(faceIm,'WhitePoint','d50');

    faceTemplateResize = imresize(faceTemplate,[size(faceIm_lab,1) size(faceIm_lab,2)]);
    faceIm_lab(~faceTemplateResize) = NaN;

    % only select saturation and hue channels
    ab = faceIm_lab(:,:,2:3);

    % reshape pixels to a single column per channel
    nrows = size(ab,1);
    ncols = size(ab,2);
    ab = reshape(ab,nrows*ncols,2);

    % search for four colors (works best to detect skin)
    warning('off')
    [cluster_idx, cluster_center] = kmeans(ab,skinDetection.nColorClusters,'distance','sqEuclidean','Replicates',1);
    warning('on')
    skinDetection.pixel_labels = reshape(cluster_idx,nrows,ncols);

    % Which of the clustered colors is skin? Use center part of face as reference for skin color
    pixel_labels_CenterFace = skinDetection.pixel_labels(floor(faceSize/2)-floor(faceSize/4):floor(faceSize/2)+floor(faceSize/4),floor(faceSize/2)-floor(faceSize/4):floor(faceSize/2)+floor(faceSize/4));

    % select color with most pixels in center face
    [hx,hy] = hist(pixel_labels_CenterFace(:),unique(pixel_labels_CenterFace(:)));
    [~,selectKIdx] = max(hx);
    skinDetection.selectK = hy(selectKIdx);

    % figure(2);,imagesc(pixel_labels_CenterFace)
    % figure(3);,imagesc(skinDetection.pixel_labels)
    % figure();,
    
    skinDetection.skinImB = zeros(size(skinDetection.pixel_labels));
    for k = skinDetection.selectK
        skinDetection.skinImB = skinDetection.skinImB+(skinDetection.pixel_labels==k);
    end
    skinDetection.skinImB     = logical(skinDetection.skinImB);

end
% skinDetection.f


if size(skinDetection.skinImB,1)~=size(faceGridVal,1) | size(skinDetection.skinImB,2)~=size(faceGridVal,2)
    skinDetection.skinImB = imresize(skinDetection.skinImB,[size(faceGridVal,1) size(faceGridVal,2)]);
end

for gridIdx = 1:skinDetection.nFaceSectors^2 % check which grid sectors have enough skin detected
    if sum(faceGridVal(:)==gridIdx & skinDetection.skinImB(:))/sum(faceGridVal(:)==gridIdx) > skinDetection.fracPixelsPresent;
        skinDetection.skinImB(faceGridVal(:)==gridIdx) = 1;
    else
        skinDetection.skinImB(faceGridVal(:)==gridIdx) = 0;
    end
end

skinImB_3   = repmat(skinDetection.skinImB,1,1,3);
    