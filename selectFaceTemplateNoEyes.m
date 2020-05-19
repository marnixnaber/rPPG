function [faceIm,skinImB_3] = selectFaceTemplateNoEyes(skinDetection, Im) 


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

% Divide face region in face sectors and rotate depending on the angle of
% the face
boxSize         = round(sqrt((skinDetection.bboxPolygon(4)-skinDetection.bboxPolygon(2))^2+(skinDetection.bboxPolygon(3)-skinDetection.bboxPolygon(1))^2))-1;
faceGridVal     = imrotate(imresize(faceGridVal,[boxSize boxSize],'Method','Nearest'),atand((skinDetection.bboxPolygon(4)-skinDetection.bboxPolygon(2))/(skinDetection.bboxPolygon(1)-skinDetection.bboxPolygon(3))));
if size(faceIm,1) > size(faceGridVal,1) | size(faceIm,2) > size(faceGridVal,2)
    if length(size(Im)) == 3
        faceIm      = faceIm(1:size(faceGridVal,1),1:size(faceGridVal,2),:);
        faceGridVal = repmat(faceGridVal(1:size(faceIm,1),1:size(faceIm,2)),1,1,3);
    else
        faceIm      = faceIm(1:size(faceGridVal,1),1:size(faceGridVal,2));
        faceGridVal = faceGridVal(1:size(faceIm,1),1:size(faceIm,2));
    end
else
    if length(size(Im)) == 3
        faceGridVal = repmat(faceGridVal(1:size(faceIm,1),1:size(faceIm,2)),1,1,3);
    else
        faceGridVal = faceGridVal(1:size(faceIm,1),1:size(faceIm,2));
    end
end
faceIm = double(faceIm);

% 10 percent of each side horizontally
% 20 percent of each side vertically
shrinkFactorH = ceil(skinDetection.nFaceSectors*0.1);
shrinkFactorV = ceil(skinDetection.nFaceSectors*0.2);

gridDeselIdx = 0:skinDetection.nFaceSectors*(shrinkFactorH); % left part
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors^2-skinDetection.nFaceSectors*(shrinkFactorH)+1:skinDetection.nFaceSectors^2]; % right part
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors*(shrinkFactorH)+1:skinDetection.nFaceSectors*(shrinkFactorH)+shrinkFactorV]; % top-left
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors*(shrinkFactorH+1)+1]; % top-left corner
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors^2-skinDetection.nFaceSectors*(shrinkFactorH)+1-shrinkFactorV:skinDetection.nFaceSectors^2-skinDetection.nFaceSectors*(shrinkFactorH)]; % bottom-right
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors^2-skinDetection.nFaceSectors*(shrinkFactorH+1)]; % bottom-right corner
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors*(shrinkFactorH+1)-shrinkFactorV+1:skinDetection.nFaceSectors*(shrinkFactorH+1)]; % bottom-left
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors*(shrinkFactorH+2)]; % bottom-left corner
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors^2-skinDetection.nFaceSectors*(shrinkFactorH+1)+1:skinDetection.nFaceSectors^2-skinDetection.nFaceSectors*(shrinkFactorH+1)+shrinkFactorV]; % top-right
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors^2-skinDetection.nFaceSectors*(shrinkFactorH+2)+1]; % top-right corner
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors+ceil(skinDetection.nFaceSectors*0.4):skinDetection.nFaceSectors:skinDetection.nFaceSectors^2]; % eyes
gridDeselIdx = [gridDeselIdx skinDetection.nFaceSectors+ceil(skinDetection.nFaceSectors*0.5):skinDetection.nFaceSectors:skinDetection.nFaceSectors^2]; % eyes



% remove grid sectors that are not skin
if length(size(Im)) == 3
    skinImB_3 = zeros(size(faceIm,1),size(faceIm,2),3);
    
    for k = gridDeselIdx
        skinImB_3(:) = skinImB_3(:)+(faceGridVal(:)==k);
    end
    skinImB_3 = ~skinImB_3;
else
    skinImB_3 = zeros(size(faceIm,1),size(faceIm,2));
    
    for k = gridDeselIdx
        skinImB_3 = skinImB_3+(faceGridVal==k);
    end
    skinImB_3 = ~skinImB_3;
end

skinImB_3 = logical(skinImB_3);