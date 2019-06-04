function [faceIm,skinImB_3] = selectFaceTemplate(skinDetection, Im) 


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

% remove grid sectors that are not skin
if length(size(Im)) == 3
    skinImB_3 = zeros(size(faceIm,1),size(faceIm,2),3);
    skinImB_3(:) = ~(faceGridVal(:)<11 | faceGridVal(:)> 53 | faceGridVal(:)== 14 | faceGridVal(:)== 15 | faceGridVal(:)== 16 | faceGridVal(:)== 49 | faceGridVal(:)== 50);
else
    skinImB_3 = zeros(size(faceIm,1),size(faceIm,2));
    skinImB_3 = ~(faceGridVal(:)<11 | faceGridVal(:)> 53 | faceGridVal(:)== 14 | faceGridVal(:)== 15 | faceGridVal(:)== 16 | faceGridVal(:)== 49 | faceGridVal(:)== 50);
end

skinImB_3 = logical(skinImB_3);