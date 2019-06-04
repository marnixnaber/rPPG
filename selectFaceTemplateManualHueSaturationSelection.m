function [skinDetection,faceIm,skinImB_3] = selectFaceTemplateManualHueSaturationSelection(skinDetection, Im) 

hueSatMapResolution     = 500;

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
    else
        faceIm      = faceIm(1:size(faceGridVal,1),1:size(faceGridVal,2));
    end
    faceGridVal = faceGridVal(1:size(faceIm,1),1:size(faceIm,2));
else
    faceGridVal = faceGridVal(1:size(faceIm,1),1:size(faceIm,2));
end
faceIm = double(faceIm);


if skinDetection.firstTimeDetected==1 & skinDetection.nTimesDetected ==1 ;  % only first time when face is detected

    % Select HSV skin color
    hsvChannels = rgb2hsv(faceIm);
    huemap = hsvChannels(:,:,1);
    satmap = hsvChannels(:,:,2);
    intmap = hsvChannels(:,:,3);

    figure();
    subplot(2,2,1)
    imshow(uint8(faceIm));
    
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
        imshow(uint8(colorspaceIm));
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
        
%         skinDetection.skinImB = huemap>TH(1) & satmap>R(1) & satmap<R(2) | huemap<TH(2) & satmap>R(1) & satmap<R(2);

        if TH(1) > TH(2)
            skinDetection.skinImB = huemap>TH(1) & huemap<=1 & satmap>R(1) & satmap<R(2) | huemap<TH(2) & huemap>=0 & satmap>R(1) & satmap<R(2);
        elseif TH(1) < TH(2)
            skinDetection.skinImB = huemap>TH(1) & huemap<TH(2) & satmap>R(1) & satmap<R(2);
        else
            skinDetection.skinImB = satmap>R(1) & satmap<R(2);
        end

        selectIm = faceIm;
        selectIm(reshape(repmat(~skinDetection.skinImB,1,3),size(faceIm,1),size(faceIm,2),3)) = NaN;
        
        h = figure(998);
        imshow(uint8(selectIm))
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
    
    
    skinDetection.TH = TH;
    skinDetection.R = R;

end

if mod(skinDetection.f,skinDetection.calcFrameRate) & skinDetection.calcFrameRate > 1 | skinDetection.calcFrameRate==1
    
    % Select HSV skin color
    hsvChannels = rgb2hsv(faceIm);
    huemap = hsvChannels(:,:,1);
    satmap = hsvChannels(:,:,2);
    intmap = hsvChannels(:,:,3);

    if skinDetection.TH(1) > skinDetection.TH(2)
        skinDetection.skinImB = huemap>skinDetection.TH(1) & huemap<=1 & satmap>skinDetection.R(1) & satmap<skinDetection.R(2) | huemap<skinDetection.TH(2) & huemap>=0 & satmap>skinDetection.R(1) & satmap<skinDetection.R(2);
    elseif skinDetection.TH(1) < skinDetection.TH(2)
        skinDetection.skinImB = huemap>skinDetection.TH(1) & huemap<skinDetection.TH(2) & satmap>skinDetection.R(1) & satmap<skinDetection.R(2);
    else
        skinDetection.skinImB = satmap>skinDetection.R(1) & satmap<skinDetection.R(2);
    end
end
        
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
    