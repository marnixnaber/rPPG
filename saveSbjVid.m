if ~(exist('plots'))
    mkdir('plots');
end

cd 'plots'

if ~(exist(sbjName, 'dir'))
    mkdir(sbjName);
end

cd(sbjName)
%         cd(['plots' filesep sbjName])

videoName = erase(files(thisVideo).name, '.mp4');

video = VideoWriter([videoName '.avi']); %create the video object
open(video); %open the file for writing


%% Make all frames same size

for thisFrame = 1:length(imToVid)
    
    allFrameWidth(thisFrame,1) = size(imToVid(thisFrame).thisImage, 2);
    
end

refFrameWidth = min(allFrameWidth);


for thisFrame = 1:length(imToVid)
    
    if size(imToVid(thisFrame).thisImage, 2) >= refFrameWidth
        
        widthDifference = size(imToVid(thisFrame).thisImage, 2) - refFrameWidth;
        
        rowToStartTrim = size(imToVid(thisFrame).thisImage, 2) - widthDifference-1;
        
        imToVid(thisFrame).thisImage(:, rowToStartTrim:end, :) = [];
        
    end
    
end

%%

for thisFrame = 1:length(imToVid)
    
    writeVideo(video, imToVid(thisFrame).thisImage);
    
end

close(video);


cd('../../')