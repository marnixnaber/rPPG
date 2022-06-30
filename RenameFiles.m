clear; close all; clc;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename)); clear tmp
%%
% Let's find all the subfolders containing videos
% Get a list of all files and folders in this folder.
cd 'C:\Gluce\onlineSchandry\extractFaceFromVideo\cohface\Videos';
videosFolder = 'C:\Gluce\onlineSchandry\extractFaceFromVideo\cohface\Videos';
files = dir(pwd);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
videosSubFolders = files(dirFlags);
% Print folder names to command window.
for k = 1 : length(videosSubFolders)
    fprintf('Sub folder #%d = %s\n', k, videosSubFolders(k).name);
end

%% Let's reorder them
[~, reindex] = sort( str2double( regexp( {videosSubFolders.name}, '\d+', 'match', 'once' )))
videosSubFolders = videosSubFolders(reindex);

for thisSbj = 1:length(videosSubFolders)
    
    cd([videosFolder filesep videosSubFolders(thisSbj).name]);
    
    disp(['This is Sbj ' videosSubFolders(thisSbj).name])
    

    %% From here on, we are in the folder with the 4 videos
    videos = dir(pwd);
    % Get a logical vector that tells which is a directory.
    dirFlags = [videos.isdir];
    % Extract only those that are directories.
    videosSubRepetitions = videos(dirFlags);
    % Print folder names to command window.
    for k = 1 : length(dirFlags)
        fprintf('Sub folder #%d = %s\n', k, videosSubRepetitions(k).name);
    end
    
    videosSubRepetitions(1:2) = [];
    
    for thisVideo = 1:length(videosSubRepetitions)
        
        disp(['This is Video ' videosSubRepetitions(thisVideo).name])
        cd(videosSubRepetitions(thisVideo).name)
        
        sbjFiles = dir(pwd);
        sbjFiles(1:2) = [];
        
        for thisFile = 1:length(sbjFiles)
            
            [filepath, name, ext] = fileparts(sbjFiles(thisFile).name);
            disp(['This is File ' num2str(thisFile) '_' ext])
            
            
            oldName = sbjFiles(thisFile).name;
            
            thisVideoName = [num2str(thisSbj) '_' num2str(thisVideo) ext];
            
            movefile(fullfile(pwd, sbjFiles(thisFile).name), ...
                fullfile('C:\Gluce\onlineSchandry\extractFaceFromVideo\cohface\renamedVideos', ...
                thisVideoName));
            
        end
        
        cd ..
        
        
    end

end










