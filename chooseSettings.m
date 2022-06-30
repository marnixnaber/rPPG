

% Choose if you want to set new parameters or load previous ones
setParameters =  menu('Do you want to set basic parameters?', 'Set Parameters', ...
    'Load Parameters', 'Use Default');
switch setParameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% User has chosen to select parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        % Choose if you want to make plots
        plotYesNo = questdlg('Do you want to make plots?', 'PlotsYesNo','Yes','No', 'Cancel', 'No');
        if strcmp(plotYesNo, 'Cancel')
            return
        end
        
        % Choose if you want to show a video of the face live, while
        % processing
        showVid = questdlg('Do you want to see a video of the face while processing?', ...
            'plotFaceFrames','Yes','No', 'Cancel', 'No');
        if strcmp(plotYesNo, 'Cancel')
            return
        end
        
        switch showVid
            case        'Yes'
                showVid = 1;
            case 'No'
                showVid = 0;
                
            case 'Cancel'
                return
        end
        
        % Choose if you want to plot each frame of the face video
        makeFaceVid = questdlg('Do you want to make a video with face and detected points?', ...
            'makeFaceFramesVideo','Yes','No', 'Cancel', 'No');
        
        switch makeFaceVid
            case        'Yes'
                makeFaceVid = 1;
            case 'No'
                makeFaceVid = 0;
                
            case 'Cancel'
                return
        end
        
        % Select the folder with all the videos
        videosFolder = uigetdir(scriptFolder, 'Select the subfolders containing the videos and Raw Data');
        if videosFolder == 0
            return
        end
        % Chose whether to overwrite previously extrated RGB values
        overwriteRGB = questdlg('Do you want to overwrite already processed raw data?', 'Overwrite Data', ...
            'true','false', 'Cancel', 'false');
        if strcmp(overwriteRGB, 'Cancel')
            return
        else
            overwriteRGB= string2boolean(overwriteRGB);
        end
        
        % Enter a guard for the peaks (in Seconds)
        
        userAnswer = inputdlg({'Insert a winsize for the Time Frequency Analysis:'}, ...
            'WinSize',[1 35],{'10'});
        if isempty(userAnswer)
            warning('You must provide a guard in seconds (e.g. 0.45) and a winSize (e.g. 10)')
            return
        end
        
        winSize = str2double(userAnswer);
        
        clear userAnswer
        
        % Choose a name for a file to save
        resultsFilename = char(inputdlg({'Choose a name for  the results file:'}, ...
            'Results Filename',[1 35],{'allSbjResults'}));
        
        % Choose a name for the folder containing all the plots
        expPlotsFolder = char(inputdlg({'Choose a name for  the folder containing all the plots:'}, ...
            'Plots folder name',[1 35],{'thisExpPlots'}));
        
        %         peakGuardInSecs = str2double(inputdlg({'Insert a guard in secs between peaks detected:'}, 'Peaks Guard', ...
        %             [1 35],{'0.45'}));
        %         if isempty(peakGuardInSecs)
        %             return
        %         end
        %
        %         % Enter a Window Size for the Time Frequency Analysis
        %         signalProcessing.tfa.winSize = str2double(inputdlg({'Insert a winsize for the Time Frequency Analysis:'},...
        %             'winSize', [1 35],{'10'}));
        %
        %         if strcmp(signalProcessing.tfa.winSize , 'Cancel')
        %             return
        %         end
        
        % Create a table with all the settings and save it for the future
        scriptSettings = table(string(plotYesNo), string(showVid), string(makeFaceVid), string(videosFolder),...
            string(overwriteRGB), winSize, string(resultsFilename), string(expPlotsFolder));
        scriptSettings.Properties.VariableNames = {'plotYesNo', 'plotFaceFrames', 'makeFaceVid', ...
            'videosFolder', 'overwriteRGB', 'winSize', 'resultsFilename', 'expPlotsFolder'};
        
        saveSettings = questdlg('Do you want to save these settings in a txt file?', 'saveSettings','Yes','No', 'Cancel', 'No');
        if strcmp(saveSettings, 'Cancel')
            return
        end
        
        
        if strcmp(saveSettings, 'Yes')
            writetable(scriptSettings, 'scriptSettings' ,'Delimiter', '\t');
            disp(['Settings saved in ''scriptSettings.txt'' in the folder: ' pwd])
        elseif strcmp(saveSettings, 'Cancel')
            return
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% User has chosen to select parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 2
        % Load a previously saved txt file containing settings
        scriptSettings = readtable([scriptFolder, filesep, 'scriptSettings']);
        plotYesNo = string(scriptSettings.plotYesNo);
        showVid = scriptSettings.plotFaceFrames;
        makeFaceVid = scriptSettings.makeFaceVid;
        videosFolder = char(scriptSettings.videosFolder);
        overwriteRGB = string2boolean(scriptSettings.overwriteRGB);
        winSize = scriptSettings.winSize;
        resultsFilename = char(scriptSettings.resultsFilename);
        expPlotsFolder = char(scriptSettings.expPlotsFolder);
        
    case 3
        % Use default settings
        disp(sprintf(['WARNING:\nif you chose default, the sript will look for videos in the subfolder named "Video"\n', ...
            'It will not make plots, it will not show the face video, it will save the face video, it will not overwrite RGB\n', ...
            'it will save the results as excel with the name allSbjResults,\nit will export the plots in a folder named: thisExpPlots']));
        plotYesNo = 'No';
        showVid = 0;
        makeFaceVid = 1;
        videosFolder = [pwd filesep 'Video'];
        overwriteRGB = false;
        winSize = 10;
        resultsFilename = 'allSbjResults';
        expPlotsFolder = char('thisExpPlots');
        
        defaultParamText = ['Plot Figures = No', newline, 'Videos Folder Path : ' videosFolder, ...
            newline, 'OverwriteRGB = false', newline ...
            'Signal Procesing tfa window size = 10', newline ...
            'results filename = allSbjResults'];
        
        disp(defaultParamText);
end



scriptSettings = readtable([scriptFolder, filesep, 'scriptSettings']);
plotYesNo = string(scriptSettings.plotYesNo);
showVid = scriptSettings.plotFaceFrames;
makeFaceVid = scriptSettings.makeFaceVid;
videosFolder = char(scriptSettings.videosFolder);
overwriteRGB = string2boolean(scriptSettings.overwriteRGB);
winSize = scriptSettings.winSize;
resultsFilename = char(scriptSettings.resultsFilename);
expPlotsFolder = char(scriptSettings.expPlotsFolder);



if ~(exist('Raw data', 'dir'))
    mkdir('Raw data');
end

if ~(exist(['Raw data' filesep expPlotsFolder], 'dir'))
    mkdir(['Raw data' filesep expPlotsFolder]);
end

if ~(exist('plots', 'dir'))
    mkdir('plots');
end
cd('plots')

if ~(exist(expPlotsFolder, 'dir'))
    mkdir(expPlotsFolder);
end

% and prepare the results folder
if  ~(exist([scriptFolder filesep 'results'], 'dir'))
    mkdir('results');
end
cd(scriptFolder)