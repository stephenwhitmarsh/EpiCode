%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deltamed2Brainvision
%
% Convert Deltamed EEG files to Brainvision EEG files.
% Convert only the channels of interest, and rename it, according to the
% config set in dtx_setparams_EEGvideo.m function.
%
% MUST BE LAUNCHED WITH 32BIT MATLAB R2015b (because of Deltamed
% Coherence5LE library).
%
% /!\ ATTENTION /!\ : the fieldtrip fonction write_brainvision_eeg has been
% modified for this script, in order to allow the reading of the header 
% and the data by Muse :
% - add of ',µV' at the end of each 'Channel Info' line in the vhdr file
% - data are written in int16 and not in float32 
% Modified lines in the fieldtrip fonction write_brainvision_eeg are
% commented with '% katia'
% Replace write_brainvision_eeg in the fieldtrip folder  by  the script 
% in this folder : write_brainvision_eeg_modified
%
% Paul Baudin and Katia Lehongre 01.2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath C:\Users\paul.baudin\Documents\MATLAB\fieldtrip;
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\epilepsy'));
addpath C:\Users\paul.baudin\Documents\MATLAB\DTX;
addpath(genpath('C:\Users\paul.baudin\Documents\MATLAB\Deltamed_to_Brainvision'));
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\Deltamed_to_Brainvision\Coherence5LE ICM\'));
addpath C:\Users\paul.baudin\Documents\MATLAB\Deltamed_to_Brainvision\Function_Coherence5LE

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for conversion

config = dtx_setparams_EEGvideo([]);
MrkTemplate = '\\lexport\iss01.epimicro\rodents\raw\TemplateEvents.mrk';


for irat = 2
    
    OrigDir     = fullfile('\\lexport\iss01.epimicro\rodents\raw\DTX-EEGRodents-Deltamed',config{irat}.foldername);
    OutputDir   = config{irat}.rawdir ;
    filelist = dir(OrigDir);
    iconversion=0;
    
    if ~(exist (OutputDir)==7)
        mkdir(OutputDir);
    else
        warning('Output folder already exists\n%s',OutputDir);
    end
    
    
    for ifile = 1:length(filelist) 
        if length(filelist(ifile).name)>3
            fileExtension = filelist(ifile).name(1,length(filelist(ifile).name)-3:length(filelist(ifile).name));
            if strncmp(fileExtension,'.EEG',4) ||  strncmp(fileExtension,'.eeg',4)
                iconversion=iconversion+1;
                fprintf('\n***************************************************************\n');
                fprintf('*********** Converting file n°%d : %s **********',iconversion,filelist(ifile).name);
                fprintf('\n***************************************************************\n\n');
                
                outputpath = fullfile(OutputDir,filelist(ifile).name(1,1:length(filelist(ifile).name)-4));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Get Deltamed data
                fprintf('\nLoading Deltamed data\n');
                
                datapath = fullfile(OrigDir,filelist(ifile).name);
                z = Coherence5LE();
                z = OpenFile(z,datapath); 
                data = GetData(z);
                data = data';
                FileInfo = GetFileInfo(z);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Create header
                %rajouter une partie pour les EMG
                fprintf('\nCreating header : %s\n', [outputpath, '_header']);
                %Search for electrodes which correspond to the chosen electrodes in the
                %config parameters
                i_outputchannel = 0;
                for ichannel=1:length(data(:,1))
                    ichannel_name = FileInfo.name(ichannel,:);
                    
                    for irawlabel=1:length(config{irat}.rawlabels)
                        
                        %find channel corresponding to config
                        if strfind(ichannel_name,config{irat}.rawlabels{irawlabel})==1
                            i_outputchannel = i_outputchannel+1;
                            channel_index(i_outputchannel)=ichannel;
                            
                            
                            %keep only letter or number of label and unit (remove "empty")
                            i_temp=1;
                            for icharacter = 1:length(FileInfo.name(i_outputchannel,:))
                                if isletter(FileInfo.name(i_outputchannel,(icharacter)))||ismember(FileInfo.name(i_outputchannel,(icharacter)),['0' '1' '2' '3' '4' '5' '6' '7' '8' '9'])
                                    header.label{i_outputchannel,1}(i_temp)=FileInfo.name(i_outputchannel,(icharacter));
                                    i_temp = i_temp+1;
                                end
                            end
                            header.label(i_outputchannel,1) = {config{irat}.labels.macro{irawlabel}}; %FileInfo.name(i_outputchannel,:);
                            i_temp=1;
                            for icharacter = 1:length(FileInfo.unit(i_outputchannel,:))
                                if isletter(FileInfo.unit(i_outputchannel,(icharacter)))
                                    header.chanunit{i_outputchannel,1}(i_temp)=FileInfo.unit(i_outputchannel,(icharacter));
                                    i_temp = i_temp+1;
                                end
                            end
                            
                            header.chantype(i_outputchannel,1) = {'eeg'}; %eeg
                            header.rawlabel(i_outputchannel,1) = {config{irat}.rawlabels{irawlabel}}; 
                        end
                    end 
                end
                
                header.nSamples = length(data(1,:));
                header.nSamplesPre = 0;
                header.nTrials = 1;
                header.Fs = FileInfo.frequency;
                header.nChans = i_outputchannel;
                header.originfilename = datapath;
                header.filename = outputpath;
                
                Y = str2num(FileInfo.date(1,16:19));
                M = str2num(FileInfo.date(1,13:14));
                D = str2num(FileInfo.date(1,10:11));
                H = str2num(FileInfo.date(1,1:2));
                MI = str2num(FileInfo.date(1,4:5));
                S = str2num(FileInfo.date(1,7:8));
                
                header.date = datetime(Y,M,D,H,MI,S);
                header.dateEnd = header.date+seconds(header.nSamples/header.Fs);
                
                %Save header (because the header created by fieldtrip loose
                %information)
                save([outputpath, '_header'],'header');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Conversion to Brainvision with fieldtrip
                fprintf('\nWriting data in the new format : %s \n',[outputpath, '.vhdr']);
                %convert only the channels of interest
                data = data(channel_index,:);
                
                ft_write_data(outputpath,data,'header',header,'dataformat','brainvision_eeg');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Create Muse marker file
                if ~exist([outputpath, '.mrk']) %to avoid overdrawing of already-placed markers
                    copyfile(MrkTemplate, [outputpath, '.mrk']);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Create timing array to check if there is missing data between 2 files
                File_name(iconversion)= {filelist(ifile).name(1,1:length(filelist(ifile).name)-4)};
                File_beginning(iconversion)={datestr(header.date)};
                File_ending(iconversion)={datestr(header.dateEnd)};
                
                
                clear z data FileInfo header channel_index
                
                
            end
        end
    end %ifile
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Create timing file to check if there is missing data between 2 files
    fprintf('\n\n Create timing file to check if there is missing data between 2 files\n\n');
    
    tablepath = fullfile(OutputDir, [config{1}.prefix 'Timings.csv']);
    tableid = fopen(tablepath,'wt');
    
    fprintf(tableid,'File;Tstart;Tend;Time difference from previous file\n');
    fprintf(tableid,'%s;%s;%s\n',File_name{1},File_beginning{1},File_ending{1});
    
    if iconversion>1
        for i=2:iconversion
            file_ending = datetime(File_ending{i-1},'InputFormat','d-MMM-y HH:mm:ss');
            nextfile_beginning = datetime(File_beginning{i},'InputFormat','d-MMM-y HH:mm:ss');
            timediff = nextfile_beginning - file_ending;
            if timediff >=0  
                fprintf(tableid,'%s;%s;%s;%s\n',File_name{i},File_beginning{i},File_ending{i},datestr(timediff,'HH:MM:SS'));
            else
                timediff = timediff*(-1);
                fprintf(tableid,'%s;%s;%s;-%s\n',File_name{i},File_beginning{i},File_ending{i},datestr(timediff,'HH:MM:SS'));
            end
        end
    end
    
    
    fclose('all');
    
end %irat
