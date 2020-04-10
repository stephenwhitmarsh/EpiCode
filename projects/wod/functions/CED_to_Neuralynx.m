%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script  CED_to_Neuralynx
% REVOIR LE TEXTE
%
% All channels are converted (list of chan types available)
%
% Convert data from Spike2 to MATLAB, in a format compatible with
% Fieldtrip. It allows to use Epicode's functions ((c) Stephen Whitmarsh).
% Call this function before any other use of Epicode.
%
% The data are converted and wrote to disk channel by channel, to avoid
% excess of RAM use. Only continuous data can be converted (adc channel or
% realwave channel). To convert events timing, use readCEDMarkers.m.
%
% ## INPUT :
% config{isubject}.prefix                    : name of the subject
% config{isubject}.CEDrawdir                 : where are the Spik2 data
% config{isubject}.directorylist{part}       : list of all data files, for each part.
%                                 Extension of file name must not be written.
% config{isubject}.rawdir                    : where to store the converted data
% Names of channels to convert.
%     - if there is a white space in Spike2 channel name, replace it by '_'
%     - if channel has no name, call it 'chan%d', where %d is the channel
%       number
% force                         : if force == true, force converting again
%                               the channel.
%
% Need of CEDS64ML interface library (loaded with CEDS64LoadLib.m), and of
% Spike2 software. Can only be ran on Windows (because of the library).
%
% Paul Baudin
% paul.baudin@live.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% To be modified by the user
% Paths needed : Fieldtrip, EpiCode external and shared, CEDS64ML
addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\wod'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML
CEDS64LoadLib('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML');

config                      = wod_setparams([]); %optionnal script with parameters of all subjects

for isubject = 1
    % if no setparams script, set parameters here :
    %     config{isubject}.prefix             = ;
    %     config{isubject}.CEDrawdir          = ;
    %     config{isubject}.directorylist      = ;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ipart = 1 : size(config{isubject}.directorylist,2)
        
        % Go through directory list
        for idir = 1 : size(config{isubject}.directorylist{ipart},2)
            
            %% Convert all continuous channels
                
            datapath =  fullfile(config{isubject}.CEDrawdir, config{isubject}.directorylist{ipart}{idir});
            fprintf('Loading data from %s \nPart %d, dir %d\n', datapath, ipart, idir);
            
            %open the file in Spike2.
            fid = CEDS64Open(datapath);
            if fid<0
                error('error while opening file %s. \nCheck that the path is correct, and that this file is not already opened in Spike2.',datapath);
            end
            
            % Go through each Spike2 channel
            n_loaded_chan = 0;
            
            for ichan = 1:CEDS64MaxChan(fid)
                
                % get channel name and rename it if needed
                chantitle = [];
                [~, chantitle] = CEDS64ChanTitle(fid, ichan);
                % if chan has no name
                if isempty(chantitle)
                    chantitle = sprintf('chan%d', ichan);
                end
                % remove white spaces
                if any(ismember(' ',chantitle))
                    chantitle = strrep(chantitle,' ','_');
                end
                
                
                
                %check channel type
                [iType] = CEDS64ChanType(fid, ichan);
                if ismember(iType, [1, 9]) %Waveform or RealWave
                    
                    n_loaded_chan = n_loaded_chan + 1;
                    
                    fname = fullfile(config{isubject}.rawdir,sprintf('%srawdata-p%d-d%d-%s.m',config{isubject}.prefix,ipart,idir,chantitle));
                    
                    
                    
                    
                    %load data
                    chanFs              = 1 / (CEDS64ChanDiv(fid, ichan) * CEDS64TicksToSecs(fid,1)); %inverse de : nb de tick par sample * durée du tick
                    maxtime             = CEDS64TicksToSecs(fid, CEDS64MaxTime(fid));
                    n_samples           = round(maxtime*chanFs);
                    
                    [n_samples, rawdata.trial{1}, ~] = CEDS64ReadWaveS(fid, ichan, n_samples, 0);
                    rawdata.trial{1}                 = rawdata.trial{1}';
                    rawdata.fsample                  = chanFs;
                    rawdata.sampleinfo               = [1,n_samples];
                    
                    %recover time in seconds
                    rawdata.time{1} = 0 : 1/chanFs : (n_samples-1)/chanFs; %-1 because begins at zero
                    
                    rawdata.label{1} = chantitle;
                    
                    % check if rawdir directory exists, if not create
                    if ~isfolder(config{isubject}.rawdir)
                        ft_notice('creating directory %s', config{isubject}.rawdir);
                        mkdir(config{isubject}.rawdir);
                    end
                    save(fname,'rawdata');
                    clear rawdata;
                    
                end
                
                %Convertir données en µV
                % Ecrire fichier comme Katia
                % un dossier par fichier spike 2. Toutes les électrodes et
                % les events dans ce dossier
                ft_write_data();
                
            end %ichan
            
            %% read events
            
            
        end %idir
    end %ipart
end %isubject


