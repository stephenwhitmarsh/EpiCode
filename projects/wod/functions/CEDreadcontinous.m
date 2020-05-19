function chandata = CEDreadcontinous(cfg,channame,ipart,idir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script  CED_to_Neuralynx
% REVOIR LE TEXTE
%
% All channels are converted (list of chan types available)
%
% Convert data from Spike2 to MATLAB, in a format compatible with
% Fieldtrip.
%
% The data are converted and wrote to disk channel by channel, to avoid
% excess of RAM use. Only continuous data can be converted (adc channel or
% realwave channel). To convert events timing, use readCEDMarkers.m.
% NOT WROTE TO DISK. CHAN BY CHAN BECAUSE OF DIFFERENT FS
%
% ## INPUT :
% cfg.prefix                    : name of the subject
% cfg.CEDrawdir                 : where are the Spik2 data
% cfg.directorylist{part}       : list of all data files, for each part.
%                                 Extension of file name must not be written.
% cfg.rawdir                    : where to store the converted data
% Names of channels to convert.
%     - if there is a white space in Spike2 channel name, replace it by '_'
%     - if channel has no name, call it 'chan%d', where %d is the channel
%       number
% force                         : if force == true, force converting again
%                               the channel.
%
% Need of CEDS64ML interface library (loaded with CEDS64LoadLib.m), and of
% Spike2 software. Can only be ran on Windows (because of the library).
% FORCE
%
% Paul Baudin
% paul.baudin@live.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


chandata = [];
n_loaded_chan = 0;

temp = dir(fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir}, '.smr*']));%.smr* because some data are .smr and other .smrx
datapath = fullfile(cfg.rawdir, temp.name);
fprintf('Loading data from %s \nPart %d, dir %d, channel %s\n', datapath, ipart, idir, channame);

%open the file in Spike2.
fid = CEDS64Open(datapath);
if fid<0
    error('error while opening file %s. \nCheck that the path is correct, and that this file is not already opened in Spike2.',datapath);
end

% Go through each Spike2 channel


for ichan = 1:CEDS64MaxChan(fid)
    
    % get channel name and rename it if needed
    chantitle = [];
    [~, chantitle] = CEDS64ChanTitle(fid, ichan);
    chantitle  = CEDrename_chan(chantitle,ichan,[],false);
       
    if strcmp(chantitle, channame)
        
        %check channel type
        [iType] = CEDS64ChanType(fid, ichan);
        if ismember(iType, [1, 9]) %Waveform or RealWave
            chandata = [];
            n_loaded_chan = n_loaded_chan+1;
            
            %load data
            chanFs              = 1 / (CEDS64ChanDiv(fid, ichan) * CEDS64TicksToSecs(fid,1)); %inverse de : nb de tick par sample * durée du tick
            maxtime             = CEDS64TicksToSecs(fid, CEDS64MaxTime(fid));
            n_samples           = round(maxtime*chanFs);
            
            [n_samples, chandata.trial{1}, ~] = CEDS64ReadWaveF(fid, ichan, n_samples, 0);
            chandata.trial{1}                 = chandata.trial{1}';
            chandata.fsample                  = chanFs;
            chandata.sampleinfo               = [1,n_samples];
            
            %recover time in seconds
            chandata.time{1} = 0 : 1/chanFs : (n_samples-1)/chanFs; %-1 because begins at zero
            chandata.label{1} = chantitle;
            [~,chandata.chanunit{1}] = CEDS64ChanUnits(fid,ichan);
            
        end
                
    end %strcmp channame chantitle
    
end %ichan

if n_loaded_chan == 0
    warning('Data of channel %s was not found in %s', channame, datapath);
    warning('Be aware that some channel names are modified by the function CEDrename_chan.m');
elseif n_loaded_chan > 1
    warning('%d channels have the name %s. Only the last one (with the higher Spike2 chan nr) was loaded', n_loaded_chan, channame);
end