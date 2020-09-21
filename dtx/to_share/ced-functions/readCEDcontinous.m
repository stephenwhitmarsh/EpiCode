function chandata = readCEDcontinous(cfg)

% Need of CEDS64ML interface library (loaded with CEDS64LoadLib.m), and of
% Spike2 software. Can only be ran on Windows (because of the library).
% FORCE

% cfg.channame = name of the channel to load

% cfg.datapath : optional field to find the data path. If it is not
% specified, this scripts uses the 'EpiCode' way of finding files : 
% cfg.rawdir and cfg.directorylist

% cfg.ipart : part number (only one part), used only if cfg.datapath is not specified
% cfg.idir  : dir number (only one dir), used if cfg.datapath is not specified

chandata = [];
n_loaded_chan = 0;

if ~isfield(cfg, 'datapath')
    temp = dir(fullfile(cfg.rawdir, [cfg.directorylist{cfg.ipart}{cfg.idir}, '.smr*']));%.smr* because some data are .smr and other .smrx
    datapath = fullfile(cfg.rawdir, temp.name);
    fprintf('Loading data from %s \nPart %d, dir %d, channel %s\n', datapath, cfg.ipart, cfg.idir, cfg.channame);
else
    datapath = cfg.datapath;
    fprintf('Loading data from %s\n', datapath);
end


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
%     chantitle  = renamechan_CED(chantitle,ichan,[],false);
% test{ichan}=chantitle;
    if strcmp(chantitle, cfg.channame)
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
%             [n_samples, chandata.trial{1}, ~] = CEDS64ReadWaveS(fid, ichan, n_samples, 0);
            chandata.trial{1}                 = chandata.trial{1}';
            chandata.fsample                  = chanFs;
            chandata.sampleinfo               = [1,n_samples];
            
            %recover time in seconds
            chandata.time{1} = 0 : 1/chanFs : (n_samples-1)/chanFs; %-1 because begins at zero
            chandata.label{1} = chantitle;
            [~,chandata.chanunit{1}] = CEDS64ChanUnits(fid,ichan);
            
        end
                
    end %strcmp cfg.channame chantitle
    
end %ichan

if n_loaded_chan == 0
    warning('Data of channel %s was not found in %s', cfg.channame, datapath);
%     warning('Be aware that some channel names are modified by the function CEDrename_chan.m');
elseif n_loaded_chan > 1
    warning('%d channels have the name %s. Only the last one (with the higher Spike2 chan nr) was loaded', n_loaded_chan, cfg.channame);
end

CEDS64Close(fid);