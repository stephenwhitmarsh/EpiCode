function [isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg)
% Return 1 in the respective variable if the format corresponds.
% Otherwise returns 0.
% Format used : Neuralynx (.ncs), Micromed (.TRC), Brainvision (.eeg)
% Error if format is different


% List all the files in the input path. 2 conditions : if dir is a file with
% all channels or if dir is a folder with a file per channel
listing = dir(fullfile(cfg.rawdir,cfg.directorylist{1}{1})); 
if isempty(listing) 
    listing = dir(cfg.rawdir);
end
if isempty(listing) 
    error('No file detected in %s or in %s',fullfile(cfg.rawdir,cfg.directorylist{1}{1}), cfg.rawdir);
    halt;
end

% remove "." and ".." files
inds = [];
k    = 1;

while k <= length(listing) 
    if strcmp(listing(k).name(1), '.')
        inds(end + 1) = k;
    end
    k = k + 1;
end

listing(inds) = [];

% remove all files exept .ncs .eeg .TRC
inds = [];
k    = 1;

while k <= length(listing) 
    if ~any(strcmp(listing(k).name(end-3:end), {'.ncs', '.TRC', '.eeg'})) 
        inds(end + 1) = k;
    end
    k = k + 1;
end

listing(inds) = [];

% check format according to file extension
if ~isempty(listing)
    datafile_extension  = listing(1).name(end-3:end);
    
    isNeuralynx = 0;
    isMicromed = 0;
    isBrainvision = 0;
    
    if strcmp(datafile_extension, '.ncs') %is neuralynx
        isNeuralynx = 1;
    elseif strcmp(datafile_extension, '.TRC')
        isMicromed = 1;
    elseif strcmp(datafile_extension, '.eeg')
        isBrainvision = 1;
        
    end
    
else
    error('Cannot detect good data format in datapath = %s \nData have to be Neuralynx (.ncs), Micromed (.TRC) or Brainvision (.eeg)', datapath);
    halt;
    
end % ~isempty(listing)

end

