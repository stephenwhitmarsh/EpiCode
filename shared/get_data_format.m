function [isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg)
% Return 1 in the respective variable if the format corresponds.
% Otherwise returns 0.
% Format used : Neuralynx (.ncs), Micromed (.TRC), Brainvision (.eeg)
% Error if format is different

datapath = fullfile(cfg.rawdir,cfg.directorylist{1}{1});
listing = dir(datapath);
% in some cases cfg.directorylist is list of files and not folders. To fix that :
if isempty(listing)
    datapath = cfg.rawdir;
    listing = dir(datapath);
end
if isempty(listing)
    error('No file detected in %s or in %s',fullfile(cfg.rawdir,cfg.directorylist{1}{1}), cfg.rawdir);
end

isNeuralynx = 0;
isMicromed = 0;
isBrainvision = 0;

for ifile = 1:length(listing)
    [~,~,file_extension] = fileparts(listing(ifile).name);
    if strcmp(file_extension, '.ncs')
        isNeuralynx        = 1;
    elseif strcmp(file_extension, '.TRC')
        isMicromed         = 1;
    elseif strcmp(file_extension, '.eeg')
        isBrainvision      = 1;
    end
end

if isNeuralynx
    fprintf('Data is Neuralynx\n');
elseif isMicromed
    fprintf('Data is Micromed\n');
elseif isBrainvision
    fprintf('Data is Brainvision\n');
end

if isNeuralynx + isMicromed + isBrainvision == 0
    error('Cannot detect good data format in datapath = %s \nData have to be Neuralynx (.ncs), Micromed (.TRC) or Brainvision (.eeg)\n', datapath);
elseif isNeuralynx + isMicromed + isBrainvision >1
    error('Several data formats are detected in datapath = %s \nData are Neuralynx (.ncs), Micromed (.TRC) or Brainvision (.eeg)\n', datapath);
end


end

