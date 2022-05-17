function data_format = get_data_format2(cfg)

% get_data_format returns 1 in the respective variable if the format corresponds
% Otherwise returns 0.
% Use as
%   [isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg)
%
% Format used : Neuralynx (.ncs), Micromed (.TRC), Brainvision (.eeg), CED
% (.smr or .smrx)
% Error if format is different

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%    EpiCode is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    EpiCode is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

datapath = fullfile(cfg.rawdir, cfg.directorylist{1}{1});
listing = dir(datapath);

% in some cases cfg.directorylist is list of files and not folders. To fix that :
if isempty(listing)
    datapath = cfg.rawdir;
    listing = dir(datapath);
end
if isempty(listing)
    error('No file detected in %s or in %s', fullfile(cfg.rawdir,cfg.directorylist{1}{1}), cfg.rawdir);
end

isNeuralynx   = 0;
isMicromed    = 0;
isBrainvision = 0;
isCED         = 0;

for ifile = 1:length(listing)
    [~,~,file_extension] = fileparts(listing(ifile).name);
    if strcmp(file_extension, '.ncs')
        isNeuralynx     = 1;
    elseif strcmp(file_extension, '.TRC')
        isMicromed      = 1;
    elseif strcmp(file_extension, '.eeg')
        isBrainvision   = 1;
    elseif strcmp(file_extension, '.smr') || strcmp(file_extension, '.smrx')
        isCED           = 1;
    end
end

if isNeuralynx
    ft_info('Data is Neuralynx\n');
    data_format = 'neuralynx';
elseif isMicromed
    ft_info('Data is Micromed\n');
    data_format = 'micromed';
elseif isBrainvision
    ft_info('Data is Brainvision\n');
    data_format = 'brainvision';
elseif isCED
    ft_info('Data is CED\n');
    data_format = 'ced';
end

if isNeuralynx + isMicromed + isBrainvision + isCED == 0
    error('Cannot detect good data format in datapath = %s \nData have to be Neuralynx (.ncs), Micromed (.TRC) or Brainvision (.eeg)\n', datapath);
elseif isNeuralynx + isMicromed + isBrainvision + isCED >1
    error('Several data formats are detected in datapath = %s \nData are Neuralynx (.ncs), Micromed (.TRC) or Brainvision (.eeg)\n', datapath);
end


end
