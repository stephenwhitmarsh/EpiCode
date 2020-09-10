function writeSpykingCircusParameters(cfg)

% WRITESPYKINGCIRCUSPARAMETERS writes .params and .prb file for Spyking Circus
%
% use as
%   writeSpykingCircusParameters(cfg)
%
% The .params file is based on the one in the template folder of EpiCode.
% Change this file to have new default values.
%
% Mandatory inputs
% cfg.prefix            % name of the data to analyse, will be appended at
%                         the begining of each data file
% cfg.directorylist     % list of folders with the neuralynx raw files (one
%                         file per electrode)
% cfg.circus.channel    % list of channels to process with Spyking-Circus
%                         (the first one is used to name the .param file)
% cfg.datasavedir       % where to save the data. The folder with Spyking-
%
% Optional cfg fields
% cfg.circus.part_list  % Array of integers with the parts numbers to
%                         analyze. Can be 'all'. Default = 'all'

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%   EpiCode is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   EpiCode is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

% get the default options
cfg.circus.part_list                = ft_getopt(cfg.circus, 'part_list', 'all');

if strcmp(cfg.circus.part_list, 'all')
    cfg.circus.part_list = 1:size(cfg.directorylist,2);
end

[p, ~, ~]               = fileparts(mfilename('fullpath'));
[p, ~, ~]               = fileparts(p);
fname_params_default    = fullfile(p,'templates','SpykingCircus.params');

for ipart = 1 : cfg.circus.part_list

    subjdir         = cfg.prefix(1:end-1);
    partdir         = ['p',num2str(ipart)];
    filename        = [cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1},'.params'];
    fname_params    = fullfile(cfg.datasavedir,subjdir,partdir,filename);
    nb_channels     = size(cfg.circus.channel,2);
    fname_prb       = ['Adtech_', num2str(nb_channels), 'chan.prb'];

    % read Spyking-Circus params file
    ini = IniConfig();
    ini.ReadFile(fname_params_default);

    % remove inline comments
    [sections, count_sections] = ini.GetSections();
    for sectioni = 1 : count_sections
        [keys, count_keys] = ini.GetKeys(sections{sectioni});
        for keysi = 1 : count_keys
            old = ini.GetValues(sections{sectioni}, keys{keysi});
            temp = split(old,{'#'});
            ini.SetValues(sections{sectioni}, keys{keysi}, temp{1});
        end
    end

    % adjust parameters
    h1 = ini.SetValues('data', {'file_format','stream_mode','mapping','suffix','overwrite','output_dir'}, {'neuralynx','None',fname_prb,'','False','SpykingCircus'});
    h2 = ini.SetValues('noedits', {'filter_done','artefacts_done','ground_done','median_done'}, {'False','False','False','False'});
    h3 = ini.SetValues('triggers', {'dead_file','dead_unit','ignore_times'}, {'SpykingCircus_artefacts_samples.dead','timestep','True'});
    if any([h1; h2; h3] ~= 1), error('Something went wrong with adjusting parameters'); end

    status = ini.WriteFile(fname_params);
    if status == false
        error('Couldn''t write file');
    end
    ini.ToString()

    % write params file
    writeProbeFile(nb_channels,fullfile(cfg.datasavedir,subjdir,partdir,fname_prb));

end
