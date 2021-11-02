function [config] = addparts(config)

temp        = dir(config.rawdir);
temp        = temp([temp.isdir]);
temp(1:2)   = []; % remove . and ..
clear flist
for i = 1 : size(temp, 1)
    flist(i, :) = temp(i).name;
end

[~, locb]       = ismember(config.directorylist{3}{end}, flist, 'rows');
flist           = flist(locb+1:end, :);
flist_datetime  = datetime(flist(:, 7:end), 'inputformat', 'yyyy-MM-dd_HH-mm', 'Format', 'd-MMM-y HH:mm:ss');
ipart_list      = floor(hours(flist_datetime - flist_datetime(1)) / 24) + 4;

% get file format
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(config);

for ipart = unique(ipart_list)'
    config.directorylist{ipart} = cellstr(flist(ipart_list == ipart, :))';
end

% check parts for data
for ipart = 1 : size(config.directorylist, 2)
    
    % loop over directories
    hasdata_dir{ipart} = false(size(config.directorylist{ipart}, 2), 1);
    
    for idir = 1 : size(config.directorylist{ipart}, 2)
        
        if isNeuralynx
            % check for all files
            hasdata_file = false(size(config.cluster.channel, 2), 1);
            for ifile = 1 : size(config.cluster.channel, 2) % one file per channel
                if isNeuralynx
                    temp = dir(fullfile(config.rawdir, config.directorylist{ipart}{idir}, ['*', config.cluster.channel{ifile},'.ncs']));
                    if ~isempty(temp)
                        hasdata_file(ifile) = true;
                    end
                end
            end
            if all(hasdata_file)
                hasdata_dir{ipart}(idir) = true;
                fprintf('Part %d, %s contains requested data\n', ipart, fullfile(config.rawdir, config.directorylist{ipart}{idir}));                                
            else
                fprintf('Part %d, %s does not contain requested data\n', ipart, fullfile(config.rawdir, config.directorylist{ipart}{idir}));
            end
        elseif isMicromed
            temp = dir(fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '.TRC']));
            if ~isempty(temp)
                hasdata_dir{ipart}(idir) = true;
                fprintf('Part %d, %s contains requested data\n', ipart, fullfile(config.rawdir, config.directorylist{ipart}{idir}));                                
            else
                fprintf('Part %d, %s does not contain requested data\n', ipart, fullfile(config.rawdir, config.directorylist{ipart}{idir}));
            end
        elseif isBrainvision
            temp = dir(fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '.eeg']));
            if ~isempty(temp)
                hasdata_dir{ipart}(idir) = true;
                fprintf('Part %d, %s contains requested data\n', ipart, fullfile(config.rawdir, config.directorylist{ipart}{idir}));                
            else
                fprintf('Part %d, %s does not contain requested data\n', ipart, fullfile(config.rawdir, config.directorylist{ipart}{idir}));
            end
        end
        
        
    end  
end
% 
for ipart = 1 : size(config.directorylist, 2)
    config.directorylist{ipart} = config.directorylist{ipart}(hasdata_dir{ipart});
    hasdata_part(ipart) = any(hasdata_dir{ipart})    ; 
end
config.directorylist = config.directorylist(hasdata_part);

