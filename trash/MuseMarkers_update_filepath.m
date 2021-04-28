
function [MuseStruct_micro, MuseStruct_macro]  = MuseMarkers_update_filepath(cfg,MuseStruct_micro, MuseStruct_macro)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function readMuseMarkers
%
% Crawls and searches through patient- and recording-directories to
% extracts marker timings (time and samples) created by Muse. These can be
% used to create an overview of events (plotmarkers.m), segment data with
% FieldTrip (e.g. plotpeaks.m) and create artefact files for Spyking-Circus
% (writeSpykingCircusDeadFile.m). The resultant MuseStruct can be edited
% and written into a new markerfile for Muse with writeMuseMarkers.m
%
% Code by Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
% with help from Jean-Didier Lemarechal & Craig Richter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all');

%% Micro electrodes, after repeat same for macro electrodes

for type = ["micro","macro"]
    
    % Go through different parts
    for ipart = 1 : size(cfg.directorylist,2)
        
        % Go through directory list
        for idir = 1 : size(cfg.directorylist{ipart},2)
            
            % check if directory has all requested files
            hasfiles = 0;
            channr = 1;
            for ichan = 1 : size(cfg.labels.(type),2)
                    d = dir2(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.labels.(type){ichan},'*.ncs']));

                    if ~isempty(d)
                        
                        % update filename and path
                        eval(sprintf('MuseStruct_%s{ipart}{idir}.directory = d.folder;',type))
                        eval(sprintf('MuseStruct_%s{ipart}{idir}.filenames{channr} = d.name;',type))
                        fprintf('Extracting file names of %s \n',fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}));
                        
                        channr = channr + 1;
                        hasfiles = 1;
                    else
                        fprintf('*** WARNING: cannot find channel %s data in %s\n',cfg.labels.(type){ichan},fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}));
                    end
            end
        end
    end
end