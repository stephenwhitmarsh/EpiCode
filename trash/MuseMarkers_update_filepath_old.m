
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

for ipart = 1 : size(cfg.directorylist,2)   
    for type = ["micro","macro"]       
        for idir = 1 : size(cfg.directorylist{ipart},2)
            eval(sprintf('MuseStruct_%s{ipart}{idir}.directory = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir});',type))
        end
    end
end
