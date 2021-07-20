function LFP = dtx_readLFP(cfg, MuseStruct, force)

%read LFP data, remove artefacted trials, flip if required, correct
%baseline, and align interictal data to the following TDW

%used in dtx_spikes_project

LFP = readLFP(cfg, MuseStruct, force);
LFP = removeArtefactedTrials(cfg, LFP);
if ~isempty(LFP)
    for ipart = 1:size(LFP)
        for markername = string(fieldnames(LFP{ipart}))'
            if isempty(LFP{ipart}.(markername))
                continue
            end
            if istrue(cfg.LFP.flip)
                for itrial = 1:size(LFP{ipart}.(markername).trial,2)
                    LFP{ipart}.(markername).trial{itrial} = LFP{ipart}.(markername).trial{itrial} .* -1;
                end
            end
            cfgtemp                           = [];
            cfgtemp.demean                    = cfg.LFP.baseline;
            cfgtemp.baselinewindow            = cfg.LFP.baselinewindow.(markername);
            LFP{ipart}.(markername)           = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));
        end
    end
    for ipart = 1:size(LFP)
        for markername = "Interictal"
            if isfield(LFP{ipart}, markername)
                for itrial = 1:size(LFP{ipart}.(markername).trialinfo,1)
                    LFP{ipart}.(markername).trialinfo.offset(itrial) = LFP{ipart}.(markername).trialinfo.offset(itrial) - LFP{ipart}.(markername).time{itrial}(end) * LFP{ipart}.(markername).fsample;
                    LFP{ipart}.(markername).time{itrial} = LFP{ipart}.(markername).time{itrial} - LFP{ipart}.(markername).time{itrial}(end);
                end
            end
        end
    end
end