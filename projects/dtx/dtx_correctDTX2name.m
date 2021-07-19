function [cfg,data] = dtx_correctDTX2name(cfg,data)
% Error during acquisition for rat DTX2 : need to rename channels ECoGS1 to
% ECoGM1G and ECoGM1 to ECoGM1D.

if strcmp(cfg.prefix,'DTX2-')
    fprintf('Rat DTX2, correct cfg and LFP channel names\n');
    for ipart = 1:length(data)
        for markername = string(fieldnames(data{ipart}))'
            if isfield(cfg.LFP.electrodetoplot, markername)
                if strcmp(cfg.LFP.electrodetoplot.(markername),'ECoGS1')
                    cfg.LFP.electrodetoplot.(markername) = 'ECoGM1G';
                end
            end
            for ilabel = 1:length(data{ipart}.(markername).label)
                if strcmp(data{ipart}.(markername).label{ilabel}, 'ECoGM1')
                    data{ipart}.(markername).label{ilabel} = 'ECoGM1D';
                end
                if strcmp(data{ipart}.(markername).label{ilabel}, 'ECoGS1')
                    data{ipart}.(markername).label{ilabel} = 'ECoGM1G';
                end
                if strcmp(cfg.labels.macro{ilabel}, 'ECoGM1')
                    cfg.labels.macro{ilabel} = 'ECoGM1D';
                end
                if strcmp(cfg.labels.macro{ilabel}, 'ECoGS1')
                    cfg.labels.macro{ilabel} = 'ECoGM1G';
                end
                if strcmp(cfg.LFP.channel{ilabel}, 'ECoGS1')
                    cfg.LFP.channel{ilabel} = 'ECoGM1G';
                end
                if strcmp(cfg.LFP.channel{ilabel}, 'ECoGM1')
                    cfg.LFP.channel{ilabel} = 'ECoGM1D';
                end
            end
        end
    end
end


