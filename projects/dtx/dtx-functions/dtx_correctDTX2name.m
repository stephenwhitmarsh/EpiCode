function [cfg,data] = dtx_correctDTX2name(cfg,data)
% Error during acquisition for rat DTX2 : need to rename channels ECoGS1 in
% ECoGM1G and ECoGM1 in ECoGM1D.
% DTX2, after the data are read.
% For all other rats than DTX2, do nothing.

if strcmp(cfg.prefix,'DTX2-')
    fprintf('Rat is DTX2, correct cfg and LFP channel names\n');
    
    for ipart = 1:length(data)
        for imarker = 1:length(data{ipart})
            
            if strcmp(cfg.LFP.electrodetoplot{imarker},'ECoGS1')
                cfg.LFP.electrodetoplot{imarker} = 'ECoGM1G';
            end
            
            for ilabel = 1:length(data{ipart}{imarker}.label)
                
                if strcmp(data{ipart}{imarker}.label{ilabel}, 'ECoGM1')
                    data{ipart}{imarker}.label{ilabel} = 'ECoGM1D';
                end
                
                if strcmp(data{ipart}{imarker}.label{ilabel}, 'ECoGS1')
                    data{ipart}{imarker}.label{ilabel} = 'ECoGM1G';
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
end

