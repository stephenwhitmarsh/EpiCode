function [SpikeTrialsPSG] = rectifywindownr(cfg,SpikeTrialsPSG)

fprintf('sdf\n');

for ipart = 1 : size(SpikeTrialsPSG,2)
    
    for itemp = 1 : size(SpikeTrialsPSG{1}.trial,2)
        
        temptrials = unique(SpikeTrialsPSG{ipart}.trial{itemp});
        
        for stage = 0 : 5
            
            trialnr = sum(find(
            
        end
        
    end
    
end

