function data = removeArtefactedTrials(cfg, data)
%remove artefacted trials,, defined by a value >0 in trialinfo.artefact

if isempty(data)
    return
end

for ipart = 1:size(data,2)
    
    fprintf('For part %d\n', ipart);
    
    for markername = string(fieldnames(data{ipart}))'
        
        fprintf('*** For %s ***\n', markername);
        
        cfg.minbadtime.(markername) = ft_getopt(cfg.minbadtime, char(markername), 0);
        
        if isempty(data{ipart}.(markername))
            continue
        end
        
        type = ft_datatype(data{ipart}.(markername));
       
        sel = data{ipart}.(markername).trialinfo.artefact_length <= cfg.minbadtime.(markername);
        cfgtemp        = [];
        cfgtemp.trials = sel;
        
        if strcmp(type, 'raw')
            data{ipart}.(markername) = ft_selectdata(cfgtemp, data{ipart}.(markername));
            
        elseif strcmp(type, 'spike')
            data{ipart}.(markername) = ft_spike_select_rmfulltrials(cfgtemp, data{ipart}.(markername));
            
        else
            error('datatype %s is not supported', type);
        end
       
    end %markername
end %ipart

end