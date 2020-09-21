function units_infos = dtx_read_unit_table(cfg, spikedata)

% SpikeTrials is used to get the labels of the units

unit_table = readtable(cfg.unit_table);
unit_table = table2struct(unit_table);

rat_idx = strcmp({unit_table.ratID}, cfg.prefix(1:end-1))';
rat_table = unit_table(rat_idx,:);
for i_unit = 1:size(spikedata.label, 2)
    unit_idx = strcmp(split(sprintf('cluster_%d,', rat_table.clusterID), ','), spikedata.label{i_unit});
    % find each element of the unit
    if sum(unit_idx) == 1
        units_infos.label{i_unit}             = spikedata.label{i_unit};
        units_infos.group{i_unit}             = rat_table(unit_idx).group;
        if strcmp(cfg.type, 'dtx')
            units_infos.code_slowwave{i_unit} = rat_table(unit_idx).code_slowwave_spikerate; 
        end
        units_infos.code_spikerate{i_unit}    = rat_table(unit_idx).code_interictal_spikerate;
        units_infos.code_cv2{i_unit}          = rat_table(unit_idx).code_interictal_cv2;
        units_infos.maxchan{i_unit}           = rat_table(unit_idx).Emax;
    end
end


%replace empty cells by nans
for i_unit = 1:size(spikedata.label, 2)
    if isempty(units_infos.maxchan{i_unit})
        units_infos.maxchan{i_unit} = NaN;
    end
    if isempty(units_infos.group{i_unit})
        units_infos.group{i_unit} = 'noise';
    end
end

