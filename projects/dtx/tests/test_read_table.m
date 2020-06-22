test = readmatrix('C:\Users\paul.baudin\Documents\Analyse Neuralynx\Analyse-marqueurs.xlsx');
test2 = readtable('Y:\rodents\raw\DTX-EEGRodents-Brainvision\DTX-EEGEMG-12\DTX-EEGEMG-12-Timings.csv');


% test = tdfread('Z:\analyses\lgi1\DTX-PROBE\Classification_units.xlsx');
unit_table = readtable('Z:\analyses\lgi1\DTX-PROBE\Classification_units.xlsx');
unit_table = table2struct(unit_table);

% for irat
rat_idx = strcmp({unit_table.ratID}, config{irat}.prefix(1:end-1))';
rat_table = unit_table(rat_idx,:);
% for i_unit = 3
unit_idx = strcmp(split(sprintf('cluster_%d,', rat_table.clusterID), ','), SpikeTrials{1}{3}.label{i_unit});
% find each element of the unit
group{i_unit} = rat_table(unit_idx).group;
% save SpikeTrials annotated