addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML
CEDS64LoadLib('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML');

%read CED markers
cfgtemp                     = [];
cfgtemp.prefix              = 'Intra_All_SlowWaves';
cfgtemp.rawdir              = 'Z:\analyses\lgi1\DTX-INTRA\Crises-spont-intra';          
cfgtemp.directorylist{1}    = {'MergeCrises'};
cfgtemp.datasavedir         = 'Z:\analyses\lgi1\DTX-INTRA\data_scripts';
CEDStruct = readCEDmarkers(cfgtemp, false);

%convert into a fieldtrip Structure
spikedata.label = fieldnames(CEDStruct{1}{1}.markers)';
spikedata.label = spikedata.label(contains(spikedata.label, 'PASW'));
spikedata.label = sort(spikedata.label);
for ichan = 1:size(spikedata.label,2)
    spikedata.time{ichan} = CEDStruct{1}{1}.markers.(spikedata.label{ichan}).synctime - CEDStruct{1}{1}.markers.SlowWave.synctime(ichan);
    spikedata.trial{ichan} = ones(size(spikedata.time{ichan}));
end
spikedata.trialtime = [-60, 60];

cfgtemp                         = [];
cfgtemp.fsample                 = 1000;   % 
cfgtemp.keeptrials              = 'yes';
cfgtemp.timwin                  = [-0.1 0.1];
cfgtemp.latency                  = [-2 2];
[sdfavg, sdfdata]               = ft_spikedensity(cfgtemp,spikedata);


%Detection automatique ne fonctionne pas car pas assez de trials, pic pas
%clair
% for i_unit = 1:8
%     [peak_sdf_intra{i_unit}, loc_sdf_intra{i_unit}] = findpeaks(sdfavg.avg(i_unit,:), sdfavg.time,'NPeaks',1,'SortStr','descend'); %Npeaks : max nr of peaks/ SortStr : peak sorting : descend = from largest to smallest
%     if isempty(peak_sdf_intra{i_unit}) 
%         peak_sdf_intra{i_unit}  = 0; 
%         loc_sdf_intra{i_unit}   = NaN; 
%     end
%     %                 plot(sdf.time, sdf.avg(i_unit,:), 'k');
%     %                 scatter(loc(i_unit), peak(i_unit), 'xr');
% end
% 
figure;hold;
plot(sdfavg.time, sdfavg.avg);
%mesure des pics à la main sur la figure sdf :
max_ol_intra = [17.51, 15.22, 40.71, 107.2, 91.79, 0,0,0,0];
% scatter([loc_sdf_intra{:}], [peak_sdf_intra{:}], 'xr');

%% compare intra extra : max décharge pendant SW
%peak pour extraellulaire : voir test_plot_sdf.m
figure;hold;
for irat = 1:5
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        if ~contains(stats{irat}{ipart}.Interictal.group{i_unit}, 'noise')
            if strcmp(stats{irat}{ipart}.Interictal.celltype{i_unit}, 'pn')
                scattertype = '^k';
            else
                scattertype = 'ok';
            end
            if  ~contains(stats{irat}{ipart}.Interictal.group{i_unit}, 'mua')
                scatter(1+rand*0.3, peak{irat}(i_unit), scattertype,'filled');
            else
                scatter(1+rand*0.3, peak{irat}(i_unit), scattertype);
            end
        end
    end
end
scatter(2+rand(size(max_ol_intra))*0.3, max_ol_intra, '^k','filled');
xlim([0 3]);
ylim([0 250]);

%essayer le plot avec sdf en image
for i=1:size(sdfavg.avg,1)
%     sdf_avg_filt(i,:) = imgaussfilt(sdf_avg(i,:),1); %smooth the spike density function values
    sdf_blcorrect(i,:) = (sdfavg.avg(i,:)) - nanmean((sdfavg.avg(i,sdfavg.time>-1.9 & sdfavg.time<-1.4))); %correct baseline
%     sdf_norm(i,:) = normalize(sdf_avg(i,:),'zscore');
end

imagesc(log10(abs(sdf_blcorrect)));
    