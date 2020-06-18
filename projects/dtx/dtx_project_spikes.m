function dtx_project_spikes(slurm_task_id)

if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    
end

ft_defaults


feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx


config = dtx_setparams_probe_spikes([]);

%% prepare data for spyking circus analysis
% if ismember(slurm_task_id, [6, 7])
%
%     irat = slurm_task_id;
%     [MuseStruct]                     = readMuseMarkers(config{irat}, false);
% 
%     %remove seizures for non-control experiments
%     if strcmp(config{irat}.type, 'dtx')
%         [MuseStruct]                     = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true, true);
%         %save deadfile without removing seizures
%         writeSpykingCircusDeadFile(config{irat},MuseStruct);
%         %remove all seizures
%         [MuseStruct]    = addMuseBAD(config{irat},MuseStruct);
%     end
%
%     %create multifile and dead file, and slurm file
%     writeSpykingCircus(config{irat}, MuseStruct, true, true);
%     %create param and prb file
%     writeSpykingCircusParameters(config{irat})
%
%     return
%
% end


%% analyse spyking circus output

for irat = slurm_task_id
    %FIXME pour enregistrements controle : définir un unique trial entre
    %Baseline_Start et Analysis_End
    
    %align markers and remove wrong seizures
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    [MuseStruct]                    = alignMuseMarkers(config{irat},MuseStruct, false);
    if strcmp(config{irat}.type, 'dtx'), MuseStruct = dtx_remove_wrong_seizure(config{irat}, MuseStruct,false, false); end
    %[mrk_clock mark_synctime] = concatenateMuseMarkers(config{irat},MuseStruct, ipart, 'SlowWave');
    
    %read LFP data
    dat_LFP                         = readLFP(config{irat}, MuseStruct, false, false);
    %for itrial=1:size(dat_LFP{1}{2}.trial,2), plot(dat_LFP{1}{2}.time{itrial}, dat_LFP{1}{2}.trial{itrial}(strcmp(dat_LFP{1}{2}.label, 'E12LFP'),:)+itrial*10000,'k'); end
    [config{irat},dat_LFP]          = dtx_correctDTX2name(config{irat},dat_LFP);
    
    %read spike data
    [SpikeRaw]                      = readSpikeRaw_Phy(config{irat},true,'all');
    %     test = readSpikeRaw_SpykingCircus(config{irat},true,'all');
    if strcmp(config{irat}.type, 'dtx')
        %make trials based on Muse Markers
        [SpikeTrials]                   = readSpikeTrials_MuseMarkers(config{irat}, MuseStruct,SpikeRaw, true, 'all');
        %FIXME rename readSpikeTrials_Muse
    elseif strcmp(config{irat}.type, 'ctrl')
        %make 10 minuts continuous trials on all the data, with same
        [SpikeTrials]                   = readSpikeTrials_continuous(config{irat}, MuseStruct,SpikeRaw, false, 'all'); %FIXME rename readSpikeTrials_Muse
        %FIXME faire le texte d'intro de cette fonction
    else
        error('config{%d}.type ''%s'' is not an available value', irat, config{irat}.type);
    end
    %FIXME pull request ft_spike_maketrials because I added the selection
    %FIXE pull request stephen avec bug si un dossier sans trials
    %of some other fields.      % and for ft_spike_select
    %FIXME rename readSpikeTrials
    
    return %REMOVEME
    clear SpikeRaw
    
%     FIXME FIND A WAY OF IMPROVE THIS FUNCTION - CORRECT BUGS
%     avec DTX2, 6 trials à rejeter
    %remove BAD LFP/Spike trials
    try %REMOVEME
    cfgtemp                         = [];
    cfgtemp                         = config{irat};
    if strcmp(config{irat}.type, 'dtx')
        cfgtemp.keeptrials               = 'yes';
        cfgtemp.LFP.electrodetoplot     = config{irat}.align.channel(1);
        cfgtemp.plotdata                = 'yes';
    else
        cfgtemp.keeptrials               = 'no';
        cfgtemp.LFP.electrodetoplot     = [];
        cfgtemp.plotdata                = 'no';
    end
    
    cfgtemp.method                  = 'remove';
    cfgtemp.markerstart             = 'BAD__START__';
    cfgtemp.markerend               = 'BAD__END__';
    cfgtemp.indexstart              = 0;
    cfgtemp.indexend                = 0;
    cfgtemp.timefrombegin           = 0;
    cfgtemp.timefromend             = 0;
    
    [dat_LFP, ~]                    = removetrials_MuseMarkers(cfgtemp, dat_LFP, MuseStruct, 'all', 'all');
    [SpikeTrials, ~]                = removetrials_MuseMarkers(cfgtemp, SpikeTrials, MuseStruct, 'all', 'all');
    end %try, REMOVEME
    
    %read spike waveforms
    [SpikeWaveforms]                = readSpikeWaveforms(config{irat}, SpikeTrials, true, 'all');
    %FIXME add _1000 for loading precomputed data quickly
    %FIXME : ctrl are set to 1000 : do again with all
    
%     for i_imagesavedir = [1 2] %REMOVEME
        
%         if i_imagesavedir == 2, config{irat}.imagesavedir = fullfile(config{irat}.imagesavedir,'Same_analysis_an_other_time'); end %FIXME REMOVEME
        %stats per unit
        stats                           = spikeratestats_Events_Baseline(config{irat},SpikeTrials,SpikeWaveforms,dat_LFP,'all',true);
        %FIXME voir les try/end dans plot_morpho
        %FIXME : trouver un meilleur moyen de raccourcir les trials trop long
        %avec la cfg input
        %stats = spikeratestats_Events_Baseline(config{irat},[],[],[],[],false);
        %stats = spikeratestats_Events_Baseline(config{irat},SpikeTrials,[],[],'all',true);
        %FIXME : ajouter doublets de PA avec spikestatsOverTime
        %FIXME : pull request Fieldtrip
        
%     end %i_imagesavedir %REMOVEME
    
end %irat
return

%% Pool over rats

%load all spiketrials{irat]
%fill spikeTrials with hand-annotated infos, based on previous plots
%plot all-units morpho, and select in/pn according to that
%plot the rest : freq, cv2, comportement OL, comportement interictal : avec
%infos sua mua, pn in, proche ou non du site d'injection
% une ligne par neurone avec représentation en couleur de l'amplitude de la
% fonction spike density. Pour 1 rat. Pour plusieurs rats, voir si merge
% avec profondeur calculée par l'histo
% corrélation spikerate et LFP


%cross-correlogram pour voir si connexions monosynaptiques. Voir si
%évolue au cours du temps interictal, et pendant SW
%SUA et MUA, sur une grande feuille A2

%Spiky pour voir évolution de la synchronie au cours du temps.
%edge effects, utiliser des trials incluant un padding
%Pour chaque trial
%Uniquement pour les SUA, et avec SUA+MUA

%Représentation spatiale avec une ligne par unit

%% Load precomputed data
for irat = 1:7
    stats{irat} = spikeratestats_Events_Baseline(config{irat},[],[],[],[],false);
%     if strcmp(config{irat}.type, 'dtx'), SpikeTrials{irat} = readSpikeTrials_MuseMarkers(config{irat}, [],[], false, 'all');end
%     if strcmp(config{irat}.type, 'ctrl'), SpikeTrials{irat} = readSpikeTrials_continuous(config{irat}, [],[], false, 'all');end
end

%REMOVEME
% for irat = 1:5
%     for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth,2)
%         if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth{i_unit})
%            stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth{i_unit} = NaN;
%         end
%         if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.peaktrough{i_unit})
%            stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.peaktrough{i_unit} = NaN;
%         end
%         if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.troughpeak{i_unit})
%            stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.troughpeak{i_unit} = NaN;
%         end
%     end
% end

%do it after filtering and keeping only sua
unit_table = readtable('Z:\analyses\lgi1\DTX-PROBE\classification_units.xlsx');
unit_table = table2struct(unit_table);
for irat = 1:7
    rat_idx = strcmp({unit_table.ratID}, config{irat}.prefix(1:end-1))';
    rat_table = unit_table(rat_idx,:);
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        unit_idx = strcmp(split(sprintf('cluster_%d,', rat_table.clusterID), ','), stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit});
        % find each element of the unit
        if sum(unit_idx) == 1
            stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}              = rat_table(unit_idx).group;
            if strcmp(config{irat}.type, 'dtx'), stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).code_slowwave{i_unit} = rat_table(unit_idx).code_slowwave_spikerate; end
            stats{irat}{ipart}.(config{irat}.spike.baselinename).code_spikerate{i_unit}    = rat_table(unit_idx).code_interictal_spikerate; 
            stats{irat}{ipart}.(config{irat}.spike.baselinename).code_cv2{i_unit}          = rat_table(unit_idx).code_interictal_cv2; 
            stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit}            = rat_table(unit_idx).Emax;
        end
    end
end

%replace empty cells by nans
for irat = 1:7
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit})
            stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit} = NaN;
        end
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit} = 'noise';
        end
    end
end

%% IN PN classification
figure;hold;
check_empty = [];
for irat = 1:7
    ipart = 1;
    stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype = cell(1,size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2));
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                x = [stats{irat}{ipart}.(config{irat}.spike.baselinename).template.halfwidth{i_unit}];
                y = [stats{irat}{ipart}.(config{irat}.spike.baselinename).template.peaktrough{i_unit}];
                z = [stats{irat}{ipart}.(config{irat}.spike.baselinename).template.troughpeak{i_unit}];
                %best clustering with halfwidth and throughpeak
                
                %values to separate in and pn determined empirically
%                 if x>2.25*10^-4 && z>6*10^-4 %PN
                if x>3*10^-4 || (x>2*10^-4&& z>4.5*10^-4 ) || z>6*10^-4%PN
                    stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit} = 'pn';
                    if strcmp(config{irat}.type, 'dtx'), plottype = '^b';end
                    if strcmp(config{irat}.type, 'ctrl'), plottype = 'ob';end
                else %IN
                    stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit} = 'in';
                    if strcmp(config{irat}.type, 'dtx'), plottype = '^r';end
                    if strcmp(config{irat}.type, 'ctrl'), plottype = 'or';end
                end
                
                if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
                    scatter(x,z,plottype, 'filled');
                else
                     scatter(x,z,plottype);
                end
            end
        end
    end
end
ax = axis;
leg1 = scatter(1,1,'^k'); %mua
leg2 = scatter(1,1,'^k','filled'); %sua
leg3 = scatter(1,1,'^b', 'filled'); %pn
leg4 = scatter(1,1,'^r', 'filled'); %in
leg5 = scatter(1,1,'^k','filled'); %dtx
leg6 = scatter(1,1,'ok','filled'); %ctrl
legend([leg1 leg2 leg3 leg4 leg5 leg6],'MUA','SUA', 'PN','IN','DTX','CTRL');
xlim([ax(1) ax(2)]); ylim([ax(3) ax(4)]);
xticklabels(xticks.*1000); %convert to ms
yticklabels(yticks.*1000); %convert to ms
set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 10);
xlabel('halfwidth(ms)');
ylabel('peak-trough (ms)');

%% norm morpho selon type
figure;hold;
check_empty = [];
for irat = 1:7
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
%             if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                    clear temp
                    for i=1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).template.values{i_unit},1)
                        temp(i,:) = stats{irat}{ipart}.(config{irat}.spike.baselinename).template.values{i_unit}{i,:}; 
                    end
                    plot(stats{irat}{ipart}.(config{irat}.spike.baselinename).template.time{1}{1}, nanmean(temp,1)/max(nanmean(temp,1)), 'b');
                elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                    clear temp
                    for i=1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).template.values{i_unit},1)
                        temp(i,:) = stats{irat}{ipart}.(config{irat}.spike.baselinename).template.values{i_unit}{i,:}; 
                    end
                    x = stats{irat}{ipart}.(config{irat}.spike.baselinename).template.time{1}{1};
%                     plot(x+x(end)+0.0005, nanmean(temp,1)/max(nanmean(temp,1)), 'r');
                    plot(x, nanmean(temp,1)/max(nanmean(temp,1)), 'r');
                end

            end
        end
    end
end

%% FREQ
figure;hold;
check_empty = [];
for irat = 1:7
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                    plottype = 'ob';
                elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                    plottype = 'or';
                end
                if strcmp(config{irat}.type, 'dtx'),   x = 1; 
                elseif strcmp(config{irat}.type, 'ctrl'), x=2; 
                end
                
                x = x+(rand-0.5)*0.2;
                y = [stats{irat}{ipart}.(config{irat}.spike.baselinename).discharge.meanfreq{i_unit}];
                
                scatter(x,y,plottype, 'filled');
            end
        end
    end
end
xlim([0 3]);

%% CV2
figure;hold;
check_empty = [];
for irat = 1:7
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                    plottype = 'ob';
                elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                    plottype = 'or';
                end
                
                if strcmp (config{irat}.type, 'dtx'),   x = 1; else, x=2; end
                x = x+(rand-0.5)*0.2;
                y = [stats{irat}{ipart}.(config{irat}.spike.baselinename).discharge.meancv2{i_unit}];
                
                scatter(x,y,plottype, 'filled');
            end
        end
    end
end
xlim([0 3]);
ylim([0 2]);

%% behavior during slowwave
figure;hold;
for irat = 1:7
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp (config{irat}.type, 'dtx')
                    x =  stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).code_slowwave{i_unit} + (rand-0.5)*0.2;
                    y = rand;
                    
                    if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                        plottype = 'ob';
                    elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                        plottype = 'or';
                    end
                    scatter(x,y,plottype, 'filled');
                end
            end
        end
    end
end

%% max discharge during slowwave, interictal spikerate and cv2, compared to deep
% à faire en fonction de la vraie profondeur
%Interessant pour le code SW
figure;hold;
for irat = 1:7
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp (config{irat}.type, 'dtx')
                    x = stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit} + (rand-0.5)*0.2;
%                     y = stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).maxfreq.max{i_unit};
                    y = stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).code_slowwave{i_unit} + (rand-0.5)*0.2;
%                     y = stats{irat}{ipart}.(config{irat}.spike.baselinename).code_spikerate{i_unit}+ (rand-0.5)*0.2;
%                     y = stats{irat}{ipart}.(config{irat}.spike.baselinename).code_cv2{i_unit}+ (rand-0.5)*0.2;
                    
                    
                    if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                        plottype = 'ob';
                    elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                        plottype = 'or';
                    end
                    scatter(x,y,plottype, 'filled');
                end
            end
        end
    end
end

%% interical behavior code
%FREQ
figure;hold;
for irat = 1:7
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp (config{irat}.type, 'dtx')
                    x =  stats{irat}{ipart}.(config{irat}.spike.baselinename).code_spikerate{i_unit} + (rand-0.5)*0.2;
                    y = rand;
                    
                    if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                        plottype = 'ob';
                    elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                        plottype = 'or';
                    end
                    scatter(x,y,plottype, 'filled');
                end
            end
        end
    end
end

%CV2
figure;hold;
for irat = 1:7
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp (config{irat}.type, 'dtx')
                    x =  stats{irat}{ipart}.(config{irat}.spike.baselinename).code_cv2{i_unit} + (rand-0.5)*0.2;
                    y = rand;
                    
                    if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                        plottype = 'ob';
                    elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                        plottype = 'or';
                    end
                    scatter(x,y,plottype, 'filled');
                end
            end
        end
    end
end

%% colorline for each neuron
figure;hold;
emax_list = nan;
sdf_avg = nan(1,length(stats{1}{ipart}.(config{1}.spike.eventsname{1}).sdf.avg));
sdf_time = stats{1}{ipart}.(config{1}.spike.eventsname{1}).sdf.time;
celltype = nan;
group = nan;
for irat = 1:5
    %for only sua
%     idx = ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'mua') & ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'noise');
    %for sua and mua
    idx = ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'noise');
    emax_list   = [emax_list, stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{idx}];
    sdf_avg     = vertcat(sdf_avg, stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).sdf.avg(idx,:));
    celltype    = [celltype, stats{irat}{ipart}.(config{1}.spike.baselinename).celltype(idx)];
    group       = [group, stats{irat}{ipart}.(config{1}.spike.baselinename).group(idx)];
end
[emax_list_sorted, deep_idx] = sort(emax_list, 'descend');

deep_idx = deep_idx(~isnan(emax_list_sorted));
emax_list_sorted = emax_list_sorted(~isnan(emax_list_sorted));
sdf_avg = sdf_avg(deep_idx,:);
celltype = celltype(deep_idx);
group = group(deep_idx);

for i=1:size(sdf_avg,1)
    sdf_avg_filt(i,:) = imgaussfilt(sdf_avg(i,:),1); %smooth the spike density function values
    sdf_filt_blcorrect(i,:) = (sdf_avg_filt(i,:)) - nanmean((sdf_avg_filt(i,1:20))); %correct baseline
end

% imagesc(log10(abs(sdf_avg_filt)));
% imagesc(log10(abs(sdf_filt_blcorrect)));
% c = caxis;
% caxis([0 c(2)]);
% clb = colorbar;
% axis tight
% xticklabels(sdf_time(xticks+1));
% set(gca,'TickDir','out');
% func_temp = @(v) sprintf('10\\^%s',v);
% clb.TickLabels = cellfun(func_temp,clb.TickLabels,'UniformOutput',false);
% 
% % plot baseline patch
% x = [0 20 20 0];
% ax = axis;
% y = [ax(3) ax(3) ax(4) ax(4)];
% patch('XData',x,'YData',y,'facecolor',[1 1 1],'edgecolor','none','facealpha',0.2);
% 
% % plot spike infos
% C = linspecer(sum(~isnan(unique(emax_list_sorted))),'qualitative');
% for ineuron = 1:size(celltype,2)
% %     if strcmp(celltype{ineuron},'pn'), plottype='^w';end
% %     if strcmp(celltype{ineuron},'in'), plottype='or';end
%     idx=find(unique(emax_list_sorted)==emax_list_sorted(ineuron));
%     scatter(size(sdf_avg,2)-1,ineuron,'s','filled','MarkerEdgeColor', C(idx,:), 'MarkerFaceColor', C(idx,:));
% end

% plot pn and in separately
figure;
for i_celltype = ["pn", "in"]
    
    if strcmp(i_celltype, "pn"), subplot(2,1,1);hold; else, subplot(2,1,2);hold; end
    imagesc(log10(abs(sdf_filt_blcorrect(strcmp(celltype,i_celltype),:))));
    c = caxis;
    caxis([0 c(2)]);
    clb = colorbar;
    axis tight
    xticklabels(sdf_time(xticks+1));
    func_temp = @(v) sprintf('10\\^%s',v);
    clb.TickLabels = cellfun(func_temp,clb.TickLabels,'UniformOutput',false);
    % plot baseline patch and zero line
    x = [0 20 20 0];
    ax = axis;
    y = [ax(3) ax(3) ax(4) ax(4)];
    patch('XData',x,'YData',y,'facecolor',[1 1 1],'edgecolor','none','facealpha',0.2);
    plot([100 100], [ax(3:4)], '--r');
    %plot electrode infos
    C = linspecer(sum(~isnan(unique(emax_list_sorted(strcmp(celltype,i_celltype))))),'qualitative');
    emax_pn = emax_list_sorted(strcmp(celltype,i_celltype));
    group_temp = group(strcmp(celltype,i_celltype));
    for ineuron = 1:size(celltype(strcmp(celltype,i_celltype)),2)
        chan_idx=find(unique(emax_pn)==emax_pn(ineuron));
        scatter(size(sdf_avg,2)+2,ineuron,'s','filled','MarkerEdgeColor', C(chan_idx,:), 'MarkerFaceColor', C(chan_idx,:));
        if ~contains(group_temp(ineuron), 'mua')%noise already removed
            scatter(size(sdf_avg,2)-1,ineuron,'<r','filled');
        end
    end
    set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 10);
    xlabel('Time to begin of SlowWave (s)');
end

%% plot of cv2 with bursts, without bursts and freq over time : avg for each neuron

%normalized at t=-60
for method = ["cv2_withoutbursts","freq"]%"cv2";"cv2_withoutbursts";
    figure;hold;
    for irat = 1:7
        if strcmp(config{irat}.type, 'dtx')
            ipart = 1;
            for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
                if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
                    check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
                else
                    if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
%                         if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                        
                        x = stats{irat}{ipart}.(config{irat}.spike.baselinename).stats_over_time.(method){i_unit}.time;
                        y = stats{irat}{ipart}.(config{irat}.spike.baselinename).stats_over_time.(method){i_unit}.avg;
                        
                        if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                            color = 'b';
                        elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                            color = 'r';
                        end
                        y = y./y(x==-60); %normalize
                        if max(y(x>-60)) <10 %avoir outliers
                            plot(x,y,'Color', color);
                        end
                    end
                end
            end
        end
    end
end
xlim([-60 0]);
ylim([0.5 1.8]);
ylim([0 3]);

%
% figure;hold;
% for iunit = 1:23
%     for itrial = 1:17
%         plot(stats{1}.Interictal.stats_over_time.freq{iunit}.time, stats{1}.Interictal.stats_over_time.freq{iunit}.values(itrial,:), 'k');
%         plot(stats{1}.Interictal.stats_over_time.freq{iunit}.time, movmean(stats{1}.Interictal.stats_over_time.freq{iunit}.values(itrial,:),4), 'r');
%     end
% end
% axis tight

%plot avg of freq over time
% figure;hold;
% for i_unit = 1:size(stats{irat}{1}.Interictal.stats_over_time.freq,2)
%     plot(stats{irat}{1}.Interictal.stats_over_time.freq{i_unit}.time,stats{irat}{1}.Interictal.stats_over_time.freq{i_unit}.avg);
% end
figure;hold;
%FREQ TIMENORM SERA REMPLACE PAR CV2WITHOUTBURSTS
for i_unit = 1:size(stats{irat}{1}.Interictal.stats_over_time.freq,2)
    plot(stats{irat}{1}.Interictal.stats_over_time.freq_timenorm{i_unit}.time,stats{irat}{1}.Interictal.stats_over_time.freq_timenorm{i_unit}.avg);
end
end
%
