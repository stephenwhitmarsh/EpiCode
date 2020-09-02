function dtx_project_patientslgi1(slurm_task_id)


%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
end


ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

config = dtx_setparams_eegvideo;
% config_origin = config;

for ipatient = [1 3]
    
    [MuseStruct]                    = readMuseMarkers(config{ipatient}, true);
    %take the last part
    ipart = size(MuseStruct,2); 
    
    %% compute eeg-emg delays and emg duration
    slowwave_begin      = concatenateMuseMarkers(MuseStruct,ipart,'SlowWave_begin');
    emg_begin           = concatenateMuseMarkers(MuseStruct,ipart,'SlowWave_EMG__START__');
    emg_end             = concatenateMuseMarkers(MuseStruct,ipart,'SlowWave_EMG__END__');
    
    if isempty(slowwave_begin.synctime)
        delays{ipatient} = NaN;
        emg_duration{ipatient} = NaN;
        continue
    end
    
    delays{ipatient}        = emg_begin.synctime - slowwave_begin.synctime;
    emg_duration{ipatient}  = emg_end.synctime - emg_begin.synctime;
    %     [bins, edges]       = histcounts(delays{ipatient},'BinWidth',0.01);
    %     bins_centers        = (edges(1:end-1)+edges(2:end))/2; %limiteinf + limitesup / 2
    %     bar(bins_centers,bins);
    %
end

%eeg emg delay
figure;hold
for ipatient = [1 3]
    scatter(rand(size(delays{ipatient}))*0.2+ipatient-0.1, delays{ipatient}, '.', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    errorbar(ipatient, mean(delays{ipatient}), std(delays{ipatient}),'--rx');
    %     errorbar(ipatient, mean(delays{ipatient}), mean(delays{ipatient})/sqrt(size(delays{ipatient},2)),'--rx'); %errbar : sem
end
xlim([0 4]);
setfig();
ylabel('eeg-emg delay (s)');

%emg duration
figure;hold
for ipatient = [1 3]
    scatter(rand(size(emg_duration{ipatient}))*0.2+ipatient-0.1, emg_duration{ipatient}, '.', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    errorbar(ipatient, mean(emg_duration{ipatient}), std(emg_duration{ipatient}),'--rx');
    %     errorbar(ipatient, mean(emg_duration{ipatient}), mean(emg_duration{ipatient})/sqrt(size(emg_duration{ipatient},2)),'--rx'); %errbar : sem
end
xlim([0 4]);
setfig();
ylabel('emg duration (s)');


    
    %% STOP ICI
    [MuseStruct]                    = alignMuseMarkers(config{ipatient},MuseStruct, false);
%     [MuseStruct]                    = alignMuseMarkersXcorr(config{ipatient},MuseStruct, true);
    [dat_LFP]                       = readLFP(config{ipatient}, MuseStruct, false); %dat_LFP{ipart}{imarker}
    
    % flip data
    if ft_getopt(config{ipatient}.LFP, 'flip', false) == true
        %inverse the data according to standard clinical visualisation
        for ipart = 1:size(dat_LFP,2)
            for imarker = 1:size(dat_LFP{ipart},2)
                if ~isempty(dat_LFP{ipart}{imarker})
                    for itrial = 1:size(dat_LFP{ipart}{imarker}.trial,2)
                        dat_LFP{ipart}{imarker}.trial{itrial} = -dat_LFP{ipart}{imarker}.trial{itrial};
                    end
                end
            end
        end
    end
    
    % correct baseline
    for ipart = 1:size(dat_LFP,2)
        for imarker = 1:size(dat_LFP{ipart},2)
            if ~isempty(dat_LFP{ipart}{imarker})
                cfgtemp                 = [];
                cfgtemp.demean          = ft_getopt(config{ipatient}.LFP, 'baseline', 'no');
                cfgtemp.baselinewindow  = ft_getopt(config{ipatient}.LFP, 'baselinewindow', []);
                dat_LFP{ipart}{imarker} = ft_preprocessing(cfgtemp, dat_LFP{ipart}{imarker});
            end
        end
    end
    
    % remove trials which intersest BAD markers
    
    [dat_LFP, ~]           = removetrials_MuseMarkers(config{ipatient}, dat_LFP, MuseStruct);
    
    % Compute CSD
    %     neigbours found here :https://github.com/fieldtrip/fieldtrip/blob/master/template/neighbours/elec1020_neighb.mat
    %     .elec was already in the template folder
    load('elec1020_neighb.mat','neighbours');
    for ipart = 1:size(dat_LFP,2)
        for imarker = 1:size(dat_LFP{ipart},2)
            if ~isempty(dat_LFP{ipart}{imarker})
                cfgtemp                     = [];
                cfgtemp.method              = 'spline';
                cfgtemp.elec                = ft_read_sens('standard_1020.elc');
                cfgtemp.neighbours          = neighbours;
                dat_CSD{ipart}{imarker}     = ft_scalpcurrentdensity(cfgtemp, dat_LFP{ipart}{imarker});
            end
        end
    end
    
    %% Topography
    % Pas de stats pour le moment
    config{ipatient}.topoplot.suffix = 'EEG';
    dtx_plot_SlowWaveTopography(config{ipatient}, dat_LFP);
    dtx_plot_SlowWaveTopographyTimecouse(config{ipatient}, dat_LFP);
    
    config{ipatient}.topoplot.suffix = 'CSD';
    dtx_plot_SlowWaveTopography(config{ipatient}, dat_CSD);
    dtx_plot_SlowWaveTopographyTimecouse(config{ipatient}, dat_CSD);
    
    
    %% hw, amplitude, for both EEG and CSD
    ipart = 1;
    cfgmorpho = config{ipatient};
    
    
    %for those patients : one marker is SW_R, and the other is SW_L
    for imarker = 1:size(dat_LFP{ipart},2)
        
        %select channel of interest
        cfgmorpho.morpho.channame = config{ipatient}.morpho.channame{imarker};
        
        %repeat the same analysis for eeg and csd
        i_analysis=1; %index for suffix
        analysis_type = {'eeg', 'csd'}; 
        for data = [dat_LFP, dat_CSD]
            
            fig = figure;hold;
            %compute data and plot results, for each trial
            for itrial = 1:size(data{ipart}{imarker}.trial,2)
                cfgtemp = [];
                cfgtemp.trials = itrial;
                data_temp = ft_selectdata(cfgtemp, data{ipart}{imarker});
                
                [stats{imarker}.(analysis_type{i_analysis}).halfwidth(itrial), ~, ~,...
                    stats{imarker}.(analysis_type{i_analysis}).amplitude(itrial)] = ...
                    plot_morpho(cfgmorpho, data_temp);
            end
            %remove text to make the figure readable
            delete(findall(gcf,'type','text'));
            
            %save fig
            set(gca,'Fontsize',15);
            if ~(exist(config{ipatient}.imagesavedir)==7)
                mkdir(config{ipatient}.imagesavedir);
                fprintf('Create forlder %s',config{ipatient}.imagesavedir);
            end
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            
            fname = [config{ipatient}.prefix,'p',num2str(ipart),'-Morpho-',config{ipatient}.name{imarker},'_',cfgmorpho.morpho.channame,'_scale',strrep(num2str(config{ipatient}.morpho.toiplot),'  ','_')];
            print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,[fname,'_',analysis_type{i_analysis},'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,[fname,'_',analysis_type{i_analysis},'.png']),'-r600');
            close all
            
            i_analysis=i_analysis+1; %index for suffix
        end
    end
    
    fname = fullfile(config{ipatient}.datasavedir, [config{ipatient}.prefix, 'stats_EEG.mat']);
    save(fname, stats, '-v7.3');
    return %STOP HERE FOR NOW
    %% count markers : time between seizures, eeg/emg delays, emg duration
    
    %% Propagation : plot canaux positifs, ou tous les canaux, normaliser, et compter délai à la main.
    
    %% scatter plot stats pooled entre les patients
    % hw, amplitudes, temps entre 2 crises, délai EEG EMG, EMG
    % duration. Faire une fonction pour scatter
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% OLD %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    %% plot morpho SW sans EMG
    %reset save path
    config = config_origin;
    config{ipatient}.imagesavedir = fullfile(config{ipatient}.imagesavedir,'..',sprintf('measurement_%s',config{ipatient}.LFP.name{1}), 'morpho_sw_maxchan');
    % check if images directory exists, if not create
    if ~isfolder(config{ipatient}.imagesavedir)
        ft_notice('creating directory %s', config{ipatient}.imagesavedir);
        mkdir(config{ipatient}.imagesavedir);
    end
    for imarker = 1:length(config{ipatient}.LFP.name)
        if ~isempty(dat_LFP{ipart}{imarker})
            %             toi.toiplot = [-2 2];
            %             toi.toibl = cfg.align.toibaseline{imarker};
            %             toi.toiac = [-1 1];
            %             plot_morpho(config{ipatient},dat_LFP,ipart,imarker,config{ipatient}.LFP.electrodetoplot{imarker}, toi,true, false,true);
            cfgtemp                     = [];
            cfgtemp.channame            = config{ipatient}.LFP.electrodetoplot{imarker};
            cfgtemp.plotstd             = 'yes';
            cfgtemp.removeoutliers      = 'no';
            cfgtemp.toiplot             = [-2 2];
            cfgtemp.toibl               = config{ipatient}.align.toibaseline{imarker};
            cfgtemp.toiac               = [-1 1];
            cfgtemp.mesurehalfwidth     = 'yes';
            cfgtemp.halfwidthmethod     = 'bl';
            cfgtemp.mesurepeaktrough    = 'yes';
            cfgtemp.name                = config{ipatient}.LFP.name{imarker};
            cfgtemp.saveplot            = 'yes';
            cfgtemp.imagesavedir        = config{ipatient}.imagesavedir;
            cfgtemp.prefix              = config{ipatient}.prefix;
            plot_morpho(cfgtemp,dat_LFP{ipart}{imarker});
        end
    end
    
    
      
    %% seizures infos
    load(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos'); %previously computed
    
    %reset save path
    config = config_origin;
    config{1}.imagesavedir = fullfile(config{1}.imagesavedir,'..',sprintf('measurement_%s',config{1}.LFP.name{1}), 'seizures_infos');
    
    % check if images directory exists, if not create
    if ~isfolder(config{1}.imagesavedir)
        ft_notice('creating directory %s', config{1}.imagesavedir);
        mkdir(config{1}.imagesavedir);
    end
    
    dtx_stats_seizures_infos(config{1},Seizure_Infos);
    
end

%svg stats en output

end





