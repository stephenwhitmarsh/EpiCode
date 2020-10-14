function dtx_LFP_stats(slurm_task_id)


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


%% batch for peak or begin


config_origin = config;


for ipatient = slurm_task_id
    
        %% TEMPORAIRE PAUL
%     config{ipatient}.imagesavedir = fullfile(config{ipatient}.imagesavedir,config{ipatient}.name{1},'test_CSD'); %FIXME remove test_CSD
    config{ipatient}.datasavedir = fullfile(config{ipatient}.datasavedir, config{ipatient}.name{1});
    
    
    %% load and align
    
    [MuseStruct]                    = readMuseMarkers(config{ipatient}, true);
%     [MuseStruct]                    = alignMuseMarkers_peak(config{ipatient},MuseStruct, true);
    [MuseStruct]                    = alignMuseMarkers_begin(config{ipatient},MuseStruct,true);
    [MuseStruct]                    = alignMuseMarkers_EMG(config{ipatient},MuseStruct,true);
    [dat_LFP]                       = readLFP(config{ipatient}, MuseStruct, true, true); %dat_LFP{ipart}{imarker}
    if strcmp(config{ipatient}.prefix, 'DTX2-')  %correct name of DTX2 (error during acquisition)
        [config{ipatient},dat_LFP]      = dtx_correctDTX2name(config{ipatient},dat_LFP);
        config_origin = config;
    end
    

    %% flip data
    if config{ipatient}.LFP.flip == true
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
    
    %% correct baseline
    for ipart = 1:size(dat_LFP,2)
        for imarker = 1:size(dat_LFP{ipart},2)
            if ~isempty(dat_LFP{ipart}{imarker})
                cfgtemp = [];
                cfgtemp.demean = 'yes';
                cfgtemp.baselinewindow = config{ipatient}.align.toibaseline{1};
                dat_LFP{ipart}{imarker} = ft_preprocessing(cfgtemp, dat_LFP{ipart}{imarker});
            end
        end
    end
    
    %% remove trials with artefacts
    %remove BAD LFP trials
    cfgtemp                         = [];
    cfgtemp                         = config{ipatient};
    cfgtemp.LFP.electrodetoplot     = [];
    cfgtemp.method                  = 'remove';
    cfgtemp.markerstart             = 'BAD__START__';
    cfgtemp.markerend               = 'BAD__END__';
    cfgtemp.indexstart              = 0;
    cfgtemp.indexend                = 0;
    cfgtemp.timefrombegin           = 0;
    cfgtemp.timefromend             = 0;
    cfgtemp.plotdata                = 'no';
    [dat_LFP, MuseStruct]           = removetrials_MuseMarkers(cfgtemp, dat_LFP, MuseStruct, 'all', 'all');
    
    %% create new dir for those analysis
    config{ipatient}.imagesavedir = fullfile(config{ipatient}.imagesavedir,'..', sprintf('measurement_%s',config{ipatient}.LFP.name{1}));
    if ~isfolder(config{ipatient}.imagesavedir)
        ft_notice('creating directory %s', config{ipatient}.imagesavedir);
        mkdir(config{ipatient}.imagesavedir);
    end
    
    %% TEMP - compute CSD
%     neigbours found here :https://github.com/fieldtrip/fieldtrip/blob/master/template/neighbours/elec1020_neighb.mat
%     .elec was already in the template folder
load('elec1020_neighb.mat','neighbours');
    for ipart = 1:size(dat_LFP,2)
        for imarker = 1:size(dat_LFP{ipart},2)
            if ~isempty(dat_LFP{ipart}{imarker})
                cfgtemp = [];
                cfgtemp.method       = 'spline';
                cfgtemp.elec         = 'standard_1020.elc';
                cfgtemp.neighbours       = neighbours;
                dat_LFP{ipart}{imarker} = ft_scalpcurrentdensity(cfgtemp, dat_LFP{ipart}{imarker});
            end
        end
    end
    
%     ipart = 3; imarker = 1;
    %% plots
    ipart = length(dat_LFP);
    
    %% plot morpho SW sans EMG
%     CORRIGER NOM PATIENTS POUR MERGE, plus fait dans la fonction
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
            cfgtemp.measurehalfwidth     = 'yes';
            cfgtemp.halfwidthmethod     = 'bl'; 
            cfgtemp.measurepeaktrough    = 'yes';
            cfgtemp.name                = config{ipatient}.LFP.name{imarker};
            cfgtemp.saveplot            = 'yes';
            cfgtemp.imagesavedir        = config{ipatient}.imagesavedir;
            cfgtemp.prefix              = config{ipatient}.prefix;
            plot_morpho(cfgtemp,dat_LFP{ipart}{imarker});
        end
    end
    
    
    
    %% morpho SW, EMG, and comparison
    % choisir si alignement d�but ou pic, � voir avec les figures d'ajd
    % ajouter savefig aux arguments input quand script fini de tourner
    %reset save path
    try
    config = config_origin;
    config{ipatient}.imagesavedir = fullfile(config{ipatient}.imagesavedir,'..',sprintf('measurement_%s',config{ipatient}.LFP.name{1}), 'morpho_sw_emg_delay');
    
    % check if images directory exists, if not create
    if ~isfolder(config{ipatient}.imagesavedir)
        ft_notice('creating directory %s', config{ipatient}.imagesavedir);
        mkdir(config{ipatient}.imagesavedir);
    end
    for imarker = 1:length(config{ipatient}.LFP.name)
        if ~isempty(dat_LFP{ipart}{imarker})
            if any(strfind(config{ipatient}.LFP.name{imarker}, 'EMG'))
                if ~strcmp(config{ipatient}.LFP.emg{imarker}, 'no')
                    dtx_plot_comparison_eeg_emg(config{ipatient},dat_LFP,ipart,imarker,true,config{ipatient}.LFP.motorcortex{imarker},config{ipatient}.LFP.emg{imarker});
                end
            end
        end
    end
    end
    
    %% SlowWave topography & statistics
    config = config_origin;
    config{ipatient}.imagesavedir = fullfile(config{ipatient}.imagesavedir,'..',sprintf('measurement_%s',config{ipatient}.LFP.name{1}), 'topo_sw_stats');
    % check if images directory exists, if not create
    if ~isfolder(config{ipatient}.imagesavedir)
        ft_notice('creating directory %s', config{ipatient}.imagesavedir);
        mkdir(config{ipatient}.imagesavedir);
    end
    
    for imarker = 1:length(config{ipatient}.LFP.name)
        if ~isempty(dat_LFP{ipart}{imarker})
            dtx_stats_SlowWaveTopography(config{ipatient},isPatient,dat_LFP,ipart,imarker,true);
        end
    end
    
    
    
end %ipatient

if ipatient == length(config)
    
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





