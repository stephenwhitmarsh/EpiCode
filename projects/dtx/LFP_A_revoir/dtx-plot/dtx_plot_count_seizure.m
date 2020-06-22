function Seizure_Infos = dtx_plot_count_seizure(cfg, MuseStruct, ipart, plotdata)

%Descriptive statistics of seizures from Muse markers
%BaselineStart = Analysis Start
%Modifs pour avake : pas d'injection. Ou alors la rentrer dans cfg ?
%pour appliquer à tous : rentrer cfg entier. préciser ipatient/irat dans
%argument à part

%rename prefix in case of "merge" data
if isfield(cfg, 'merge')
    if cfg.merge == true
        if ipart > 1 && ipart == length(cfg.directorylist) %last part = merge (except if only one part, nothing to merge)
            cfg.prefix = [cfg.prefix, 'MERGED-'];
        else
            cfg.prefix = [cfg.prefix, cfg.directorylist{ipart}{:}, '-'];
        end
    end
end

%get all marker times
[SlowWave_R_clock, ~] = concatenateMuseMarkers(cfg,MuseStruct,ipart, 'SlowWave_R');
[SlowWave_L_clock, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'SlowWave_L');
[SlowWave_R_EMG_clock, ~] = concatenateMuseMarkers(cfg,MuseStruct,ipart, 'SlowWave_R_EMG');
[SlowWave_L_EMG_clock, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'SlowWave_L_EMG');
[SlowWave_clock, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'SlowWave');
[Crise_Start_clock, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'Crise_Start');
[Crise_End_clock, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'Crise_End');
[Baseline_Start_clock, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'Baseline_Start');
[Analysis_End_clock, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'Analysis_End');
[Injection_clock, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'Injection');
[BAD, ~] = concatenateMuseMarkers(cfg, MuseStruct,ipart,'BAD__START__');

%recognize rat or patient data
isRat = 0;
isPatient = 0;
if ~isempty(SlowWave_R_clock) || ~isempty(SlowWave_L_clock)
    isPatient = 1;
elseif ~isempty(SlowWave_clock)
    isRat = 1;
end

%load header
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);
for idir = 1:length(MuseStruct{ipart})
    if isNeuralynx
        filelist = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.ncs'));
        hdr{idir} = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},filelist(1).name));
    elseif isMicromed
        hdr{idir} = ft_read_header(fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']));
    elseif isBrainvision
        hdr{idir} = ft_read_header(fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']));
    end
end

%avoir error in crise annotation
if isRat
    if size(Crise_Start_clock,2) ~= size(Crise_End_clock,2)
        error('%s : there are not the same number of Crise_Start than Crise_End', cfg.prefix(1:end-1));
    elseif isempty(Baseline_Start_clock) || isempty(Analysis_End_clock)
        error('%s : miss Baseline_Start or Analysis_End', cfg.prefix(1:end-1));
    end
end


%Total time
if isRat
    BaselineDuration = cfg.injectiontime-Baseline_Start_clock{1,1};
    RecordDuration = Analysis_End_clock{1,1}-Baseline_Start_clock{1,1};
    PostInjDuration = RecordDuration - BaselineDuration;
    length_previous_dirs(1) = seconds(0);
    for idir = 1:length(MuseStruct{ipart})
        length_previous_dirs(idir+1) = length_previous_dirs(idir) + seconds(hdr{idir}.nSamples/hdr{idir}.Fs);
    end
    
elseif isPatient
    length_previous_dirs = seconds(0);
    RecordDuration = seconds(0);
    for idir = 1:length(MuseStruct{ipart})
        length_previous_dirs(idir+1) = length_previous_dirs(idir) + seconds(hdr{idir}.nSamples/hdr{idir}.Fs);
        RecordDuration = RecordDuration + seconds(hdr{idir}.nSamples/hdr{idir}.Fs); %en secondes
    end
end


%Nb of seizures
if isRat
    nbSlowWaves = size(SlowWave_clock,2);
elseif isPatient
    nbSlowWaves = size(SlowWave_R_clock,2)+size(SlowWave_L_clock,2);
end

nbSeizures = size(Crise_Start_clock,2);

%nb of artefacts : TO ADAPT TO NEW ARTEFACT DETECTION METHOD

% for imarker = 1:length(cfg.LFP.name)
%     n_artefact.(cfg.LFP.name{imarker}) = [];
%     for idir = 1:length(MuseStruct{ipart})
%         if isfield(MuseStruct{ipart}{idir}, cfg.LFP.name{imarker})
%             if isfield(MuseStruct{ipart}{idir}.(cfg.LFP.name{imarker}), 'hasartefact')
%                 n_artefact.(cfg.LFP.name{imarker}) = n_artefact.(cfg.LFP.name{imarker}) + sum(MuseStruct{ipart}{idir}.(cfg.LFP.name{imarker}).hasartefact);
%             end
%         end
%     end
% end


%Time between 2 SlowWaves
timeBetween2SlowWaves     = duration.empty;
timeBetween2SlowWaves_R   = duration.empty;
timeBetween2SlowWaves_L   = duration.empty;
seizureLength             = duration.empty;
x_timeBetween2SlowWaves   = duration.empty;
x_timeBetween2SlowWaves_R = duration.empty;
x_timeBetween2SlowWaves_L = duration.empty;
x_seizureLength           = duration.empty;
if isPatient
    iresult = 1;
    for iSlowWave = 1:size(SlowWave_R_clock,2)-1
        if SlowWave_R_clock{2,iSlowWave} == SlowWave_R_clock{2,iSlowWave+1} %ignore if they are not in the same file
            timeBetween2SlowWaves_R(iresult) = SlowWave_R_clock{1,iSlowWave+1}-SlowWave_R_clock{1,iSlowWave}; %Time between two slow waves
            x_timeBetween2SlowWaves_R(iresult) = SlowWave_R_clock{1,iSlowWave}-MuseStruct{ipart}{SlowWave_R_clock{2,iSlowWave}}.starttime+length_previous_dirs(SlowWave_R_clock{2,iSlowWave});
            iresult = iresult +1;
        end
    end
    iresult = 1;
    for iSlowWave = 1:size(SlowWave_L_clock,2)-1
        if SlowWave_L_clock{2,iSlowWave} == SlowWave_L_clock{2,iSlowWave+1} %ignore if they are not in the same file
            timeBetween2SlowWaves_L(iresult) = SlowWave_L_clock{1,iSlowWave+1}-SlowWave_L_clock{1,iSlowWave}; %Time between two slow waves
            x_timeBetween2SlowWaves_L(iresult) = SlowWave_L_clock{1,iSlowWave}-MuseStruct{ipart}{SlowWave_L_clock{2,iSlowWave}}.starttime+length_previous_dirs(SlowWave_L_clock{2,iSlowWave});
            iresult = iresult +1;
        end
    end
elseif isRat
    iresult = 1;
    for iSlowWave = 1:size(SlowWave_clock,2)-1
        if SlowWave_clock{2,iSlowWave} == SlowWave_clock{2,iSlowWave+1} %ignore if they are not in the same file
            timeBetween2SlowWaves(iresult) = SlowWave_clock{1,iSlowWave+1}-SlowWave_clock{1,iSlowWave}; %Time between two slow waves
            x_timeBetween2SlowWaves(iresult) = SlowWave_clock{1,iSlowWave}-MuseStruct{ipart}{SlowWave_clock{2,iSlowWave}}.starttime+length_previous_dirs(SlowWave_clock{2,iSlowWave});
            iresult = iresult +1;
        end
    end
end

%Seizure length
if isRat
    iresult = 1;
    for iSeizure = 1:size(Crise_Start_clock,2)-1
        if Crise_Start_clock{2,iSeizure} == Crise_End_clock{2,iSeizure} %ignore if they are not in the same file
            seizureLength(iresult) = Crise_End_clock{1,iSeizure}-Crise_Start_clock{1,iSeizure}; %Time between two slow waves
            x_seizureLength(iresult) = Crise_Start_clock{1,iSeizure}-MuseStruct{ipart}{Crise_Start_clock{2,iSeizure}}.starttime+length_previous_dirs(Crise_Start_clock{2,iSeizure});
            iresult = iresult +1;
        end
    end
end


%CV CV2 FF
%Mean	SD	SEM	Min	Med	Max	CV2	Fano Factor
%% Matlab result variable creation

Seizure_Infos.ID                              = cfg.prefix(1:end-1);
Seizure_Infos.nbSeizures                      = nbSeizures;
Seizure_Infos.nbSlowWaves                     = nbSlowWaves;
Seizure_Infos.RecordDuration                  = RecordDuration;
if isPatient
    Seizure_Infos.timeBetween2SlowWaves_L     = timeBetween2SlowWaves_L;
    Seizure_Infos.timeBetween2SlowWaves_R     = timeBetween2SlowWaves_R;
    Seizure_Infos.x_timeBetween2SlowWaves_L   = x_timeBetween2SlowWaves_L;
    Seizure_Infos.x_timeBetween2SlowWaves_R   = x_timeBetween2SlowWaves_R;
elseif isRat
    Seizure_Infos.seizureLength               = seizureLength;
    Seizure_Infos.timeBetween2SlowWaves       = timeBetween2SlowWaves;
    Seizure_Infos.x_seizureLength             = x_seizureLength;
    Seizure_Infos.x_timeBetween2SlowWaves     = x_timeBetween2SlowWaves;
    Seizure_Infos.BaselineDuration            = BaselineDuration;
    Seizure_Infos.PostInjDuration             = PostInjDuration;
end


%% Table creation

if ~(exist (cfg.imagesavedir)==7)
    mkdir(cfg.imagesavedir);
    warning('%s did not exist for saving images, create now',cfg.imagesavedir);
end

tablepath = fullfile(cfg.imagesavedir,[cfg.prefix, 'seizure_info.csv']);
tableid = fopen(tablepath,'wt');
fprintf(tableid,'%s\n',cfg.prefix(1:end-1));
fprintf(tableid,'Number of total Slow Waves :;%g\n',nbSlowWaves);
% if isRat
%     fprintf(tableid,'Number of artefacted Slow Waves :;%g\n',sum(n_artefact.SlowWave));
% elseif isPatient
%     fprintf(tableid,'Number of artefacted Slow Waves :;%g\n',sum(n_artefact.SlowWave_L)+sum(n_artefact.SlowWave_R));
%     fprintf(tableid,'Number of right Slow Waves :;%g\n',length(SlowWave_R_clock));
%     fprintf(tableid,'Number of left Slow Waves :;%g\n',length(SlowWave_L_clock));
% end
fprintf(tableid,'Number of seizures :;%g\n\n',nbSeizures);
fprintf(tableid,';Minutes;Hours\n');
fprintf(tableid,'Record duration :;%g;%g\n',minutes(RecordDuration),hours(RecordDuration));
if isRat
    fprintf(tableid,'Baseline duration :;%g;%g\n',minutes(BaselineDuration),hours(BaselineDuration));
    fprintf(tableid,'Post injection duration :;%g;%g\n',minutes(PostInjDuration),hours(PostInjDuration));
end
fprintf(tableid,'\n');

%minutes
fprintf(tableid,'MINUTES :\r');
fprintf(tableid,';Mean;SD;SEM;Min;Med;Max;CV;CV2;Fano Factor\r');
if isRat
    [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(minutes(timeBetween2SlowWaves));
    fprintf(tableid,'Time between 2 SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(minutes(seizureLength));
    fprintf(tableid,'Seizure length;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    
elseif isPatient
    if ~isempty(SlowWave_R_clock)
        [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(minutes(timeBetween2SlowWaves_R));
        fprintf(tableid,'Time between 2 right SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    else
        fprintf(tableid,'No right SlowWaves\n');
    end
    if ~isempty(SlowWave_L_clock)
        [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(minutes(timeBetween2SlowWaves_L));
        fprintf(tableid,'Time between 2 left SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    else
        fprintf(tableid,'No left SlowWaves\n');
    end
    
end
fprintf(tableid,'\n');

%seconds
fprintf(tableid,'SECONDS :\n');
fprintf(tableid,';Mean;SD;SEM;Min;Med;Max;CV;CV2;Fano Factor\r');

if isRat
    [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(seconds(timeBetween2SlowWaves));
    fprintf(tableid,'Time between 2 SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(seconds(seizureLength));
    fprintf(tableid,'Seizure length;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    
elseif isPatient
    if ~isempty(SlowWave_R_clock)
        [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(seconds(timeBetween2SlowWaves_R));
        fprintf(tableid,'Time between 2 right SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    else
        fprintf(tableid,'No right SlowWaves\n');
    end
    if ~isempty(SlowWave_L_clock)
        [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(seconds(timeBetween2SlowWaves_L));
        fprintf(tableid,'Time between 2 left SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    else
        fprintf(tableid,'No left SlowWaves\n');
    end
    
end

fclose('all');

%% Plot
%ispatient : time between 2 slow waves only
%israt : add seizure length
if plotdata
    
    fig=figure;
    
    %% slow waves over time
    if isPatient
        subplot(3,3,1:3);
        hold;
        for iR=1:size(SlowWave_R_clock,2)
            t = minutes(SlowWave_R_clock{1,iR}-MuseStruct{ipart}{SlowWave_R_clock{2,iR}}.starttime)+ minutes(length_previous_dirs(SlowWave_R_clock{2,iR}));
            plot_R{iR} = plot([t,t],[-1,1],'-b');
        end
        for iL=1:size(SlowWave_L_clock,2)
            t = minutes(SlowWave_L_clock{1,iL}-MuseStruct{ipart}{SlowWave_L_clock{2,iL}}.starttime)+ minutes(length_previous_dirs(SlowWave_L_clock{2,iL}));
            plot_L{iL} = plot([t,t],[-1,1],'-r');
        end
        
        %plot limit betweeen 2 files
        for idir = 1:length(MuseStruct{ipart})
            if idir >1
                t = minutes(length_previous_dirs(idir));
                plot_dirlim = plot([t,t],[-2,2],'Color',[0.5 0.5 0.5],'Linewidth', 2, 'LineStyle', '--');
            end
        end
        
    elseif isRat
        subplot(5,3,1:3)
        hold;
        for i=1:size(SlowWave_clock,2)
            t = minutes(SlowWave_clock{1,i}-cfg.injectiontime);
            plot_SW{i} = plot([t,t],[-1,1],'b');
        end
        %Analysis_start, injection time and analysis_end
        t_start = minutes(Baseline_Start_clock{1,1}-cfg.injectiontime);
        plot_start = plot([t_start,t_start],[-2,2],'g','Linewidth', 2);
        
        plot_inj = plot([0,0],[-2,2],'r','Linewidth', 2);

        t_end = minutes(Analysis_End_clock{1,1}-cfg.injectiontime);
        plot_end = plot([t_end,t_end],[-2,2],'g','Linewidth', 2);
        
    end
    
    
    if isPatient
        title(sprintf('%s : %d SlowWave_R and %d SlowWave_L',cfg.prefix(1:end-1), size(SlowWave_R_clock,2),size(SlowWave_L_clock,2)),'Fontsize',18,'Interpreter','none');
        xlabel('Time (minutes)');
        if idir >1
            legend([plot_R{iR}, plot_L{iL}, plot_dirlim], cfg.LFP.name{1:2}, 'Change of EEG file','Location','eastoutside', 'Interpreter','none');
        else
            legend([plot_R{iR}, plot_L{iL}], cfg.LFP.name{1:2}, 'Interpreter','none');
        end
        
    elseif isRat
        title(sprintf('%s : %d SlowWaves',cfg.prefix(1:end-1), size(SlowWave_clock,2)),'Fontsize',18,'Interpreter','none');
        legend([plot_SW{i},plot_start,plot_inj, plot_end], cfg.LFP.name{1}, 'Analysis_Start', 'Injection', 'Analysis_End','Location','eastoutside', 'Interpreter','none');
        axis tight
        xlabel('Time from injection (minutes)');
    end

    ylabel(sprintf('Slow deflexions \noccurences'));
    set(gca,'TickDir','out','FontWeight','bold');
    yticklabels([]);
    yticks(10);
    ylim([-2 2]);
    
    
    
    %% time between 2 slow waves
    if isPatient
        subplot(3,3,[4,7])
        hold;
        if ~isempty(timeBetween2SlowWaves_R)
            plot(minutes(x_timeBetween2SlowWaves_R),minutes(timeBetween2SlowWaves_R),'-ob');
        end
        if ~isempty(timeBetween2SlowWaves_L)
            plot(minutes(x_timeBetween2SlowWaves_L),minutes(timeBetween2SlowWaves_L),'-or');
        end
        
        %plot limit betweeen 2 eeg
        ylim_dirlim = get(gca,'YLim');
        for idir = 1:length(MuseStruct{ipart})
            if idir >1
                t = minutes(length_previous_dirs(idir));
                plot([t,t],ylim_dirlim,'Color',[0.5 0.5 0.5],'Linewidth', 2, 'LineStyle', '--');
            end
        end
        
    elseif isRat
        subplot(5,3,[4,5,7,8])
        hold;
        delayinj = cfg.injectiontime - MuseStruct{ipart}{1}.starttime;
        plot(minutes(x_timeBetween2SlowWaves)-minutes(delayinj),minutes(timeBetween2SlowWaves),'-ob');
        ylim_dirlim = get(gca,'YLim');
        plot_start = plot([t_start,t_start],ylim_dirlim,'g','Linewidth', 2);
        plot_inj = plot([0,0],ylim_dirlim,'r','Linewidth', 2);
        plot_end = plot([t_end,t_end],ylim_dirlim,'g','Linewidth', 2);
    end
    

    
    %title('Time between 2 slow deflexions','Fontsize',15);
    xlabel('Time at recording');
    ylabel(sprintf('Time between \n2 slow deflexions (minutes)'));
    set(gca,'TickDir','out','FontWeight','bold');
    
    if isPatient
        if idir >1
            legend(cfg.LFP.name{1:(2)}, 'New File', 'Interpreter','none');
        else
            legend(cfg.LFP.name{1:(2)}, 'Interpreter','none');
        end
    elseif isRat
        legend(cfg.LFP.name{1}, 'Interpreter','none');
        axis tight
    end
    
    %% Slowwaves distrib
    
    if isPatient
        if ~isempty(timeBetween2SlowWaves_L)
            subplot(3,3,[5,8])
            hold;
            histogram(minutes(timeBetween2SlowWaves_L),'BinWidth',3,'FaceColor','r');
            %title('Distribution of time between 2 LEFT slow deflexions','Fontsize',15);
            xlabel(sprintf('Time between \n2 LEFT slow deflexions (minutes)'));
            ylabel('Nb of occurences');
            set(gca,'TickDir','out','FontWeight','bold');
        end
        if ~isempty(timeBetween2SlowWaves_R)
            subplot(3,3,[6,9])
            hold;
            histogram(minutes(timeBetween2SlowWaves_R),'BinWidth',3,'FaceColor','b');
            %title('Distribution of time between 2 RIGHT slow deflexions','Fontsize',15);
            xlabel(sprintf('Time between \n2 RIGHT slow deflexions (minutes)'));
            ylabel('Nb of occurences');
            set(gca,'TickDir','out','FontWeight','bold');
        end
    elseif isRat
        subplot(5,3,[6,9])
        hold;
        histogram(minutes(timeBetween2SlowWaves),'BinWidth',1,'FaceColor','b');
        %title('Distribution of time between 2 slow deflexions','Fontsize',15);
        xlabel(sprintf('Time between 2 slow deflexions (minutes)'));
        ylabel('Nb of occurences');
        set(gca,'TickDir','out','FontWeight','bold');
    end
    
    
    
    %% Seizure length over time
    if isRat
        subplot(5,3,[10,11,13,14])
        hold;
        plot(minutes(x_seizureLength)-minutes(delayinj),minutes(seizureLength),'-ob');
        ylim_dirlim = get(gca,'YLim');
        plot_start = plot([t_start,t_start],ylim_dirlim,'g','Linewidth', 2);
        plot_inj = plot([0,0],ylim_dirlim,'r','Linewidth', 2);
        plot_end = plot([t_end,t_end],ylim_dirlim,'g','Linewidth', 2);
               
        %title('Time between 2 slow deflexions','Fontsize',15);
        xlabel('Time at recording (minutes)');
        ylabel('Seizure length (minutes)');
        set(gca,'TickDir','out','FontWeight','bold');
        
        legend('Seizure length', 'Interpreter','none');
        axis tight
    end
    
    
    %% Seizure length distrib
    if isRat
        
        subplot(5,3,[12,15])
        hold;
        histogram(minutes(seizureLength),'BinWidth',0.5,'FaceColor','b');
        %title('Distribution of time between 2 slow deflexions','Fontsize',15);
        xlabel('Seizure length (minutes)');
        ylabel('Nb of occurences');
        set(gca,'TickDir','out','FontWeight','bold');
    end
    
    
    %% Print to file
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('%s did not exist for saving images, create now',cfg.imagesavedir);
    end
    
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix, 'seizures_infos']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix, 'seizures_infos']),'-r600');
    close all
    
    
end


% Save result variable
save(fullfile(cfg.datasavedir,[cfg.prefix, 'seizure_infos.mat']),'Seizure_Infos');


end

