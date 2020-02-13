function dtx_plot_count_seizure(cfg, MuseStruct, savedata)

%Descriptive statistics of seizures from Muse markers
%BaselineStart = Analysis Start
%Modifs pour avake : pas d'injection. Ou alors la rentrer dans cfg ?
%pour appliquer à tous : rentrer cfg entier. préciser ipatient/irat dans
%argument à part

%get all marker times
SlowWave_R_total = dtx_concatenateMuseMarkers(cfg,MuseStruct, 'SlowWave_R');
SlowWave_L_total = dtx_concatenateMuseMarkers(cfg, MuseStruct,'SlowWave_L');
SlowWave_total = dtx_concatenateMuseMarkers(cfg, MuseStruct,'SlowWave');
Crise_Start_total = dtx_concatenateMuseMarkers(cfg, MuseStruct,'Crise_Start');
Crise_End_total = dtx_concatenateMuseMarkers(cfg, MuseStruct,'Crise_End');
Baseline_Start = dtx_concatenateMuseMarkers(cfg, MuseStruct,'Baseline_Start');
Analysis_End = dtx_concatenateMuseMarkers(cfg, MuseStruct,'Analysis_End');
Injection = dtx_concatenateMuseMarkers(cfg, MuseStruct,'Injection');

%recognize rat or patient data
isRat = 0;
isPatient = 0;
if ~isempty(SlowWave_R_total) || ~isempty(SlowWave_L_total)
    isPatient = 1;
elseif ~isempty(SlowWave_total)
    isRat = 1;
end

%avoir error in crise annotation
if isRat
    if size(Crise_Start_total,2) ~= size(Crise_End_total,2)
        error('%s : there are not the same number of Crise_Start than Crise_End', cfg.prefix(1:end-1));
    elseif isempty(Baseline_Start) || isempty(Analysis_End)
        error('%s : miss Baseline_Start or Analysis_Start', cfg.prefix(1:end-1));
    end
end


%Nb of seizures
if isRat
    nbSlowWaves = size(SlowWave_total,2);
elseif isPatient
    nbSlowWaves = size(SlowWave_R_total,2)+size(SlowWave_L_total,2);
end

nbSeizures = size(Crise_Start_total,2);

%Time between 2 SlowWaves
timeBetween2SlowWaves     = duration.empty;
timeBetween2SlowWaves_R   = duration.empty;
timeBetween2SlowWaves_L   = duration.empty;
seizureLength             = duration.empty;
x_timeBetween2SlowWaves   = datetime.empty;
x_timeBetween2SlowWaves_R = datetime.empty;
x_timeBetween2SlowWaves_L = datetime.empty;
x_seizureLength           = datetime.empty;
if isPatient
    iresult = 1;
    for iSlowWave = 1:size(SlowWave_R_total,2)-1
        if SlowWave_R_total{2,iSlowWave} == SlowWave_R_total{2,iSlowWave+1} %ignore if they are not in the same file
            timeBetween2SlowWaves_R(iresult) = SlowWave_R_total{1,iSlowWave+1}-SlowWave_R_total{1,iSlowWave}; %Time between two slow waves
            x_timeBetween2SlowWaves_R(iresult) = SlowWave_R_total{1,iSlowWave};
            iresult = iresult +1;
        end
    end
    iresult = 1;
    for iSlowWave = 1:size(SlowWave_L_total,2)-1
        if SlowWave_L_total{2,iSlowWave} == SlowWave_L_total{2,iSlowWave+1} %ignore if they are not in the same file
            timeBetween2SlowWaves_L(iresult) = SlowWave_L_total{1,iSlowWave+1}-SlowWave_L_total{1,iSlowWave}; %Time between two slow waves
            x_timeBetween2SlowWaves_L(iresult) = SlowWave_L_total{1,iSlowWave};
            iresult = iresult +1;
        end
    end
elseif isRat
    iresult = 1;
    for iSlowWave = 1:size(SlowWave_total,2)-1
        if SlowWave_total{2,iSlowWave} == SlowWave_total{2,iSlowWave+1} %ignore if they are not in the same file
            timeBetween2SlowWaves(iresult) = SlowWave_total{1,iSlowWave+1}-SlowWave_total{1,iSlowWave}; %Time between two slow waves
            x_timeBetween2SlowWaves(iresult) = SlowWave_total{1,iSlowWave};
            iresult = iresult +1;
        end
    end
end

%Seizure length
if isRat
    iresult = 1;
    for iSeizure = 1:size(Crise_Start_total,2)-1
        if Crise_Start_total{2,iSeizure} == Crise_End_total{2,iSeizure} %ignore if they are not in the same file
            seizureLength(iresult) = Crise_End_total{1,iSeizure}-Crise_Start_total{1,iSeizure}; %Time between two slow waves
            x_seizureLength(iresult) = Crise_Start(1,iSeizure);
            iresult = iresult +1;
        end
    end
end

%Total time
if isRat
    BaselineDuration = Injection{1,1}-Baseline_Start{1,1};
    RecordDuration = Analysis_End{1,1}-Baseline_Start{1,1};
    PostInjDuration = RecordDuration - BaselineDuration;
elseif isPatient
    RecordDuration = 0;
    for idir = 1:length(cfg.directorylist{1})
        hdr = ft_read_header(fullfile(cfg.rawdir,[cfg.directorylist{1}{idir} '.TRC']));
        RecordDuration = seconds(RecordDuration + hdr.nSamples/hdr.Fs); %en secondes
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

tablepath = fullfile(cfg.imagesavedir,'seizures_info',[cfg.prefix, 'seizure_info.csv']);
tableid = fopen(tablepath,'wt');
fprintf(tableid,'%s\n',cfg.prefix(1:end-1));
fprintf(tableid,'Number of Slow Waves :;%g\n',nbSlowWaves);
fprintf(tableid,'Number of seizures :;%g\n\n',nbSeizures);
fprintf(tableid,';Minutes;Hours\n');
fprintf(tableid,'Record duration :;%g;%g\n',minutes(RecordDuration),hours(RecordDuration));
if isRat
    fprintf(tableid,'Baseline duration :;%g;%g\n',minutes(BaselineDuration),hours(BaselineDuration));
    fprintf(tableid,'Post injection duration :;%g;%g\n',minutes(PostInjDuration),hours(PostInjDuration));
end
fprintf(tableid,'\n');

%minuts
fprintf(tableid,'MINUTES :\r');
fprintf(tableid,';Mean;SD;SEM;Min;Med;Max;CV;CV2;Fano Factor\r');
if isRat
    [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(minutes(timeBetween2SlowWaves));
    fprintf(tableid,'Time between 2 SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(minutes(seizureLength));
    fprintf(tableid,'Seizure length;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);

elseif isPatient
    if ~isempty(SlowWave_R_total)
        [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(minutes(timeBetween2SlowWaves_R));
        fprintf(tableid,'Time between 2 right SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    else
        fprintf(tableid,'No right SlowWaves\n');
    end
    if ~isempty(SlowWave_L_total)
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
    if ~isempty(SlowWave_R_total)
        [Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor] = dtx_descr_stats(seconds(timeBetween2SlowWaves_R));
        fprintf(tableid,'Time between 2 right SlowWaves;%g;%g;%g;%g;%g;%g;%g;%g;%g\n',Mean,SD,SEM,Min,Med,Max,cv,cv2,FanoFactor);
    else
        fprintf(tableid,'No right SlowWaves\n');
    end
    if ~isempty(SlowWave_L_total)
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

fig=figure;

%slow waves over time
if isPatient
    subplot(3,3,1:3);
    hold;
    for i=1:size(SlowWave_R_total,2)
        t = SlowWave_R_total{1,i};
        line([t,t],[-1,1],'color','b');
    end
    for i=1:size(SlowWave_L_total,2)
        t = SlowWave_L_total{1,i};
        line([t,t],[-1,1],'color','r');
    end
    
elseif isRat
    subplot(5,3,1:3)
    hold;
    for i=1:size(SlowWave_total,2)
        t = SlowWave_total{1,i};
        line([t,t],[-1,1],'color','b'); 
    end
end

title(sprintf('%s : seizures quantification',cfg.prefix(1:end-1)),'Fontsize',18,'Interpreter','none');
xlabel('Time at recording');
ylabel('Slow deflexions occurences');
set(gca,'TickDir','out','FontWeight','bold');
yticklabels([]);
yticks(10);
ylim([-2 2]);
legend(cfg.LFP.name,'Interpreter','none');


%time between 2 slow waves
if isPatient
    subplot(3,3,[4,7])
    hold;
    if ~isempty(timeBetween2SlowWaves_R)
        plot(x_timeBetween2SlowWaves_R,timeBetween2SlowWaves_R,'-ob');
    end
    if ~isempty(timeBetween2SlowWaves_L)
        plot(x_timeBetween2SlowWaves_L,timeBetween2SlowWaves_L,'-or');
    end
    
elseif isRat
    subplot(5,3,[4,5,7,8])
    plot(x_timeBetween2SlowWaves,timeBetween2SlowWaves,'-ob');
end


%title('Time between 2 slow deflexions','Fontsize',15);
xlabel('Time at recording');
ylabel('Time between 2 slow deflexions');
set(gca,'TickDir','out','FontWeight','bold');
legend(cfg.LFP.name,'Interpreter','none');

%Slowwaves distrib

if isPatient
    if ~isempty(timeBetween2SlowWaves_L)
        subplot(3,3,[5,8])
        hold;
        histogram(timeBetween2SlowWaves_L,'BinWidth',seconds(60),'FaceColor','r');
        %title('Distribution of time between 2 LEFT slow deflexions','Fontsize',15);
        xlabel('Time between 2 LEFT slow deflexions');
        ylabel('Nb of occurences');
        set(gca,'TickDir','out','FontWeight','bold');
    end
    if ~isempty(timeBetween2SlowWaves_R)
        subplot(3,3,[6,9])
        hold;
        histogram(timeBetween2SlowWaves_R,'BinWidth',seconds(60),'FaceColor','b');
        %title('Distribution of time between 2 RIGHT slow deflexions','Fontsize',15);
        xlabel('Time between 2 RIGHT slow deflexions');
        ylabel('Nb of occurences');
        set(gca,'TickDir','out','FontWeight','bold');
    end
elseif isRat
    subplot(5,3,[6,9])
    hold;
    histogram(timeBetween2SlowWaves,'BinWidth',seconds(10),'FaceColor','b');
    %title('Distribution of time between 2 slow deflexions','Fontsize',15);
    xlabel('Time between 2 slow deflexions');
    ylabel('Nb of occurences');
    set(gca,'TickDir','out','FontWeight','bold');
end



% Seizure length
if isRat
    subplot(5,3,[10,11,13,14])
    plot(x_seizureLength,seizureLength,'-ob');
    %title('Seizure length','Fontsize',15);
    xlabel('Time at recording');
    ylabel('Seizure length');
    set(gca,'TickDir','out','FontWeight','bold');
    
    subplot(5,3,[12,15])
    hold;
    histogram(x_seizureLength,'BinWidth',seconds(10),'FaceColor','r');
    %title('Distribution of time between 2 seizures','Fontsize',15);
    xlabel('Time between 2 seizures');
    ylabel('Nb of occurences');
    set(gca,'TickDir','out','FontWeight','bold');
end

%% Print to file
if savedata
    
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,'seizures_info',[cfg.prefix, 'seizures_infos']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,'seizures_info',[cfg.prefix, 'seizures_infos']),'-r600');
    close all
    
    
    % Save result variable
    save(fullfile(cfg.datasavedir,[cfg.prefix, 'seizure_infos.mat']),'Seizure_Infos');
    
end

end

