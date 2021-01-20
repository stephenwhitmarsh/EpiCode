%% extract all markers

addpath('~/fieldtrip');
ft_defaults

patient_directory       = '/Volumes/epimicro/Donnees-analyses/Stephen/2230 files with markers/';
directory_searchstring  = '02230*';
data_searchstring       = '*m1*.ncs';

eventtype_start         = 'BAD__START__';
eventtype_end           = 'BAD__END__';
[bad.clocktime,bad.clocktime_append,bad.sample,bad.data,bad.ldir,bad.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'P';
eventtype_end           = 'P';
[P.clocktime,P.clocktime_append,P.sample,P.data,P.ldir,P.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'RR__START__';
eventtype_end           = 'RR__END__';
[RR.clocktime,RR.clocktime_append,RR.sample,RR.data,RR.ldir,RR.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'VF4__START__';
eventtype_end           = 'VF4__END__';
[VF4.clocktime,VF4.clocktime_append,VF4.sample,VF4.data,VF4.ldir,VF4.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'VF1__START__';
eventtype_end           = 'VF1__END__';
[VF1.clocktime,VF1.clocktime_append,VF1.sample,VF1.data,VF1.ldir,VF1.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'PP__START__';
eventtype_end           = 'PP__END__';
[PP.clocktime,PP.clocktime_append,PP.sample,PP.data,PP.ldir,PP.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'BAD';
eventtype_end           = 'BAD';
[bad2.clocktime,bad2.clocktime_append,bad2.sample,bad2.data,bad2.ldir,bad2.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'VF3__START__';
eventtype_end           = 'VF3__END__';
[VF3.clocktime,VF3.clocktime_append,VF3.sample,VF3.data,VF3.ldir,VF3.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'StartRecord';
eventtype_end           = 'StartRecord';
[start.clocktime,start.clocktime_append,start.sample,start.data,start.ldir,start.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'StopRecord';
eventtype_end           = 'StopRecord';
[stop.clocktime,stop.clocktime_append,stop.sample,stop.data,stop.ldir,stop.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'StartRecord';
eventtype_end           = 'StopRecord';
[rec.clocktime,rec.clocktime_append,rec.sample,rec.data,rec.ldir,rec.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

eventtype_start         = 'PP__START__';
eventtype_end           = 'PP__END__';
[PP.clocktime,PP.clocktime_append,PP.sample,PP.data,PP.ldir,PP.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

%% Plot overview

h = figure; 

subplot(8,1,1); hold;
fill([rec.clocktime_append(1,1),rec.clocktime_append(end,2),rec.clocktime_append(end,2),rec.clocktime_append(1,1)],[0 0 1 1],[1 0 0],'EdgeColor','none','facealpha',1);
for i = 1 : size(rec.clocktime_append,1)
    fill([rec.clocktime_append(i,1),rec.clocktime_append(i,2),rec.clocktime_append(i,2),rec.clocktime_append(i,1)],[0 0 1 1],[0 1 0],'EdgeColor',[0 0 1],'facealpha',1);
end
title('Data');
axis tight
axx = xlim;
xlim(axx)

subplot(8,1,2); hold;
for i = 1 : size(bad.clocktime_append,1)
    fill([bad.clocktime_append(i,1),bad.clocktime_append(i,2),bad.clocktime_append(i,2),bad.clocktime_append(i,1)],[0 0 1 1],[0 0 0],'EdgeColor','none','facealpha',1);
end
title('Artifacts');
xlim(axx)

subplot(7,1,3); hold;
for i = 1 : size(RR.clocktime_append,1)
    fill([RR.clocktime_append(i,1),RR.clocktime_append(i,2),RR.clocktime_append(i,2),RR.clocktime_append(i,1)],[0 0 1 1],[1 0 0],'EdgeColor',[0 0 0],'facealpha',1);
end
title('Fast ripples');
xlim(axx)

subplot(8,1,4); hold;
for i = 1 : size(VF1.clocktime_append,1)
    fill([VF1.clocktime_append(i,1),VF1.clocktime_append(i,2),VF1.clocktime_append(i,2),VF1.clocktime_append(i,1)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',1);
end
title('VF1');
xlim(axx)

% subplot(8,1,5); hold;
% for i = 1 : size(VF3.clocktime_append,1)
%     fill([VF3.clocktime_append(i,1),VF3.clocktime_append(i,2),VF3.clocktime_append(i,2),VF3.clocktime_append(i,1)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',1);
% end
% title('VF3');
% xlim(axx)

subplot(8,1,6); hold;
for i = 1 : size(VF4.clocktime_append,1)
    fill([VF4.clocktime_append(i,1),VF4.clocktime_append(i,2),VF4.clocktime_append(i,2),VF4.clocktime_append(i,1)],[0 0 1 1],[0 0 1],'EdgeColor',[0 0 0],'facealpha',1);
end
title('VF4');
xlim(axx)

subplot(8,1,7); hold;
for i = 1 : size(PP.clocktime_append,1)
    fill([PP.clocktime_append(i,1),PP.clocktime_append(i,2),PP.clocktime_append(i,2),PP.clocktime_append(i,1)],[0 0 1 1],[0 0 1],'EdgeColor',[0 0 0],'facealpha',1);
end
title('PP');
xlim(axx)




subplot(8,1,8); hold;

[P.minutes, P.edges]    = discretize(P.clocktime_append(:,1),'minute');
[P.histcount,edges,bin] = histcounts(P.minutes,-0.5:1:P.minutes(end)+0.5);

bar(P.edges,P.histcount,1,'facecolor','k')
ylabel('Spikes per minute');
axis tight;
xlim(axx)

c = smooth(P.histcount);
plot(P.edges,c,'r','linewidth',1)

% 
% % set(gca,'Xtick',linspace(axx(1),axx(2),12),'XTickLabelRotation',45)
% set(gca,'Xtick',linspace(axx(1),axx(2),12),'XTickLabelRotation',45)
% set(gca,'Xtick',linspace(axx(1),axx(2),12),'XTickLabelRotation',45)

set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'overview_visual.pdf','-r600');
% 
% % polar
% figure;
% 
% [P.hours, P.edges]    = discretize(P.clocktime_append(:,1),'hour');
% [P.histcount,edges,bin] = histcounts(P.hours,-0.5:1:P.hours(end)+0.5);
% 
% figure
% h = polarhistogram(P.minutes / 1440*pi*2,linspace(0,pi*2,1440));
% 
% figure
% h = polarhistogram(P.hours / 24*pi*2,linspace(0,pi*2,24));
% 
% h.DisplayStyle = 'stairs';
% 
% addpath ~/Documents/MATLAB/CircStat2012a/
% figure; circ_plot(P.hours / 24*pi*2,'density','r.',true,'linewidth',2,'color','r')
% figure; circ_plot(P.minutes / 1440*pi*2)

%% Text files and figures for overview markers

% Overview textfile artefacts %%%%%%%%%%%%%%%%%%%%%

Starttime       = [];
Endtime         = [];
Startsample     = [];
Endsample       = [];
Duration        = [];
Filename        = [];
Folder          = [];
Event           = [];

for idir = 1 : size(bad.sample,2)
    for ievent = 1 : size(bad.sample{idir},1)
        Event       = [Event; 'bad'];
        Starttime   = [Starttime;   bad.clocktime{idir}(ievent,1)];
        Endtime     = [Endtime;     bad.clocktime{idir}(ievent,2)];
        Startsample = [Startsample; bad.sample{idir}(ievent,1)];
        Endsample   = [Endsample;   bad.sample{idir}(ievent,2)];  
        Duration    = [Duration;    (bad.sample{idir}(ievent,2) - bad.sample{idir}(ievent,1))/32000;];
        Filename    = [Filename;    bad.ldir(idir).name];
        Folder      = [Folder;      bad.ldir(idir).folder];
    end  
end

overview_bad = table(Event,Folder,Filename,Starttime,Endtime,Duration,Startsample,Endsample);
writetable(overview_bad);

h = figure;
subplot(5,1,1);
histogram(Duration,50);
title(sprintf('Artefacts \n Median: %.2f s Min: %.2f s, Max: %.2f s',median(Duration),min(Duration),max(Duration)));
xlabel('Duration (s)');
ylabel('Nr. of events');

% Overview textfile VF1 %%%%%%%%%%%%%%%%%%%%%%%%%%

Starttime       = [];
Endtime         = [];
Startsample     = [];
Endsample       = [];
Duration        = [];
Filename        = [];
Folder          = [];
Event           = [];

for idir = 1 : size(VF1.sample,2) 
    for ievent = 1 : size(VF1.sample{idir},1)
        Event       = [Event; 'VF1'];
        Starttime   = [Starttime;   VF1.clocktime{idir}(ievent,1)];
        Endtime     = [Endtime;     VF1.clocktime{idir}(ievent,2)];
        Startsample = [Startsample; VF1.sample{idir}(ievent,1)];
        Endsample   = [Endsample;   VF1.sample{idir}(ievent,2)];  
        Duration    = [Duration;    (VF1.sample{idir}(ievent,2) - VF1.sample{idir}(ievent,1))/32000;];
        Filename    = [Filename;    VF1.ldir(idir).name];
        Folder      = [Folder;      VF1.ldir(idir).folder];
    end  
end

overview_VF1 = table(Event,Folder,Filename,Starttime,Endtime,Duration,Startsample,Endsample);
writetable(overview_VF1);

subplot(5,1,2);
histogram(Duration,50);
title(sprintf('VF1 \n Median: %.2f s Min: %.2f s, Max: %.2f s',median(Duration),min(Duration),max(Duration)));
xlabel('Duration (s)');
ylabel('Nr. of events');

% Overview textfile VF3 %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Starttime       = [];
% Endtime         = [];
% Startsample     = [];
% Endsample       = [];
% Duration        = [];
% Filename        = [];
% Folder          = [];
% Event           = [];
% 
% for idir = 1 : size(VF3.sample,2)
%     for ievent = 1 : size(VF3.sample{idir},1)
%         Event       = [Event; 'VF3'];
%         Starttime   = [Starttime;   VF3.clocktime{idir}(ievent,1)];
%         Endtime     = [Endtime;     VF3.clocktime{idir}(ievent,2)];
%         Startsample = [Startsample; VF3.sample{idir}(ievent,1)];
%         Endsample   = [Endsample;   VF3.sample{idir}(ievent,2)];
%         Duration    = [Duration;    (VF3.sample{idir}(ievent,2) - VF3.sample{idir}(ievent,1))/32000;];
%         Filename    = [Filename;    VF3.ldir(idir).name];
%         Folder      = [Folder;      VF3.ldir(idir).folder];
%     end
% end
% 
% overview_VF3 = table(Event,Folder,Filename,Starttime,Endtime,Duration,Startsample,Endsample);
% writetable(overview_VF3);
% 
% subplot(5,1,3);
% histogram(Duration,50);
% title(sprintf('VF3 \n Median: %.2f s Min: %.2f s, Max: %.2f s',median(Duration),min(Duration),max(Duration)));
% xlabel('Duration (s)');
% ylabel('Nr. of events');

% Overview textfile VF4 %%%%%%%%%%%%%%%%%%%%%%%%%%

Starttime       = [];
Endtime         = [];
Startsample     = [];
Endsample       = [];
Duration        = [];
Filename        = [];
Folder          = [];
Event           = [];

for idir = 1 : size(VF4.sample,2) 
    for ievent = 1 : size(VF4.sample{idir},1)
        Event       = [Event; 'VF1'];
        Starttime   = [Starttime;   VF4.clocktime{idir}(ievent,1)];
        Endtime     = [Endtime;     VF4.clocktime{idir}(ievent,2)];
        Startsample = [Startsample; VF4.sample{idir}(ievent,1)];
        Endsample   = [Endsample;   VF4.sample{idir}(ievent,2)];  
        Duration    = [Duration;    (VF4.sample{idir}(ievent,2) - VF4.sample{idir}(ievent,1))/32000;];
        Filename    = [Filename;    VF4.ldir(idir).name];
        Folder      = [Folder;      VF4.ldir(idir).folder];
    end  
end

overview_VF4 = table(Event,Folder,Filename,Starttime,Endtime,Duration,Startsample,Endsample);
writetable(overview_VF4);

subplot(5,1,4);
histogram(Duration,50);
title(sprintf('VF4 \n Median: %.2f s Min: %.2f s, Max: %.2f s',median(Duration),min(Duration),max(Duration)));
xlabel('Duration (s)');
ylabel('Nr. of events');

% Overview textfile fast ripples %%%%%%%%%%%%%%%%%%%%%%%%%%

Starttime       = [];
Endtime         = [];
Startsample     = [];
Endsample       = [];
Duration        = [];
Filename        = [];
Folder          = [];
Event           = [];

for idir = 1 : size(RR.sample,2)
    for ievent = 1 : size(RR.sample{idir},1)
        Event       = [Event; 'RR'];
        Starttime   = [Starttime;   RR.clocktime{idir}(ievent,1)];
        Endtime     = [Endtime;     RR.clocktime{idir}(ievent,2)];
        Startsample = [Startsample; RR.sample{idir}(ievent,1)];
        Endsample   = [Endsample;   RR.sample{idir}(ievent,2)];  
        Duration    = [Duration;    (RR.sample{idir}(ievent,2) - RR.sample{idir}(ievent,1))/32000;];
        Filename    = [Filename;    RR.ldir(idir).name];
        Folder      = [Folder;      RR.ldir(idir).folder];
    end  
end

overview_RR = table(Event,Folder,Filename,Starttime,Endtime,Duration,Startsample,Endsample);
writetable(overview_RR);

subplot(5,1,5);
histogram(Duration,50);
title(sprintf('Fast Ripples\n Median: %.2f s Min: %.2f s, Max: %.2f s',median(Duration),min(Duration),max(Duration)));
xlabel('Duration (s)');
ylabel('Nr. of events');

set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'overview_distribution.pdf','-r600');

% Overview textfile peaks %%%%%%%%%%%%%%%%%%%%%%%%%%
        % append data over files, and add real-time to timings
        clocktime_append = [];
        for i = 1 : size(clocktime,2)
            clocktime_append = [clocktime_append; dt{i} + seconds(clocktime{i})];
        end
        
Starttime       = [];
Endtime         = [];
Startsample     = [];
Endsample       = [];
Duration        = [];
Filename        = [];
Folder          = [];
Event           = [];

for idir = 1 : size(P.sample,2)
    for ievent = 1 : size(P.sample{idir},1)
        Event       = [Event; 'P'];
        Starttime   = [Starttime;   P.clocktime{idir}(ievent,1)];
        Endtime     = [Endtime;     P.clocktime{idir}(ievent,2)];
        Startsample = [Startsample; P.sample{idir}(ievent,1)];
        Endsample   = [Endsample;   P.sample{idir}(ievent,2)];  
        Duration    = [Duration;    (P.sample{idir}(ievent,2) - P.sample{idir}(ievent,1))/32000;];
        Filename    = [Filename;    P.ldir(idir).name];
        Folder      = [Folder;      P.ldir(idir).folder];
    end  
end

overview_P = table(Event,Folder,Filename,Starttime,Endtime,Duration,Startsample,Endsample);
writetable(overview_P);

%% Place artefact file for Spyking-Circus in each data directory

for idir = 1 : size(bad.clocktime,2)
    filename = fullfile(bad.ldir(idir).folder,bad.ldir(idir).name,'artifact.txt');
    fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
    % has to be in ms instead of seconds
    dlmwrite(filename,bad.clocktime{idir}*1000,'delimiter',' ','precision','%.0f'); % high precision, i.e sub-millisecond, does not make sense
end
