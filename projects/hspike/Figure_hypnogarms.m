function Figure2

restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/cbrewer/cbrewer/sigstar-master    
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));    
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epishare-master'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/SPIKY_apr_2021'))    
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\git\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\cbrewer\cbrewer    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
end

ft_defaults
MuseStruct{8}           = [];
LFP_cluster{8}          = [];
LFP_cluster_detected{8} = [];
clusterindx{8}          = [];
config                 = hspike_setparams;

for ipatient = 1 : 8
    MuseStruct{ipatient}                                            = readMuseMarkers(config{ipatient}, false);
    MuseStruct{ipatient}                                            = padHypnogram(MuseStruct{ipatient});
    [marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}]  = hypnogramMuseStats(config{ipatient}, MuseStruct{ipatient}, true);
end

t = table;
t_position = table;

for ipatient = 1 : 8
    for ipart = 1 : 3
        disp(ipart)
        hypnogram{ipatient}.epochtime   = convertTo(hypnogram{ipatient}.starttime ,'epochtime','Epoch','2001-01-01','TicksPerSecond',1);
        indx                            = hypnogram{ipatient}.part == ipart;        
        hypnogram{ipatient}.x(indx)     = hypnogram{ipatient}.epochtime(indx) - hypnogram{ipatient}.epochtime(find(indx, 1));
        [dt, X, Y, label]               = plot_hyp_lines(hypnogram{ipatient}(hypnogram{ipatient}.part == ipart, :));
        
        if X(1).Hour < 12
            X0 = datetime(X.Year,X.Month,X.Day-1,0,0,0, 'Format', 'dd/MM/yyyy HH:mm');
        else
            X0 = datetime(X.Year,X.Month,X.Day,0,0,0, 'Format', 'dd/MM/yyyy HH:mm');
        end

        t_temp          = table;
        t_temp.X        = X';
        t_temp.Y        = Y';
        t_temp.dt       = dt';
        t_temp.label    = label';
        t_temp.part     = ones(height(t_temp), 1) * ipart;
        t_temp.patient  = ones(height(t_temp), 1) * ipatient;
        t_temp.hour     = (X.Hour/24 + X.Minute/24/60 + X.Second/24/60/60)';
        t_temp.hour0    = double(convertTo(X','epochtime','Epoch',X0(1),'TicksPerSecond',1))/24/60/60-1;
        t_temp.radian   = t_temp.hour * pi * 2;
        t_temp.degree   = t_temp.hour * 360;

        t_temp_position             = table;
        t_temp_position.offset      = X.Hour(1)/24 + X.Minute(1)/24/60 + X.Second(1)/24/60/60;
        t_temp_position.duration    = double(convertTo(X(end),'epochtime','Epoch',X(1),'TicksPerSecond',1))/24/60/60;
        t_temp_position.part        = ones(height(t_temp_position), 1) * ipart;
        t_temp_position.patient     = ones(height(t_temp_position), 1) * ipatient;

        t           = [t; t_temp];
        t_position  = [t_position; t_temp_position];
    end
end




% parameters for subplots
nrows       = 24;
w           = 0.92;
hratio      = 1/3;
vratio      = 1/nrows * 0.95;
htop        = 1/nrows * 0.8;
hbottom     = 1/nrows * 0.8;
spacehor    = 0.1;
rshift      = 0.02;
upshift     = -0.00;
imarker     = 1;


fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'position', [0 0 1000/sqrt(2) 1000]);
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');

imarker = 1;
row  = 1;
cmap = cbrewer('qual', 'Set2', 8);
range = [-6, 12];
for ipatient = 1 : 8
    for ipart = 1 : 3
        
        s1 = axes('Position', [hratio*(imarker-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
        sel = t(t.patient == ipatient & t.part == ipart & t.label ~= "NO_SCORE", :);
        plot(sel.hour0*24, sel.Y, 'color', cmap(ipatient, :));
        hold on
        sel_REM = t(t.patient == ipatient & t.part == ipart & t.label == "REM", :);
        for i = 1 : 2 : size(sel_REM)
            plot([sel_REM.hour0(i)*24, sel_REM.hour0(i+1)*24], [sel_REM.Y(i), sel_REM.Y(i+1)], 'linewidth', 2, 'color', cmap(ipatient, :));
        end
        
        ylim([1,5]);
        xlim(range);
        
        % set xlabels
        if ipatient == 8 && ipart == 3
            xt = mod(get(gca,'Xtick') + 24, 24);
            xticklabels(xt);
            set(gca,'TickLabelInterpreter', 'none', 'XColor', 'k', 'YColor', 'none', 'box', 'off', 'TickDir', 'out');
            xlabel('Time of day');
        else
            set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'YColor', 'none', 'box', 'off', 'ytick', [], 'yticklabel', [], 'TickDir', 'out');
        end
        
        % set ylabels 
        if ipatient == 1 && ipart == 1 
            set(gca,'YColor', 'k');
            yticks(1:5);
            yticklabels(["S3","S2","S1","REM","AWAKE"]);
            for i = 1 : 5
                lh = plot(range, [i, i], 'k');
                lh.Color    = [lh.Color 0.1];
            end
        end    
%         title(sprintf('Patient%d Night %d', ipatient, ipart));
        row = row + 1;
    end
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 6)

% write to figure
fname = fullfile('D:\Dropbox\Apps\Overleaf\Hspike', 'hypnograms');
exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
exportgraphics(fig, strcat(fname, '.pdf'));

