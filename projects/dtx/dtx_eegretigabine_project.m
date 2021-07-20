%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\development
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\natsort
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/development
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/natsort
end

ft_defaults

%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info off
ft_debug off
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

configretig = dtx_eegretigabine_setparams;
configdtx   = dtx_eegvideo_setparams;
config = [configdtx(:)', configretig(:)'];

%% #DATA# read Muse marker data
for ipatient = 1:size(config,2)
    fprintf('\nRead marker data for %s\n',config{ipatient}.prefix(1:end-1));
    MuseStruct{ipatient} = readMuseMarkers(config{ipatient},false);
    MuseStruct_concat{ipatient} = concatenateMuseMarkers(config{ipatient},MuseStruct{ipatient},false);    %count seizure infos
end

%% Nb. of seizures over time
firstwin = -3;
lastwin = 7;
winsize = 1;
imeasure = "nb_seizures";
clear nb_crises crises_cumlength data* pval* param

for ipatient = 1:size(config, 2)
    i = 0;
    if strcmp(config{ipatient}.group, 'dtx') % no treatment with retigabine or DMSO
        config{ipatient}.injectionretig = MuseStruct_concat{ipatient}{1}.markers.Analysis_Start.clock(1) + hours(3);
    end
    temps_crises = hours(MuseStruct_concat{ipatient}{1}.markers.Crise_Start.clock - config{ipatient}.injectionretig);
    duree_crises = seconds(MuseStruct_concat{ipatient}{1}.markers.Crise_End.clock - MuseStruct_concat{ipatient}{1}.markers.Crise_Start.clock);
    
    for iwin = firstwin:winsize:lastwin-winsize
        i = i+1;
        startrecord = hours(MuseStruct_concat{ipatient}{1}.markers.Analysis_Start.clock - config{ipatient}.injectionretig);
        endrecord = hours(MuseStruct_concat{ipatient}{1}.markers.Analysis_End.clock - config{ipatient}.injectionretig);
        t_start = iwin;
        t_end = iwin + winsize;
        sel = temps_crises>t_start & temps_crises<t_end;
        param(ipatient, i) = sum(sel) * 1/winsize;
        group(1, ipatient) = string(config{ipatient}.group);
        
        if t_start >= startrecord && t_end < endrecord
            continue
        elseif t_start < startrecord && t_end > startrecord
            perc = (t_end - startrecord)/winsize;
            param(ipatient, i) = param(ipatient, i) * 1/perc;
        elseif t_end > endrecord && t_start < endrecord
            perc = (endrecord - t_start)/winsize;
            param(ipatient, i) = param(ipatient, i) * 1/perc;
        else
            param(ipatient, i) = nan;
        end
    end
end
param(7, 3:end) = nan; %2 rats with less than 10 hours of recording, replace 0 by nans 
param(11,5:end) = nan;

pvalretig_group = [];
pvaldmso_group = [];
datadtx = param(strcmp(group, 'dtx')', :);
datadmso = param(strcmp(group, 'dmso')', :);
dataretig = param(strcmp(group, 'retigabine')', :);

X = [];
y = [];
for irat = 1:size(param, 1)
    y = [y; param(irat, :)'];
    for j = 1:size(param, 2)
        X(end+1, 1) = j;
        if strcmp(config{irat}.group, 'dtx')
            X(end, 2) = 1;
        elseif strcmp(config{irat}.group, 'dmso')
            X(end, 2) = 2;
        elseif strcmp(config{irat}.group, 'retigabine')
            X(end, 2) = 3;
        end
    end
end
statstable       = table.empty;
statstable.data  = y;
statstable.time  = X(:,1);
statstable.group = X(:,2);

%search for significant differences explained by time
mdl   = fitlm(statstable, 'data ~ time + group + group:time'); 
stats = anova(mdl,'component');

statsctrl = statstable(statstable.group == 1, :);
mdl       = fitlm(statsctrl, 'data ~ time'); 
stats1    = anova(mdl,'component');

statsdmso = statstable(statstable.group == 2, :);
mdl       = fitlm(statsdmso, 'data ~ time');
stats2    = anova(mdl,'component');

statsretig = statstable(statstable.group == 3, :);
mdl        = fitlm(statsretig, 'data ~ time');
stats3     = anova(mdl,'component');

%compare each group to its baseline
[pvalctrl, pvalretig, pvaldmso] = deal([]);
i = 0;
for iwin = firstwin:winsize:lastwin-winsize
    i = i+1;
    pvalctrl(i) = ranksum(reshape(datadtx(:,1:3), 1, []), datadtx(:,i));
    pvaldmso(i) = ranksum(reshape(datadmso(:,1:3), 1, []), datadmso(:,i));
    pvalretig(i) = ranksum(reshape(dataretig(:,1:3), 1, []), dataretig(:,i));
end
[h, crit_p, adj_ci_cvrg, pvalretig_corr] = fdr_bh(pvalretig);
[h, crit_p, adj_ci_cvrg, pvaldmso_corr]  = fdr_bh(pvaldmso);
[h, crit_p, adj_ci_cvrg, pvalctrl_corr]  = fdr_bh(pvalctrl);

% plot results
fig = figure; hold on;
for igroup = unique(group)
    rat_list = strcmp(group, igroup)';
    data = param(rat_list, :);
    switch igroup
        case "dtx"
            scattype = '-ok';
            color = 'k';
            shift = -winsize/10;
        case "dmso"
            scattype = '-^k';
            color = [0.6 0.6 0.6];
            shift = 0;
        case "retigabine"
            scattype = '-sk';
            color = [1 0.6 0.2];
            shift = winsize/10;
    end
    errorbar((firstwin+winsize:winsize:lastwin)+shift, mean(data,1,'omitnan'), std(data,0,1,'omitnan'), scattype, 'markerfacecolor', color, 'markersize', 10, 'capsize', 10);
end
xlim([firstwin - winsize, lastwin + winsize]);
set(gca, 'tickdir', 'out', 'fontsize', 15)
plot([firstwin lastwin], [0 0], 'k');

% plot stats
posretig = pvalretig_corr < 0.05;
posdmso  = pvaldmso_corr < 0.05;
t = firstwin+winsize:winsize:lastwin;
y = ylim;
plot(t(posretig), ones(size(t(posretig)))*y(1), '*', 'markeredgecolor', [1 0.6 0.2], 'markersize', 10, 'linewidth', 2);
plot(t(posdmso), ones(size(t(posdmso)))*(y(1) - diff(y)/10), '*', 'markeredgecolor', [0.6 0.6 0.6], 'markersize', 10, 'linewidth', 2);

figname = fullfile(config{end}.imagesavedir, '..', 'condition_comparison', sprintf('%s_win%s', imeasure, strrep(num2str(winsize), '.','_')));
dtx_savefigure(fig, figname, 'png', 'pdf', 'fig', 'close');