%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\wod'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML
    CEDS64LoadLib('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML');

elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip;
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/wod'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'));
end
ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neuralynx

config = dtx_intra_setparams;

%problème merge crises : "fillgaps", puis créer virtual channels, puis svg
%virtual channels

for irat = 1
    ipart = 1;
    
    %% Get LFP data
    
% 	For readCEDMarkers and CED2mat : Can only be done on windows because 
%   of CEDS64 library. For Linux, set force argument to false, it will 
%   load the results.


%si besoin : ajouter ipart aux noms d'images et de data enregistrées

%read all events
%directorylist CED
%datasavedir pour data coupées et event
    [CEDStruct]               = readCEDmarkers(config{irat}, true);
    dat_LFP                   = CED_maketrials(config{irat},CEDStruct,true,true);
    
    %demean
    for ipart = 1:length(dat_LFP)
        for imarker = 1:length(dat_LFP{ipart})
            cfgtemp                 = [];
            cfgtemp.demean          = config{irat}.LFP.demean{imarker};
            cfgtemp.baselinewindow  = config{irat}.LFP.baselinewindow{imarker};
            dat_LFP{ipart}{imarker} = ft_preprocessing(cfgtemp, dat_LFP{ipart}{imarker});
        end
    end
    
    %compute derivative of Vm
    cfgtemp             = [];
    cfgtemp.derivative  = 'yes';
    der                 = ft_preprocessing(cfgtemp, dat_LFP{1}{2});
   
    % plot data, ECoG-M1G
    figure;
    subplot(2,1,1);hold;
    for itrial = 1:size(dat_LFP{1}{1}.trial,2)
        plot(dat_LFP{1}{1}.time{itrial},-dat_LFP{1}{1}.trial{itrial});
    end
    xlim([-1 1]);
  

    % plot data, enveloppe of Vm
    subplot(2,1,2);hold;
    for itrial = 1:size(dat_LFP{1}{1}.trial,2)
%        y = envelope(dat_LFP{1}{2}.trial{itrial},1000,'rms');
       
       y =  dat_LFP{1}{2}.trial{itrial};  
       t = dat_LFP{1}{2}.time{itrial};
       [AP_vals, AP_locs] = findpeaks(y,t,'MinPeakProminence',20);
       toremove = false(size(y));
       for i = 1:size(AP_locs,2)
           toremove = (toremove | (t>AP_locs(i)-0.001 & t<AP_locs(i)+0.001) );
       end
       y(toremove) = NaN;
       toremove2 = der.trial{itrial}>0.2 | der.trial{itrial}<-0.2;
       y(toremove2) = NaN;
       
       plot(dat_LFP{1}{2}.time{itrial},y);
    end
    xlim([-1 1]);
    ylim([-10 40]);
    
    % sum intra
    all_intra_summed = dat_LFP{1}{2}.trial{1};
    for itrial = 2:size(dat_LFP{1}{2}.trial,2)
        all_intra_summed = all_intra_summed + dat_LFP{1}{2}.trial{itrial};
    end
    plot(dat_LFP{1}{2}.time{1},all_intra_summed);
    
    %% Rm dyn
    rmdyn = readtable('Z:\analyses\lgi1\DTX-INTRA\neurones-DTX\test2 RmDyn\Rm_dyn_vals.txt'); %output xy view from Spike2
    rmdyn = table2array(rmdyn)';
    %1st line : time
    %2nd line : rm values
    figure;
    scatter(rmdyn(1,:), rmdyn(2,:),'.k');
    rmdyn_movmean = movmean(rmdyn(2,:),20); %sur 200ms
    plot(rmdyn(1,:), rmdyn_movmean,'LineWidth',2);
%     plot(rmdyn(1,:), rmdyn(2,:));
    xlim([2416.004 2418.598])
    
    %% FI 
    FI = readtable('Z:\analyses\lgi1\DTX-INTRA\temporaire_fi.csv'); %output xy view from Spike2
    %FI.DTX_ith(end-1) = nan; FI.DTX_gain(end-1) = nan; FI.Ctrl_ith(4) = nan; FI.Ctrl_gain(4) = nan;
    p1 = ranksum(FI.DTX_ith, FI.Ctrl_ith);
    p2 = ranksum(FI.DTX_gain, FI.Ctrl_gain);

    fig = figure; hold on;
    bar(1, nanmean(FI.DTX_ith),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    bar(2, nanmean(FI.Ctrl_ith),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    scatter(rand(size(FI.DTX_ith))*0.2+1-0.1, FI.DTX_ith, 100, 'o','filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'white');
    scatter(rand(size(FI.Ctrl_ith))*0.2+2-0.1, FI.Ctrl_ith, 100, 'o','filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'white');
    errorbar(1, nanmean(FI.DTX_ith),0, nanstd(FI.DTX_ith),'k','CapSize',10);
    errorbar(2, nanmean(FI.Ctrl_ith),0, nanstd(FI.Ctrl_ith),'k','CapSize',10);
    setfig();
    xticks([]);
    title(sprintf('p=%g', p1));
    fname = 'Z:\analyses\lgi1\DTX-INTRA\temporaire_fi_ith';
    dtx_savefigure(fig, fname, 'png', 'pdf', 'close');
    
    fig = figure; hold on;
    bar(1, nanmean(FI.DTX_gain),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    bar(2, nanmean(FI.Ctrl_gain),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    scatter(rand(size(FI.DTX_gain))*0.2+1-0.1, FI.DTX_gain, 100, 'o','filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'white');
    scatter(rand(size(FI.Ctrl_gain))*0.2+2-0.1, FI.Ctrl_gain, 100, 'o','filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'white');
    errorbar(1, nanmean(FI.DTX_gain),0, nanstd(FI.DTX_gain),'k','CapSize',10);
    errorbar(2, nanmean(FI.Ctrl_gain),0, nanstd(FI.Ctrl_gain),'k','CapSize',10);
    setfig();
    xticks([]);
    title(sprintf('p=%g', p2));
    fname = 'Z:\analyses\lgi1\DTX-INTRA\temporaire_fi_gain';
    dtx_savefigure(fig, fname, 'png', 'pdf', 'close');
end
