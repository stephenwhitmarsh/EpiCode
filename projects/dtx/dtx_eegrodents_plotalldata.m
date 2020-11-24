function dtx_eegrodents_plotalldata(ipatient, config_script)

%config_script : name of the configuration script which will be called with
%'eval'. So this script can be used with several parameters scripts


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

feature('DefaultCharacterSet', 'CP1252'); % To fix bug for weird character problems in reading neurlynx

config = eval(config_script);%dtx_setparams_eegvideo;

ipart = 1;

%% load precomputed data

fprintf('Reading %s\n', fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']));
temp = load(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP');
LFP{ipatient} = temp.LFP;
clear temp

%% plot overdraw of each channel
for do_hpfilter = ["not_refiltered", "hpfilt_0_15", "hpfilt_1"]
    for idata = ["EEG", "EMG"]
        for markername = string(config{ipatient}.name)
            if ~isfield(LFP{ipatient}{ipart},markername)
                continue
            end
            if isempty(LFP{ipatient}{ipart}.(markername))
                continue
            end
            clear data*
            
            %prepare eeg data to plot
            switch idata
                case "EEG"
                    %remove EMG
                    cfgtemp = [];
                    cfgtemp.channel = {'all', '-EMG*'};
                    data_plot = ft_selectdata(cfgtemp, LFP{ipatient}{ipart}.(markername));
                case "EMG"
                    %select M1G
                    cfgtemp         = [];
                    cfgtemp.channel = config{ipatient}.LFP.motorcortex;
                    data_plot   = ft_selectdata(cfgtemp, LFP{ipatient}{ipart}.(markername));
            end
            
            fprintf('Low pass filtering before plotting %s\n', idata);
            cfgtemp                             = [];
            cfgtemp.lpfilter                    = 'yes';
            cfgtemp.lpfreq                      = 30;
            cfgtemp.lpfilttype                  = 'but';
            cfgtemp.lpinstabilityfix            = 'reduce';
            data_plot                           = ft_preprocessing(cfgtemp,data_plot);
            
            if strcmp(do_hpfilter, "hpfilt_0_15")
                fprintf('High pass filtering (0.15Hz) before plotting %s\n', idata);
                cfgtemp                             = [];
                cfgtemp.hpfilter                    = 'yes';
                cfgtemp.hpfreq                      = 0.15;
                cfgtemp.hpfilttype                  = 'but';
                cfgtemp.hpinstabilityfix            = 'reduce';
                data_plot                           = ft_preprocessing(cfgtemp,data_plot);
            elseif strcmp(do_hpfilter, "hpfilt_1")
                fprintf('High pass filtering (1Hz) before plotting %s\n', idata);
                cfgtemp                             = [];
                cfgtemp.hpfilter                    = 'yes';
                cfgtemp.hpfreq                      = 1;
                cfgtemp.hpfilttype                  = 'but';
                cfgtemp.hpinstabilityfix            = 'reduce';
                data_plot                           = ft_preprocessing(cfgtemp,data_plot);
            end
            
            for toi = [1,2]
                fig = figure;hold on
                
                %remove EMG
                %data_plot = LFP{ipatient}{ipart}.(markername);
                
                %data_plot_avg = ft_timelockanalysis([],data_plot);
                
                color = 'k';%color of trials when plotted
                h = config{ipatient}.plotseizure.h /2;%5000;% FIXME à ajuster

                %plot each trial
                for itrial = 1:size(data_plot.trialinfo,1)
                    waitbar(itrial/size(data_plot.trialinfo,1));
                    %plot EEG
                    ichan = 0;
                    for channame = string(config{ipatient}.LFP.electrodetoplot)% ["PtA","M1D","M1G"]
                        if ~any(strcmp(data_plot.label, channame))
                            continue
                        end
                        ichan = ichan+1;
                        %select channel
                        cfgtemp = [];
                        cfgtemp.channel = convertStringsToChars(channame);
                        data_1chan = ft_selectdata(cfgtemp,data_plot);
                        if isempty(data_1chan.label)
                            continue
                        end
                        p = plot(data_1chan.time{itrial},data_1chan.trial{itrial}+(numel(data_plot.label)+1)*h-h*ichan,'Color', color); %first on top
                        p.Color(4) = 0.2;
                        %                     %plot average after the last trial
                        %                     if itrial == size(data_plot.trialinfo,1)
                        %                         chan_idx = strcmp(channame, data_plot_avg.label);
                        %                         plot(data_plot_avg.time, data_plot_avg.avg(chan_idx,:)+(numel(data_plot.label)+1)*h-h*ichan, 'b', 'LineWidth', 2);
                        %                     end
                        
                    end
                    
                    %plot EMG if needed
                    if strcmp(idata,'EMG')

                        cfgtemp = [];
                        cfgtemp.channel = ft_getopt(config{ipatient}, 'EMG.(markername)', []);
                        EMG = ft_selectdata(cfgtemp,LFP{ipatient}{ipart}.(markername));
                        if ~isempty(EMG.label)
                            ichan = ichan+2;
                            rescale = 2 * h / (max(EMG.trial{itrial}) - min(EMG.trial{itrial}));
                            p = plot(EMG.time{itrial}, EMG.trial{itrial}.*rescale + (numel(data_plot.label)+1)*h-h*ichan, 'Color',color); %first on top
                            p.Color(4) = 0.2;
                            %                     %plot average of rectified emg after the last trial
                            %                     if itrial == size(data_plot.trialinfo,1)
                            %                         cfgtemp = [];
                            %                         cfgtemp.channel = convertStringsToChars(channame);
                            %                         cfgtemp.rectify = 'yes';
                            %                         emg_rectified = ft_preprocessing(cfgtemp, data_plot);
                            %                         emg_rectified_avg = ft_timelockanalysis([],emg_rectified);
                            %                         cfgtemp           = [];
                            %                         cfgtemp.baseline  = config{ipatient}.LFP.baselinewindow.(markername);
                            %                         emg_rectified_avg = ft_timelockbaseline(cfgtemp, emg_rectified_avg);
                            %
                            %                         rescale = h / (max(emg_rectified_avg.avg) - min(emg_rectified_avg.avg));
                            %                         plot(emg_rectified_avg.time, emg_rectified_avg.avg.*rescale+ (numel(data_plot.label)+1)*h-h*ichan, 'b', 'LineWidth', 2);
                            %                     end
                            %                 elseif length(EMG.label)>1
                            %                     error('several EMG channel. It should have only one');
                            %
                        end
                    end
                end
                
                
                %set figure display
                axis tight
                ylim([-h (numel(data_plot.label)+2)*h]);
                xlim([-toi toi]);
                xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
                ylabel('Channel name', 'Fontsize',15);
                tick = h;
                yticks(h : tick : numel(data_plot.label)*h);
                %chan_name_all = [string(config{ipatient}.LFP.channel), string(EMG.label)];
                chan_name_all = ["PtA","M1D","M1G"];%[string(data_plot.label)'];
                set(gca, 'YTickLabel',flip(chan_name_all));
                set(gca, 'FontWeight','bold', 'Fontsize',15);
                set(gca,'TickDir','out');
                
                title(sprintf('%s : %s %s',config{ipatient}.prefix(1:end-1),markername, do_hpfilter),'Interpreter','none','Fontsize',18);
                
                %print to file
                fname = fullfile(config{ipatient}.imagesavedir, '..', 'plot_each_seizure_overdraw',convertStringsToChars(markername),convertStringsToChars(idata), sprintf('toi_%ds',toi),convertStringsToChars(do_hpfilter),sprintf('%soverdraw_all_%s_toi_%ds_%s',config{ipatient}.prefix,markername,toi,do_hpfilter));
                dtx_savefigure(fig,fname,'pdf','png','close');
            end%toi
        end%markername
    end%idata
end %do_hpfilter

return

%% plot average of each channel
for markername = string(config{ipatient}.name)
    if ~isfield(LFP{ipatient}{ipart},markername)
        continue
    end
    if isempty(LFP{ipatient}{ipart}.(markername))
        continue
    end
    iplot = 0;
    fig = figure;hold on
    iplot = iplot+1;
    
    %compute avg data : already done in the previous loop
%     data_plot = LFP{ipatient}{ipart}.(markername);
%     data_plot_avg = ft_timelockanalysis([],data_plot);
%     %compute avg rectified emg
%     cfgtemp = [];
%     cfgtemp.channel = config{ipatient}.EMG.(markername);
%     cfgtemp.rectify = 'yes';
%     emg_rectified = ft_preprocessing(cfgtemp, data_plot);
%     emg_rectified_avg = ft_timelockanalysis([],emg_rectified);
%     cfgtemp           = [];
%     cfgtemp.baseline  = config{ipatient}.LFP.baselinewindow.(markername);
%     emg_rectified_avg = ft_timelockbaseline(cfgtemp, emg_rectified_avg);
    
    toi = [-2 2];
    
    clear leg
    for ichan = 1:size(data_plot_avg.label,1)
        if strfind(data_plot_avg.label{ichan}, 'EMG')
            continue
%             leg{ichan} = plot(data_plot_avg.time, emg_rectified_avg.avg.*100 + 1000);
        else
            leg{ichan} = plot(data_plot_avg.time, data_plot_avg.avg(ichan,:));
        end
    end
    
    %add patch of selected period
    ax  = axis;
    toi_temp = config{ipatient}.morpho.toiac;
    x   = [toi_temp(1) toi_temp(2) toi_temp(2) toi_temp(1)];
    y   = [ax(3) ax(3) ax(4) ax(4)];
    patch('XData',x,'YData',y,'facecolor',[0 0 0],'edgecolor','none','facealpha',0.1);
    
    %set figure display
    legend([leg{:}],data_plot_avg.label,'Location','northeastoutside');
    axis tight
    xlim(toi);
    xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
    %ylabel('ï¿½V', 'Fontsize',15);
    set(gca, 'FontWeight','bold', 'Fontsize',15);
    set(gca,'TickDir','out');
    
    title(sprintf('%s : %s',config{ipatient}.prefix(1:end-1),convertStringsToChars(markername)),'Interpreter','none','Fontsize',18);
    
    %print to file
    fname = fullfile(config{ipatient}.imagesavedir, '..', 'plot_each_seizure_avg',sprintf('%s%s_avg',config{ipatient}.prefix,markername));
    dtx_savefigure(fig,fname,'pdf','png','fig','close');
    
end


end
