function data_cleaned = select_trials_visual(cfg,data,ipart,markername,force)

%get default parameters
cfg.rejectvisual = ft_getopt(cfg,'rejectvisual',[]);
suffix           = ft_getopt(cfg.rejectvisual, 'suffix', []);
method           = ft_getopt(cfg.rejectvisual, 'method', 'seltrials');%trials, periods
write            = ft_getopt(cfg.rejectvisual, 'write', true);

fname = fullfile(cfg.datasavedir, sprintf('%sp%d_%s_VisualSelection_%s.mat', cfg.prefix,ipart,convertStringsToChars(markername), suffix));

%load preselected data if required
if exist(fname,'file') && force == false
    fprintf('*********************************\n');
    fprintf('** loading pre-selected trials **\n');
    fprintf('*********************************\n\n');
    load(fname,'data_cleaned');
    return
else
    fprintf('**********************\n');
    fprintf('** selecting trials **\n');
    fprintf('**********************\n\n');
end

switch method
    case 'trials'
        %select data for slow wave morpho (less selective, check only the channel to measure)
        %for each trial, can reject trial or channel 
        waitfor(msgbox(sprintf('%s, p%d, %s, %s, reject trials and/or channels',cfg.prefix(1:end-1),ipart,markername, suffix)));
        cfgtemp             = [];
        cfgtemp.method      = 'trial';
        cfgtemp.keepchannel = 'nan';
        cfgtemp.box         = 'yes';
        cfgtemp.axis        = 'yes';
        data_cleaned        = ft_rejectvisual(cfgtemp,data);
    case 'periods'
        
        %prepare data for fieldtrip's GUI
%         h = 120;
%         chanlabels = [cfg.LFP.channel, ft_getopt(cfg.EMG,convertStringsToChars(markername),[])];
%         [~,chan_plot_order] = match_str(chanlabels,data.label);
%         for itrial = 1:size(data.trial,2)
%             for ichan = chan_plot_order
%                 h_temp = h*(length(chan_plot_order)+1) - h * ichan;
%                 data.trial{itrial}(ichan,:) = data.trial{itrial}(ichan,:) + h_temp;
%             end
%         end

        waitfor(msgbox(sprintf('%s, p%d, %s, %s, remove artefacted trials and channels before selecting periods',cfg.prefix(1:end-1),ipart,markername, suffix)));
        %check channels
        cfgtemp             = [];
        cfgtemp.method      = 'trial';
        cfgtemp.keepchannel = 'nan';
        cfgtemp.box         = 'yes';
        data_cleaned        = ft_rejectvisual(cfgtemp,data);

        waitfor(msgbox(sprintf('%s, p%d, %s, %s, select periods',cfg.prefix(1:end-1),ipart,markername, suffix)));
        cfg_artefacts                   = ft_databrowser([], data_cleaned);
        cfg_artefacts.artfctdef.reject 	= 'partial';
        cfg_artefacts.artfctdef.invert	= 'yes';
        
        data_cleaned = ft_rejectartifact(cfg_artefacts,data_cleaned);
        
       
        
end

if write
    fprintf('saving data to %s\n',fname);
    save(fname,'data_cleaned', '-v7.3');
end

end


%                 for idata = ["lfp csd"]
%                     %for xcorr propagation : select only the slowwave period
%                     cfg_sel         = [];
%                     cfg_sel.channel = config{ipatient}.(markername).channel; %all frontopolar unilateral channels
%                     datatemp        = ft_selectdata(cfg_sel, data_sel);
%                     config{ipatient}.rejectvisual.suffix          	= 'propagation';
%                     config{ipatient}.rejectvisual.method            = 'periods';
%                     data.(idata).propagation{ipart}.(markername)    = select_trials_visual(config{ipatient},datatemp,ipart,markername, true);
%                     
%                     %for topoplot : remove artefacted trials or channels :
%                     %between -2 and 2
%                     cfg_sel         = [];
%                     cfg_sel.channel = config{ipatient}.LFP.channel; %all lfp channels
%                     datatemp        = ft_selectdata(cfg_sel, data_sel);
%                     config{ipatient}.rejectvisual.suffix           	= 'topoplot';
%                     config{ipatient}.rejectvisual.method            = 'trials';
%                     config{ipatient}.rejectvisual.channel       	= config{ipatient}.LFP.channel;
%                     data.(idata).topoplot{ipart}.(markername)       = select_trials_visual(config{ipatient}, datatemp,ipart,markername, true);
%                     
%                     %for slowwave morphology : focus on selected channel and
%                     %removed artefacted ones
%                     cfg_sel         = [];
%                     cfg_sel.channel = config{ipatient}.align.channel.(markername); %only slowwave maximum channel
%                     datatemp        = ft_selectdata(cfg_sel, data_sel);
%                     config{ipatient}.rejectvisual.suffix         	= 'slowwavemorpho';
%                     config{ipatient}.rejectvisual.method            = 'trials';
%                     data.(idata).slowwavemorpho{ipart}.(markername) = select_trials_visual(config{ipatient}, datatemp,ipart,markername, true);
%                 end
