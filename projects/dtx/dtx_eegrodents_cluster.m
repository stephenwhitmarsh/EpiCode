function dtx_eegrodents_cluster(ipatient, config_script)

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
do_align = true;
do_morpho = true;

%% align LFP, remove bad trials, correct baseline
if do_align
    MuseStruct = readMuseMarkers(config{ipatient},true);
    MuseStruct_concat = concatenateMuseMarkers(config{ipatient},MuseStruct,true); %for sliding time window of morpho
    
    %read data according to not re-aligned SlowWave, EMG begin, and Crise_End
    config{ipatient}.LFP.name   = {'SlowWave_EMG_begin','Crise_End','SlowWave_not_aligned'};
    LFP_others                   = readLFP(config{ipatient},MuseStruct,true);
    
    %align slowwave to peak and read slowwave
    config{ipatient}.align.name = {'SlowWave'};
    config{ipatient}.LFP.name   = {'SlowWave'};
    MuseStruct                  = alignMuseMarkers(config{ipatient},MuseStruct,true);
    LFP_SlowWave                = readLFP(config{ipatient},MuseStruct,true);
    
    %align slowwave to begin and read slowwave
    config{ipatient}.align.name = {'SlowWave_begin'};
    config{ipatient}.LFP.name   = {'SlowWave_begin'};
    MuseStruct                  = alignMuseMarkers(config{ipatient},MuseStruct,true);
    LFP_SlowWave_begin          = readLFP(config{ipatient},MuseStruct,true);
    
    %append LFP data
    LFP{1}.SlowWave             = LFP_SlowWave{1}.SlowWave;
    LFP{1}.SlowWave_begin       = LFP_SlowWave_begin{1}.SlowWave_begin;
    LFP{1}.SlowWave_EMG_begin   = LFP_others{1}.SlowWave_EMG_begin;
    LFP{1}.Crise_End            = LFP_others{1}.Crise_End;
    LFP{1}.SlowWave_not_aligned = LFP_others{1}.SlowWave_not_aligned; 
    clear LFP_SlowWave* LFP_others
    
    cfgtemp = config{ipatient};
    cfgtemp.rmtrials.plotdata = 'yes';
    cfgtemp.rmtrials.electrodetoplot.SlowWave               = config{ipatient}.align.channel.SlowWave;
    cfgtemp.rmtrials.electrodetoplot.SlowWave_begin         = config{ipatient}.align.channel.SlowWave;
    cfgtemp.rmtrials.electrodetoplot.SlowWave_EMG_begin     = config{ipatient}.align.channel.SlowWave;
    cfgtemp.rmtrials.electrodetoplot.Crise_End              = config{ipatient}.align.channel.SlowWave;
    cfgtemp.rmtrials.electrodetoplot.SlowWave_not_aligned   = config{ipatient}.align.channel.SlowWave;
    LFP = removetrials_MuseMarkers(cfgtemp,LFP,MuseStruct,true);
    
    %remove ref EMG, flip data if required, and correct baseline
    for markername = string(config{ipatient}.name)
        
        if isempty(LFP{ipart}.(markername))
            continue
        end
        
        cfgtemp = [];
        cfgtemp.channel = {'all', '-EMG2'};
        LFP{ipart}.(markername) = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
        
        if istrue(config{ipatient}.LFP.flip)
            for itrial = 1:size(LFP{ipart}.(markername).trial,2)
                LFP{ipart}.(markername).trial{itrial} = LFP{ipart}.(markername).trial{itrial} .* -1;
            end
        end
        
        cfgtemp                           = [];
        cfgtemp.demean                    = config{ipatient}.LFP.baseline;
        cfgtemp.baselinewindow            = config{ipatient}.LFP.baselinewindow.(markername);
        LFP{ipart}.(markername) = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));
        
    end
    
    save(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP', '-v7.3');
else
    MuseStruct = readMuseMarkers(config{ipatient},false);
    MuseStruct = alignMuseMarkers(config{ipatient},MuseStruct,false);
    MuseStruct_concat = concatenateMuseMarkers(config{ipatient},MuseStruct,false); %for sliding time window of morpho
    load(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP');
end

%% compute slow wave morpho 
if do_morpho
    datasavedir = config{ipatient}.datasavedir;
    for do_hpfilter = [false 0.15 1]
        for markername = "SlowWave"
            if ~isfield(LFP{ipart},markername)
                continue
            end
            if isempty(LFP{ipart}.(markername))
                continue
            end
            
            %lp filter
            fprintf('Low pass filtering before computing SW morpho\n');
            cfgtemp                             = [];
            cfgtemp.lpfilter                    = 'yes';
            cfgtemp.lpfreq                      = 30;
            cfgtemp.lpfilttype                  = 'fir';
            LFP{ipart}.(markername)   = ft_preprocessing(cfgtemp,LFP{ipart}.(markername));
            
            i_filt = "raw";
            if do_hpfilter >0
                fprintf('High pass filtering before computing SW morpho\n');
                cfgtemp                 = [];
                cfgtemp.hpfilter        = 'yes';
                cfgtemp.hpfilttype      = 'but';
                cfgtemp.hpinstabilityfix= 'reduce';
                
                if do_hpfilter == 0.15
                    i_filt = "hpfilt_0_15";
                    cfgtemp.hpfreq = 0.15;%0.15 comme Micromed
                elseif do_hpfilter == 1
                    i_filt = "hpfilt_1";
                    cfgtemp.hpfreq = 1;%1 comme filtres anesth
                end
                
                LFP{ipart}.(markername) = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));
            end
            
            %pas besoin de s�lectionner channel, d�j� un seul channel dans data
            %cfgmorpho.morpho.channame = config{ipatient}.morpho.channame.(markername);
            
            fig = figure;hold;
            %compute data and plot results, for each trial
            for itrial = 1:size(LFP{ipart}.(markername).trial,2)
                cfgtemp = [];
                cfgtemp.trials = itrial;
                LFP_onetrial = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
                
                config{ipatient}.morpho.channame = config{ipatient}.align.channel.(markername);
                try
                    [hw, ~, ~,amp] = plot_morpho(config{ipatient}, LFP_onetrial);
                catch
                    warning('cannot find slowwave in trial %d', itrial);
                    hw = nan;
                    amp = nan;
                end
                
                morpho.(i_filt).halfwidth(itrial) = hw;
                morpho.(i_filt).amplitude(itrial) = amp;
                startdir                              = MuseStruct{1}{LFP_onetrial.trialinfo.idir}.starttime;
                morpho.(i_filt).time(itrial)      = startdir + seconds((LFP_onetrial.trialinfo.begsample - LFP_onetrial.trialinfo.offset) / LFP_onetrial.fsample);
            end
            %remove text to make the figure readable
            delete(findall(gcf,'type','text'));
            
            %print to file
            %fname = fullfile(config{ipatient}.imagesavedir,'sw_morpho',sprintf('%smorpho_%s',config{ipatient}.prefix,markername));
            fname = fullfile(config{ipatient}.imagesavedir,'..','morpho_eachrat',convertStringsToChars(i_filt),sprintf('%smorpho_%s_%s',config{ipatient}.prefix,markername,i_filt));
            dtx_savefigure(fig,fname,'pdf','png','close');
            
            
            %% average over a sliding time window
            morpho.(i_filt).statsovertime.winsize = config{ipatient}.morpho.winsize;
            morpho.(i_filt).statsovertime.winstep = config{ipatient}.morpho.winstep;
            morpho.(i_filt).statsovertime.time_unit = 'seconds';
            
            %first window
            twin_start = MuseStruct_concat{1}.markers.Analysis_Start.clock(1);
            twin_end   = twin_start + seconds(config{ipatient}.morpho.winsize);
            i_window   = 0;
            
            %go trough each window
            while twin_end < MuseStruct_concat{1}.markers.Analysis_End.clock
                
                i_window   = i_window+1;
                
                %keep window time
                if isfield(config{ipatient}, 'injectiontime')
                    %time relative to injection
                    morpho.(i_filt).statsovertime.starttime_orig(i_window) = twin_start;
                    morpho.(i_filt).statsovertime.endtime_orig(i_window)   = twin_end;
                    morpho.(i_filt).statsovertime.starttime(i_window)      = twin_start - config{ipatient}.injectiontime;
                    morpho.(i_filt).statsovertime.endtime(i_window)        = twin_end - config{ipatient}.injectiontime;
                else
                    morpho.(i_filt).statsovertime.starttime(i_window)      = twin_start;
                    morpho.(i_filt).statsovertime.endtime(i_window)        = twin_end;
                end
                
                for iparam= ["halfwidth", "amplitude"]
                    %mean and std of param
                    idx = morpho.(i_filt).time > twin_start & morpho.(i_filt).time < twin_end;
                    morpho.(i_filt).statsovertime.(sprintf('%s_nb_sw',iparam))(i_window) = sum(idx);
                    morpho.(i_filt).statsovertime.(sprintf('%s_mean',iparam))(i_window) = nanmean(morpho.(i_filt).(iparam)(idx));
                    morpho.(i_filt).statsovertime.(sprintf('%s_std',iparam))(i_window)  = nanstd(morpho.(i_filt).(iparam)(idx));
                end
                
                %remove values if this is part of missing file
                keepwindow = true;
                if ~isempty(config{ipatient}.missingdata)
                    if twin_start > config{ipatient}.missingdata(1) && twin_start > config{ipatient}.missingdata(2)
                        keepwindow = true;
                    elseif twin_end < config{ipatient}.missingdata(1) && twin_end < config{ipatient}.missingdata(2)
                        keepwindow = true;
                    else
                        keepwindow = false;
                    end
                end
                if ~keepwindow
                    for iparam= ["halfwidth", "amplitude"]
                        %mean and std of param
                        morpho.(i_filt).statsovertime.(sprintf('%s_nb_sw',iparam))(i_window) = nan;
                        morpho.(i_filt).statsovertime.(sprintf('%s_mean',iparam))(i_window) = nan;
                        morpho.(i_filt).statsovertime.(sprintf('%s_std',iparam))(i_window)  = nan;
                    end
                end
                
                %go to the next window
                twin_start = twin_start + seconds(config{ipatient}.morpho.winstep);
                twin_end   = twin_end   + seconds(config{ipatient}.morpho.winstep);
                
            end %while
        end
    end %i_hpfilter
    
    %save computed morpho
    fprintf('save morpho values to %s\n',fullfile(datasavedir,sprintf('%sslowwave_morpho.mat',config{ipatient}.prefix)));
    save(fullfile(datasavedir,sprintf('%sslowwave_morpho.mat',config{ipatient}.prefix)), 'morpho', '-v7.3');
end


end
