%% Create filename list
clear RLfilename fout_table fout_data

% loop over patient
ifname = 1;
electrodetype = [];

for ipatient = 1:length(vect_pat)
    
    % loop over sleep vs wake
    for iperiod = 1:2
        
        % loop over macro vs micro
        for electype = 1:2
            
            for elec = 1 : nrelectrodes(electype,iperiod,ipatient)
                
                % if in sleep
                if iperiod == 1
                    
                    % if different patterns
                    if diff_pat(ipatient)
                        for ipart = 1:2
                            RLfilename{ifname}      = join([num2str(vect_pat(ipatient)),'-',p_som{ipart},'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'.rhfe'],'');
                            fout_table{ifname}      = join([num2str(vect_pat(ipatient)),'-',p_som{ipart},'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'-table'],'');
                            fout_data{ifname}       = join([num2str(vect_pat(ipatient)),'-',p_som{ipart},'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'-data'],'');
                            electrodetype(ifname)   = electype;
                            ifname                  = ifname + 1;
                        end
                    % if NOT different patterns        
                    else
                        RLfilename{ifname}      = join([num2str(vect_pat(ipatient)),'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'.rhfe'],'');
                        fout_table{ifname}      = join([num2str(vect_pat(ipatient)),'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'-table'],'');
                        fout_data{ifname}       = join([num2str(vect_pat(ipatient)),'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'-data'],'');
                        electrodetype(ifname)   = electype;
                        ifname                  = ifname + 1;
                    end
                
                % if awake
                else
                    RLfilename{ifname}      = join([num2str(vect_pat(ipatient)),'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'.rhfe'],'');
                    fout_table{ifname}      = join([num2str(vect_pat(ipatient)),'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'-table'],'');
                    fout_data{ifname}       = join([num2str(vect_pat(ipatient)),'-',period{iperiod},'-',electrode{electype},'-',num2str(elec),'-data'],'');
                    electrodetype(ifname)   = electype;
                    ifname                  = ifname + 1;
                end
            end
        end
    end
end


%% Descriptives and creation of fieldtrip signal file
%% Already calculated, go down to load data

for ifile = 1:length(RLfilename)
    
    fname = fullfile(datadir,RLfilename{ifile});
    
    if exist(fname,'file')
        
        % Load data
        dat_rl      = load(fname,'-mat'); %
        
        % Extract header
        f           = fields(dat_rl);                               % structure du fichier
        mont        = f{find(~strcmp(f,'st_FileData'))};            % nom du montage
        temp        = dat_rl.(f{find(~strcmp(f,'st_FileData'))});   % donnes sur les evenements
        trials      = temp.v_Intervals;                             % partie du fichier
        
        Fs          = dat_rl.st_FileData.v_SampleRate;

        label       = dat_rl.st_FileData.v_Labels;                  %name of the montage
        
        % extraction par evenement
        evenements      = temp.st_HFOInfo.m_EvtLims;                            % absolute limit of event (number of sample)
        evenements_rel  = temp.st_HFOInfo.m_Rel2IntLims;                        % relative limit of event (number of sample)
        debut           = temp.s_StartSample;                                   %
        Fs              = temp.st_HFOInfo.s_Sampling;                           % sampling rate
        
        % Create fieldtrip structure
        dat                     = [];
        dat.label               = label;
        dat.trialinfo           = temp.st_HFOInfo.v_EvType; % name of the event (Ripples,FastRipple,..)
        dat.trialinfo(:,2:3)    = evenements_rel;
        
        for itrial = 1 : size(trials,1)
            dat.trial{itrial}   = trials{itrial}';
            dat.time{itrial}    = linspace(-length(trials{itrial}) / Fs / 2, length(trials{itrial}) / Fs / 2, length(trials{itrial}) );
        end
        
        fout = fullfile(analysisdatadir,fout_data{ifile});
        fprintf('Saving FieldTrip data to: %s\n',fout)
        save(fout,'dat')
        
        % create table with results
        tbl             = table;
        tbl.Eventnr     = (1:length(evenements(:,1)))';
        tbl.Eventtype   = temp.st_HFOInfo.v_EvType;
        tbl.EventStart  = ((debut/Fs)+(evenements(:,1)/Fs));
        tbl.Duration    = (evenements(:,2)-evenements(:,1))/(Fs/1000);    % duratin of the event
        
        % max min analysis
        duree_select    = 0.5;   % duration to calculate max and min (in second)
        
        % part of signal for extremum analysis
        for i=1:size(trials,1)
            
            if evenements_rel(i,1)+duree_select*Fs <= length(trials{i})
                trial_select    = trials{i}(evenements_rel(i,1) : (evenements_rel(i,1)+(duree_select*Fs)));
            else
                trial_select    = trials{i}(evenements_rel(i,1) : end);
            end
            
            [M, index]          = max(trial_select);
            index_s             = index/Fs;
            tbl.Max_mV(i)       = M;
            tbl.Max_s(i)        = index_s;
            
            [m, index]          = min(trial_select);
            index_s             = index/Fs;
            tbl.Min_mV(i)       = M;
            tbl.Min_s(i)        = index_s;
            tbl.ampl_debut(i)   = trials{i}(evenements_rel(i,1));
            
            tbl.typename = categorical (tbl.Eventtype,[1:24],{'None','Gamma','Ripple','FastRipple','Spike','Artifact','Ripple_Multi','Ripple_HFA','Fast_on_Ripple','Plat_iso','Plat_OL','HFA_iso','HFA_OL','cSharp_iso','cSharp_OL','eSpike_iso','eSpike_OL','HighRipple','pattern17','pattern18','pattern19','pattern20','pattern21','pattern22'});
        end
        
        %% save data
        fout = fullfile(analysisdatadir,fout_table{ifile});
        fprintf('Saving table to:          %s\n',fout)      
        save(fout,'tbl')
    else
        fprintf('File does not exist:      %s\n',fname)
    end
    
end
disp('Done!');


% collect data for plotting
eventcodes = [7, 3, 4, 9];
for ifile = 1 : length(RLfilename)
    
    temp            = load([fullfile(analysisdatadir,fout_data{ifile}),'.mat']);
    dat             = temp.dat;
    
    cfg             = [];
    cfg.resamplefs  = 1000;
    dat             = ft_resampledata(cfg,dat);
    

    for ievent = 1:4
        cfg                 = [];
        cfg.trials          = find(dat.trialinfo == eventcodes(ievent));
        if ~isempty(cfg.trials)
            eventdat{ifile,ievent}    = ft_selectdata(cfg,dat);
            eventdat{ifile,ievent}.label{1} = 'cleared'; % to allow averaging over differnt electrodes
        end
    end
end

% save data
fout = fullfile(analysisdatadir,'eventdat');
fprintf('Saving data to:  %s\n','eventdat');
save(fout,'eventdat','-v7.3')


%% load data

fout = fullfile(analysisdatadir,'eventdat');
load(fout,'eventdat')

hasdat = ind2sub(size(eventdat),~cellfun('isempty', eventdat));


% combine all data
cfg = [];
cfg.keepsampleinfo = 'no';

ER_ripple_multi_macro        = ft_appenddata(cfg,eventdat{hasdat(:,1) & electrodetype' == 1,1});
ER_ripple_single_macro       = ft_appenddata(cfg,eventdat{hasdat(:,2) & electrodetype' == 1,2});
ER_fastripple_multi_macro    = ft_appenddata(cfg,eventdat{hasdat(:,4) & electrodetype' == 1,4});
ER_fastripple_single_macro   = ft_appenddata(cfg,eventdat{hasdat(:,3) & electrodetype' == 1,3});
ER_ripple_multi_micro        = ft_appenddata(cfg,eventdat{hasdat(:,1) & electrodetype' == 2,1});
ER_ripple_single_micro       = ft_appenddata(cfg,eventdat{hasdat(:,2) & electrodetype' == 2,2});
ER_fastripple_multi_micro    = ft_appenddata(cfg,eventdat{hasdat(:,4) & electrodetype' == 2,4});
ER_fastripple_single_micro   = ft_appenddata(cfg,eventdat{hasdat(:,3) & electrodetype' == 2,3});

