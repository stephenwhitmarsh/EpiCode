addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));

addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
ft_defaults

addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML
CEDS64LoadLib('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML');

%% paramètres :
datapath{1}    = 'Z:\echanges\paul.baudin\Neurones Manon\P_21082020.smrx';
datapath{2}    = 'Z:\echanges\paul.baudin\Neurones Manon\P_25082020.smrx';
channame    = {'Vm','Vm'};
new_chan_number = {15,15};
new_chan_name   = 'Vm_detrend';
neurone_startend = {[0 2435], [0 1325]};

for ineuron = 1:size(datapath,2)
%     fprintf('%s\n',datapath{ineuron}); %
% end
    %% lire le canal d'intérêt
    cfgtemp             = [];
    cfgtemp.datapath    = datapath{ineuron};
    cfgtemp.channame    = channame{ineuron};
    data                = readCEDcontinous(cfgtemp);
    
    %ouvrir le fichier Spike2
    fid         = CEDS64Open(datapath,0);
    if fid <=0
        error('Error while loading %s. It can be : \n- file already loaded in MATLAB \n- file already opened in Spike2 \n- wrong file path', datapath{ineuron});
    end
    
    %voir les données chargées
    % plot(data.time{1}, data.trial{1});
    
    %% corriger la tendance linéaire avec Fieldtrip
    
    %sélectionner la période choisie
    cfgtemp                         = [];
    cfgtemp.latency                 = neurone_startend{ineuron};%[1500 data.time{1}(end)];
    data_neurone                    = ft_selectdata(cfgtemp,data);
    % plot(data_neurone.time{1}, data_neurone.trial{1});
    
    %remove linear trend
    cfgtemp                         = [];
    cfgtemp.detrend                 = 'yes';
    data_detrend                    = ft_preprocessing(cfgtemp,data_neurone);
    
    %vérifier le résultat
    % plot(data_detrend.time{1}, data_detrend.trial{1});
    
    %convert time from seconds to ticks
    start_time = CEDS64SecsToTicks(fid,data_detrend.time{1}(1));
    waveform   = int16(data_detrend.trial{1}*100); %*100 pour garder une précision des integer.
    
    %écrire le nouveau canal dans le fichier Spike2, et fermer le fichier
    i64Div = 1.0/(data.fsample*CEDS64TimeBase(fid)); %trouver le nombre de divisions entre 2 points
    check = CEDS64SetWaveChan(fid,new_chan_number,i64Div,1,data.fsample); %créer un nouveau canal
    check = CEDS64ChanTitle(fid, new_chan_number, new_chan_name); %renommer le nouveau canal
    check = CEDS64WriteWave(fid,new_chan_number,waveform,start_time); %écrire les données corrigées dans ce canal
    CEDS64Close(fid); %ferme le fichier spike2
    
    unloadlibrary ceds64int; 
end
