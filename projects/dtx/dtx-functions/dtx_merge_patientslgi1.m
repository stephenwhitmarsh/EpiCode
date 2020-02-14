function [cfg_patient,MuseStruct_merged, data_LFP_EMG, data_LFP] = dtx_merge_patientslgi1(merge_index,ipatient, MuseStruct_temp, data)
%ipatient : all the patients to merge ie {1, 2, 3}
%data : a structure with all the data{ipart}{imarker} :
%data{ipatient}{ipart}{imarker}
%Musestruct{ipatient}{ipart}{idir}
%'align' electrode is the one of the first patient/eeg
%Script not yet adapted for others markers than SlowWave_R and SlowWave_L


merge_index = merge_index{ipatient};
ipart = 1;
[config, ~] = dtx_setparams_patients_lgi1([]);


%% Merge MuseStruct : 1 idir per EEG
dir_index = 0;
for i_pat = 1:length(merge_index)
    for idir = 1:length(MuseStruct_temp{i_pat}{ipart})
        MuseStruct_merged{ipart}{idir + dir_index} = MuseStruct_temp{i_pat}{ipart}{idir};
        dir_index = dir_index+1;
    end
end


%% Merge data
if size(data,2)~=size(merge_index,2)
    error('There must be one data field per patient to merge');
end

idata=1;

count_SlowWave_R = 0;
count_SlowWave_L = 0;
count_SlowWave_R_EMG = 0;
count_SlowWave_L_EMG = 0;

for i_pat = merge_index
    for imarker = 1:length(config{i_pat}.LFP.name)
        hasEMG = false;
        
        %search and rename EMG (so if not same EMG name, it can still be
        %appended)
        for ichannel=1:length(data{i_pat}{ipart}{imarker}.label)
            if strcmp(data{i_pat}{ipart}{imarker}.label{ichannel},config{i_pat}.LFP.emg{imarker})
                data{i_pat}{ipart}{imarker}.label{ichannel} = 'EMG';
                hasEMG = true;
            end
        end
        
        %find which marker it is
       
        if strcmp(config{i_pat}.LFP.name{imarker},'SlowWave_R')
            
            if hasEMG
                count_SlowWave_R_EMG = count_SlowWave_R_EMG + 1;
                eval(sprintf('dataEEGEMG_R_%d = data{idata}{ipart}{imarker};',count_SlowWave_R_EMG));
            end
            
            count_SlowWave_R = count_SlowWave_R + 1;
            cfgtemp = [];
            cfgtemp.channel = {'all','-EMG*'};
            eval(sprintf('dataEEG_R_%d = ft_selectdata(cfgtemp,data{idata}{ipart}{imarker});',count_SlowWave_R));
            
        elseif strcmp(config{i_pat}.LFP.name{imarker},'SlowWave_L')
            
            if hasEMG
                count_SlowWave_L_EMG = count_SlowWave_L_EMG +1 ;
                eval(sprintf('dataEEGEMG_L_%d = data{idata}{ipart}{imarker};',count_SlowWave_L_EMG));
            end
            
            count_SlowWave_L = count_SlowWave_L + 1;
            cfgtemp = [];
            cfgtemp.channel = {'all','-EMG*'};
            eval(sprintf('dataEEG_L_%d = ft_selectdata(cfgtemp,data{idata}{ipart}{imarker});',count_SlowWave_L));
        
        end 
        
    end %imarker
    idata=idata+1;
end %ipat


%string for appenddata
dataEEG_R_string = [];
dataEEGEMG_R_string = [];
dataEEG_L_string = [];
dataEEGEMG_L_string = [];

for i_R = 1:count_SlowWave_R
    dataEEG_R_string = [dataEEG_R_string, sprintf(',dataEEG_R_%d',i_R)];
end

for i_R = 1:count_SlowWave_R_EMG
    dataEEGEMG_R_string = [dataEEGEMG_R_string, sprintf(',dataEEGEMG_R_%d',i_R)];
end

for i_L = 1:count_SlowWave_L
    dataEEG_L_string = [dataEEG_L_string, sprintf(',dataEEG_L_%d',i_L)];
end

for i_L = 1:count_SlowWave_L_EMG
    dataEEGEMG_L_string = [dataEEGEMG_L_string, sprintf(',dataEEGEMG_L_%d',i_L)];
end


%Appenddata, and set cfg.LFP.name and cfg.LFP.emg
config{ipatient}.LFP.emg      = {'EMG', 'EMG'};

if count_SlowWave_R > 0 && count_SlowWave_L > 0
    config{ipatient}.LFP.name     = {'SlowWave_R','SlowWave_L'};
    eval(sprintf('data_LFP{ipart}{1} = ft_appenddata([]%s);',dataEEG_R_string));
    eval(sprintf('data_LFP{ipart}{2} = ft_appenddata([]%s);',dataEEG_L_string));
    if count_SlowWave_R_EMG > 0
        eval(sprintf('data_LFP_EMG{ipart}{1} = ft_appenddata([]%s);',dataEEGEMG_R_string));
    else
        config{ipatient}.LFP.emg{1} = 'no';
    end
    if count_SlowWave_L_EMG > 0
        eval(sprintf('data_LFP_EMG{ipart}{2} = ft_appenddata([]%s);',dataEEGEMG_L_string));
    else
        config{ipatient}.LFP.emg{2} = 'no';
    end
    
elseif count_SlowWave_R > 0 && count_SlowWave_L == 0
    config{ipatient}.LFP.name     = {'SlowWave_R'};
    eval(sprintf('data_merged{ipart}{1} = ft_appenddata([]%s);',dataEEG_R_string));
    if count_SlowWave_R_EMG > 0
        eval(sprintf('data_LFP_EMG{ipart}{1} = ft_appenddata([]%s);',dataEEGEMG_R_string));
    else
        config{ipatient}.LFP.emg{1} = 'no';
    end
    
elseif count_SlowWave_R == 0 && count_SlowWave_L > 0
    config{ipatient}.LFP.name     = {'SlowWave_L'};
    eval(sprintf('data_merged{ipart}{1} = ft_appenddata([]%s);',dataEEG_L_string));
    if count_SlowWave_L_EMG > 0
        eval(sprintf('data_LFP_EMG{ipart}{2} = ft_appenddata([]%s);',dataEEGEMG_L_string));
    else
        config{ipatient}.LFP.emg{2} = 'no';
    end
    
else
    error('Unable to detect good markers SlowWave_R or SlowWave_L');
end


%% create new config
config{ipatient}.prefix         = [config{ipatient}.prefix(1:12),'-MERGED-'];
config{ipatient}.ismerged       = true;
cfg_patient = config{ipatient};

end

