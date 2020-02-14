function [cfg_patient, data_merged] = dtx_merge_patientslgi1(merge_index,ipatient, MuseStruct, data)
%ipatient : all the patients to merge ie {1, 2, 3}
%data : a structure with all the data{ipart}{imarker} :
%data{ipatient}{ipart}{imarker}
%Musestruct{ipatient}{ipart}{idir}
%'align' electrode is the one of the first patient/eeg
%Script not yet adapted for others markers than SlowWave_R and SlowWave_L

merge_index = merge_index{ipatient};
ipart = 1;
[config, ~] = dtx_setparams_patients_lgi1([]);

if size(data,2)~=size(merge_index,2)
    error('There must be one data field per patient to merge');
end

idata=1;

count_SlowWave_R = 0;
count_SlowWave_L = 0;

for i_pat = merge_index
    for imarker = 1:length(config{i_pat}.LFP.name)
        
        %rename EMG
        for ichannel=1:length(data{i_pat}{ipart}{imarker}.label)
            if strcmp(data{i_pat}{ipart}{imarker}.label{ichannel},config{i_pat}.LFP.emg{imarker})
                data{i_pat}{ipart}{imarker}.label{ichannel} = 'EMG';
            end
        end
        
        %find which marker it is
        for markername = {'SlowWave_R', 'SlowWave_L'}
            
            if strcmp(config{i_pat}.LFP.name{imarker},markername)
                
                if strcmp('SlowWave_R',markername)
                    count_SlowWave_R = count_SlowWave_R+1;
                    %temp = ft_selectdata(data{idata}{ipart}{imarker} 
                    eval(sprintf('data_R_%d = data{idata}{ipart}{imarker};',count_SlowWave_R));
                elseif strcmp('SlowWave_L',markername)
                    count_SlowWave_L = count_SlowWave_L+1;
                    eval(sprintf('data_L_%d = data{idata}{ipart}{imarker};',count_SlowWave_L));
                else
                    error('Marker must be SlowWave_R or SlowWave_L');
                end

            end
            
        end %markername
    end %imarker
    idata=idata+1;
end %ipat


%string for appendata
data_R_string = [];
data_L_string = [];
for i_R = 1:count_SlowWave_R
    data_R_string = [data_R_string, sprintf(',data_R_%d',i_R)];
end

for i_L = 1:count_SlowWave_L
    data_L_string = [data_L_string, sprintf(',data_L_%d',i_L)];
end


%Set marker LFP.name and dat

if count_SlowWave_R > 0 && count_SlowWave_L > 0
    config{ipatient}.LFP.name     = {'SlowWave_R','SlowWave_L'};
    eval(sprintf('data_merged{ipart}{1} = ft_appenddata([]%s);',data_R_string));
    eval(sprintf('data_merged{ipart}{2} = ft_appenddata([]%s);',data_L_string));
    
elseif count_SlowWave_R > 0 && count_SlowWave_L == 0
    config{ipatient}.LFP.name     = {'SlowWave_R'};
    eval(sprintf('data_merged{ipart}{1} = ft_appenddata([]%s);',data_R_string));
    
    
elseif count_SlowWave_R == 0 && count_SlowWave_L > 0
    config{ipatient}.LFP.name     = {'SlowWave_L'};
    eval(sprintf('data_merged{ipart}{1} = ft_appenddata([]%s);',data_L_string));
    
else
    error('Unable to detect good markers SlowWave_R or SlowWave_L');
end


%adapt new config
config{ipatient}.prefix       = [config{ipatient}.prefix(1:12),'-MERGED-'];
config{ipatient}.LFP.emg      = {'EMG', 'EMG'};

cfg_patient = config{ipatient};

end

