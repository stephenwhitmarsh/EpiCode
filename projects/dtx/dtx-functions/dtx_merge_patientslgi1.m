function [config, data] = dtx_merge_patientslgi1(ipatient_to_merge,data)
%ipatient : all the patients to merge ie {1, 2, 3}
%data : a structure with all the data{ipart}{imarker} :
%data{ipatient}{ipart}{imarker}
%'align' electrode is the one of the first patient/eeg
%Script not yet adapted for others markers than SlowWave_R and SlowWave_L

ipart = 1;
config = dtx_setparams_patients_lgi1([]);

if size(data,2)~=size(patient,2)
    error('There must be one data field per patient to merge');
end

idata=1;

count_SlowWave_R = 0;
count_SlowWave_L = 0;

for i_pat = ipatient_to_merge
    for imarker = length(config{i_pat}.LFP.name)
        for markername = ['SlowWave_R', 'SlowWave_L']
            
            if strcmp(config{i_pat}.LFP.name{imarker},markername)
                
                %get which marker has this eeg
                if strcmp('SlowWave_R',markername)
                    count_SlowWave_R = count_SlowWave_R+1;
                    eval(sprintf('data_R_%d = data{idata}{ipart}{imarker}',count_SlowWave_R);
                elseif strcmp('SlowWave_R',markername)
                    has_SlowWave_L = 1;
                else
                    error('Marker must be SlowWave_R or SlowWave_L');
                end
                
                %rename EMG
                for ichannel=1:length(data{i_pat}{ipart}{imarker}.label)
                    if strcmp(data{i_pat}{ipart}{imarker}.label{ichannel},config{i_pat}.LFP.emg{imarker})
                        data{i_pat}{ipart}{imarker}.label{ichannel} = 'EMG';
                    end
                end
                                
            else
                error('Error : %s in config.LFP.name (%s) : \nthis script is not yet adapted for other markers than SlowWave_R or SlowWave_L.',...
                    config{i_pat}.prefix(1:end-1),config{i_pat}.LFP.name{imarker})
            end
            
        end %markername
    end %imarker
    idata=idata+1;
end %ipat


%appendata
for i_R = 1:count_SlowWave_R
    data_R_string = [];
    data_R_string = [data_R_string, sprintf(',data_R_%d',i_R);
end

for i_L = 1:count_SlowWave_L
    data_L_string = [];
    data_L_string = [data_L_string, sprintf(',data_R_%d',i_L);
end


%Set marker LFP.name and data
if count_SlowWave_R > 0 && count_SlowWave_L > 0
    config.LFP.name     = {'SlowWave_R','SlowWave_L'};
    eval(sprintf('data{ipart}{1} = ft_appenddata([]%s);',data_R_string)
    eval(sprintf('data{ipart}{2} = ft_appenddata([]%s);',data_L_string)
    
elseif count_SlowWave_R > 0 && count_SlowWave_L = 0
    config.LFP.name     = {'SlowWave_R'};
    eval(sprintf('data{ipart}{1} = ft_appenddata([]%s);',data_R_string)
    
    
elseif count_SlowWave_R = 0 && count_SlowWave_L > 0
    config.LFP.name     = {'SlowWave_L'};
    eval(sprintf('data{ipart}{1} = ft_appenddata([]%s);',data_L_string)
    
else
    error('Unable to detect good markers SlowWave_R or SlowWave_L');
end


%adapt new config
config              = config{ipatient_to_merge{ipart}};
config.prefix       = [config.prefix(1:12),'-MERGED-'];
config.LFP.emg      = {'EMG', 'EMG'};


end

