if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    addpath \\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
    addpath /network/lustre/iss01/charpier/echanges/scripts-paul/Spike2_vers_MATLAB
    
end

ft_defaults

config = DC_setparams;

detect_images= fullfile(config{1}.imagesavedir,'Start_detect');

if ~isfolder(detect_images)
    mkdir(detect_images);
end

temp=load(fullfile(config{1}.datasavedir,'DC_timings.mat'));
DC_timings=temp.DC_timings;
clear temp

chanlist=["DC_sup" "DC_dep"];

%% Separate data from depth to superficial for raw signal
for ifilt= ["raw" "filt"]
    for channame= chanlist
        for idata= ["start" "min_slope" "max_slope" "peak" "area"]
            
            if idata== "start"
                data.(ifilt).(channame).(idata).baseline= DC_timings.(ifilt).(idata).baseline.(channame);
                data.(ifilt).(channame).(idata).max_slope= DC_timings.(ifilt).(idata).max_slope.(channame);
            else
                data.(ifilt).(channame).(idata)= DC_timings.(ifilt).(idata).(channame);
            end
            
        end %idata
    end %channame
    
    Sup_timings.(ifilt)=data.(ifilt).DC_sup;
    Dep_timings.(ifilt)=data.(ifilt).DC_dep;
end %ifilt
clear data DC_timings

%% Test between depths

%Test start times for both thresholding methods

for ifilt= ["raw" "filt"]
    for ithr= ["baseline" "max_slope"]
        [p,h,stats]=ranksum(Sup_timings.(ifilt).start.(ithr),Dep_timings.(ifilt).start.(ithr));
        Tests.start.(ifilt).(ithr).p=p;
        Tests.start.(ifilt).(ithr).h=h;
        Tests.start.(ifilt).(ithr).stats=stats;
        clear p h stats
    end %ithr
    
    %Test other data
    
    for idata= ["peak" "min_slope" "max_slope" "area"]
        
        [p,h,stats]=ranksum(Sup_timings.(ifilt).(idata),Dep_timings.(ifilt).(idata));
        Tests.(idata).(ifilt).p=p;
        Tests.(idata).(ifilt).h=h;
        Tests.(idata).(ifilt).stats=stats;
        clear p h stats
        
    end %idata
end %ifilt