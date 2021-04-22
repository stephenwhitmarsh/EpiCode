if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    addpath \\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB
    addpath \\lexport\iss01.charpier\analyses\wod\fdr_bh
  

  
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
    addpath /network/lustre/iss01/charpier/echanges/scripts-paul/Spike2_vers_MATLAB
    addpath /network/lustre/iss01/charpier/analyses/wod/fdr_bh

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

DC_stats_datadir=fullfile(config{1}.datasavedir,'statistics');

if ~isfolder(DC_stats_datadir)
    mkdir(DC_stats_datadir);
end

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
a=0;
        
for ifilt= ["raw" "filt"]
    for ithr= ["baseline" "max_slope"]
        a=a+1;
        [p,h,stats]=ranksum(Sup_timings.(ifilt).start.(ithr),Dep_timings.(ifilt).start.(ithr));
        Tests.start.(ifilt).(ithr).p=p;
        Tests.start.(ifilt).(ithr).h=h;
        Tests.start.(ifilt).(ithr).stats=stats;
        P_val(a,1)=p;
        clear p h stats
    end %ithr
    
    %Test other data
    for idata= ["peak" "min_slope" "max_slope" "area"]
        [p,h,stats]=ranksum(Sup_timings.(ifilt).(idata),Dep_timings.(ifilt).(idata));
        Tests.(idata).(ifilt).p=p;
        Tests.(idata).(ifilt).h=h;
        Tests.(idata).(ifilt).stats=stats;
        P_val(end+1,1)=p;
        clear p h stats
        
        
    end %idata
    
    % Test between two thresholding methods for same electrode
    A=Sup_timings.raw.start.baseline;
    B=Sup_timings.raw.start.max_slope;
    
    [p h stats]= signrank(A,B);
    clear A B
    
    Tests.start_thr.Sup=p;
    
    A=Dep_timings.raw.start.baseline;
    B=Dep_timings.raw.start.max_slope;
    
    [p h stats]= signrank(A,B);
    clear A B
    
    Tests.start_thr.Dep=p;
end %ifilt

save(fullfile(DC_stats_datadir,'test_betw_DC.mat'),'Tests');




%% Test values with WoD and WoR peak at same 

configscript_lfp='wod_setparams';
config_lfp=eval(configscript_lfp);
detectiondatapath= fullfile(config_lfp{4}.datasavedir,'Detection');

%load LFP timings data

temp= load(fullfile(detectiondatapath,'WoD_data.mat'));
WoD_data=temp.WOD_data;
clear temp
temp= load(fullfile(detectiondatapath,'WoR_data.mat'));
WoR_data=temp.WOR_data;
clear temp
temp= load(fullfile(detectiondatapath,'Depth_electrode.mat'));
Depth=temp.Electrode_depth;
clear temp

%Determine electrode of same depth




