function Intra_overdraw(configscript)

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(fileparts(scriptpath)))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects',filesep, 'wod']))
addpath (genpath([epicodepath,'projects',filesep, 'dtx']))
addpath (genpath([epicodepath,'projects',filesep, 'dtx',filesep,'dtx_functions']))
addpath (genpath([epicodepath,'projects',filesep, 'wod',filesep,'wod_functions']))
addpath (genpath([epicodepath,'projects',filesep, 'wod',filesep,'Intra']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    addpath \\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
    addpath /network/lustre/iss01/charpier/echanges/scripts-paul/Spike2_vers_MATLAB
    
end

ft_defaults

config = eval(configscript);


chanlist= ["Vm","EEG-S1-L"]; %chanlist Intra


%% detect peak WOD on ECoG and make it t0
fig_wod=figure;
for iprot=1:size(config,2)
    
    %loading EEG channel
    matlabpath=fullfile(config{iprot}.rawdir,'matlab_structures');
    temp=load(fullfile(matlabpath,sprintf('%s_%s.mat',config{iprot}.directorylist{1}{1},chanlist(2))));
    data=temp.data;
    clear temp
    temp=load(fullfile(matlabpath,sprintf('%s_events.mat',config{iprot}.directorylist{1}{1})));
    CEDStruct=temp.CEDStruct;
    clear temp
    
    
    %search time window
    t1=CEDStruct.markers.WoD.synctime+config{iprot}.Intra.wod_toisearch(1);
    t2=CEDStruct.markers.WoD.synctime+config{iprot}.Intra.wod_toisearch(2);
    t_sel= [t1 t2];
    
    %cut data to search for peak
    cfgtemp=[];
    cfgtemp.latency=t_sel;
    WOD_cut=ft_selectdata(cfgtemp,data);
    
    %find peak
    [v_peak_wod t_peak_wod]= findpeaks(-WOD_cut.trial{1},WOD_cut.time{1},'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
    time_peak.(sprintf('prot_%i',iprot))=t_peak_wod;
    clear data t1 t2 t_sel
    
    %% Overdraw channels aligned on WoD peak
    
    
    for channame= chanlist
        %loading channels
        temp=load(fullfile(matlabpath,sprintf('%s_%s.mat',config{iprot}.directorylist{1}{1},channame)));
        data=temp.data;
        clear temp
        
        % plotting window
        t1= CEDStruct.markers.WoD.synctime-130;
        t2= CEDStruct.markers.WoD.synctime+80;
        t_sel=[t1 t2];
        
        
        %cut data to search for plotting
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        data_plot=ft_selectdata(cfgtemp,data);
        
        %plotting
        if channame=="EEG-S1-L"
            data_plot.trial{1}=(data_plot.trial{1}*10)+20;
        end
        
        color={[0.5 0.5 0.5 0.5], [0 0 0 0.5],[0.8500 0.3250 0.0980]};
        
        x_wod=data_plot.time{1}-time_peak.(sprintf('prot_%i',iprot));
%         y_wod=ones(1,size(data_plot.trial{1},2))*config{iprot}.Intra.dep;
%         z_wod=data_plot.trial{1};
        y_wod=data_plot.trial{1};

        
        
        
        if channame=="EEG-S1-L"
            plot(x_wod,y_wod,'Color',color{2});
        elseif channame=="Vm" & config{iprot}.Intra.rename{1}=='Intra_sup'
            plot(x_wod,y_wod,'Color',color{3});
        else
            plot(x_wod,y_wod,'Color',color{1});
        end
        
        
        
        hold on
        xlim([-130 40]);
        xline(0,'Color','r')
    end %channame
    
    clear x_wod y_wod z_wod
    
end %iprot
clear data_plot data
overdraw_imagesdir=fullfile(config{1}.imagesavedir,'overdraw');

if ~isfolder(overdraw_imagesdir)
    mkdir(overdraw_imagesdir);
end

fname=fullfile(overdraw_imagesdir,'depth_intra_wod_aligned');
dtx_savefigure(fig_wod,fname,'png','pdf','close');



%% Overdraw chans aligned on 50% on depolarization
fig_thr=figure;
for iprot=1:size(config,2)
    
    %loading Vm channel
    matlabpath=fullfile(config{iprot}.rawdir,'matlab_structures');
    temp=load(fullfile(matlabpath,sprintf('%s_%s.mat',config{iprot}.directorylist{1}{1},chanlist(1))));
    data=temp.data;
    clear temp
    temp=load(fullfile(matlabpath,sprintf('%s_events.mat',config{iprot}.directorylist{1}{1})));
    CEDStruct=temp.CEDStruct;
    clear temp
    
    %search time window
    t1=CEDStruct.markers.WoD.synctime+config{iprot}.Intra.wod_toisearch(1);
    t2=CEDStruct.markers.WoD.synctime+config{iprot}.Intra.wod_toisearch(2);
    t_sel= [t1 t2];
    
    %cut data to search for peak
    cfgtemp=[];
    cfgtemp.latency=t_sel;
    WOD_cut=ft_selectdata(cfgtemp,data);
    
    %define threshold
    idx_beg=find(WOD_cut.time{1}<=t1+1);
    idx_beg=min(idx_beg);
    idx_depo=find(WOD_cut.time{1}>=CEDStruct.markers.WoD.synctime);
    idx_depo=min(idx_depo);
    v_beg= WOD_cut.trial{1}(1,idx_beg);
    v_depo=WOD_cut.trial{1}(1,idx_depo);
    
    thr_50=v_beg+0.6*abs(v_beg-v_depo);
    
    x1= WOD_cut.time{1};
    y1= WOD_cut.trial{1};
    x2=WOD_cut.time{1};
    y2= ones(1,size(WOD_cut.trial{1},2))*thr_50;
    [x_wodintersect, y_wodintersect] = intersections(x1, y1, x2, y2);
    clear x1 y1 x2 y2
    for channame= chanlist
        %loading channels
        temp=load(fullfile(matlabpath,sprintf('%s_%s.mat',config{iprot}.directorylist{1}{1},channame)));
        data=temp.data;
        clear temp
        
        % plotting window
        t1= CEDStruct.markers.WoD.synctime-60;
        t2= CEDStruct.markers.WoD.synctime+80;
        t_sel=[t1 t2];
        
        
        %cut data to search for plotting
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        data_plot=ft_selectdata(cfgtemp,data);
        
        %plotting
        if channame=="EEG-S1-L"
            data_plot.trial{1}=data_plot.trial{1}*10;
        end
        
        color={[0.5 0.5 0.5 ], [0 0 0],[0.8500 0.3250 0.0980]};
        
        x_wod=data_plot.time{1}-x_wodintersect(1);
%         y_wod=ones(1,size(data_plot.trial{1},2))*config{iprot}.Intra.dep;
%         z_wod=data_plot.trial{1};
        y_wod=data_plot.trial{1};

        if channame=="EEG-S1-L"
            plot(x_wod,y_wod,'Color',color{2});
        elseif channame=="Vm" & config{iprot}.Intra.rename{1}=='Intra_sup'
            plot(x_wod,y_wod,'Color',color{3});
        else
            plot(x_wod,y_wod,'Color',color{1});
        end
        
        
        
        hold on
        xlim([-60 40]);
        xline(CEDStruct.markers.VentOff.synctime-x_wodintersect(1),'Color','r')
    end %channame
    
    
end % iprot

fname=fullfile(overdraw_imagesdir,'depth_intra_thr_aligned');
dtx_savefigure(fig_thr,fname,'png','pdf','close');
clear fname

%% Overdraw chans aligned at Vent Off

fig_voff=figure;
for iprot=1:size(config,2)
    for channame= chanlist
        %loading channels
        matlabpath=fullfile(config{iprot}.rawdir,'matlab_structures');
        temp=load(fullfile(matlabpath,sprintf('%s_%s.mat',config{iprot}.directorylist{1}{1},channame)));
        data=temp.data;
        clear temp
        temp=load(fullfile(matlabpath,sprintf('%s_events.mat',config{iprot}.directorylist{1}{1})));
        CEDStruct=temp.CEDStruct;
        clear temp
        
        % plotting window
        t1= CEDStruct.markers.VentOff.synctime-80;
        t2= CEDStruct.markers.VentOff.synctime+170;
        t_sel=[t1 t2];
        
        
        %cut data to search for plotting
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        data_plot=ft_selectdata(cfgtemp,data);
        
        %plotting
        if channame=="EEG-S1-L"
            data_plot.trial{1}=(data_plot.trial{1}*10)+20;
        end
        
        color={[0.5 0.5 0.5 0.5], [0 0 0 0.5],[0.8500 0.3250 0.0980]};
        
        x_wod=data_plot.time{1}-CEDStruct.markers.VentOff.synctime;
        %y_wod=ones(1,size(data_plot.trial{1},2))*config{iprot}.Intra.dep;
        %z_wod=data_plot.trial{1};
        y_wod=data_plot.trial{1};
        if channame=="EEG-S1-L"
            c=color{2};
        elseif channame=="Vm" & config{iprot}.Intra.rename{1}=='Intra_sup'
            c=color{3};
        else
            c=color{1};
        end
        
        plot(x_wod,y_wod,'Color',c);
        
        
        hold on
        %xlim([-60 40]);
        xline(0,'Color','r');
    end %channame
end % iprot

fname=fullfile(overdraw_imagesdir,'depth_intra_voff_aligned');
dtx_savefigure(fig_voff,fname,'png','pdf','close');
clear fname


%% Overdraw filtered normalized delta Vm


fig_baselinecorr=figure;
for iprot=1:size(config,2)
    
    %loading Vm channel
    matlabpath=fullfile(config{iprot}.rawdir,'matlab_structures');
    temp=load(fullfile(matlabpath,sprintf('%s_%s.mat',config{iprot}.directorylist{1}{1},chanlist(1))));
    data=temp.data;
    clear temp
    temp=load(fullfile(matlabpath,sprintf('%s_events.mat',config{iprot}.directorylist{1}{1})));
    CEDStruct=temp.CEDStruct;
    clear temp
    
    %select baseline
    t1= CEDStruct.markers.VentOff.synctime-60;
    t2= CEDStruct.markers.VentOff.synctime;
    t_sel=[t1 t2];
    
    cfgtemp=[];
    cfgtemp.latency=t_sel;
    baseline= ft_selectdata(cfgtemp,data);
    clear t1 t2 t_sel
    
    % find AP
    [~, t_ap]=findpeaks(baseline.trial{1},baseline.time{1},'MinPeakProminence',20);
    % replace neighbouring samples of AP timings by NaNs Fsample= 1000Hz
    for i= 1:size(t_ap,2)
        idx_ap=find(baseline.time{1}==t_ap(i));
        remove_ap= [idx_ap-2 :1: idx_ap+2];
        
        %replace +/- 2 samples aound AP by NaNs
        baseline.trial{1}(1,remove_ap)=NaN;
    end %ap

    %interpolate
    baseline.trial{1}=fillmissing(baseline.trial{1},'pchip');
    %average baseline Vm
    baseline_mean=nanmean(baseline.trial{1});
    
    %
    x=data.time{1}-CEDStruct.markers.VentOff.synctime;
    y=data.trial{1}-baseline_mean;
    
    color={[0.5 0.5 0.5],[0.8500 0.3250 0.0980]};
    
    if config{iprot}.Intra.rename{1}=='Intra_sup'
        c=color{2};
    else
        c=color{1};
    end
    plot(x,y,'Color',c);
    hold on
    xlim([-80 180]);
    xline(0,'r');
    ylim([-10 100]);
end % iprot

fname=fullfile(overdraw_imagesdir,'intra_baselinecorr_wod_aligned');
dtx_savefigure(fig_baselinecorr,fname,'png','pdf','close');
clear fname

end %function