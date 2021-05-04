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


for iprot= 1:size(config,2)
    mat_datapath= fullfile(config{iprot}.datasavedir,'Matlab_struct');
    %% Load data into a structure
    chanlist = ["Vm", "DC"];

    if iprot >= 13
        chanlist= "Vm";
    end
    
    for channame=1:size(chanlist,2)
        temp= load(fullfile(mat_datapath,sprintf('%s_%s.%s',config{iprot}.directorylist{1}{1},chanlist(channame),'mat')));
        dataprot.raw.(config{iprot}.DC.rename{channame})=temp.data_raw;
        clear temp
        temp= load(fullfile(mat_datapath,sprintf('%s_%s_%s.%s',config{iprot}.directorylist{1}{1},chanlist(channame),'filt','mat')));
        dataprot.filt.(config{iprot}.DC.rename{channame})=temp.data_filt;
        clear temp
        temp= load(fullfile(mat_datapath,sprintf('%s_%s_%s.%s',config{iprot}.directorylist{1}{1},chanlist(channame),'events','mat')));
        dataprot.events=temp.CEDStruct;
        clear temp
    end %channame
    
    chanlist=["DC_sup" "DC_dep"];
     if iprot >= 13
        chanlist= "DC_dep";
    end
    %% Calculate slope detect start times for WOD
    for channame= chanlist
        for isignal= ["raw" "filt"]

            %smoothing signal
        dataprot.(isignal).(channame).trial{1}=movmean(dataprot.(isignal).(channame).trial{1},1000); 
            
            
        %calculate slope for raw signal and exatrct min slope and max slope
        %timings
        dataslope.(isignal).(channame)= dataprot.(isignal).(channame);
        dataslope.(isignal).(channame).trial{1}=ft_preproc_derivative(dataslope.(isignal).(channame).trial{1});
        
        %select baseline
        t1=dataprot.events.markers.WoD.synctime-30;
        t2=dataprot.events.markers.WoD.synctime-20;
        t_sel= [t1 t2];
        
        cfgtemp=[];
        cfgtemp.latency= t_sel;
        cfgtemp.trials= 'all';
        cfgtemp.channel= 'all';
        blslope.(isignal).(channame)= ft_selectdata(cfgtemp,dataslope.(isignal).(channame));
        clear t1 t2 t_sel
        
        %calculate threshold
        thr_slopebl= -3*std(blslope.(isignal).(channame).trial{1});
        
        %cut data to center detection on AD
        t1=dataprot.events.markers.WoD.synctime-20;
        t_sel= [t1 dataprot.(isignal).(channame).time{1}(1,end)];
        
        
        cfgtemp=[];
        cfgtemp.latency= t_sel;
        cfgtemp.trials= 'all';
        cfgtemp.channel= 'all';
        AD_slope.(isignal).(channame)= ft_selectdata(cfgtemp,dataslope.(isignal).(channame));
        clear t1 t2 t_sel
            
                
        %find indices where data crossing threshold in negative
        idx_slope= find(AD_slope.(isignal).(channame).trial{1}< thr_slopebl);
        idx_start= min(idx_slope);
        time_startbl= AD_slope.(isignal).(channame).time{1}(1,idx_start);
        
        
        
        
        %in positive
        t1=dataprot.events.markers.VentOn.synctime+2;
        t_sel= [t1 dataprot.(isignal).(channame).time{1}(1,end)];

        
        
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        cfgtemp.trials='all';
        cfgtemp.channel='all';
        WOR_slope.(isignal).(channame)= ft_selectdata(cfgtemp,dataslope.(isignal).(channame));

        clear t1 t_sel
        
        
        idx_slopewor= find(WOR_slope.(isignal).(channame).trial{1}> thr_slopebl);
        idx_start= min(idx_slopewor);
        timewor_startbl= WOR_slope.(isignal).(channame).time{1}(1,idx_start);

        
%         Functions to find intersection points with threshold and/or other curve
%         x1 = AD_slope.(isignal).(channame).time{1};
%         y1 = AD_slope.(isignal).(channame).trial{1};
%         x2 = x1;
%         y2 = ones(size(y1)) .* thr_slopebl;
%         [x0,y0] = intersections(x1,y1,x2,y2);
        
        clear idx_slope idx_start
        
        %find min and max slopes
        [v_minslope t_minslope]= findpeaks(-AD_slope.(isignal).(channame).trial{1},AD_slope.(isignal).(channame).time{1},'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        [v_maxslope t_maxslope]= findpeaks(AD_slope.(isignal).(channame).trial{1},AD_slope.(isignal).(channame).time{1},'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        
        %Calculate threshold
        thr_slope= -0.2*v_minslope;
        %Get start time
        idx_slope= find(AD_slope.(isignal).(channame).trial{1}< thr_slope);
        idx_start= min(idx_slope);
        time_start= AD_slope.(isignal).(channame).time{1}(1,idx_start);
        
        %Create directory for threshold detection
        detect_thr=fullfile(detect_images,'baseline');
        
        if ~isfolder(detect_thr)
            mkdir(detect_thr);
        end
        
        
        
        %plot detection for each channel
%         fig1=figure;
%         plot(AD_slope.(isignal).(channame).time{1},AD_slope.(isignal).(channame).trial{1});
%         hold on
%         plot(blslope.(isignal).(channame).time{1},blslope.(isignal).(channame).trial{1},'Color','r');
%         yline(thr_slopebl);
%         yline(thr_slope,'Color','b');
%         
%         
%         fname1= fullfile(detect_thr,sprintf('%s_prot_%i_%s_%s_slope',config{iprot}.prefix,iprot,channame,isignal));
%         dtx_savefigure(fig1,fname1,'png','pdf','close');
%         
        
        
%         fig2=figure;
%         plot(dataprot.(isignal).(channame).time{1},dataprot.(isignal).(channame).trial{1});
%         xline(t_minslope,'Color','g');
%         xline(min(blslope.(isignal).(channame).time{1}));
%         xline(max(blslope.(isignal).(channame).time{1}));
%         xline(time_startbl,'Color','r');
%         xline(time_start,'Color','b');
%         xline(timewor_startbl,'Color','r');
%         
%         %save figures for visual verification
%         fname2= fullfile(detect_thr,sprintf('%s_prot_%i_%s_%s',config{iprot}.prefix,iprot,channame,isignal));
%         dtx_savefigure(fig2,fname2,'png','pdf','close');
        
        %Calculate Area under curve
        Area_uc=trapz(dataprot.(isignal).(channame).time{1},dataprot.(isignal).(channame).trial{1});
        
        
        
        
        
        %store data in a structure
        DC_timings.(isignal).wod_start.baseline.(channame)(iprot,1)= time_startbl-dataprot.events.markers.VentOff.synctime;
        DC_timings.(isignal).wod_start.max_slope.(channame)(iprot,1)= time_start-dataprot.events.markers.VentOff.synctime;
        DC_timings.(isignal).min_slope.(channame)(iprot,1)=t_minslope-dataprot.events.markers.VentOff.synctime;
        DC_timings.(isignal).max_slope.(channame)(iprot,1)=t_maxslope-dataprot.events.markers.VentOff.synctime;
        DC_timings.(isignal).area.(channame)(iprot,1)= Area_uc;
        
        
        %% Detection of DC shift peak
        
        t1= t_minslope-30;
        t2= t_minslope+30;
        t_sel= [t1 t2];
        
        cfgtemp=[];
        cfgtemp.latency= t_sel;
        cfgtemp.trials= 'all';
        cfgtemp.channel='all';
        dataAD.(isignal).(channame)=ft_selectdata(cfgtemp,dataprot.(isignal).(channame));
        
        [v_peak t_peak]= findpeaks(-dataAD.(isignal).(channame).trial{1},dataAD.(isignal).(channame).time{1},'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        
        %Store timings in DC_timings
        DC_timings.(isignal).wod_peak.(channame)(iprot,1)=t_peak-dataprot.events.markers.VentOff.synctime;
        
        %plot peak detection
        fig3=figure;
        plot(dataprot.(isignal).(channame).time{1},dataprot.(isignal).(channame).trial{1});
        xline(t1);
        xline(t2);
        xline(t_peak,'Color','r');
        
        %Create directory for peak detection
        detect_peak=fullfile(detect_images,'peak');
        
        if ~isfolder(detect_peak)
            mkdir(detect_peak);
        end
        
        %save figures for visual verification
        fname3= fullfile(detect_peak,sprintf('%s_prot_%i_%s_%s_peak',config{iprot}.prefix,iprot,channame,isignal));
        dtx_savefigure(fig3,fname3,'png','pdf','close');
        
        clear v_peak t_peak dataAD

        
        
        %% Detection of rising part
        
        %Detect maximum rising slope
        t= dataslope.(isignal).(channame).time{1};
        t1=t>(dataprot.events.markers.VentOn.synctime);
        t2=t<(dataprot.events.markers.VentOn.synctime+60);
        t_sel=t1 & t2;
        
        [v_peakslope t_peakslope]=findpeaks(dataslope.(isignal).(channame).trial{1}(1,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        
        clear t t1 t2 t_sel
        
        %Detect DC shift peak
        t= dataprot.(isignal).(channame).time{1};
        t1=t>(dataprot.events.markers.VentOn.synctime);
        t2=t<(dataprot.events.markers.VentOn.synctime+100);
        t_sel=t1 & t2;
        
        [v_peakrising t_peakrising]=findpeaks(dataprot.(isignal).(channame).trial{1}(1,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');

        clear t t1 t2 t_sel
        
        %store timings in structure
        DC_timings.(isignal).wor_slope.(channame)(iprot,1)=t_peakslope-dataprot.events.markers.VentOn.synctime;
        DC_timings.(isignal).wor_repol.(channame)(iprot,1)=t_peakrising-dataprot.events.markers.VentOn.synctime;
        
        %plot for visual control
        fig_wor=figure;
        
        plot(dataprot.(isignal).(channame).time{1},dataprot.(isignal).(channame).trial{1});
        hold on 
        xline(dataprot.events.markers.VentOn.synctime,'Color','k');
        xline(t_peakslope,'Color','r');
        xline(t_peakrising,'Color','b');        
        
        legend('Voltage (mV)','Vent On','Maximum slope','Peak rising phase','Location','eastoutside');
        
        detect_rising= fullfile(detect_images,'Rising');
        
        if ~isfolder(detect_rising)
            mkdir(detect_rising);
        end
        
        fname4=fullfile(detect_rising,sprintf('%s_prot_%i_%s_%s_rising',config{iprot}.prefix,iprot,channame,isignal));
        dtx_savefigure(fig_wor,fname4,'png', 'pdf','close');

        
        
        end %isgnal
        end %channame
    
end %iprot
save(fullfile(config{1}.datasavedir,'DC_timings.mat'),'DC_timings');


