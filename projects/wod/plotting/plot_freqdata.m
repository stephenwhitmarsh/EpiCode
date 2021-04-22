if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    
    
    
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
    
end

ft_defaults

config = wod_setparams;
rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\';
statsavedir=fullfile(config{4}.statsavedir);
freqstatpath= fullfile(statsavedir,'freq_data');
comp_woddir= fullfile(statsavedir,'comp_wod');

analysis_names={'timefreq_wod','timefreq_wod_blcorrected'};

temp= load(fullfile(freqstatpath,'peak_time.mat'));
all_data.time= temp.time;
clear temp
temp= load(fullfile(freqstatpath,'peak_value.mat'));
all_data.value= temp.value;
clear temp


all_data.Depth= table2array(depthelectrodewod);
all_data.Depth(1,:)=[];
all_data.Depth([9:end],:)=[];

bands={"HF", "MF", "MLF", "LF"};

%% Plot frequency times and values
for idata= 1:size(analysis_names,2)
    for ival= ["time" "value"]
        for iband= 1,size(bands,2)
            chan_list = fieldnames(all_data.(ival).(bands{iband}).(analysis_names{idata}));
            for ichan= 1: numel(chan_list)
                chan_name=chan_list{ichan};
                
                data_plot.(bands{iband}).(ival)(:,ichan)=all_data.(ival).(bands{iband}).(analysis_names{idata}).(chan_name);
                
            end %ichan
            
            
            data_plot.(bands{iband}).(ival)(:,end)=[];
            t= 1:size(data_plot.(bands{iband}).(ival),2);
            t_interp= linspace(t(1),t(end),60);
            
            A=[];
            B=[];
            for iwod= 1:size(data_plot.(bands{iband}).(ival),1)
                %remove nans before interpolation
                sel = ~isnan(data_plot.(bands{iband}).(ival)(iwod,:));
                A(iwod,:)= pchip(t(sel),data_plot.(bands{iband}).(ival)(iwod,sel),t_interp);
                
                %replace data by nan to have the same data size for each
                %wod
                i1 = find(~isnan(data_plot.(bands{iband}).(ival)(iwod,:)), 1, 'first');
                i2 = find(~isnan(data_plot.(bands{iband}).(ival)(iwod,:)), 1, 'last');
                t_sel(1) = t(i1);
                t_sel(2) = t(i2);
                
                A(iwod, t_interp<t_sel(1) | t_interp>t_sel(2)) = nan;
                
                %Interpolate depth electrode
                B(:,iwod)= pchip(t,all_data.Depth(:,iwod),t_interp);
            end %iwod
            
            A=A.';
            fig= figure;
            for ichan= 1:size(A,1)
                x_med(ichan,1)=nanmedian(A(ichan,:));
                y_med(ichan,1)=nanmedian(B(ichan,:));
            end %ichan
            
            
            for iwod= 1:size(A,2)
                sgtitle(sprintf('%s %s %s', (bands{iband}),(ival),analysis_names{idata}), 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
                
                %plot all times
                x=A(:,iwod);
                y=B(:,iwod);
                plot(x,y,'Color',[0 0 0 0.5],'LineWidth',0.5);
                hold on
                clear x y
                %plot median on top
                
                plot(x_med,y_med,'Color','k','LineWidth',2);
                
                unit=[];
                
                if ival== "time"
                    unit= "(s)";
                else
                    unit= "(mV²)";
                end
                
                xlabel([(ival),(unit)]);
                ylabel('Depth (µm)');
                
            end %iwod
            
            imagesavedir= fullfile(config{4}.imagesavedir,'freqdata_plots');
            
            if ~isfolder(imagesavedir)
                mkdir(imagesavedir);
            end
            
            fname= fullfile(imagesavedir,sprintf('%s_%s_%s',analysis_names{idata},ival,bands{iband}));
            dtx_savefigure(fig,fname,'pdf','png','close');
        end %iband
        clear A B y_med x_med x y
    end %ival
end %idata

%% Plot normalized frequency times

for idata= 1:size(analysis_names,2)
    for ival= ["time"]
        for iband= 1:size(bands,2)
            chan_list = fieldnames(all_data.(ival).(bands{iband}).(analysis_names{idata}));
            for ichan= 1: numel(chan_list)
                chan_name=chan_list{ichan};
                
                data_plot.(bands{iband}).(ival)(:,ichan)=all_data.(ival).(bands{iband}).(analysis_names{idata}).(chan_name);
                
            end %ichan
            data_plot.(bands{iband}).(ival)(:,end)=[];
            
            for iwod= 1:size(data_plot.(bands{iband}).(ival),1)
                data_plot.(bands{iband}).time_norm(iwod,:)= data_plot.(bands{iband}).(ival)(iwod,:)-min(data_plot.(bands{iband}).(ival)(iwod,:));
            end
            
            %all_data.Depth=table2array(all_data.Depth);
            t= 1:size(data_plot.(bands{iband}).time_norm,2);
            t_interp= linspace(t(1),t(end),60);
            
            A=[];
            B=[];
            for iwod= 1:size(data_plot.(bands{iband}).time_norm,1)
                A(iwod,:)= abs(pchip(t,data_plot.(bands{iband}).time_norm(iwod,:),t_interp));
                %Interpolate depth electrode
                B(:,iwod)= abs(pchip(t,all_data.Depth(:,iwod),t_interp));
            end %iwod
            
            A=A.';
            fig= figure;
            for ichan= 1:size(A,1)
                x_med(ichan,1)=nanmedian(A(ichan,:));
                y_med(ichan,1)=nanmedian(B(ichan,:));
            end %ichan
            
            
            for iwod= 1:size(A,2)
                sgtitle([(bands{iband}),(ival),analysis_names{idata}], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
                
                %plot all times
                x=A(:,iwod);
                y=B(:,iwod);
                plot(x,y,'Color',[0 0 0 0.5],'LineWidth',0.5);
                hold on
                clear x y
                %plot median on top
                
                plot(x_med,y_med,'Color','k','LineWidth',2);
                
                unit=[];
                
                if ival== "time"
                    unit= "(s)";
                else
                    unit= "(mV²)";
                end
                
                xlabel([(ival),(unit)]);
                ylabel('Depth (µm)');
                
            end %iwod
            imagesavedir= fullfile(config{4}.imagesavedir,'freqdata_plots');
            
            if ~isfolder(imagesavedir)
                mkdir(imagesavedir);
            end
            
            fname= fullfile(imagesavedir,sprintf('%s_%s_%s',analysis_names{idata},'time_norm',bands{iband}));
            dtx_savefigure(fig,fname,'pdf','png','close');
        end %iband
    end %ival
end %idata



%% Plot frequency times vs. WoD delays

%load Comp structure
temp= load(fullfile(comp_woddir,'comp_wod_toplot.mat'));
Comp_wod= temp.Comp_wod;
clear temp


for idata= 1:size(analysis_names,2)
    for iband= 1:size(bands,2)
        for ival= ["time" "time_norm"]
            fig= figure;
            x_wod= Comp_wod.(bands{iband}).(analysis_names{idata}).(ival)(:,1);
            y_freq=Comp_wod.(bands{iband}).(analysis_names{idata}).(ival)(:,2);
            
            %delete rows with NaN values
            A= [x_wod y_freq];
            A(any(isnan(A),2),:)=[];
            x_wod= A(:,1);
            y_freq= A(:,2);
            clear A
            
            % make y=x curve
            x=[0:1:max(x_wod)];
            y=[0:1:max(x_wod)];
            
            %plot data
            scatter(x_wod,y_freq,'Marker','.');
            xlabel('WoD normalized delay (s)');
            
            %make linear regression of data
            X_wod=[ones(length(x_wod),1),x_wod];
            coef_wod= X_wod\y_freq;
            linregr=X_wod*coef_wod;
            
            unit=[];
            if ival== "time_norm"
                unit="Peak normalized delay (s)";
            else
                unit="Peak delay (s)"
            end
            ylabel(unit)
            hold on
            plot(x,y,'Color','r','LineStyle','--');
            hold on
            
            plot(x_wod,linregr,'Color',[0.8500 0.3250 0.0980],'LineStyle','--')
            legend('Data','y=x',sprintf('%s=%f%s+%f','y',coef_wod(2),'x',coef_wod(1)),'Location','best');
            axis tight
            
            imagesavedir= fullfile(config{4}.imagesavedir,'freqdata_plots','compare_wod');
            
            if ~isfolder(imagesavedir)
                mkdir(imagesavedir);
            end
            
            %save plots
            fname= fullfile(imagesavedir,sprintf('%s_%s_%s',analysis_names{idata},ival,bands{iband}));
            dtx_savefigure(fig,fname,'pdf','png','close');
            
            %store fitting data in structure
            Coef_dir.(analysis_names{idata}).(bands{iband}).(ival)=coef_wod(2);
            save(fullfile(freqstatpath,'coeff_linreg_wod.mat'),'Coef_dir');
        end %ival
    end %iband
end %idata



%% Plot all frequency bands + SD on same fig
bands={"HF", "MF", "MLF", "LF"};
for idata= 1:size(analysis_names,2)
    for ival= "time"
        for iband= 1:size(bands,2)
            chan_list= fieldnames(all_data.(ival).(bands{iband}).(analysis_names{idata}));                                  
            for ichan= 1: numel(chan_list)
                chan_name=chan_list{ichan};
                
                data_plot.(bands{iband}).(ival)(:,ichan)=all_data.(ival).(bands{iband}).(analysis_names{idata}).(chan_name);
                
            end %ichan
            
            
            data_plot.(bands{iband}).(ival)(:,end)=[];
            t= 1:size(data_plot.(bands{iband}).(ival),2);
            t_interp= linspace(t(1),t(end),60);
            
            A=[];
            B=[];
            for iwod= 1:size(data_plot.(bands{iband}).(ival),1)
                %remove nans before interpolation
                sel = ~isnan(data_plot.(bands{iband}).(ival)(iwod,:));
                A(iwod,:)= pchip(t(sel),data_plot.(bands{iband}).(ival)(iwod,sel),t_interp);
                
                %replace data by nan to have the same data size for each
                %wod
                i1 = find(~isnan(data_plot.(bands{iband}).(ival)(iwod,:)), 1, 'first');
                i2 = find(~isnan(data_plot.(bands{iband}).(ival)(iwod,:)), 1, 'last');
                t_sel(1) = t(i1);
                t_sel(2) = t(i2);
                
                A(iwod, t_interp<t_sel(1) | t_interp>t_sel(2)) = nan;
                
                %Interpolate depth electrode
                B(:,iwod)= pchip(t,all_data.Depth(:,iwod),t_interp);
            end %iwod
            
            A=A.';
            
            
            if iband== 1
            fig= figure;
            end
            for ichan= 1:size(A,1)
                x_med(ichan,1)=nanmedian(A(ichan,:));
                y_med(ichan,1)=nanmedian(B(ichan,:));
                std_time(ichan,1)=std(A(ichan,:),'omitnan');
                SD_neg(ichan,1)= x_med(ichan,1)-std_time(ichan,1);
                SD_plus(ichan,1)= x_med(ichan,1)+std_time(ichan,1);
                
            end %ichan
            
            
            sgtitle(sprintf('%s %s %s', bands{iband},(ival),analysis_names{idata}), 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
            
            Color_med={'k','r','g','b'};
            Color_std={[0.5 0.5 0.5],[1 0.5 0.5],[0.5 1 0.6],[0.2 0.5 0.9]};
                
            %plot median
            plot(x_med,y_med,'Color',Color_med{iband},'LineWidth',0.5);
            
            clear x y
            hold on
%             %plot std on top
%             plot(SD_neg,y_med,'Color',Color_std{iband},'LineWidth',0.5);
%             plot(SD_plus,y_med,'Color',Color_std{iband},'LineWidth',0.5);
           
            
            unit="(s)";
            
            xlabel([(ival),(unit)]);
            ylabel('Depth (µm)');
            
            legend(bands{1},bands{2},bands{3},bands{4});
            
            
        end %iband
    clear A B y_med x_med x y
    hold on
    end %ival
    imagesavedir= fullfile(config{4}.imagesavedir,'freqdata_plots','all_bands');
    
    if ~isfolder(imagesavedir)
        mkdir(imagesavedir);
    end
    
    fname= fullfile(imagesavedir,sprintf('%s_%s',analysis_names{idata},ival));
    dtx_savefigure(fig,fname,'pdf','png','close');
end %idata