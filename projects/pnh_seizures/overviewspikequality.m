
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
datapath = '/network/lustre/iss01/epimicro/patients/raw';
outputdir = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/qualitycheck';



dlist = dir(fullfile(outputdir,'*.mat'));

for ipatient = 1 : size(dlist,1)
    fprintf('Loading %s...',fullfile(dlist(ipatient).folder,dlist(ipatient).name));   
    dat{ipatient} = load(fullfile(dlist(ipatient).folder,dlist(ipatient).name));
    fprintf('OK \n');
end

%% ORGANIZE DATA

clear q m
mad = 6;
all = [];
for ipatient = 1 : size(dlist,1)
    i = 1;
    for idir = 1 : size(dat{ipatient}.quality.dat_mad_cnt,2)
        for imad = 1 : 5
            for ifile = 1 : size(dat{ipatient}.quality.dat_mad_cnt{idir},2)
                q{ipatient}(imad,ifile,idir)    = dat{ipatient}.quality.dat_mad_cnt{idir}{ifile}(imad);
                m{ipatient}(ifile,idir)         = dat{ipatient}.quality.dat_mad{idir}{ifile};
            end
        end
       
    end
    all = [all m{ipatient}(:)'];
end

%% 
% % 
% figure; hist(all,1000)
% size(find(all > 300),2)
% figure; hist(log(all),100)

thresh = 12;

% thresh = 10;
% remove those that cross threshold
for ipatient = 1 : size(dlist,1)
    mask{ipatient} = m{ipatient} > thresh;
end

cmax = -inf;
cmin = inf;
mmax = -inf;
mmin = inf;
mad = 6;
% calculate max/min etc.
for ipatient = 1 : size(dlist,1)
    

        qsel = squeeze(q{ipatient}(mad-3,:,:));
        qsel(mask{ipatient}) = nan;
    
        if max(max(max(qsel))) > cmax
            cmax = max(max(max(qsel)));
        end
        if min(min(min(qsel))) < cmin
            cmin = min(min(min(qsel)));
        end
        
        if max(max(max(m{ipatient}))) > mmax
            mmax = max(max(max(m{ipatient})));
        end
        if min(min(min(m{ipatient}))) < mmin
            mmin = min(min(min(m{ipatient})));
        end    
               
end


w = 3;
h = 13;

%% PLOT SPIKE NUMBER
fig = figure;
for ipatient = 1 : size(dlist,1)  
    subplot(h,w,ipatient);    
    qsel = (squeeze(q{ipatient}(mad-3,:,:)) / cmax) * 1000;
    qsel(mask{ipatient}) = nan;    
    imagesc(qsel);
    set(gca,'ytick',[]);
    temp = regexp(dat{ipatient}.quality.fname{i}{1},filesep,'split');
    temp = strrep(temp,'_',' ');
    temp = temp{end-1}(1:5);
    title(temp);
    set(gca,'fontsize',6);

end
colormap parula(1000)

fig = figure;
for ipatient = 1 : size(dlist,1)  
    subplot(h,w,ipatient);    
      imagesc(mask{ipatient});
    set(gca,'ytick',[]);
    temp = regexp(dat{ipatient}.quality.fname{i}{1},filesep,'split');
    temp = strrep(temp,'_',' ');
    temp = temp{end-1}(1:5);
    title(temp);
    set(gca,'fontsize',6);

end
colormap parula(1000)



fname_out = fullfile(outputdir,'images',sprintf('qualitycheck_spikes_mad%d.png',mad));
print(fig, '-dpng', fullfile(fname_out));
        
%% PLOT SPIKE NUMBER (LOG)
fig = figure;
for ipatient = 1 : size(dlist,1)  
    subplot(h,w,ipatient);    
    qsel = (log(squeeze(q{ipatient}(mad-3,:,:))) / log(cmax)) * 1000;
    image(qsel);
    set(gca,'ytick',[]);   
    temp = regexp(dat{ipatient}.quality.fname{i}{1},filesep,'split');
    temp = strrep(temp,'_',' ');
    temp = temp{end-1}(1:5);
    title(temp);
    set(gca,'fontsize',6);

end
colormap parula(1000)
set(gca,'fontsize',6);
fname_out = fullfile(outputdir,'images',sprintf('qualitycheck_spikes_log_mad%d.png',mad));
print(fig, '-dpng', fullfile(fname_out));

%% PLOT MAD
fig = figure;
for ipatient = 1 : size(dlist,1)  
    subplot(h,w,ipatient);
    z = (m{ipatient} ./ mmax) * 1000;
    image(z);
    set(gca,'ytick',[]);
    temp = regexp(dat{ipatient}.quality.fname{i}{1},filesep,'split');
    temp = strrep(temp,'_',' ');
    temp = temp{end-1}(1:5);
    title(temp); 
set(gca,'fontsize',6);

end
colormap parula(200)
set(gca,'fontsize',6);
fname_out = fullfile(outputdir,'images',sprintf('qualitycheck_mad%d.png',mad));
print(fig, '-dpng', fullfile(fname_out));

%% PLOT MAD (LOG)
fig = figure;
for ipatient = 1 : size(dlist,1)  
    subplot(h,w,ipatient);
    z = (log(m{ipatient}) ./ log(mmax));
    z(z == -inf) = 0;
    z = z * 1000;
    image(z);
    set(gca,'ytick',[]);
    temp = regexp(dat{ipatient}.quality.fname{i}{1},filesep,'split');
    temp = strrep(temp,'_',' ');
    temp = temp{end-1}(1:5);
    title(temp);
set(gca,'fontsize',6);    
end
colormap parula(1000)
fname_out = fullfile(outputdir,'images',sprintf('qualitycheck_log_mad%d.png',mad));
print(fig, '-dpng', fullfile(fname_out));






for mad = 4 : 8
    
    for ipatient = 1 : size(dlist,1)
        
        %     q_norm = q{ipatient} ./ max(q{ipatient},[],1);
        %     q_norm = q ;
        
        labels = dat{ipatient}.quality.fname{1};
        for i = 1 : size(labels,2)
            indx = findstr(labels{i},'_m');
            labels{i} = labels{i}(indx+1:end-4);
            labels{i} = strrep(labels{i},'_','');
        end
        
        for i = 1 : size(dat{ipatient}.quality.fname,2)
            temp = regexp(dat{ipatient}.quality.fname{i}{1},filesep,'split');
            temp = strrep(temp,'_',' ');
            fnames{i} = temp{end-1};
        end
        
        % plotting
        
        fig = figure('Numbertitle','off','Name',sprintf('Patient %d',ipatient));
        set(fig,'PaperOrientation','portrait');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        
        title(num2str(ipatient));
               
        % mad
        subplot(4,1,1);
        imagesc(m{ipatient}./max(m{ipatient},[],1))
        yticks(1:size(m{ipatient},1))
        yticklabels(labels);
        xticklabels([]);
        axis tight
        set(gca,'fontsize',4);
        title('MAD, scaled over channels');

        
        subplot(4,1,2);
        imagesc(m{ipatient}./max(m{ipatient},[],2))
        yticks(1:size(m{ipatient},1))        
        yticklabels(labels);
        xticklabels([]);
        axis tight
        set(gca,'fontsize',4);
        title('MAD, scaled over time');

        
        subplot(4,1,3);
        qsel = squeeze(q{ipatient}(mad-3,:,:));      
        
        imagesc(qsel./max(qsel,[],1))
        yticks(1:size(qsel,1))
        yticklabels(labels);
        xticklabels([]);
        axis tight
        set(gca,'fontsize',4);
        title('Scaled over channels');
        
        subplot(4,1,4);
        imagesc(qsel./max(qsel,[],2))
        yticks(1:size(qsel,1))
        yticklabels(labels);
        xticks(1:size(qsel,2))
        xticklabels(fnames);
        xtickangle(90)
        axis tight
        set(gca,'fontsize',4);
        title('Scaled over time');
        
        
        fname_out = fullfile(outputdir,'images',sprintf('qualitycheck_%d_mad_%d.png',ipatient,mad));

        print(fig, '-dpng', fullfile(fname_out));
        
        close all
    end
end
