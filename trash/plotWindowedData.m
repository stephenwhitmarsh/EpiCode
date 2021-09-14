function plotWindowedData(cfg, MuseStruct, dat)

ncols   = 1;
nrows   = size(dat, 2);
cm      = cool(5);
ipart   = ft_getopt(cfg, 'ipart', 1);


i = 1;
clear bad_t1 bad_t2
for idir = 1 : size(MuseStruct, 2)
    for ibad = 1 : size(MuseStruct{idir}.markers.BAD__START__.clock, 2)
        bad_t1(i) = MuseStruct{idir}.markers.BAD__START__.clock(ibad);
        bad_t2(i) = MuseStruct{idir}.markers.BAD__END__.clock(ibad);
        i = i + 1;
    end
end

figure;

for irow = 1 : nrows
    
    dat{irow}.plotart = ft_getopt(dat{irow}, 'plotart', false);
    dat{irow}.offset  = ft_getopt(dat{irow}, 'offset', []);
    
    subplot(nrows, ncols, irow); hold on;
    title(dat{irow}.title);
    
    
    switch dat{irow}.type
        case 'hypnogram'
            
            plot_hyp_lines(dat{irow}.data);
            axis tight;
            plot_hyp_colors(dat{irow}.data, cm, ylim);
            set(gca,'children',flipud(get(gca,'children')))
            time = xlim;
                
        case 'power'
            
            cfg = [];
            cfg.frequency = dat{irow}.frequency;
            cfg.channel = dat{irow}.channel;
            cfg.avgoverfreq = 'yes';
            cfg.avgoverchan = 'yes';
            data = ft_selectdata(cfg, dat{irow}.data);
            time = data.trialinfo.starttime + (data.trialinfo.endtime - data.trialinfo.starttime) / 2;  
            if ~isempty(dat{irow}.offset)
                time = second(time - offset);
            end
            plot(time, log(data.powspctrm), 'k');
            
            if dat{irow}.plotart
                y = ylim;
                for ibad = 1 : size(bad_t1, 2)
                    fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y(1) y(1) y(2) y(2)], 'r', 'edgecolor', 'none', 'facealpha', 0.4);
                end
            end
            
        case 'spike'
            
            time = dat{irow}.data.trialinfo.starttime + (dat{irow}.data.trialinfo.endtime - dat{irow}.data.trialinfo.starttime) / 2;
            if ~isempty(dat{irow}.offset)
                time = second(time - offset);
            end            
            plot(time, dat{irow}.data.(dat{irow}.metric), 'k');
            if dat{irow}.plotart
                y = ylim;
                for ibad = 1 : size(bad_t1, 2)
                    fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y(1) y(1) y(2) y(2)], 'r', 'edgecolor', 'none', 'facealpha', 0.4);
                end
            end

        otherwise
            sprintf('Can not recognized type %s', dat{irow}.type);
            continue
            
    end
    
    if irow == 1
        xleft  = time(1);
        xright = time(2);
    else
        xleft   = max(time(1), xleft);
        xright  = min(time(end), xright);
    end
    xlim([xleft, xright]);
    
end



function plot_hyp_colors(h, cm, y)
% plot hypnogram in colors
for im = 1 : height(h)
    switch h.hyplabel{im}
        case 'NO_SCORE'
            ci = 1;
        case 'REM'
            ci = 1;
        case 'AWAKE'
            ci = 2;
        case 'PHASE_1'
            ci = 3;
        case 'PHASE_2'
            ci = 4;
        case 'PHASE_3'
            ci = 5;
    end
    fill([h.starttime(im), h.endtime(im), h.endtime(im), h.starttime(im)],[y(1) y(1) y(2) y(2)], cm(ci, :) ,'edgecolor', 'none', 'facealpha',0.4);
end

function plot_hyp_lines(h)

% plot hypnogram lines
X = [];
Y = [];
for im = 1 : height(h)
    if ~isempty(X)
        % if there's a gap, 'fill' with 0
        if h.starttime(im) ~= X(end)
            X = [X, X(end) h.starttime(im)];
            Y = [Y, 0,  0];
        end
    end
    X = [X, h.starttime(im), h.endtime(im)];

    % height in hypnogram is based on order of cfg.hyp.contains
    switch cell2mat(h.hyplabel(im))
        case 'NO_SCORE'
            y = 5;
        case 'AWAKE'
            y = 5;
        case 'REM'
            y = 4;
        case 'PHASE_1'
            y = 3;
        case 'PHASE_2'
            y = 2;
        case 'PHASE_3'
            y = 1;
    end
    Y = [Y, y, y];
end

for i = 1 : length(X)-1
    if Y(i) ~= 0 && Y(i+1) ~= 0
        if Y(i) == 4 && Y(i+1) == 4 % REM gets thicker line
            plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k','LineWidth',3);
        else
            plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k');
        end
    end
end
set(gca,'Layer','top');
set(gca,'Ytick', 1 : 5, 'Yticklabels',{'STAGE 3','STAGE 2','STAGE 1','REM','WAKE'});