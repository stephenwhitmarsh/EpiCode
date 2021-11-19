function [dt, X, Y, label] = plot_hyp_lines(h)

% plot hypnogram lines
X = [];
Y = [];
label = [];
for im = 1 : height(h)

    X = [X, h.starttime(im), h.endtime(im)];
    
    % height in hypnogram is based on order of cfg.hyp.contains
    switch cell2mat(h.hyplabel(im))
        case 'NO_SCORE'
            y = 5;
        case 'PRESLEEP'
            y = 5;            
        case 'POSTSLEEP'
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
    
    if ~isempty(X)
        % if there's a gap, 'fill' with 0
        if h.endtime(im) ~= X(end)
            X = [X, X(end) h.starttime(im)];
            Y = [Y, y,  y];
            label = [label, h.hyplabel(im)];
            label = [label, h.hyplabel(im)];
        end
    end
    
    Y = [Y, y, y];
    label = [label, h.hyplabel(im)];
    label = [label, h.hyplabel(im)];
end

dt = X;
t = convertTo(X,'epochtime','Epoch','2001-01-01','TicksPerSecond', 1);
dt = t - t(1);

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
set(gca,'Ytick', 1 : 5, 'Yticklabels',{'Stage 3','Stage 2','Stage 1','REM','Awake'});