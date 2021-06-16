

function testaxes

w           = 1/3;
htop        = 1/16 * 0.6;
hbottom     = 1/16 * 0.6;
spacehor    = 0.1;
w = w * (1-spacehor);

% first print average LFP and TFR
figure
for ipatient = 1 : 4
    
    for imarker = 1 : 3
    
        for itemp = 1 : 2
            
            row = (itemp-1) * 2 + (ipatient - 1) * 4 + 1;
            s1 = axes('Position', [1/3*(imarker-1) + w*(spacehor/2), 1-(1/16*row), w, htop]);
            title(sprintf('ipatient:%d, imarker:%d, itemp:%d', ipatient, imarker, itemp));

            row = (itemp-1) * 2 + (ipatient - 1) * 4 + 1 + 1;
            s2 = axes('Position', [1/3*(imarker-1) + w*(spacehor/2), 1-(1/16*row) + (1/16-hbottom), w, hbottom]);
            
            disp('');
             
        end
    end
end


figure

% first print average LFP and TFR
for ipatient = 1 : 4
    
    imarker = 1;
    
    for markername = string(cfg_fig.plot.name)
        
        if ~isfield(SpikeDensity{ipatient}{ipart}.sdf_bar, markername)
            imarker = imarker + 1;
            continue
        end
        
        for itemp = 1 : 2
            
                        
            row = (itemp-1) * 2 + (ipatient - 1) * 4 + 1;
            s1 = axes('Position', [1/3*(imarker-1) + w*(spacehor/2), 1-(1/16*row), w, htop]);
            title(sprintf('ipatient:%d, imarker:%d, itemp:%d', ipatient, imarker, itemp));

            row = (itemp-1) * 2 + (ipatient - 1) * 4 + 1 + 1;
            s2 = axes('Position', [1/3*(imarker-1) + w*(spacehor/2), 1-(1/16*row) + (1/16-hbottom), w, hbottom]);
            
            disp('');
             
            
            
            
        end
        imarker = imarker + 1;
    end
end
