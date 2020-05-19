% subplot(3,3,9); hold;
%             temp        = dir(fullfile(cfg.datasavedir,[cfg.prefix,'all_data_',cfg.circus.channel{1}(1:end-2),'_*.ncs']));
%             hdr_fname   = fullfile(temp(1).folder,temp(1).name);
%             hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
% 
% remove MinPeakDistance

for itemp = 1:17
tempsel = squeeze(SpikeRaw_npy{ipart}.template{itemp}(:,SpikeRaw_npy{ipart}.template_maxchan(itemp)+1,:));%+1 because electrodes are zero-based
temptime = ((0:size(SpikeRaw_npy{ipart}.template{itemp},3)-1)/hdr.Fs*1000)';
figure; hold;
% interpolate template
temptime_int = linspace(temptime(1),temptime(end),10000);
tempsel_int = pchip(temptime,tempsel,temptime_int);

if size(tempsel_int,1) > 1
    plot(temptime_int,tempsel_int,'Color', [0.6 0.6 0.6]);
    tempsel_int = nanmean(tempsel_int,1);
end
    
plot(temptime_int,tempsel_int,'k', 'LineWidth', 2);

axis tight
% Find the higher peak :
[Ypos,Xpos] = findpeaks(tempsel_int, temptime_int,'NPeaks',1,'SortStr','descend'); %Npeaks : max nr of peaks/ SortStr : peak sorting : descend = from largest to smallest
% Search first peak near XPos, before and after
[Yneg,Xneg_temp, w, p] = findpeaks(-tempsel_int,temptime_int);%,'NPeaks',2,'SortStr','descend');%,'MinPeakDistance',0.5);
Xneg_temp = Xneg_temp(p>1|w>0.5); %0.5 and 1 empiric. To remove bad peak detection

if length(Xneg_temp) >= 2
    Xneg(1) = max(Xneg_temp(Xneg_temp-Xpos < 0));
    Xneg(2) = min(Xneg_temp(Xneg_temp-Xpos > 0));
    
    plot([Xpos,Xneg(1)], [Yneg(1)*0.3, Yneg(1)*0.3],'-o','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
    plot([Xpos,Xneg(2)],-[Yneg(2)*0.3, Yneg(2)*0.3],'-o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
    
    x = double((Xpos + Xneg(1))/2);
    y = double(Yneg(1)*0.1);
    text(x,y,sprintf('%.0fus',abs(Xpos-Xneg(1))*1000),'HorizontalAlignment','center','FontWeight', 'bold');
    stats.template_pt(itemp) = abs(Xpos-Xneg(1))*1000;
    
    x = double((Xpos + Xneg(2))/2);
    y = double(-Yneg(2)*0.1);
    text(x,y,sprintf('%.0fus',abs(Xpos-Xneg(2))*1000),'HorizontalAlignment','center','FontWeight', 'bold');
    stats.template_tp(itemp)  = abs(Xpos-Xneg(2))*1000;
else
    stats.template_pt(itemp) = NaN;
    stats.template_tp(itemp)  = NaN;
end

xlabel('time')

midline = min(tempsel_int) + (max(tempsel_int) - min(tempsel_int)) / 2 ;
zci  = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
indx = zci(tempsel_int - midline);
indx = indx([ceil(length(indx)/2-0.5), ceil(length(indx)/2+0.5)]);
plot(temptime_int(indx),ones(1,length(indx))*midline,'-o','Color',[0 1 0],'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);

x = double(sum(temptime_int(indx))/length(indx));
y = double(midline*1.1);
text(x,y,sprintf('%.0fus',(temptime_int(indx(2))-temptime_int(indx(1)))*1000),'HorizontalAlignment','center','FontWeight', 'bold');
stats.template_width(itemp)  = (temptime_int(indx(2))-temptime_int(indx(1)))*1000;

title(SpikeRaw_npy{ipart}.label{itemp}, 'Interpreter', 'none');
end