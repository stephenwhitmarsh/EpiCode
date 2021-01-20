%% script to test different ftypes in msdf on real and simulated psths
clear all,close all

% determine whether to use real data (1) or simulated data (0)
realSpikes = 1;

% specify w for each ftype
w_boxcar  = 100;      % note that if mod(w_boxcar,2)=0, msdf (but not smooth.m) will make w -> w+1
w_Gauss   = 100;
w_exp     = 100;
w_exGauss = [100 100]; % SD of Gaussian and tau of exponential distribution
%% load or generate psth
if realSpikes
  load unitForMLIBTesting
  psth = mpsth(spx.timings,spx.eventtimings(spx.eventmarkers==1),'fr',1);
else
  win  = [-1,+1];                   % psth duration in seconds
  nspx = 50;                        % number of spikes
  t    = (win(1)*1000:win(2)*1000)';
  
  % generate psth, generate random spike train
  psth(:,1) = t;
  psth(:,2) = [ones(nspx,1);zeros(numel(t)-nspx,1)];
  r         = randperm(size(psth,1));
  psth(:,2) = psth(r,2);
end
%% plot sdfs and kernels
figure('units','normalized','position',[.1 .5 .8 .4])

subplot(121),title('SDFs'),hold on
bar(psth(psth(:,2)>0.1,1),psth(psth(:,2)>0,2)/max(psth(:,2)),'FaceColor',ones(1,3)*0.7,'EdgeColor',ones(1,3)*0.7)
% plot(psth(psth(:,2)>0.1,1),psth(psth(:,2)>0,2)/max(psth(:,2)),'k.','MarkerSize',10)
plot([0 0],[0 1],'k:')

[sdf kernel] = msdf(psth,'boxcar',w_boxcar);
[~,peak(1)] = max(sdf(:,2));
subplot(121)
plot(sdf(:,1),sdf(:,2),'m')
subplot(122),title('normalized kernels'),hold on
plot(kernel(:,1),kernel(:,2)/max(kernel(:,2)),'m')

% msdf(x,'boxcar',101) and smooth(x,101) yield the same curves - with the exception of the first w/2 ticks!!!
sdf(:,2) = smooth(psth(:,2),w_boxcar);
[~,peak(2)] = max(sdf(:,2));
subplot(121)
plot(sdf(:,1),sdf(:,2),'r')

[sdf kernel] = msdf(psth,'Gauss',w_Gauss);
[~,peak(3)] = max(sdf(:,2));
subplot(121)
plot(sdf(:,1),sdf(:,2),'b')
subplot(122)
plot(kernel(:,1),kernel(:,2)/max(kernel(:,2)),'b')

[sdf kernel] = msdf(psth,'exp',w_exp);
[~,peak(4)] = max(sdf(:,2));
subplot(121)
plot(sdf(:,1),sdf(:,2),'g')
subplot(122)
plot(kernel(:,1),kernel(:,2)/max(kernel(:,2)),'g')

[sdf kernel] = msdf(psth,'exGauss',w_exGauss);
[~,peak(5)] = max(sdf(:,2));
subplot(121)
plot(sdf(:,1),sdf(:,2),'k')
subplot(122)
plot(kernel(:,1),kernel(:,2)/max(kernel(:,2)),'k')

subplot(121),legend('psth','','boxcar','smooth','Gauss','exp','exGauss')
subplot(122),legend('boxcar','Gauss','exp','exGauss')

disp(' ')
disp('peak indices')
disp(peak)