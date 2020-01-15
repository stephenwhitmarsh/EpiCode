%% axes only
figure('name','Normal');
subplot(3,2,1);plot(randn(9,3));set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
subplot(3,2,[2 4]);plot(randn(9,3));set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
subplot(3,2,3);plot(randn(9,3));set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
subplot(3,2,[5 6]);plot(randn(9,3));set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
figure('name','subplot_er','Visible','off');
subplot_er(3,2,1);plot(randn(9,3));set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
subplot_er(3,2,[2 4]);plot(randn(9,3));set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
subplot_er(3,2,3);plot(randn(9,3));set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
subplot_er(3,2,[5 6]);plot(randn(9,3));set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
set(gcf,'Visible','on')

%% subplot_er
% Create data set
t = (-5:0.1:5)';
f = sin(t);
% Normal subplot
figure('name','Normal');
subplot(2,2,1);
plot(t,f);xlabel('Time [s]');title('f(t) = sin(t)');
subplot(2,2,[2 4]);
plot(t,f);xlabel('Time [s]');legend('f(t) = sin(t)');
subplot(2,2,3);
plot(t,f);ylabel('f(t)');
% New subplot
figure('name','subplot_er','Visible','off');
subplot_er(2,2,1);
plot(t,f);xlabel('Time [s]');title('f(t) = sin(t)');
subplot_er(2,2,[2 4]);
plot(t,f);xlabel('Time [s]');legend('f(t) = sin(t)');
subplot_er(2,2,3);
plot(t,f);ylabel('f(t)');
set(gcf,'Visible','on')

%% subplotXmanyY_er
% Create data set
t = (-10:0.1:10)';
f = sin(t) + 0.2*randn(numel(t),1);
dfdt = [0;diff(f)];
% Normal subplot
figure('name','Normal');
subplot(2,1,1);
plot(t,f);ylabel('Dispalcement [m]');title('f(t) = sin(t) + noise');grid on
subplot(2,1,2);
plot(t,dfdt);legend('f''');ylabel('Velocity [m/s]');grid on
xlabel('Time [s]')
% New subplot
figure('name','subplotXmanyY_er','Visible','off');
ax(1) = subplotXmanyY_er(2,1);
plot(t,f);ylabel('Dispalcement [m]');title('f(t) = sin(t) + noise');grid on
set(gca,'XTickLabel',{})
ax(2) = subplotXmanyY_er(2,2);
plot(t,dfdt);legend('f''');ylabel('Velocity [m/s]');grid on
xlabel('Time [s]');grid on
set(gcf,'Visible','on')
linkaxes(ax,'x');   % Link axes for zooming

%% -> Manually resize figures and enjoy auto-update

