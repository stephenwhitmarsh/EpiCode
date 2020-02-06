function z = display(z)
%returns String representation of z
if z.unlocked
    unlockstr = '';
else
    unlockstr = ' (LOCKED)';
end
    
disp(['Coherence5LE' unlockstr ' ' num2str(z.Version.major) '.' num2str(z.Version.minor)]);
if (isempty(z.FileInfo))
    disp('No File loaded');
else
    %disp([z.Filename ' recorded on ' z.FileInfo.date]);
    fprintf(['\tDuration  : ' num2str(z.FileInfo.duration) ' seconds\n']);
    fprintf(['\tFrequency : ' num2str(z.FileInfo.frequency) ' Hz\n']);
    fprintf(['\tDuration  : ' num2str(z.FileInfo.electrodes) ' electrodes\n']);
end