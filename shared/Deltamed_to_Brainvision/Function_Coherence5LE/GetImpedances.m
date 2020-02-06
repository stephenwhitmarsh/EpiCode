function impedances = GetImpedances(z, startpos )
%impedances = GetImpedances(z, startpos )
% returns impedances
z = ensureload(z);

if ~exist('startpos','var')
    startpos = z.FileInfo.duration*z.FileInfo.frequency;
end

% Get Impedances
imps = libpointer('TImpedances');
imps.Value.pos=0; % allocate space for imps
[res, i] =calllib(z.libname,'Eeg3_GetImpedances',startpos,imps);
impedances = i;
impedances.text = char(impedances.text);
if (res~=0)
    impedances = []; % on error return empty
end
clear imps;
