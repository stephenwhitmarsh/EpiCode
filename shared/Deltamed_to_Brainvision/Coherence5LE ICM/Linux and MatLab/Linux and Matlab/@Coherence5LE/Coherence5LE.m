function z = Coherence5LE( unlockstruct )
% base class for Coherence5LE 
% Coherence5LE(unlockcode)
% If no unlock code is supplied the one hardcoded in @Coherence5LE/unlock.m is used
%   INPUT: unlockcode(optional) i.e. [ 12; 23 ;34 ;45]
%   OUTPUT: Coherence5LE object

global Coherence_lib_isunlocked; % global bool to save if lib is allready unlocked



public.unlockstruct=[];
public.unlocked = Coherence_lib_isunlocked;
public.libheader = [];
public.libfullname = '';
public.libname = '';
public.Filename='';
public.FileInfo = [];
public.Version = '';
public.maxchannels = 1024;
public.Impedances = [];

% create class
z = class(public, mfilename);

% load library
os = getenv('OS');
if (strfind(os,'Windows')>0)
    z.libfullname = 'Coherence5LE.dll';
    z.libheader = 'MatlabCoherence5LEdll.h';
    z.libname = 'Coherence5LE';
else
    z.libfullname = 'libCoherence5LEWrapper.so';
    z.libheader = 'MatlabCoherence5LEWrapper.h';
    z.libname = 'libCoherence5LEWrapper';
end

if (exist('unlockstruct','var'))
    z.unlockstruct = unlockstruct;
end



z = ensureload(z);

% Version

v = libpointer('TVersion');
v.Value.major=0; % allocate space for v
C=calllib(z.libname,'Eeg3_Version',v);
z.Version = v.Value;
clear v;



