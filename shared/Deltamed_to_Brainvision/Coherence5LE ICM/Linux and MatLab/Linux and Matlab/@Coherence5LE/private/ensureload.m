function z = ensureload(z)
%ENSURELOAD 
% makes sure that library is loaded and the right file is opened
global Coherence_lib_filename;

if ~libisloaded(z.libname) % check if lib is loaded
    warning off; % skip warnings during loadlibrary
    [notfound,loadlibwarnings]=loadlibrary(z.libfullname,z.libheader,'alias',z.libname);
    warning on;
    
    C=calllib(z.libname,'Eeg3_DebugFileSwitch',10);
    
    % init
    res=calllib(z.libname,'Eeg3_Initialisation');
    %unlock
    z = unlock(z);
end

% check for right file

% open file
% if we have filename and it's different to the one currently used by lib
if length(z.Filename)>0 && ~strcmp(z.Filename,Coherence_lib_filename)
    res = calllib(z.libname,'Eeg3_CloseFile');
    Coherence_lib_filename = z.Filename;
    p = libpointer('TCoh3');
    p.Value.duration=0; % allocate space for p
    res = calllib(z.libname,'Eeg3_OpenFile',z.Filename,p);
    if res==-102
        error(['File not found. Make sure filename is correct and you are in the right working directoy']);
    elseif res~=0
        error(['Opening File with code: ' num2str(res)]);
    end
    clear p;
end