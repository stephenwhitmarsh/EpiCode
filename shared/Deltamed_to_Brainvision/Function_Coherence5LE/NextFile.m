function [z,res] = NextFile( z, direction )
% [z,errorcode] = NextFile(z1, direction)

global Coherence_lib_filename;


z = ensureload(z); % this functions handles to OpenFile()

p = libpointer('TCoh3');
p.Value.duration=0; % allocate space for p

filename = libpointer('string');
filename.Value(1:261) = ' ';

[res, filename]=calllib(z.libname,'Eeg3_NextFile',direction,  filename,p);
if res==-102
    error(['File not found. Make sure filename is correct and you are in the right working directoy']);
elseif res~=0
    error(['Opening File with code: ' num2str(res)]);
end

Coherence_lib_filename = filename;
z.Filename = filename;
z.FileInfo = p.Value;
z.FileInfo.name= char(reshape(z.FileInfo.name,8,128)');
z.FileInfo.data = char(z.FileInfo.date);
z.FileInfo.unit= char(reshape(z.FileInfo.unit,4,128)');



clear p;

