function z = OpenFile( z, filename )
%z = OpenFile( z, filename )
% OpenFile
global Coherence_lib_filename;


%filename = which(filename); 

if ~exist(filename,'file')
    error(['File ''' filename ''' does not exist' ]);
end


z = ensureload(z); % this functions handles to OpenFile()

%st=libstruct('TCoh3Le')
%tmp=get(st)   
%p = libpointer('TCoh3',st);
p = libpointer('TCoh3');

p.Value.duration=0; % allocate space for p
res=calllib(z.libname,'Eeg3_OpenFile', filename,p);

if res==-102
    error(['File not found. Make sure filename is correct and you are in the right working directoy']);
elseif res~=0
    error(['Opening File with code: ' num2str(res)]);
end

Coherence_lib_filename = filename;
z.Filename = filename;
z.FileInfo = p.Value;

%z.FileInfo.name= char(reshape(z.FileInfo.name,8,1024)');
%z.FileInfo.date = char(z.FileInfo.date);
%z.FileInfo.unit= char(reshape(z.FileInfo.unit,4,1024)');

reshape(z.FileInfo.date,24,256); 
z.FileInfo.date = char(z.FileInfo.date.value)';

reshape(z.FileInfo.name,8,1024);
z.FileInfo.name= char(z.FileInfo.name.value)';

reshape(z.FileInfo.unit,4,1024)
z.FileInfo.unit= char(z.FileInfo.unit.value');

reshape(z.FileInfo.type,4,1024)
z.FileInfo.type= char(z.FileInfo.type.value');

clear p;

