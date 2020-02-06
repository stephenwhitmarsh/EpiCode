function [data, errorcode] = GetEeg( z, start, duration)
%[data, errorcode] = GetEeg( z, start, duration)
% returns next 'duration' samples from eeg-data beginning at start

if length(z.Filename)==0
    data=[];
    disp('No File loaded, can''t get data');
    errorcode=-1;
    return;
end

z = ensureload(z);  % ensure Library is loaded



% different Function Names and allocating problems in Linux and Windows
os = getenv('OS');
if (strfind(os,'Windows')>0)
    GetEegFuncName = 'Eeg3_GetEeg2';
    multiplier = 2;
else
    GetEegFuncName = 'Eeg3_GetEeg';
    multiplier = 1;
end
% create buffer
b = libpointer('int16Ptr',1:z.maxchannels);
b.Value = int16(zeros(1,z.maxchannels*duration)); % allocate


errorcode = calllib(z.libname,GetEegFuncName,start,duration,b);
if errorcode>0
    data = reshape(b.Value(1:z.FileInfo.electrodes* errorcode),z.FileInfo.electrodes, errorcode);
    data = data';
else
    data = [];
end

clear b;
