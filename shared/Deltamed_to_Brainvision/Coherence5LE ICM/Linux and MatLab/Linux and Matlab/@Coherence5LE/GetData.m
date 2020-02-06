function Data  = GetData(z)
% Data = GetData(z)
% returns all eeg data from z

z = ensureload(z);  % ensure Library is loaded

Data = [];

res =0;
blocksize=z.FileInfo.frequency*2;
start = 0;
while (res>=0)    
   [d, res] = GetEeg(z,start,blocksize);
   start = start+blocksize;
   Data = [Data ; d];
end







