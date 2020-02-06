[notfound,warning]=loadlibrary('libCoherence5LEWrapper.so','Coherence5LEWrapper.h');

libisloaded('libCoherence5LEWrapper')
m = libfunctions('libCoherence5LEWrapper')

res=calllib('libCoherence5LEWrapper','Eeg3_Initialisation')

clef.int1=1;
clef.int2=2;
clef.int3=3;
clef.int4=4;
B=calllib('libCoherence5LEWrapper','Eeg3_Unlock',clef)


v = libpointer('TVersion');
v.Value.major=0; % allocate space for v
C=calllib('libCoherence5LEWrapper','Eeg3_Version',v)
get(v,'Value')
clear v;

p = libpointer('TCoh3');
p.Value.duration=0; % allocate space for p

C=calllib('libCoherence5LEWrapper','Eeg3_OpenFile','test_0001.Eeg',p)


b = libpointer('int16Ptr');
b.Value = int16(zeros(p.Value.electrodes,200));

clear p;
%C=calllib('libCoherence5LEWrapper','Eeg3_GetEeg',1,200,b);
clear b;

tmarker = libpointer('TMarker');
tmarker.Value.pos = 1;
tmarker.Value.pos = 1;
tmarker.Value.text = int8('test');
tmarker.Value.evttype = 1;
returncode = calllib('libCoherence5LEWrapper','Eeg3_PutMarker', tmarker);
clear tmarker;


res =calllib('libCoherence5LEWrapper','Eeg3_GetMarkersNumber', 0,100)

clear p;
unloadlibrary('libCoherence5LEWrapper');