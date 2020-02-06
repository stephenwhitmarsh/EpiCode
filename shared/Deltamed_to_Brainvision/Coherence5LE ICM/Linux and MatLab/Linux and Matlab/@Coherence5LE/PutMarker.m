function [z, returncode ] = PutMarker( z,pos, duration, text , eventtype)
% PutMarker(z, pos , duration, text, eventype)
% puts a marker at given position

z = ensureload(z);

tmarker = libpointer('TMarker');
tmarker.Value.pos = pos;
tmarker.Value.duration = duration;
tmarker.Value.text = int8(text);
tmarker.Value.evttype = eventtype;
returncode = calllib(z.libname,'Eeg3_PutMarker', tmarker);

