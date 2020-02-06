function z = unlock( z, unlockcode)
% unlock(z,unlockcode);
%  if called unlock(z) hardcoded unlock code is used
%  if called unlock(z,code) with code=[1;2;3;4] code is used

%% UNLOCK CODE
% insert your code here and you don't need to care about unlock anymore
% set values to [] and unlock every time manually
clef.int1=[654125];
clef.int2=[894265];
clef.int3=[512743];
clef.int4=[554190];
%%

global Coherence_lib_isunlocked;

z = ensureload(z);

if (exist('unlockcode','var') && length(unlockcode)==4);
    clef.int1 = unlockcode(1);
    clef.int2 = unlockcode(2);
    clef.int3 = unlockcode(3);
    clef.int4 = unlockcode(4);
elseif ~isempty(z.unlockstruct);
    clef = z.unlockstruct;
end


res = calllib(z.libname,'Eeg3_Unlock',clef);
if res==0
    z.unlocked = 1;
    Coherence_lib_isunlocked = 1;
else
    error(['Unlock failed with code: ' num2str(res)]);
end

disp Coherence_lib_isunlocked