winegcc Coherence5LE.c -o Coherence5LE
gcc -fPIC -c Coherence5LEWrapper.c
ld -shared -soname libCoherence5LEWrapper.so.1.0 -o libCoherence5LEWrapper.so.1.0 -lc Coherence5LEWrapper.o
ln -sf libCoherence5LEWrapper.so.1.0 libCoherence5LEWrapper.so 
#   Uncomment this for testing
gcc -o testwrapper testwrapper.c -L. -lCoherence5LEWrapper
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH

