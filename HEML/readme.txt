in folder Debug/ run:

make clean
make all

run:
./HEML filename(string) isYfirst(0,1) iter(long) learnPortion(double) is3approx(0,1) isEncrypted(0,1)

example:
./HEML "../data/data103x1579.txt" 1 7 0.9 0 1

filename - path to file
isYfirst - is y first OR last
iter - number of iterations
learnPortion - portion of data for learning
is3approx - used 3 OR 5 degree polynomial aprroximation
isEncrypted - encrypted OR unencrypted