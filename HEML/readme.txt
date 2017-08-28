in folder Debug/ run:

make clean
make all
./HEML

in file src/HEML.cpp in main method you can change parameters:


filename - path to file (currently on train file)
isYfirst - is y first or last (in train file y is first)
iter - number of iterations (currently iter = 7)
learnPortion - portion of data for learning (currently learnPortion = 0.9)
is3approx - used 3 degree or 5 degree polynomial aprroximation (currently is3approx = false means 5 degree polynomial approximation is used)