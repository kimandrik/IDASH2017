new:
	make clean
	make object
	make exe

exe:
	g++ src/HEML.cpp -o HEML ../../HEAANBOOT/HEAANBOOT/lib/libHEAANBOOT.a -lntl -lgmp -pthread -lm -O2 -std=c++11 -I ../../HEAANBOOT/HEAANBOOT/src src/CipherGD.o src/GD.o src/TestGD.o

object:
	g++ -g -c src/CipherGD.cpp -o src/CipherGD.o -lntl -lgmp -lm -std=c++11 -I ../../HEAANBOOT/HEAANBOOT/src
	g++ -g -c src/GD.cpp -o src/GD.o -lntl -lgmp -lm -std=c++11 -I ../../HEAANBOOT/HEAANBOOT/src
	g++ -g -c src/TestGD.cpp -o src/TestGD.o -lntl -lgmp -lm -std=c++11 -I ../../HEAANBOOT/HEAANBOOT/src

clean:
	rm -f src/*.o
	rm -f HEML
