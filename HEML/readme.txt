Dear IDASH organizers

Sending description of our final version of program on server.

We installed gmp, and ntl libraries in default folders.

In folder idash/HEAAN there is a FHE library that we used in our implementation, more information here:
https://eprint.iacr.org/2016/421
https://github.com/kimandrik/HEAAN

Code for the task 3 is in folder idash/HEML/src. We used a variant of Nesterov Logistic Gradient Descent (NLGD) in implementation.
In folder idash/HEML/data there are some datasets that we used for testing, including testing dataset from IDASH (data103x1579.txt)
makefile is in folder idash/HEML/Debug. To compile the program(it already compiled) in folder idash/HEML/Debug run:
make clean
make all

To run the program in folder idash/HEML/Debug run:
./HEML filename(string) isYfirst(bool) iter(long) learnPortion(double) approx(long) isEncrypted(bool)

example: ./HEML "../data/data103x1579.txt" 1 7 0.9 7 1

parameters:
filename - path to file
isYfirst - {0,1} Y parameter is first OR last
iter - number of iterations
learnPortion - portion of data used for learning (randomly chosen from sample set)
approx - {3,5,7} polynomial approximation degree used (better use degree 5 or 7) 
isEncrypted - {0,1} encrypted OR unencrypted learning

current files that in data folder (filename isYfirst):
"../data/data5x500.txt" false
"../data/data9x1253.txt" false
"../data/data15x1500.txt" false
"../data/data16x101.txt" false
"../data/data27x148.txt" false
"../data/data51x653.txt" false
"../data/data67x216.txt" false
"../data/data103x1579.txt" true

if approx = 3 suggested number of iterations: 4, 9, 18, 36, ...
if approx = 5 suggested number of iterations: 3, 7, 14, 28, ...
if approx = 7 suggested number of iterations: 3, 7, 14, 28, ...

Explanation: increasing number of iterations increases ciphertext modulus space. To obtain 80 bits of security we have to increase dimension N of RLWE scheme. Suggested number of iterations means that adding one iterations increases N, making program run slower.

In file idash/HEML/src/HEML.cpp we also have some other parameters that could be changed for optimizations and better results, but we tried to hardcore them so they work for any dataset (for example we have xyBits = 37, wBits=37, but if the dataset is not so big we can lower these parameters a little).

After each iteration we decrypt ciphers to check the correctness, and then continue with the same ciphers (it is not recryption, just decryption for checking. You can remove this check step in idash/HEML/src/HEML.cpp).

We attached files with some results on data103x1579.txt.

./HEML "../data/data103x1579.txt" 1 3 1 5 1
results_iter3_deg5.txt   total time ~4 mins (~1.3 mins for one iteration)

./HEML "../data/data103x1579.txt" 1 3 1 7 1
results_iter3_deg7.txt   total time ~5 mins (~1.5 mins for one iteration)

./HEML "../data/data103x1579.txt" 1 7 1 5 1
results_iter7_deg5.txt   total time ~25 mins (~3.5 mins for one iteration)

./HEML "../data/data103x1579.txt" 1 7 1 7 1
results_iter7_deg7.txt   total time ~30 mins (~4.2 mins for one iteration)

./HEML "../data/data103x1579.txt" 1 14 1 5 1
results_iter14_deg5.txt  total time ~3 hours (~12 mins for one iteration)

./HEML "../data/data103x1579.txt" 1 14 1 7 1
results_iter14_deg7.txt  total time ~4 hours (~15 mins for one iteration)

each iteration need less time than previous, so we just calculated average

Kind Regards, Andrey Kim