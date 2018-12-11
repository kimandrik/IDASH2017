# IDASH2017

A winning solution of Track 3 at iDASH privacy and security competition 2017 (http://www.humangenomeprivacy.org/2017/).

## License
Copyright (c) by CryptoLab inc.
This program is licensed under a
Creative Commons Attribution-NonCommercial 3.0 Unported License.
You should have received a copy of the license along with this
work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.

## Run

run: ./IDASH2017 trainfile isYfirst numIter k gammaUp gammaDown isInitZero fold isEncrypted testfile

./IDASH2017 string bool long long double double bool long bool string

example: ./IDASH2017 "../data/data103x1579.txt" 1 7 5 1 -1 1 5 1

example: ./IDASH2017 "../data/1_training_data_csv" 1 7 5 1 -1 1 0 1 "../data/1_testing_data_csv"

parameters:

trainfile - path to train file

isYfirst - {0,1} y parameter first OR last

numIter - number of iterations

kdeg - degree of sigmoid approximation function k in {3,5,7}

gammaUp - corresponds to learning rate

gammaDown - corresponds to learning rate

isInitZero - is initial weights zero or average

fold - folding method if arguments <= 8 we use folding method

isEncrypted - encrypted or plain

testfile - path to test file (checks if number of arguments > 8 then we use standard method



current files that in data folder (filename isYfirst):

"../data/data5x500.txt" false

"../data/data9x1253.txt" false

"../data/data15x1500.txt" false

"../data/data16x101.txt" false

"../data/data27x148.txt" false

"../data/data51x653.txt" false

"../data/data67x216.txt" false

"../data/data103x1579.txt" true

"../data/1_training_data.csv" true

"../data/1_testing_data.csv" true



FYI: approx 3 suggested iter: 4, 9, 18, 36, ...

FYI: approx 5 suggested iter: 3, 7, 14, 28, ...

FYI: approx 7 suggested iter: 3, 7, 14, 28, ...

