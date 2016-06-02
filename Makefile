CC=g++
CFLAGS = -std=c++0x -O3 -pthread -g 

CSI-MAC: cMACprf.cpp cPRFCluster.cpp base.cpp
	$(CC) $(CFLAGS) $^ -o bin/$@


