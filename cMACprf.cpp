

#include "cPRFCluster.h"

int main(int argc, const char* argv[]) {
	
	cPRFCluster kk;

	try {		
		kk.Run(argc, argv);
		//kk.Run('Attacin-C_DmDs_pol.fas','Attacin-C_DmDs_div.fas');
	}
	catch (...) {
	}

	return 1;
}


 
