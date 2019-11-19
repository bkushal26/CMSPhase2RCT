#include "algo_top_parameters.h"
#include "algo_top.h"
#include <algorithm>
#include <utility>

#include "objects.h"
#include "makeTower.h"

using namespace std;
using namespace algo;



// Each input link carries the information of a 5x5 region
// Last 16-bits are reserved for CRC
bool unpackInputLink(hls::stream<axiword> &link, Crystal crystals[5][5]) {
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME   
#pragma HLS INTERFACE axis port=link
#pragma HLS INLINE

#pragma HLS ARRAY_PARTITION variable=crystals complete dim=0

   uint8_t carry = 0;

   for (size_t i = 0; i < N_WORDS_PER_FRAME; i++) {
#pragma LOOP UNROLL
#ifndef __SYNTHESIS__
      // Avoid simulation warnings
      if (link.empty()) continue;
#endif

      uint64_t data = link.read().data;

      switch (i) {
	 case 0: {
		    for (size_t k = 0; k < 4; k++) {
#pragma LOOP UNROLL
		       crystals[0][k] = Crystal(data >> (k * 14));
		    }
		 } break;
	 case 1: {
		    crystals[0][4] = Crystal((data << 8) | carry);
		    for (size_t k = 0; k < 4; k++) {
#pragma LOOP UNROLL
		       crystals[1][k] = Crystal(data >> (k * 14 + 6));
		    }
		 } break;
	 case 2: {
		    crystals[1][4] = Crystal(data);
		    for (size_t k = 0; k < 3; k++) {
#pragma LOOP UNROLL
		       crystals[2][k] = Crystal(data >> (k * 14 + 14));
		    }
		 } break;
	 case 3: {
		    crystals[2][3] = Crystal((data << 8) | carry);
		    crystals[2][4] = Crystal(data >> 6);
		    for (size_t k = 0; k < 3; k++) {
#pragma LOOP UNROLL
		       crystals[3][k] = Crystal(data >> (k * 14 + 20));
		    }
		 } break;
	 case 4: {
		    for (size_t k = 0; k < 2; k++) {
#pragma LOOP UNROLL
		       crystals[3][k+3] = Crystal(data >> (k * 14));
		    }
		    for (size_t k = 0; k < 2; k++) {
#pragma LOOP UNROLL
		       crystals[4][k] = Crystal(data >> (k * 14 + 28));
		    }
		 } break;
	 case 5: {
		    crystals[4][2] = Crystal((data << 8) | carry);
		    for (size_t k = 0; k < 2; k++) {
#pragma LOOP UNROLL
		       crystals[4][k+3] = Crystal(data >> (k * 14 + 6));
		    }
		 } break;
	 default: break;
      }

      // Remaining data to be used on the next cycle
      carry = data >> 56;
   }

   return true;
}

void algo_top(hls::stream<axiword> link_in[N_INPUT_LINKS], hls::stream<axiword> link_out[N_OUTPUT_LINKS]) {
#pragma HLS INTERFACE axis port=link_in
#pragma HLS INTERFACE axis port=link_out
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0

#ifndef ALGO_PASSTHROUGH

   // Step 1: Unpack links
   Tower ecalTowers[TOWERS_IN_PHI][TOWERS_IN_ETA];
#pragma HLS ARRAY_PARTITION variable=ecalTowers complete dim=0

   for (size_t i = 0; i < TOWERS_IN_ETA * TOWERS_IN_PHI; i++) {
#pragma LOOP UNROLL

      Crystal crystal[5][5];
      bool unpackLinks = unpackInputLink(link_in[i], crystal);

      size_t  ieta = i / TOWERS_IN_PHI;  
      size_t  iphi = i % TOWERS_IN_PHI;

      ecalTowers[iphi][ieta] = makeTower(crystal);

//- #ifndef __SYNTHESIS__
//-       cout<<"ecalTowers["<<iphi<<"]["<<ieta<<"]: "<<ecalTowers[iphi][ieta].toString()<<std::endl;
//- #endif
   }

#ifndef __SYNTHESIS__
   cout <<endl<<"-->> ECAL Towers: *** pre-stitching *** "<<std::endl;
   for (size_t i = 0; i < TOWERS_IN_ETA ; i++) {
      cout << "[0]["<<i<<"]:"<< ecalTowers[0][i].toString() ;
      cout << "  ||  "; 
      cout << "[1]["<<i<<"]:"<< ecalTowers[1][i].toString() ;
      cout<< endl;
   }
#endif

   // Step 2: stitching towers in eta
   Tower etaStitchedTower_phi1[TOWERS_IN_ETA];
   Tower etaStitchedTower_phi2[TOWERS_IN_ETA];
#pragma HLS ARRAY_PARTITION variable=etaStitchedTower_phi1 complete dim=0
#pragma HLS ARRAY_PARTITION variable=etaStitchedTower_phi2 complete dim=0

   bool etaStitch_phi1 = stitchInEta( ecalTowers[0], etaStitchedTower_phi1);
   bool etaStitch_phi2 = stitchInEta( ecalTowers[1], etaStitchedTower_phi2);

#ifndef __SYNTHESIS__
   cout <<endl<<"-->> ECAL Towers: *** ETA-stitched *** "<<std::endl;
   for (size_t i = 0; i < TOWERS_IN_ETA ; i++) {
      cout << "[0]["<<i<<"]:"<< etaStitchedTower_phi1[i].toString() ;
      cout << "  ||  "; 
      cout << "[1]["<<i<<"]:"<< etaStitchedTower_phi2[i].toString() ;
      cout<< endl;
   }
#endif

   // Step 3: stitching towers in phi
   Tower stitchedTower[TOWERS_IN_PHI * TOWERS_IN_ETA];
#pragma HLS ARRAY_PARTITION variable=stitchedTower complete dim=0 

 bool phiStitch = stitchInPhi(etaStitchedTower_phi1, etaStitchedTower_phi2, stitchedTower); 

#ifndef __SYNTHESIS__
   cout <<endl<<"-->> ECAL Towers: *** PHI-stitched *** "<<std::endl;
   for (size_t i = 0; i < TOWERS_IN_ETA ; i++) {
      cout << "[0]["<<i<<"]:"<< stitchedTower[i].toString() ;
      cout << "  ||  "; 
      cout << "[1]["<<i<<"]:"<< stitchedTower[i+TOWERS_IN_ETA].toString() ;
      cout<< endl;
   }
#endif

   // Step 4: Pack the outputs
   for (size_t i = 0; i < N_OUTPUT_LINKS; i++) {
#pragma LOOP UNROLL
      for (size_t j = 0; j < N_OUTPUT_WORDS_PER_FRAME-1; j++) {
#pragma LOOP UNROLL

	 const size_t phi = i/2; // Two links carry information for each eta
	 const size_t eta =  j*2 + (i%2)*10;
	 const size_t ilink= eta+(phi*TOWERS_IN_ETA) ;

	 axiword r; r.last = 0; r.user = 0;


	 if(TOWERS_IN_ETA%2 && eta < TOWERS_IN_ETA-1){ // it eta has odd geometry
	    r.data = ((uint64_t)stitchedTower[ilink+1] << 32) |
	       ((uint64_t)stitchedTower[ilink]);
	 }
	 else if(TOWERS_IN_ETA%2 && eta == TOWERS_IN_ETA-1){
	    r.data = ((uint64_t)stitchedTower[ilink]);
	 }
	 else if(TOWERS_IN_ETA%2 && eta > TOWERS_IN_ETA-1)
	    r.data = 0;

	 link_out[i].write(r);

      }

      // Last word is CRC
      axiword r; r.last = 0; r.user = 0; r.data = 0;
      link_out[i].write(r);
   }

#else
   // Algo passthrough (for testing)
   for (size_t i = 0; i < N_WORDS_PER_FRAME; i++) {
      axiword r[N_INPUT_LINKS];

      // Read all inputs
      for (size_t l = 0; l < N_INPUT_LINKS; l++)
	 r[l] = link_in[l].read();

      // Write inputs to outputs
      for (size_t l = 0; l < N_OUTPUT_LINKS; l++) {
	 if (l >= N_INPUT_LINKS) {
	    link_out[l].write(r[N_INPUT_LINKS-1]);
	 } else {
	    link_out[l].write(r[l]);
	 }
      }
   }


#endif /* !ALGO_PASSTHROUGH */

}
