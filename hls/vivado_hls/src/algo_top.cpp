#include "algo_top_parameters.h"
#include "algo_top.h"
#include <algorithm>
#include <utility>

#include "ClusterFinder.h"

using namespace std;
using namespace algo;

// Each input link carries the information of a 5x5 region
// Last 16-bits are reserved for CRC
Cluster unpackInputLink(hls::stream<axiword> &link) {
#pragma HLS INTERFACE axis port=link
#pragma HLS INLINE

   Crystal crystals[5][5];
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


   Cluster tower;
   //compute clusters in each tower using 5x5 crystals
   tower = computeCluster(crystals);

   return tower;
}

void algo_top(hls::stream<axiword> link_in[N_INPUT_LINKS], hls::stream<axiword> link_out[N_OUTPUT_LINKS]) {
#pragma HLS INTERFACE axis port=link_in
#pragma HLS INTERFACE axis port=link_out
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0

#ifndef ALGO_PASSTHROUGH

   // Step 1: Unpack links and compute clusters
   Cluster ecalClusters[TOWERS_IN_PHI][TOWERS_IN_ETA];
#pragma HLS ARRAY_PARTITION variable=ecalClusters complete dim=0

   for (size_t i = 0; i < TOWERS_IN_ETA * TOWERS_IN_PHI; i++) {
#pragma LOOP UNROLL

      size_t  ieta = i / TOWERS_IN_PHI;  
      size_t  iphi = i % TOWERS_IN_PHI;
      ecalClusters[iphi][ieta] = unpackInputLink(link_in[i]);

      ecalClusters[iphi][ieta].tower_eta = ieta;
      ecalClusters[iphi][ieta].tower_phi = iphi;

//- #ifndef __SYNTHESIS__
//-       cout<<"ecalClusters["<<iphi<<"]["<<ieta<<"]: "<<ecalClusters[iphi][ieta].toString()<<std::endl;
//- #endif

   }

#ifndef __SYNTHESIS__
   for (size_t i = 0; i < TOWERS_IN_ETA ; i++) {
      cout << "Clustering[0]["<<i<<"]:"<< ecalClusters[0][i].toString() ;
      cout << "  ||  "; 
      cout << "Clustering[1]["<<i<<"]:"<< ecalClusters[1][i].toString() << endl;
   }
#endif
 
   // Step 2: Merge neighbours in eta
   Cluster etaStitched_phi1[TOWERS_IN_ETA];
   Cluster etaStitched_phi2[TOWERS_IN_ETA];
#pragma HLS ARRAY_PARTITION variable=etaStitched_phi1 complete dim=0 
#pragma HLS ARRAY_PARTITION variable=etaStitched_phi2 complete dim=0 

   bool etaStitch_phi1 = stitchInEta( ecalClusters[0], etaStitched_phi1);
   bool etaStitch_phi2 = stitchInEta( ecalClusters[1], etaStitched_phi2);

#ifndef __SYNTHESIS__
   cout << " *** ETA stitched *** "<<std::endl;
   for (size_t i = 0; i < TOWERS_IN_ETA ; i++) {
      cout << "Clustering[0]["<<i<<"]:"<< etaStitched_phi1[i].toString() ;
      cout << "  ||  "; 
      cout << "Clustering[1]["<<i<<"]:"<< etaStitched_phi2[i].toString() << endl;
   }
#endif

   //Step 3: Merge neighbours in phi
   Cluster ecalClustersStitched[TOWERS_IN_PHI][TOWERS_IN_ETA];
#pragma HLS ARRAY_PARTITION variable=ecalClustersStitched complete dim=0
   
   bool phiStitch = stitchInPhi(etaStitched_phi1, etaStitched_phi2, ecalClustersStitched);

#ifndef __SYNTHESIS__
   for (size_t i = 0; i < TOWERS_IN_ETA; i++) {
      for (size_t j = 0; j < TOWERS_IN_PHI; j++) {
	 cout << "Stitched["<<j<<"]["<<i<<"]:"<< ecalClustersStitched[j][i].toString() << endl;
      }
   }
#endif

   // Step 4: Pack the outputs
   for (size_t i = 0; i < N_OUTPUT_LINKS; i++) {
#pragma LOOP UNROLL
      for (size_t j = 0; j < N_OUTPUT_WORDS_PER_FRAME-1; j++) {
#pragma LOOP UNROLL

	 const size_t phi = i/2; // Two links carry information for each eta
	 const size_t eta =  j*2 + (i%2)*10;

	 axiword r; r.last = 0; r.user = 0;

	 if( TOWERS_IN_ETA%2 == 0 && eta < TOWERS_IN_ETA){ // if eta is an even number
	    r.data = ((uint64_t)ecalClustersStitched[phi][eta+1] << 32) |
	       ((uint64_t)ecalClustersStitched[phi][eta]);
	 }
	 else if(TOWERS_IN_ETA%2 == 0 && eta >= TOWERS_IN_ETA)
	    r.data = 0;

	 if(TOWERS_IN_ETA%2 && eta < TOWERS_IN_ETA-1){ // it eta has odd geometry
	    r.data = ((uint64_t)ecalClustersStitched[phi][eta+1] << 32) |
	       ((uint64_t)ecalClustersStitched[phi][eta]);
	 }
	 else if(TOWERS_IN_ETA%2 && eta == TOWERS_IN_ETA-1){
	    r.data = ((uint64_t)ecalClustersStitched[phi][eta]);
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
