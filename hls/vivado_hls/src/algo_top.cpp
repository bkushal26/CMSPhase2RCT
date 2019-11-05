#include "algo_top_parameters.h"
#include "algo_top.h"
#include <algorithm>
#include <utility>

#include "ClusterFinder.h"

using namespace std;
using namespace algo;



// Each input link carries the information of a 5x5 region
// Last 16-bits are reserved for CRC
inline Tower unpackInputLink(hls::stream<axiword> &link) {
#pragma HLS INTERFACE axis port=link
#pragma HLS INLINE

   Tower tower;

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
		       tower.crystals[0][k] = Crystal(data >> (k * 14));
		    }
		 } break;
	 case 1: {
		    tower.crystals[0][4] = Crystal((data << 8) | carry);
		    for (size_t k = 0; k < 4; k++) {
#pragma LOOP UNROLL
		       tower.crystals[1][k] = Crystal(data >> (k * 14 + 6));
		    }
		 } break;
	 case 2: {
		    tower.crystals[1][4] = Crystal(data);
		    for (size_t k = 0; k < 3; k++) {
#pragma LOOP UNROLL
		       tower.crystals[2][k] = Crystal(data >> (k * 14 + 14));
		    }
		 } break;
	 case 3: {
		    tower.crystals[2][3] = Crystal((data << 8) | carry);
		    tower.crystals[2][4] = Crystal(data >> 6);
		    for (size_t k = 0; k < 3; k++) {
#pragma LOOP UNROLL
		       tower.crystals[3][k] = Crystal(data >> (k * 14 + 20));
		    }
		 } break;
	 case 4: {
		    for (size_t k = 0; k < 2; k++) {
#pragma LOOP UNROLL
		       tower.crystals[3][k+3] = Crystal(data >> (k * 14));
		    }
		    for (size_t k = 0; k < 2; k++) {
#pragma LOOP UNROLL
		       tower.crystals[4][k] = Crystal(data >> (k * 14 + 28));
		    }
		 } break;
	 case 5: {
		    tower.crystals[4][2] = Crystal((data << 8) | carry);
		    for (size_t k = 0; k < 2; k++) {
#pragma LOOP UNROLL
		       tower.crystals[4][k+3] = Crystal(data >> (k * 14 + 6));
		    }
		 } break;
	 default: break;
      }

      // Remaining data to be used on the next cycle
      carry = data >> 56;
   }

   return tower;
}

void algo_top(hls::stream<axiword> link_in[N_INPUT_LINKS], hls::stream<axiword> link_out[N_OUTPUT_LINKS]) {
#pragma HLS INTERFACE axis port=link_in
#pragma HLS INTERFACE axis port=link_out
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0

#ifndef ALGO_PASSTHROUGH

   // Step 1: Unpack links
   Tower ecalTowers[TOWERS_IN_ETA][TOWERS_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ecalTowers complete dim=0

   for (size_t i = 0; i < TOWERS_IN_ETA * TOWERS_IN_PHI; i++) {
#pragma LOOP UNROLL

      size_t  ieta = i / TOWERS_IN_PHI;  
      size_t  iphi = i % TOWERS_IN_PHI;
      ecalTowers[ieta][iphi] = unpackInputLink(link_in[i]);

#ifndef __SYNTHESIS__
      cout<<"ecalTowers["<<ieta<<"]["<<iphi<<"]:"<<endl<<ecalTowers[ieta][iphi].toString()<<std::endl;
#endif

   }

   // Step 2: Compute clusters
   Cluster ecalClusters[TOWERS_IN_ETA][TOWERS_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ecalClusters complete dim=0

   for (size_t i = 0; i < TOWERS_IN_ETA * TOWERS_IN_PHI; i++) {
#pragma LOOP UNROLL
      size_t  ieta = i / TOWERS_IN_PHI;  
      size_t  iphi = i % TOWERS_IN_PHI;
      uint16_t towerEt = 0;
      ecalClusters[ieta][iphi] = ecalTowers[ieta][iphi].computeCluster(ieta, iphi, towerEt);

#ifndef __SYNTHESIS__
      cout << "Clustering["<<ieta<<"]["<<iphi<<"]:"<< ecalClusters[ieta][iphi].toString() << endl;
#endif
   }

   // Step 3: Merge neighbor clusters
   Cluster ecalClustersStitched[TOWERS_IN_ETA][TOWERS_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ecalClustersStitched complete dim=0

   Cluster::stitchRegion<TOWERS_IN_ETA,TOWERS_IN_PHI>(ecalClusters, ecalClustersStitched);

#ifndef __SYNTHESIS__
   for (size_t i = 0; i < TOWERS_IN_ETA; i++) {
      for (size_t j = 0; j < TOWERS_IN_PHI; j++) {
	 cout << "Stitched["<<i<<"]["<<j<<"]:"<< ecalClustersStitched[i][j].toString() << endl;
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
	    r.data = ((uint64_t)ecalClustersStitched[eta+1][phi] << 32) |
	       ((uint64_t)ecalClustersStitched[eta][phi]);
	 }
	 else if(TOWERS_IN_ETA%2 == 0 && eta >= TOWERS_IN_ETA)
	    r.data = 0;

	 if(TOWERS_IN_ETA%2 && eta < TOWERS_IN_ETA-1){ // it eta has odd geometry
	    r.data = ((uint64_t)ecalClustersStitched[eta+1][phi] << 32) |
	       ((uint64_t)ecalClustersStitched[eta][phi]);
	 }
	 else if(TOWERS_IN_ETA%2 && eta == TOWERS_IN_ETA-1){
	    r.data = ((uint64_t)ecalClustersStitched[eta][phi]);
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