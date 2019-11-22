#ifndef __MAKETOWER_H__
#define __MAKETOWER_H__

#include "algo_top_parameters.h"
//#include <algorithm>
//#include <utility>

using namespace std;

void getCrystalInfo(ap_uint<64> words[N_WORDS_PER_FRAME], ap_uint<14> crystals[5][5]){
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME
#pragma HLS ARRAY_PARTITION variable=words complete dim=0
#pragma HLS ARRAY_PARTITION variable=crystals complete dim=0
#pragma HLS latency min=1 max=6

   crystals[0][0] = words[0].range(13,  0);
   crystals[0][1] = words[0].range(27, 14);
   crystals[0][2] = words[0].range(41, 28);
   crystals[0][3] = words[0].range(55, 42);
   ap_uint<8> leftover0 = words[0].range(63, 56);
   
   crystals[0][4] = words[1].range(13,  0);
   crystals[1][0] = words[1].range(27, 14);
   crystals[1][1] = words[1].range(41, 28);
   crystals[1][2] = words[1].range(55, 42);
   ap_uint<8> leftover1 = words[1].range(63, 56);
   
   crystals[1][3] = words[2].range(13,  0);
   crystals[1][4] = words[2].range(27, 14);
   crystals[2][0] = words[2].range(41, 28);
   crystals[2][1] = words[2].range(55, 42);
   ap_uint<8> leftover2 = words[2].range(63, 56);
   
   crystals[2][2] = words[3].range(13,  0);
   crystals[2][3] = words[3].range(27, 14);
   crystals[2][4] = words[3].range(41, 28);
   crystals[3][0] = words[3].range(55, 42);
   ap_uint<8> leftover3 = words[3].range(63, 56);
   
   crystals[3][1] = words[4].range(13,  0);
   crystals[3][2] = words[4].range(27, 14);
   crystals[3][3] = words[4].range(41, 28);
   crystals[3][4] = words[4].range(55, 42);
   ap_uint<8> leftover4 = words[4].range(63, 56);
   
   crystals[4][0] = words[5].range(13,  0);
   crystals[4][1] = words[5].range(27, 14);
   crystals[4][2] = words[5].range(41, 28);
   crystals[4][3] = words[5].range(56, 42);
   ap_uint<8> leftover5 = words[5].range(63, 56);
   
   ap_uint<14> crystalData44 = (((ap_uint<14>) (leftover1.range(5, 0) & 0x3F)) << 8)  | ((ap_uint<14>) leftover0);
   crystals[4][4] = crystalData44;
   
   ap_uint<32> crc = (((ap_uint<32>) leftover5) << 24) | (((ap_uint<32>) leftover4) << 16) | (((ap_uint<32>) leftover3) << 8) | ((ap_uint<32>) leftover2);

}

ap_uint<3> getPeakBinOf5(const ap_uint<12> et[5], const ap_uint<16> etSum) {
#pragma HLS ARRAY_PARTITION variable=et complete dim=0
#pragma HLS INLINE
   uint16_t iEtSum =
      (et[0] >> 1)                +  // 0.5xet[0]
      (et[1] >> 1) + et[1]        +  // 1.5xet[1]
      (et[2] >> 1) + (et[2] << 1) +  // 2.5xet[2]
      (et[3] << 2) - (et[3] >> 1) +  // 3.5xet[3]
      (et[4] << 2) + (et[4] >> 1) ;  // 4.5xet[4]
   uint16_t iAve = 0xBEEF;
   if(     iEtSum <= etSum) iAve = 0;
   else if(iEtSum <= (etSum << 1)) iAve = 1;
   else if(iEtSum <= (etSum + (etSum << 1))) iAve = 2;
   else if(iEtSum <= (etSum << 2)) iAve = 3;
   else iAve = 4;
   return iAve;
}


Tower makeTower(const ap_uint<14> crystals[5][5]) {
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME
#pragma HLS ARRAY_PARTITION variable=crystals complete dim=0

   ap_uint<12> phi_strip[5], eta_strip[5];
#pragma HLS ARRAY_PARTITION variable=phi_strip complete dim=0
#pragma HLS ARRAY_PARTITION variable=eta_strip complete dim=0

   // Compute strips
   for (size_t eta = 0; eta < 5; eta++) {
#pragma LOOP UNROLL
      eta_strip[eta] = 0;
      for (size_t phi = 0; phi < 5; phi++) {
#pragma LOOP UNROLL
	 eta_strip[eta] += crystals[eta][phi].range(9, 0);
      }
   }

   for (size_t phi = 0; phi < 5; phi++) {
#pragma LOOP UNROLL
      phi_strip[phi] = 0;
      for (size_t eta = 0; eta < 5; eta++) {
#pragma LOOP UNROLL
	 phi_strip[phi] += crystals[eta][phi].range(9, 0);
      }
   }

   // Compute tower energy
   ap_uint<16> tet = 0;
   for (size_t eta = 0; eta < 5; eta++) {
#pragma LOOP UNROLL
      tet += eta_strip[eta];
   }
   ap_uint<10> towerEt;
   if(tet > 0x3FF){
      towerEt = 0x3FF;
   }else{
      towerEt = tet;
   }

   // Compute peak locations
   ap_uint<3> peakEta = getPeakBinOf5(eta_strip, tet);
   ap_uint<3> peakPhi = getPeakBinOf5(phi_strip, tet);

   // Small cluster ET is just the 3x5 around the peak
   uint16_t clusterEt = 0;
   for (int dEta = -1; dEta <= 1; dEta++) {
#pragma LOOP UNROLL
      int eta = peakEta + dEta;
      if(eta >= 0 && eta < 5){
	 clusterEt =  clusterEt + eta_strip[eta];
      } else{
	 clusterEt = clusterEt;
      }
   }


   //for now setting peakTime and HoE == 0
   Tower tower(clusterEt, towerEt, peakPhi, peakEta, 0 ,0);

   return tower;
}

void stitchNeigbours( Tower Ai, Tower Bi, Tower &Ao, Tower &Bo) {
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME

   // Check that the clusters are neigbhors in eta or phi
   ap_uint<12> clustered_et = Ai.cluster_et() + Bi.cluster_et();

   ap_uint<10> cluster_et;
   if(clustered_et > 0x3FF){
      cluster_et = 0x3FF;
   }else {
      cluster_et = clustered_et;
   }

   if(Ai.cluster_et() > Bi.cluster_et()){
      // Merge 2 in to 1, and set 2 to remnant energy centered in tower
      Ao = Tower(cluster_et, Ai.tower_et(), Ai.peak_phi(), Ai.peak_eta(), Ai.peak_time(), Ai.hOe()) ;
      Bo = Tower(0, Bi.tower_et(), 2, 2, Bi.peak_time(), Bi.hOe()) ;
   }
   else{
      // Merge 1 in to 2, and set 1 to remnant energy centered in tower
      Ao = Tower(0, Ai.tower_et(), 2, 2, Ai.peak_time(), Ai.hOe()) ;
      Bo = Tower(cluster_et, Bi.tower_et(), Bi.peak_phi(), Bi.peak_eta(), Bi.peak_time(), Bi.hOe()) ;
   }

#ifndef __SYNTHESIS__
   cout << "Stitching cluster cluster [(" + to_string(Ai.peak_phi()) + ", " + to_string(Ai.peak_eta()) + ")]:"<<to_string(Ai.cluster_et())<<" with [(" <<
      to_string(Bi.peak_phi()) + ", " + to_string(Bi.peak_eta()) +")]:"<<to_string(Bi.cluster_et()) << endl;
#endif

}


bool stitchInEta(Tower in[TOWERS_IN_ETA], Tower out[TOWERS_IN_ETA]){
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_PARTITION variable=in complete dim=0
#pragma HLS ARRAY_PARTITION variable=out complete dim=0

   Tower stitch_eta_step1[TOWERS_IN_ETA];
#pragma HLS ARRAY_PARTITION variable=stitch_eta_step1 complete dim=0

   // Stitch in eta: first set of pairs (0,1), (2,3),...(14,15) 
   stitch_eta_step1[TOWERS_IN_ETA-1] = in[TOWERS_IN_ETA-1]; //needed for ODD eta geometry (UNcomment it)
   for (size_t teta = 0; teta < TOWERS_IN_ETA-1; teta += 2) {
#pragma LOOP UNROLL
      if(in[teta].peak_eta() == 4 && in[teta+1].peak_eta() == 0 && in[teta].peak_phi() == in[teta+1].peak_phi()){
	 stitchNeigbours(in[teta], in[teta+1], stitch_eta_step1[teta], stitch_eta_step1[teta+1]);
      } 
      else{
	 stitch_eta_step1[teta]   = in[teta];
	 stitch_eta_step1[teta+1] = in[teta+1];
      }
      //cout<<teta<<":"<<stitch_eta_step1[teta].toString()<<"   "<<teta+1<<":"<<stitch_eta_step1[teta+1].toString()<<std::endl;
   }

   // Stitch in eta: second set of pairs (1,2), (3,4),...(15,16) 
   out[0]   = stitch_eta_step1[0];
   for (size_t teta = 1; teta < TOWERS_IN_ETA-1; teta += 2) {
#pragma LOOP UNROLL
      if(in[teta].peak_eta() == 4 && in[teta+1].peak_eta() == 0 && in[teta].peak_phi() == in[teta+1].peak_phi()){
	 stitchNeigbours(stitch_eta_step1[teta], stitch_eta_step1[teta+1], out[teta], out[teta+1]);
      } 
      else{
	 out[teta]   = stitch_eta_step1[teta];
	 out[teta+1] = stitch_eta_step1[teta+1];
      }
   }

   return true;
}

bool stitchInPhi(Tower inPhi1[TOWERS_IN_ETA], Tower inPhi2[TOWERS_IN_ETA], Tower out[TOWERS_IN_PHI*TOWERS_IN_ETA]){
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_PARTITION variable=inPhi1 complete dim=0
#pragma HLS ARRAY_PARTITION variable=inPhi2 complete dim=0
#pragma HLS ARRAY_PARTITION variable=out complete dim=0

   // Stitch in phi
   for (size_t teta = 0; teta < TOWERS_IN_ETA; teta++) {
#pragma LOOP UNROLL
      if(inPhi1[teta].peak_phi() == 4 && inPhi2[teta].peak_phi() == 0 && inPhi1[teta].peak_eta() == inPhi2[teta].peak_eta()){
	 stitchNeigbours(inPhi1[teta], inPhi2[teta], out[teta], out[teta+TOWERS_IN_ETA]);
      }
      else{
	 out[teta] = inPhi1[teta];
	 out[teta+TOWERS_IN_ETA] = inPhi2[teta];
      }
   }

   return true;
}


#endif /*!__MAKETOWER_H__*/
