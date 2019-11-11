#ifndef __CLUSTERFINDER_H__
#define __CLUSTERFINDER_H__

#include "algo_top_parameters.h"

using namespace std;


/* ECAL crystal object definition */
struct Crystal {
   Crystal() : energy(0), timing(0), spike(false) {};

   Crystal(ap_uint<14> i) {
      this->energy = i;
      this->timing = i >> 10;
      this->spike = i >> 13;
   }

   Crystal& operator=(const Crystal& rhs) {
      energy = rhs.energy;
      timing = rhs.timing;
      spike  = rhs.spike;
      return *this;
   }

   inline operator uint16_t() {
      return  ((uint16_t)(this->spike) << 13) |
	 ((uint16_t)(this->timing) << 10) |
	 this->energy;
   }

#ifndef __SYNTHESIS__
   string toString() {
      return "energy = " + to_string(this->energy) + ", timing = " + to_string(this->timing) + ", spike = " + to_string(this->spike);
   }
#endif

   ap_uint<10> energy;
   ap_uint<3>  timing;
   ap_uint<1>  spike;
};

/* ECAL cluster object definition */
struct Cluster {
   Cluster() : et(0),  peak_phi(0), peak_eta(0), tower_et(0) {};

   Cluster(uint32_t i) {
      this->et = i;
      this->peak_phi = i >> 16;
      this->peak_eta = i >> 19;
      this->tower_et = i >> 22; 
   }

   Cluster& operator=(const Cluster& rhs) {
      et        = rhs.et;
      peak_phi  = rhs.peak_phi;
      peak_eta  = rhs.peak_eta;
      tower_et  = rhs.tower_et;

      return *this;
   }

   inline friend bool operator >(const Cluster& c1, const Cluster& c2) {
      if (c1.et > c2.et) return true;
      else return false;
   }

   inline operator uint32_t() {
#pragma HLS INLINE
      return  ((uint32_t)(this->tower_et) << 22) |
	 ((uint32_t)(this->peak_eta) << 19) |
	 ((uint32_t)(this->peak_phi) << 16) |
	 (uint32_t)this->et;
   }

#ifndef __SYNTHESIS__
   string toString() const {
      return "Cluster (" + to_string(this->peak_phi) + "," + to_string(this->peak_eta) + ")"  + ": " + to_string(this->et);
   }
#endif

   ap_uint<16> et;
   ap_uint<3>  peak_phi;
   ap_uint<3>  peak_eta;
   ap_uint<10> tower_et;
};

inline uint16_t getPeakBinOf5(const uint16_t et[5], const uint16_t etSum) {
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

Cluster computeCluster(const Crystal crystals[5][5]) {
//#pragma HLS data_pack variable=crystals
#pragma HLS ARRAY_PARTITION variable=crystals complete dim=0

   uint16_t phi_strip[5], eta_strip[5];
#pragma HLS ARRAY_PARTITION variable=phi_strip complete dim=0
#pragma HLS ARRAY_PARTITION variable=eta_strip complete dim=0

   // Compute strips
   for (size_t eta = 0; eta < 5; eta++) {
#pragma LOOP UNROLL
      eta_strip[eta] = 0;
      for (size_t phi = 0; phi < 5; phi++) {
#pragma LOOP UNROLL
	 eta_strip[eta] += crystals[eta][phi].energy;
      }
   }

   for (size_t phi = 0; phi < 5; phi++) {
#pragma LOOP UNROLL
      phi_strip[phi] = 0;
      for (size_t eta = 0; eta < 5; eta++) {
#pragma LOOP UNROLL
	 phi_strip[phi] += crystals[eta][phi].energy;
      }
   }

   // Compute tower energy
   uint16_t tet = 0;
   for (size_t eta = 0; eta < 5; eta++) {
#pragma LOOP UNROLL
      tet += eta_strip[eta];
   }

   ap_uint<10> towerEt;

   towerEt = (tet > 0x3FF) ? 0x3FF : tet;
      

   // Compute peak locations
   ap_uint<3> peakEta = getPeakBinOf5(eta_strip, tet);
   ap_uint<3> peakPhi = getPeakBinOf5(phi_strip, tet);

   // Small cluster ET is just the 3x5 around the peak
   uint16_t clusterEt = 0;
   for (int dEta = -1; dEta <= 1; dEta++) {
#pragma LOOP UNROLL
      int eta = peakEta + dEta;
      clusterEt = (eta >= 0 && eta < 5)? clusterEt + eta_strip[eta] : clusterEt;
   }

   Cluster cluster;
//#pragma HLS data_pack variable=cluster

   cluster.et = clusterEt;
   cluster.peak_eta = peakEta;
   cluster.peak_phi = peakPhi;
   cluster.tower_et = towerEt;

   return cluster;
}

void stitchNeigbours(const Cluster Ai, const Cluster Bi, Cluster &Ao, Cluster &Bo) {
   // Check that the clusters are neigbhors in eta or phi
   if (Ai.peak_eta == Bi.peak_eta || Ai.peak_phi == Bi.peak_phi ) {
      if(Ai.et > Bi.et){
	 // Merge 2 in to 1, and set 2 to remnant energy centered in tower
	 Ao.et = Ai.et + Bi.et;
	 Ao.peak_eta = Ai.peak_eta;
	 Ao.peak_phi = Ai.peak_phi;
	 Bo.et = 0 ;
	 Bo.peak_eta = 2 ;
	 Bo.peak_phi = 2 ;
      }
      else{
	 // Merge 1 in to 2, and set 1 to remnant energy centered in tower
	 Ao.et = 0;
	 Ao.peak_eta = 2;
	 Ao.peak_phi = 2;
	 Bo.et = Ai.et + Bi.et ;
	 Bo.peak_eta = Bi.peak_eta ;
	 Bo.peak_phi = Bi.peak_phi ;
      }

#ifndef __SYNTHESIS__
      cout << "Stitching cluster cluster [(" + to_string(Ai.peak_phi) + ", " + to_string(Ai.peak_eta) + ")]:"<<to_string(Ai.et)<<" with [(" <<
	 to_string(Bi.peak_phi) + ", " + to_string(Bi.peak_eta) +")]:"<<to_string(Bi.et) << endl;
#endif

   } else {
      Ao = Ai;
      Bo = Bi;
   }
}

bool stitchInEta(const Cluster in[TOWERS_IN_ETA], Cluster out[TOWERS_IN_ETA]){
//#pragma HLS data_pack variable=in
//#pragma HLS data_pack variable=out
#pragma HLS ARRAY_PARTITION variable=in complete dim=0
#pragma HLS ARRAY_PARTITION variable=out complete dim=0

   Cluster stitch_eta_step1[TOWERS_IN_ETA];
   Cluster stitch_eta_step2[TOWERS_IN_ETA];
//#pragma HLS data_pack variable=stitch_eta_step1
//#pragma HLS data_pack variable=stitch_eta_step2

   #pragma HLS ARRAY_PARTITION variable=stitch_eta_step1 complete dim=0
   #pragma HLS ARRAY_PARTITION variable=stitch_eta_step2 complete dim=0

   // Stitch in eta: first set of pairs (0,1), (2,3),...(14,15) 
   stitch_eta_step1[TOWERS_IN_ETA-1] = in[TOWERS_IN_ETA-1]; //needed for ODD eta geometry (UNcomment it)
   for (size_t teta = 0; teta < TOWERS_IN_ETA-1; teta += 2) {
//#pragma HLS dataflow
#pragma LOOP UNROLL
      if(in[teta].peak_eta == 4 && in[teta+1].peak_eta == 0 ){
	 stitchNeigbours(in[teta], in[teta+1], stitch_eta_step1[teta], stitch_eta_step1[teta+1]);
      } 
      else{
	 stitch_eta_step1[teta]   = in[teta];
	 stitch_eta_step1[teta+1] = in[teta+1];
      }
   //cout<<teta<<":"<<stitch_eta_step1[teta].toString()<<"   "<<teta+1<<":"<<stitch_eta_step1[teta+1].toString()<<std::endl;
   }
      
   // Stitch in eta: second set of pairs (1,2), (3,4),...(15,16) 
   stitch_eta_step2[0]   = stitch_eta_step1[0];
//   stitch_eta_step1[TOWERS_IN_ETA-1] = in[TOWERS_IN_ETA-1]; //needed for even geometry of ETA
   for (size_t teta = 1; teta < TOWERS_IN_ETA-1; teta += 2) {
//#pragma HLS dataflow
#pragma LOOP UNROLL
      if(in[teta].peak_eta == 4 && in[teta+1].peak_eta == 0 ){
	 stitchNeigbours(stitch_eta_step1[teta], stitch_eta_step1[teta+1], stitch_eta_step2[teta], stitch_eta_step2[teta+1]);
      } 
      else{
	 stitch_eta_step2[teta]   = stitch_eta_step1[teta];
	 stitch_eta_step2[teta+1] = stitch_eta_step1[teta+1];
      }
   }

   for (size_t teta = 0; teta < TOWERS_IN_ETA; teta++) {
#pragma LOOP UNROLL
      out[teta] = stitch_eta_step2[teta];
   }

   return true;
}

bool stitchInPhi(const Cluster inPhi1[TOWERS_IN_ETA], const Cluster inPhi2[TOWERS_IN_ETA], Cluster out[TOWERS_IN_PHI][TOWERS_IN_ETA]){
//#pragma HLS data_pack variable=inPhi1
//#pragma HLS data_pack variable=inPhi2
//#pragma HLS data_pack variable=out
#pragma HLS ARRAY_PARTITION variable=inPhi1 complete dim=0
#pragma HLS ARRAY_PARTITION variable=inPhi2 complete dim=0
#pragma HLS ARRAY_PARTITION variable=out complete dim=0

//#pragma HLS stable variable=A

   Cluster stitch_phi[TOWERS_IN_PHI][TOWERS_IN_ETA];
#pragma HLS ARRAY_PARTITION variable=stitch_phi complete dim=0

   // Stitch in phi
   for (size_t teta = 0; teta < TOWERS_IN_ETA; teta++) {
#pragma LOOP UNROLL
//#pragma HLS dataflow

      if(inPhi1[teta].peak_phi == 4 && inPhi2[teta].peak_phi == 0 ){
	 stitchNeigbours(inPhi1[teta], inPhi2[teta], stitch_phi[0][teta], stitch_phi[1][teta]);
      }
      else{
	 stitch_phi[0][teta] = inPhi1[teta];
	 stitch_phi[1][teta] = inPhi2[teta];
      }
   }

   for (size_t teta = 0; teta < TOWERS_IN_ETA; teta++) {
#pragma LOOP UNROLL
      out[0][teta] = stitch_phi[0][teta];
      out[1][teta] = stitch_phi[1][teta];
   }

   return true;
}

#endif /*!__CLUSTERFINDER_H__*/
