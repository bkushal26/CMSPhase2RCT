#include "algo_top_parameters.h"
#include "algo_top.h"
#include <algorithm>
#include <utility>

using namespace std;
using namespace algo;

// Each input link carries the information of a 5x5 region
// Last 16-bits are reserved for CRC
void processInputData(hls::stream<axiword> &link, ap_uint<64> words[N_WORDS_PER_FRAME]) {
#pragma HLS INLINE
#pragma HLS PIPELINE II=6
#pragma HLS INTERFACE axis port=link
 input_word: for (size_t i = 0; i < N_WORDS_PER_FRAME; i++) {
#ifndef __SYNTHESIS__
    // Avoid simulation warnings
    if (link.empty()) continue;
#endif
    words[i] = link.read().data;
  }
}

void processOutputData(hls::stream<axiword> &link, ap_uint<64> words[N_WORDS_PER_FRAME]) {
#pragma HLS INLINE
#pragma HLS PIPELINE II=6
#pragma HLS INTERFACE axis port=link
 output_word: for (size_t i = 0; i < N_OUTPUT_WORDS_PER_FRAME; i++) {
    axiword r;
    r.user = 0;
    if (i == (N_OUTPUT_WORDS_PER_FRAME - 1)) {
      r.last = 0;
    }
    else {
      r.last = 1;
    }
    r.data = words[i];
    link.write(r);
  }
}

void algo_top(hls::stream<axiword> link_in[N_INPUT_LINKS], hls::stream<axiword> link_out[N_OUTPUT_LINKS]) {
#pragma HLS INTERFACE axis port=link_in
#pragma HLS INTERFACE axis port=link_out

  // Step First: Unpack
  ap_uint<64> input[N_INPUT_LINKS][N_WORDS_PER_FRAME];
#pragma HLS ARRAY_PARTITION variable=input complete
 input: for (size_t i = 0; i < N_INPUT_LINKS; i++) {
#pragma HLS UNROLL
    processInputData(link_in[i], input[i]);
  }

  // Pack the outputs

  ap_uint<64> output[N_OUTPUT_LINKS][N_WORDS_PER_FRAME];
#pragma HLS ARRAY_PARTITION variable=output complete
  ap_uint<64> ecalData;
 pack_i: for (size_t i = 0; i < N_OUTPUT_LINKS; i++) {
#pragma HLS UNROLL
  pack_j: for (size_t j = 0; j < N_OUTPUT_WORDS_PER_FRAME; j++) {
#pragma HLS UNROLL
      ap_uint<64> towerData = input[i][j];
      output[i][j] = towerData;
    }
  }

  // Step Last: Write the output

 output: for (size_t i = 0; i < N_OUTPUT_LINKS; i++) {
#pragma HLS UNROLL
    processOutputData(link_out[i], output[i]);
  }

}
