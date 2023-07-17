/*
  Software defined radio wideband FM receiver

  Copyright (c) 2013, Cosmin Gorgovan
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of the project's author nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>

#include "rtl-sdr.h"

//#define DEFAULT_SAMPLE_RATE   2560000 //500000
#define DEFAULT_SAMPLE_RATE   1000000
#define DEFAULT_DEVICE_INDEX  0
#define BUFFER_SIZE           (2048*2)
#define FFT_SIZE              (BUFFER_SIZE/2)
#define DEFAULT_FREQ          97700000
#define CAPTURE_COUNT         (64*24)

FILE *record_file = NULL;
rtlsdr_dev_t *device;

typedef int (*get_samples_f)(uint8_t *buffer, int buffer_len);

// File reader
int get_samples_file(uint8_t *buffer, int buffer_len) {
  int len;
  len = fread(buffer, sizeof(uint8_t), buffer_len, record_file);
  
  return len == buffer_len;
}

// rtl-sdr reader
int get_samples_rtl_sdr(uint8_t *buffer, int buffer_len) {
  int len;
  rtlsdr_read_sync(device, buffer, buffer_len, &len);
  
  return len == buffer_len;
}


void receive(get_samples_f get_samples , int buffer_size, int capture_count, int sample_rate) {
  uint8_t buffer[BUFFER_SIZE];
  int len;
  double complex sample;
  double complex product;
  double complex prev_sample = 0 + 0 * I;
  
  int fft_size = buffer_size/2;
  int16_t output_buffer[fft_size];
  double fft_in[fft_size];
  double fft_out[fft_size];
  int cc=0;
  int inCount = 0;
  int outCount = 0;
  int overflowCount = 0;
  fftw_complex *filter_mid = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (fft_size));
  fftw_plan plan_filter_in = fftw_plan_dft_r2c_1d(fft_size, fft_in, filter_mid, FFTW_ESTIMATE);
  
  fftw_complex *filter_mid2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (fft_size));
  fftw_plan plan_filter_out = fftw_plan_dft_c2r_1d(fft_size, filter_mid2, fft_out, FFTW_ESTIMATE);
  
static float min=999999;
static float max=-999999;
  while(get_samples(buffer, buffer_size) == 1) {
    /* Compute amplitudes in the time domain
       from phase difference between succesive samples */
    for (int i = 0; i < buffer_size; i+=2) {
      sample = (((int)buffer[i]) - 127) + (((int)buffer[i+1]) - 127) * I;
      product = sample * conj(prev_sample);
      fft_in[i/2] = atan2(cimag(product), creal(product));
      prev_sample = sample;
	++inCount;
    }
    
    fftw_execute(plan_filter_in);
    filter_mid2[0] = 0;
    int cutoff = 20000*fft_size/sample_rate;//(fft_size)/28;
    for (int i = 1; i < cutoff;i++) {
      filter_mid2[i] = filter_mid[i];
      filter_mid2[fft_size-i] = conj(filter_mid2[i]);
    }
    for (int i = cutoff; i < fft_size/2;i++) {
      filter_mid2[i] = 0 + 0 *I;
      filter_mid2[fft_size-i] = conj(filter_mid2[i]);
    }
    
    fftw_execute(plan_filter_out);
    
    int b = 0;
static float sum=0.;
static int oc=0;
static int div = 0;
	float fix = 10000.0/fft_size;
    for (int i = 0; i < fft_size; i++) {
float f = fft_out[i];
if(f>max) max=f;
if(f<min) min=f;
	sum+=f;
	++oc;
	div += 44100;
	if(div<sample_rate) continue;
	div-=sample_rate;
	int v = sum/oc*fix;
	oc=0.;
	sum=0.;
	if(v>32767) {v=32767;++overflowCount;}
	if(v<-32768) {v=-32768;++overflowCount;}
      output_buffer[b++] = v;
    }
	outCount += b;

    fwrite(output_buffer, sizeof(int16_t), b, stdout);
	if(capture_count && ++cc>=capture_count) break;
  }
fprintf(stderr, "inCount=%d, outCount=%d, overflowCount=%d\n", inCount, outCount, overflowCount);
fprintf(stderr, "min=%f, max=%f\n", min, max);
}

void print_usage() {
  fprintf(stderr,
"Usage: fm_receiver [OPTIONS]\n\n\
Valid options:\n\
  -r <file>         Use recorded data from <file> instead of an rtl-sdr device\n\
  -f <frequency>    Frequency to tune to, in Hz (default: %.2f MHz)\n\
  -d <device_index> Rtl-sdr device index (default: 0)\n\
  -h                Show this\n",
  DEFAULT_FREQ/1000000.0);
}

void main(int argc, char **argv) {
  int opt;
  char *end_ptr;

  /* SDR settings */
  uint32_t device_index = DEFAULT_DEVICE_INDEX;
  uint32_t sample_rate = DEFAULT_SAMPLE_RATE;
  int buffer_size = BUFFER_SIZE;
  uint32_t freq = DEFAULT_FREQ; 

  /* Parse command line arguments */
  errno = 0;
  while ((opt = getopt(argc, argv, "f:d:r:h")) != -1) {
    switch (opt) {
      case 'r':
        record_file = fopen(optarg, "r");
        if (record_file == NULL) {
          fprintf(stderr, "Error opening record file: %s\n", strerror(errno));
          exit(EXIT_FAILURE);
        }
        break;
      case 'f':
        freq = strtol(optarg, &end_ptr, 10);
        if (errno != 0 || *end_ptr != '\0') {
          fprintf(stderr, "Invalid frequency specified.\n");
          print_usage();
          exit(EXIT_FAILURE);
        }
        break;
      case 'd':
        device_index = strtol(optarg, &end_ptr, 10);
        if (errno != 0 || *end_ptr != '\0') {
          fprintf(stderr, "Invalid device index specified.\n");
          print_usage();
          exit(EXIT_FAILURE);
        }
        break;
      case 'h':
        print_usage();
        exit(EXIT_SUCCESS);
        break;
      default:
        print_usage();
        exit(EXIT_FAILURE);
    }
  }

  if (record_file != NULL) {
    receive(get_samples_file, buffer_size, 0, sample_rate);
  } else {
    /* Open device and set it up */
    if (rtlsdr_get_device_count() - 1 < device_index) {
      fprintf(stderr, "Device %d not found\n", device_index);
      print_usage();
      exit(EXIT_FAILURE);
    }
    
    if (rtlsdr_open(&device, device_index) != 0) {
      fprintf(stderr, "Error opening device %d: %s", device_index, strerror(errno));
      print_usage();
      exit(EXIT_FAILURE);
    }
int gains[1024];
int count = rtlsdr_get_tuner_gains(device, gains);
fprintf(stderr, "rtlsdr_get_tuner_gains returned %d\n", count);
//int i;for(i=0;i<count;++i) fprintf(stderr, "%2d: %d\n", i, gains[i]);
    // Tuner (ie. E4K/R820T/etc) auto gain
    rtlsdr_set_tuner_gain_mode(device, 0); // 0 = auto, 1 = manual
//	rtlsdr_set_tuner_gain(device, gains[count>>1]);
    // RTL2832 auto gain off
    rtlsdr_set_agc_mode(device, 0+1);
    rtlsdr_set_sample_rate(device, sample_rate);
    rtlsdr_set_center_freq(device, freq);
    fprintf(stderr, "Frequency set to %d\n", rtlsdr_get_center_freq(device));
    
    // Flush the buffer
    rtlsdr_reset_buffer(device);
    receive(get_samples_rtl_sdr, buffer_size, CAPTURE_COUNT, sample_rate);
  }
}
