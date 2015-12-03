/*
*** TOOLS FOR LOSSLESS CODEC
*/

#ifndef CODEC_H_INCLUDED
#define CODEC_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define	W	512		// window size
#define	N	256		// frame shift
#define	P	32		// linear prediction order

#define MAX_CODE_SIZE 272 // 256 entropy coded residual error + 16 LPC coefficients.

typedef	struct wave
{
    char riff[4] ;
    unsigned long NCD ;
    char wavefmt[8] ;
    unsigned long N10H ;
    unsigned short int N1H ;
    unsigned short int mode ;
    unsigned long sampling_rate ;
    unsigned long output_rate ;
    unsigned short int bytepersample ;
    unsigned short int bitpersample ;
    char data[4] ;
    unsigned long num_data ;
} WAVE ;

// TOOLS FOR LINEAR PREDICTION
void ham_win(double d[], double t[], int n);
void calR(double db[], double rb[], int n, int p);
double calA(double r[], double a[], int p);
void convert_to_double(short int bIn[], double bOut[], int w);
void convert_to_short_int(double *in, short int *out, int taille);
void get_predictor(double dBuff[], double aBuff[], int w, int n, int p);
void predict(double *aBuff, double *iBuff, double *pBuff, int n, int p);
void get_residual(short int *iBuff, short int *pBuff, short int *eBuff, int n);
void quantize(double *iBuff, short int *oBuff, int p);

// TOOLS FOR FILE READING
long prepareFile(FILE *fpi, WAVE *header);
int	readdata(FILE *fp, short int ibuf[], int n);
int	readdata_stereo(FILE *fp, short int ibuf[], short int sbuf[], int n);

// TOOLS FOR DATA PACKING
void subpack(unsigned char c, int nbit, int *outbitcount, int *outbytecount, unsigned char codedata[]);
void pack(unsigned long c, int nbit, int *outbitcount, int *outbytecount, unsigned char codedata[]);

// ENCODER FUNCTIONS
void encode_mono(FILE *fpi, FILE *fpo);
void encode_stereo(FILE *fpi, FILE *fpo);
void encoder_main(char input_name[], char output_name[]);

// DECODER FUNCTIONS
void decode_mono(FILE *fpi, FILE *fpo);
void decode_stereo(FILE *fpi, FILE *fpo);
void decoder_main(char input_name[], char output_name[]);

// ADD TOOLS FOR PERFORMANCE MEASUREMENTS.
void fill_entropy_table(short int *buf, short int *table, int tab_size);
double compute_entropy_table(short int *table, int tab_size);

#endif // CODEC_H_INCLUDED
