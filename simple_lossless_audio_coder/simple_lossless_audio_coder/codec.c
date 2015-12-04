/*
*** TOOLS FOR LOSSLESS CODEC
*/
#include "stdafx.h"
#include "codec.h"
#include "a_coder.h"

typedef struct {
	uint32_t pos;
	char* buf;
} byte_buffer;

// Fill table containing each symbole and its number of occurences.
void fill_entropy_table(short int *buf, short int *table, int tab_size){
    int i;
    for(i = 0; i < tab_size; i++)
    {
        if(buf[i] < 0)
            table[abs(buf[i]) + 32767] += 1;
        else
            table[buf[i]] += 1;
    }
}

// From a table containing symboles and their occurence number,
// compute total entropy of a signal.
double compute_entropy_table(short int *table, int tab_size){
    int i, nb_symb;
    double entropy = 0.0;
    double probability = 0.0;

    nb_symb = 0;
    for(i = 0; i < tab_size; i++)
    {
        nb_symb += table[i];
    }

    for(i = 0; i < tab_size; i++)
    {
        if(table[i] != 0)
        {
            probability = ((double)table[i])/nb_symb;
            entropy += probability*log2((1/probability));
        }
    }
    return entropy;
}

// Perform Hamming windowing
void ham_win(double d[], double t[], int n){
    int i ;
    for(i=0; i<n; i++)
        t[i] = d[i] * (0.54 - 0.46 * cos(2*M_PI*i/n)) ;
}

// Calculate autocorrelation
void calR(double db[], double rb[], int n, int p){
    int	i, j ;
    double	data ;
    for(i=0; i<p; i++)
    {
        data = 0.0 ;
        for(j=0; j<n-i; j++)
        {
            data += db[j] * db[i+j] ;
            //printf("j = %d\tdata = %d\n", j, data);
        }
        rb[i] = data ;
    }
}

// Calculate LPC coefficients using Levinson's algorithm. Note p is the predictor order
double calA(double r[], double a[], int p){
    double	alpha, beta, k ;
    double	tmp[200] ;
    int	i, j ;
    alpha = r[0] ;
    beta = r[1] ;
    k = beta / alpha ;
    a[0] = k ;
    for(i=1; i<p; i++)
    {
        alpha *= (1 - k*k) ;
        beta = r[i+1] ;
        for(j=0; j<i; j++)
        {
            beta -= a[j] * r[i-j] ;
            //printf("r[%d] = %d\n", i-j, r[i-j]);
        }
        k = beta / alpha ;
        a[i] = k ;
        for(j=0; j<i; j++)
            tmp[j] = a[j] - k * a[i-1-j] ;
        for(j=0; j<i; j++)
            a[j] = tmp[j] ;
    }
    alpha *= (1 - k*k) ;
    return(alpha) ;
}

// Read WAV filer header and prepare file pointer for reading the data chunk
long prepareFile(FILE *fpi, WAVE *header){
    unsigned long filesize, totaldata ;
    long sampling_rate ;
    int	stereo ;
    char buf[5] ;
    int	freq ;
    int	i, n ;

    //read(fpi, header, sizeof(WAVE)) ;
    fread(header, sizeof(WAVE), 1, fpi);
    if(strncmp(header->wavefmt, "WAVEfmt ", 8))
    {
        printf("Format Error!\n") ;
        fflush(NULL) ;
        return(-1) ;
    }

    stereo = header->mode-1;
    freq = header->sampling_rate;

    if(strncmp(header->data, "data", 4))
    {
        buf[4] = 0;
        rewind(fpi);
        //read(fpi, buf, 4*sizeof(char)) ;
        fread(buf, sizeof(char), 4, fpi);
        n = sizeof(WAVE);
        while(strncmp(buf, "data", 4))
        {
            for(i=0; i<3; i++)
                buf[i] = buf[i+1] ;
            //read(fpi, &buf[3], sizeof(char));
            fread(&buf[3], sizeof(char), 1, fpi);
            n += 1 ;
        }
        //read(fpi, &(header->num_data), sizeof(long));
        fread(&(header->num_data), sizeof(long), 1, fpi);
    }
    filesize = header->num_data;
    totaldata = filesize>>2;
    return(totaldata) ;
}

int	readdata(FILE *fp, short int ibuf[], int n){
    int	m ;

    m = fread(ibuf, sizeof(short int), n, fp);
    if(m == -1)
        return(-1) ;
    return(m) ;
}

int	readdata_stereo(FILE *fp, short int ibuf[], short int sbuf[], int n){
    int	m, s, x;

    for(x = 0; x < n; x++){
        m = fread(ibuf, sizeof(short int), 1, fp);
        s = fread(ibuf, sizeof(short int), 1, fp);
        if(m == -1 || s == -1)
            return(-1) ;
    }
    return(m + s) ;
}

// Rountine for packing data in 1 to 16 bits long into codedata buffer
void subpack(unsigned char c, int nbit, int *outbitcount, int *outbytecount, unsigned char codedata[]){
    int	x ;

    if(nbit == 0)
        return ;
    x = 7 - (*outbitcount & 0x07) ;
    if(x < nbit-1)
    {
        codedata[*outbytecount] |= (unsigned char)(c >> (nbit - 1 - x)) ;
        codedata[*outbytecount + 1] |= (unsigned char)(c << (9 - nbit + x)) ;
    }
    else
        codedata[*outbytecount] |= (unsigned char)(c << (x - nbit + 1)) ;

    *outbitcount += nbit ;
    *outbytecount = *outbitcount >> 3 ;
}

// Rountine for packing data in 1 to 16 bit long into codedata buffer
void pack(unsigned long c, int nbit, int *outbitcount, int *outbytecount, unsigned char codedata[]){
    if(nbit > 16)
    {
        subpack((unsigned char)(c>>16), nbit - 16, outbitcount, outbytecount, codedata);
        nbit = 16 ;
    }
    if(nbit > 8)
    {
        subpack((unsigned char)((c>>8) & 0xff), nbit - 8,  outbitcount, outbytecount, codedata);
        nbit = 8 ;
    }
    subpack((unsigned char)(c & 0xff), nbit,  outbitcount, outbytecount, codedata) ;
}

void convert_to_double(short int bIn[], double bOut[], int w){
    int i=0, x= 0;
    for(i=0; i < w; i++)
    {
        bOut[i] = ((double)bIn[i])/((int)32768);
    }
}

void convert_to_short_int(double *in, short int *out, int taille){
    int i=0, x= 0;
    for(i=0; i < taille; i++)
    {
        out[i] = (short int)(in[i] * 0x7fff);
    }
}

void LPCCalculation(short int ibuf[], double nbuf[], double abuf[], short int pbuf[], double dpbuf[], short int ebuf[]){
        printf("Converting...");
        convert_to_double(ibuf, nbuf, W) ;            		// convert 16 bit data to floating point data
        printf("\tDone\n\n");

        printf("Getting predictor...\n");
        get_predictor(nbuf, abuf, W, N, P) ;				        // find linear predictor coefficients

        printf("\nPredicting future values... ");
        predict(abuf, &nbuf[(W-N)/2], dpbuf, N, P) ;		    // implement predictor and generate predicted signal
        printf("\tDone\n\n");

        convert_to_short_int(dpbuf, pbuf, N);               // convert predicted values back to 16-bit data.

        printf("Getting residual values... ");
        get_residual(&ibuf[(W-N)/2], pbuf, ebuf, N) ;	    // compute residual signal between 16bit original
        printf("\tDone\n\n");                               // and predicted signal
}

void get_predictor(double dBuff[], double aBuff[], int w, int n, int p){
    double cSig[W] = {0.0};
    double rSig[P] = {0.0};

    double refCoeff = 0.0;

    // Perform Hamming windowing on various signals.
    printf("\tPerforming Hamming Windowing... ");
    ham_win(dBuff, cSig, w);
    printf("\tDone\n");

    // Calculate autocorrelation of the windowed signal.
    printf("\tCalculating autocorrelation of the windowed signal... ");
    calR(&cSig[(W-N)/2], rSig, N, P);
    printf("\tDone\n");

    // Calculate lpc coefficients using Levinson algorithm.
    printf("\tCalculating LPC... ");
    refCoeff = calA(rSig, aBuff, P);

    //printf("  Variance : %f  ", refCoeff);
    printf("\tDone\n");
}

void predict(double *aBuff, double *iBuff, double *pBuff, int n, int p){
    int i, j, x;
    double tmpPBuff[N] = {0.0};

    for(i = 0; i < n; i++)
    {
        tmpPBuff[i] = 0;
        for(j = 0; j < P ; j++)
        {
            //printf("pBuff[%d] -= %f*%f\n", i, aBuff[j], iBuff[i-j-1]);
            tmpPBuff[i] -= (aBuff[j])*iBuff[(i - j) - 1];
        }
        pBuff[i] = tmpPBuff[i];
        //printf("\n\n\tiBuff[%d] - pBuff[%d] = %f + %f = %f\n\n", i, i, iBuff[i], tmpPBuff[i], iBuff[i] + tmpPBuff[i]);
        //for(x = 0; x < 50000000; x++);
    }
}

void get_residual(short int *iBuff, short int *pBuff, short int *eBuff, int n){
    int i, x;
    for(i = 0; i < n; i++)
    {
        //printf("eBuff[%d] = iBuff[%d] - pBuff[%d] = %d - %d = %d\n", i, (W-N)/2 + i,i, iBuff[i], pBuff[i], iBuff[i] - pBuff[i]);
        eBuff[i] = iBuff[i] - pBuff[i];
        //for(x = 0; x < 1000000000; x++);
    }
}

void update_output_buffers(char *outBuf, char *outCoeff, short int ebuf[], short int qabuf[], int indice){
    int x;
    // Update output buffer for encoder.
    outBuf = realloc(outBuf, 512*indice*sizeof(char));
    if(outBuf == NULL)
        return;
    for(x = (indice - 1)*2*N ; x < indice*2*N; x += 2){
        outBuf[x] = ebuf[x/2 - (indice - 1)*N] & 0xff;
        outBuf[x + 1] = ebuf[x/2 - (indice - 1)*N] >> 8;
    }

    outCoeff = realloc(outCoeff, 32*indice*sizeof(char));
    if(outCoeff == NULL)
        return;
    for(x = (indice -1)*2*P; x < indice*2*P; x += 2){
        outCoeff[x] = qabuf[x/2 - (indice - 1)*P] & 0xff;
        outCoeff[x+1] = qabuf[x/2 - (indice - 1)*P] >> 8;
    }
}

uint32_t predictor_mono(FILE *fpi, FILE *fpo, char *outBuf, char *outCoeff){

    *outBuf = malloc(512*sizeof(char));
    *outCoeff = malloc(32*sizeof(char));

    short int ibuf[W] = {0} ;		// input buffer for 16 bit data
    double nbuf[W] = {0.0} ;		// input buffer for floating arithmetic
    short int pbuf[N] = {0} ;		// buffer for the predicted signal
    double dpbuf[N] = {0} ;		    // buffer for the predicted signal
    short int ebuf[N] = {0} ;		// buffer for the residual (error) signal
    double abuf[P] = {0.0} ;		// predictor coefficients
    short int qabuf[P] = {0} ;		// predictor coefficients

    // Entropy calculation variables
    short int residual_entropy_table[65535] = {0};
    short int original_entropy_table[65535] = {0};
    double residual_entropy = 0.0;
    double original_entropy = 0.0;

    char c = 0xff;
    int i = 0, x = 0, indice = 0, n;
    clock_t begin, end;
    double time_spent;
    unsigned char codedata[MAX_CODE_SIZE] = {0} ;
    int	outbytecount = 0 ;
    int	outbitcount = 0 ;

    // Begin clock for the execution time calculation.
    begin = clock();

    n = readdata(fpi, ibuf, (W-N)) ;				// read W-N samples first to fill the overlapping window
    if(n != (W-N))
    {
        printf("\nn = %d and W-N = %d", n, W-N);
        printf("\nFlushing files");
        fflush(NULL) ;
        return 0;
    }

    for(;;)
    {
        for(x = 0; x < 1250000000; x++);
        system("cls");

        printf("SAMPLE %d\n\n", indice++);

        n = readdata(fpi, &ibuf[W-N], N) ;
        if(n != N)
        {
            // file end, need to finish off the remaining data here

            system("cls");

            residual_entropy = compute_entropy_table(residual_entropy_table, 65535);
            original_entropy = compute_entropy_table(original_entropy_table, 65535);

            // Print out summary.
            printf("******** SUMMARY ********\n\n");
            printf("FILTER ORDER : %d\n\n", P);
            printf("NUMBER OF SAMPLES ENCODED : %d\n\n", indice);
            printf("INPUT SIGNAL ENTROPY : %6f  bits per symbol\n", original_entropy);
            printf("RESIDUAL ERROR ENTROPY : %6f  bits per symbol\n", residual_entropy);
            printf("ENTROPY GAINED : %f", residual_entropy - original_entropy);

            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

            printf("\n\nExecution time : %f seconds\n\n", time_spent);
            printf("Flushing buffers\n\n");

            fclose(fpo); fclose(fpi);
            fflush(NULL) ;
            return indice;
        }


        LPCCalculation(ibuf, nbuf, abuf, pbuf, dpbuf, ebuf);

        // Quantize coefficients (this is just a conversion to 16-bit integer, not a proper quantization)
        convert_to_short_int(abuf, qabuf, P);

        // Update output buffers for entropy coder and output file.
       update_output_buffers(outBuf, outCoeff, ebuf, qabuf, indice);

    /**/
        printf("CHECK\n");
        for(x = 0; x < 42; x++)
            printf("\t[%2d]   Original Signal : %6d\t  Residual Signal = %6d\t[%d]   Original Signal : %6d\t  Residual Signal = %6d\t[%3d]   Original Signal : %6d\tResidual Signal = %6d\n", x, ibuf[(W-N)/2 + x], ebuf[x], x + 115, ibuf[(W-N)/2 + x + 115], ebuf[x + 115], x + 213, ibuf[(W-N)/2 + x + 213], ebuf[x + 213]);
    /**/

        // Copy last part of ibuf -> 50% Overlapping
        memcpy(ibuf, &ibuf[N], W);

        // Fill entropy tables
        fill_entropy_table(ebuf, residual_entropy_table, N);
        fill_entropy_table(ibuf, original_entropy_table, N);
    }
}

uint32_t predictor_stereo(FILE *fpi, FILE *fpo, char *outBuf, char *outCoeff,  char *soutBuf, char *soutCoeff){

    short int lbuf[W] = {0};        // input buffer for 16bit integer data - Left
    short int rbuf[W] = {0};        // input buffer for 16bit integer data - Right

    double nlbuf[W] = {0.0} ;		// input buffer for floating arithmetic - Left
    double nrbuf[W] = {0.0} ;		// input buffer for floating arithmetic - Right

    short int plbuf[N] = {0} ;		// buffer for the predicted signal - Left
    short int prbuf[N] = {0} ;		// buffer for the predicted signal - Right

    double dplbuf[N] = {0} ;		    // buffer for the predicted signal - Left
    double dprbuf[N] = {0} ;		    // buffer for the predicted signal - Right

    short int elbuf[N] = {0} ;		// buffer for the residual (error) signal - Left
    short int erbuf[N] = {0} ;		// buffer for the residual (error) signal - Right

    double albuf[P] = {0.0} ;		// predictor coefficients - Left
    double arbuf[P] = {0.0} ;		// predictor coefficients - Right

    short int qalbuf[P] = {0} ;		// predictor coefficients - Left
    short int qarbuf[P] = {0} ;		// predictor coefficients - Right

    // Entropy calculation variables
    short int residual_entropy_table[65535] = {0};
    short int original_entropy_table[65535] = {0};
    double residual_entropy = 0.0;
    double original_entropy = 0.0;

    char c = 0xff;
    int i = 0, x = 0, indice = 0, n;
    clock_t begin, end;
    double time_spent;
    unsigned char codedata[MAX_CODE_SIZE] = {0} ;
    int	outbytecount = 0 ;
    int	outbitcount = 0 ;

    // Begin clock for execution time calculation.
    begin = clock();

    n = readdata_stereo(fpi, rbuf, lbuf, W-N);
    if(n != 2*(W-N))
    {
        printf("\nn = %d and W-N = %d", n, W-N);
        printf("\nFlushing files");
        fflush(NULL) ;
        return 0;
    }

    for(;;)
    {
        //for(x = 0; x < 1250000000; x++);
        system("cls");
        printf("SAMPLE %d\n\n", indice++);

    n = readdata_stereo(fpi, &rbuf[N], &lbuf[N], W-N);
    if(n != 2*(W-N))
        {
            system("cls");

            // file end, need to finish off the remaining data here
            residual_entropy = compute_entropy_table(residual_entropy_table, 65535);
            original_entropy = compute_entropy_table(original_entropy_table, 65535);

            // Print out summary.
            printf("******** SUMMARY ********\n\n");
            printf("FILTER ORDER : %d\n\n", P);
            printf("NUMBER OF SAMPLES ENCODED : %d\n\n", indice);
            printf("INPUT SIGNAL ENTROPY : %6f  bits per symbol\n", original_entropy);
            printf("RESIDUAL ERROR ENTROPY : %6f  bits per symbol\n", residual_entropy);
            printf("ENTROPY GAINED : %f", residual_entropy - original_entropy);

            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

            printf("\n\nExecution time : %f seconds\n\n", time_spent);
            printf("Flushing buffers\n\n");

            fclose(fpo); fclose(fpi);
            fflush(NULL) ;
            return 0;
        }

        printf("Mixing stereo signals... ");                // Mix the two stereo signals.
        for(x = 0; x < W; x++){
            lbuf[x] = lbuf[x] + rbuf[x];
            rbuf[x] = lbuf[x] - 2*rbuf[x];
        }
        printf("\tDone\n\n");

        // Compute LPC coefficients, predict signal, and compute residual signal.
        LPCCalculation(lbuf, nlbuf, albuf, plbuf, dplbuf, elbuf);
        LPCCalculation(rbuf, nrbuf, arbuf, prbuf, dprbuf, erbuf);

        // Quantize coefficients (this is just a conversion to 16-bit integer, not a proper quantization)
        convert_to_short_int(albuf, qalbuf, P);
        convert_to_short_int(arbuf, qarbuf, P);

        // Update output buffers for entropy coder and output file.
        update_output_buffers(outBuf, outCoeff, elbuf, qalbuf, indice);
        update_output_buffers(outBuf, outCoeff, erbuf, qarbuf, indice);

        memcpy(lbuf, &lbuf[N], W);
        memcpy(rbuf, &rbuf[N], W);

        fill_entropy_table(elbuf, residual_entropy_table, N);
        fill_entropy_table(erbuf, residual_entropy_table, N);
        fill_entropy_table(lbuf, original_entropy_table, N);
        fill_entropy_table(rbuf, original_entropy_table, N);
    }
}

void write_to_buf(void* buf, uint32_t input) {
	((byte_buffer*)buf)->buf[((byte_buffer*)buf)->pos++] = input;
}

void main_coder(char *outBuf, char *outCoeff, uint32_t size, byte_buffer** out_buf_ptr) {
	probability_model p;
	more_zero_adaptor s;
	encoder en;
	uint32_t count_table[256];
	byte_buffer byte_buf;
	get_count(count_table, 255, outBuf, size*N*2);
	bin_p_init(&p, get_low_count(outBuf, size*N * 2), size * N * 2 * 8);
	s_zero_adap_init(&s, count_table);
	byte_buf.pos = 0;
	byte_buf.buf = malloc(size * N * 2 * 8 + s.max_index);
	for (int i = 0; i <= s.max_index; ++i)
		write_to_buf(&byte_buf, s.saved_indices[i]);
	en_init(&en, &p, &s, &byte_buf, write_to_buf);
	for (uint32_t i = 0; i < size * N * 2 * 8; ++i) {
		coder_input(&en.rc, outBuf[i]);
	}
	en_finalise(&en);
	*out_buf_ptr = &byte_buf;
}

void encoder_main(char input_name[], char output_name[]){

    FILE	*fpi, *fpo ;

	uint32_t size;
    char *predictor_output;
    char *predictor_coeff;
    char *s_predictor_output;
    char *s_predictor_coeff;

    if((fpi = fopen(input_name, "rb")) == NULL)
    {
        printf("Can't open %s\n", input_name) ;
        return ;
    }

    if((fpo = fopen(output_name, "wb")) == NULL)
    {
        printf("Can't open %s\n", output_name) ;
        return ;
    }

    WAVE header;

    long totaldata = prepareFile(fpi, &header);
    if(totaldata == -1)
    {
        fflush(NULL) ;
        return ;
    }

    // perhaps you need to copy header to the compressed file here
    fwrite(&header, sizeof(WAVE), 1, fpo);

    if(header.sampling_rate != 44100)
    {
        printf("Sampling rate not 44100 Hz\n") ;
        fflush(NULL) ;
        return ;
    }
    switch(header.mode)
    {
    case  1  :
		size = predictor_mono(fpi, fpo, predictor_output, predictor_coeff) ;
		byte_buffer** buf;
		main_coder(*predictor_output, *predictor_coeff, size, buf);
        // entropy encoder
        // Free pointers
        break ;
    case  2  :
        predictor_stereo(fpi, fpo, predictor_output, predictor_coeff, s_predictor_output, s_predictor_coeff) ;
        // entropy encoder
        // Free pointers
        break ;
    default  :
        printf("Only Mono and Stereo Modes are supported!") ;
        fflush(NULL) ;
        return ;
    }
    fflush(NULL) ;
}

void decode_mono(FILE *fpi, FILE *fpo){
    short int ibuf[W] = {0} ;		// input buffer for 16 bit data
    double nbuf[W] = {0.0} ;		// input buffer for floating arithmetic
    short int pbuf[N] = {0} ;		// buffer for the predicted signal
    double dpbuf[N] = {0} ;		    // buffer for the predicted signal
    short int ebuf[N] = {0} ;		// buffer for the residual (error) signal
    double abuf[P] = {0.0} ;		// predictor coefficients

    // Entropy calculation variables
    short int residual_entropy_table[65535] = {0};
    short int original_entropy_table[65535] = {0};
    double residual_entropy = 0.0;
    double original_entropy = 0.0;

    char c = 0xff;
    int i = 0, x = 0, indice = 0, n;
    clock_t begin, end;
    double time_spent;
    unsigned char codedata[MAX_CODE_SIZE] = {0} ;
    int	outbytecount = 0 ;
    int	outbitcount = 0 ;

    // Begin clock for the execution time calculation.
    begin = clock();

    n = readdata(fpi, ibuf, (W-N)) ;				// read W-N samples first to fill the overlapping window
    if(n != (W-N))
    {
        printf("\nn = %d and W-N = %d", n, W-N);
        printf("\nFlushing files");
        fflush(NULL) ;
        return ;
    }

    for(;;)
    {
        //for(x = 0; x < 1250000000; x++);
        //scanf("%d", &x);
        system("cls");

        printf("SAMPLE %d\n\n", indice++);

        n = readdata(fpi, &ibuf[W-N], N) ;
        if(n != N)
        {
            system("cls");

            // file end, need to finish off the remaining data here
            residual_entropy = compute_entropy_table(residual_entropy_table, 65535);
            original_entropy = compute_entropy_table(original_entropy_table, 65535);

            // Print out summary.
            printf("******** SUMMARY ********\n\n");
            printf("FILTER ORDER : %d\n\n", P);
            printf("NUMBER OF SAMPLES DECODED : %d\n\n", indice);
            printf("INPUT SIGNAL ENTROPY : %6f  bits per symbol\n", original_entropy);

            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

            printf("\n\nExecution time : %f seconds\n\n", time_spent);
            printf("Flushing buffers\n\n");

            fclose(fpo); fclose(fpi);
            fflush(NULL) ;
            return ;
        }

        printf("Converting...");
        convert_to_double(ibuf, nbuf, W) ;            		// convert 16 bit data to floating point data
        printf("\tDone\n\n");

        printf("Getting predictor...\n");
        get_predictor(nbuf, abuf, W, N, P) ;				        // find linear predictor coefficients

        printf("\nPredicting future values... ");
        predict(abuf, &nbuf[(W-N)/2], dpbuf, N, P) ;		    // implement predictor and generate predicted signal
        printf("\tDone\n\n");

        convert_to_short_int(dpbuf, pbuf, N);               // convert predicted values back to 16-bit data.

        printf("Getting residual values... ");
        get_residual(&ibuf[(W-N)/2], pbuf, ebuf, N) ;	    // compute residual signal between 16bit original
        printf("\tDone\n\n");                               // and predicted signal


        // Write binary output to file.
        /*
        for(i = 0; i < N; i++){
            c &= ebuf[i];
            fprintf(fpo, "%c", c);
            fprintf(fpo, "%c", ebuf[i] >> 8);
            c = 0xff;
        }
        */

        //fwrite(abuf, sizeof(short int), 16, fpo);
        //fwrite(ebuf, sizeof(short int), N, fpo);

        // need to quantize linear predictive coefficients and save them as side information
        // printf("Quantizing LPC... ");
        // quantize(abuf, qabuf, P);
        // printf("\tDone\n\n");

    /**/
        printf("CHECK\n");
        for(x = 0; x < 42; x++)
            printf("\t[%2d]   Original Signal : %6d\t  Residual Signal = %6d\t[%d]   Original Signal : %6d\t  Residual Signal = %6d\t[%3d]   Original Signal : %6d\tResidual Signal = %6d\n", x, ibuf[(W-N)/2 + x], ebuf[x], x + 42, ibuf[(W-N)/2 + x + 42], ebuf[x + 42], x + 84, ibuf[(W-N)/2 + x + 84], ebuf[x + 84]);
    /**/

        memcpy(ibuf, &ibuf[N], W);

        /*
            write(fpo, oBuff, p);
            entropy_code(ebuf, codedata, N) ;				// entropy coding
            write(fpo, codedata, outbytecount) ;			// write compressed data to file
            codedata[0] = codedata[outbytecount] ;			// reset packing code buffer
            outbitcount = outbitcount - 8*outbytecount ;
            outbytecount = 0 ;
        */

        fill_entropy_table(ebuf, residual_entropy_table, N);
        fill_entropy_table(ibuf, original_entropy_table, N);

    }
}

void decoder_main(char input_name[], char output_name[]){
    FILE	*fpi, *fpo ;

    if((fpi = fopen(input_name, "rb")) == NULL)
    {
        printf("Can't open %s\n", input_name) ;
        return ;
    }

    if((fpo = fopen(output_name, "wb")) == NULL)
    {
        printf("Can't open %s\n", output_name) ;
        return ;
    }

     // ADDED
    WAVE header;

    long totaldata = prepareFile(fpi, &header);
    if(totaldata == -1)
    {
        fflush(NULL) ;
        return ;
    }

     // perhaps you need to copy header to the compressed file here
    fwrite(&header, sizeof(WAVE), 1, fpo);

    if(header.sampling_rate != 44100)
    {
        printf("Sampling rate not 44100 Hz\n") ;
        fflush(NULL) ;
        return ;
    }

    switch(header.mode)
    {
    case  1  :
        decode_mono(fpi, fpo) ;
        break ;
    case  2  :
        //decode_stereo(fpi, fpo) ;
        break ;
    default  :
        printf("Only Mono and Stereo Modes are supported!") ;
        fflush(NULL) ;
        return ;
    }
    fflush(NULL) ;
}

