

#define	W	512		// window size
#define	N	256		// frame shift
#define	P	16		// linear prediction order


// WAV file header definition
typedef	struct	wave
	{
		char		riff[4] ;
		unsigned	long	NCD ;
		char		wavefmt[8] ;
		unsigned	long	N10H ;
		unsigned	short int	N1H ;
		unsigned	short int	mode ;
		unsigned	long	sampling_rate ;
		unsigned	long	output_rate ;
		unsigned	short int	bytepersample ;
		unsigned	short int	bitpersample ;
		char		data[4] ;
		unsigned	long	num_data ;
	}	WAVE ;



short int ibuf[W] = {0} ;		// input buffer for 16 bit data
double	nbuf[W] = {0.0} ;		// input buffer for floating arithmetic
short int pbuf[N] = {0} ;		// buffer for the predicted signal
short int ebuf[N] = {0} ;		// buffer for the residual (error) signal
double	abuf[P] = {0.0} ;		// predictor coefficients


unsigned	char	codedata[MAX_CODE_SIZE] = {0} ;
int	outbytecount = 0 ;
int	outbitcount = 0 ;



// tools for linear prediction

// Perform Hamming windowing 
void	Ham_win(double d[], double t[], int n)
{
	int i ;
	for(i=0; i<n; i++)
        t[i] = d[i] * (0.54 - 0.46 * cos(2*M_PI*i/n)) ;
}


// Calculate autocorrelation
void	calR(double db[], double rb[], int n, int p)
{
	int	i, j ;
	double	data ;
	for(i=0; i<p; i++)
	{
		data = 0.0 ;
		for(j=0; j<n-i; j++)
			data += db[j] * db[i+j] ;
		rb[i] = data ;
	}
}

// Calculate lpc coefficients using Levinson algorithm. Note p is the predictor order
double	calA(double r[], double a[], int p)
{
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
			beta -= a[j] * r[i-j] ;
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



// read WAV filer header and prepare file pointer for reading the data chunk
long	Prepare_File(FILE *fpi, WAVE *header)
{
	unsigned long filesize, totaldata ;
	long	sampling_rate ;
	int	stereo ;
	char	buf[5] ;
      	int	freq ;
	int	i, n ;

	read(fpi->fd, header, sizeof(WAVE)) ;
	if(strncmp(header->wavefmt, "WAVEfmt ", 8))
	{
		printf("Format Error!\n") ;
		fcloseall() ;
		return(-1) ;
	}

	stereo = header->mode-1;
	freq = header->sampling_rate ;
	if(strncmp(header->data, "data", 4))
	{
		buf[4] = 0 ;
		rewind(fpi) ;
		read(fpi->fd, buf, 4*sizeof(char)) ;
		n = sizeof(WAVE) ;
		while(strncmp(buf, "data", 4))
		{
			for(i=0; i<3; i++)
				buf[i] = buf[i+1] ;
			read(fpi->fd, &buf[3], sizeof(char)) ;
			n += 1 ;
		}
		read(fpi->fd, &(header->num_data), sizeof(long)) ;
	}
	filesize = header->num_data ;
	totaldata=filesize>>2;
	return(totaldata) ;        
}


int	readdata(FILE *fp, short int ibuf[], int n)
{
	int	m ;

	m = read(fp->fd, ibuf, n*sizeof(short int)) ;
	if(m == -1)
	      return(-1) ;
	m = m / sizeof(short int) ;
	
  	return(m) ;
}   


// rountine for packing data in 1 to 16 bits long into codedata buffer
void	subpack(unsigned char c, int nbit)
{
	int	x ;

	if(nbit==0)
		return ;
	x = 7 - (outbitcount & 0x07) ;
	if(x < nbit-1)
	{
		codedata[outbytecount] |= (unsigned char)(c >> (nbit-1-x)) ;
		codedata[outbytecount+1] |= (unsigned char)(c << (9-nbit+x)) ;
	}
	else
		codedata[outbytecount] |= (unsigned char)(c << (x - nbit + 1)) ;
	outbitcount += nbit ;
	outbytecount = outbitcount >> 3 ;
}


// rountine for packing data in 1 to 16 bit long into codedata buffer
void	pack(unsigned long c, int nbit)
{
	if(nbit > 16)
	{
		subpack((unsigned char)(c>>16), nbit-16) ;
		nbit = 16 ;
	}
	if(nbit > 8)
	{
		subpack((unsigned char)((c>>8) & 0xff), nbit-8) ;
		nbit = 8 ;
	}
	subpack((unsigned char)(c & 0xff), nbit) ;
}



void	encode_mono(FILE *fpi, FILE *fpo)
{
	int	n ;


	n = readdata(fpi, ibuf, (W-N)) ;				// read W-N samples first to fill the overlapping windoe
	if(n != (W-N))
	{
		fcloseall() ;
		return ;
	}

	for(;;)
	{
		n = readdata(fp, &ibuf[W-N], N) ;			// read N samples of 16 bit data
		if(n != N)
		{
			// file end, need to finish off the remaining data here
			fcloseall() ;
			return ;
		}

/* write your program code here

		convert(ibuf, nbuf, W) ;            			// convert 16 bit data to floating point data
		get_predictor(nbuf, W, abuf, P) ;			// find linear predictor coefficients
	        predict(abuf, &ibuf[(W-N)/2], pbuf, N) ;		// implement predictor and generate predicted signal
		get_residual(&ibuf[(W-N)/2], pbuf, ebuf, N) ;		// compute residual signal
// need to quantize linear predictove coefficeints and save them as side information
		entropy_code(ebuf, codedata, N) ;			// entropy coding
		write(fpo, codedata, outbytecount) ;			// write compressed data to file

		codedata[0] = codedata[outbytecount] ;			// reset packing code buffer
		outbitcount = outbitcount - 8*outbytecount ;
		outbytecount = 0 ;

		copy(&ibuf[N], ibuf, W-N) ;				// copy the last part of input buffer for overlapping
*/

	}
}



void	encode_stereo(FILE *fpi, FILE *fpo)
{
//	write your code here
}




void encoder_main(char input_name[], char output_name[])
{
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

	totaldata = Prepare_File(fpi, &header);
        if(totaldata == -1)
        {
	        fcloseall() ;
     		return ;
	}

// perhaps you need to copy header to the compressed file here


   if(header.sampling_rate != 44100)
   {
      printf("Sampling rate not 44100 Hz\n") ;
      fcloseall() ;
      return ;
   }


   switch(header.mode)
   {
      case  1  :  encode_mono(fpi, fpo) ;
                  break ;
      case  2  :  encode_stereo(fpi, fpo) ;
                  break ;
      default  :  printf("Only Mono and Stereo Modes are supported!") ;
                  fcloseall() ;
                  return ;
   }                     
   fcloseall() ;
}

