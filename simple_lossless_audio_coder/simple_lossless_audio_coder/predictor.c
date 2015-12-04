#include "predictor.h"
#include "stdafx.h"

// Calculate cross-correlation between two signals
void cal_CR(double l[], double r[], double rb[], int n, int p) {
	int	i, j;
	double	data;
	for (i = 0; i<p; i++)
	{
		data = 0.0;
		for (j = 0; j<n - i; j++)
		{
			data += l[j] * r[i + j];
			//printf("j = %d\tdata = %d\n", j, data);
		}
		rb[i] = data;
	}
}
