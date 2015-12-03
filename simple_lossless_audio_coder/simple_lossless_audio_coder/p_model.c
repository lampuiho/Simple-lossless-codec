#include "stdafx.h"
#include "p_model.h"
#include "math.h"

void bin_p_init(probability_model* p_model, uint32_t low_count, uint32_t total_count) {
	p_model->maxIndex = 0;
	p_model->symbol_cdf = malloc(1);
	*(p_model->symbol_cdf) = (uint32_t)(((uint64_t)low_count << shiftvalue) / total_count);
}


uint32_t bin_cdf(probability_model* p_model) {
	return *(p_model->symbol_cdf);
}

uint32_t get_entropy(uint32_t* count, uint32_t max_symbol) {
	double total = 0;
	for (uint32_t i = 0; i <= max_symbol; ++i)
		total += count[i];
	double entropy = 0;
	for (uint32_t i = 0; i <= max_symbol; ++i) {
		if (count[i] > 0)
			entropy += count[i] * log2(total / count[i]);
	}
	return (uint32_t)ceil(entropy);
}

void get_count(uint32_t* count, uint32_t max_symbol, uint32_t* input, uint32_t size)
{
	for (uint32_t i = 0; i < max_symbol; ++i)
		count[i] = 0;
	for (uint32_t i = 0; i < size; ++i)
		count[input[i]]++;
}

uint32_t get_low_count(uint8_t* input, uint32_t size)
{
	uint32_t low = 0;
	for (uint32_t i = 0; i < size; ++i)
	{
		uint8_t temp = input[i];
		for (uint32_t j = 0; j < 8; ++j)
		{
			if ((temp & 1) == 0)
				low++;
			temp >>= 1;
		}
	}
	return low;
}