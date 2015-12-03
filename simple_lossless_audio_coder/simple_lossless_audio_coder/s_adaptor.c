#include "stdafx.h"
#include "s_adaptor.h"

typedef struct {
	uint8_t index;
	uint32_t value;
} index_pair;

int cmp_indexed_desc_func(const void * a, const void * b)
{
	const index_pair *elem1 = a;
	const index_pair *elem2 = b;
	if (elem1->value < elem2->value)
		return -1;
	else if (elem1->value > elem2->value)
		return 1;
	else
		return 0;
}
uint8_t get_combinations(uint8_t** output, uint8_t one_count) {
	uint8_t buffer[256]; uint8_t index = 0;
	for (int i = 0; i <= 255; ++i)
	{
		int count = 0;
		for (uint8_t j = 0; j < 8; ++j)
		{
			int temp = (i >> j);
			if ((temp & 1) == 1)
				count++;
		}
		if (count == one_count)
			buffer[index++] = i;
	}
	*output = malloc(index);
	for (int i = 0; i <= index; ++i)
		*output[index] = buffer[i];

	return index;
}
void get_lookup(uint8_t output[256])
{
	uint8_t index = 0;
	uint8_t** combinations = malloc(1);
	for (int i = 0; i <= 8; i++)
	{
		uint8_t count = get_combinations(combinations, i);
		for (int j = 0; j < count; ++j)
			output[index++] = (*combinations)[j];
		free(*combinations);
	}
	free(combinations);
}
void s_zero_adap_init(more_zero_adaptor* s_adap, uint32_t* count_table) {
	uint8_t new_symbol_lookup[256];
	get_lookup(new_symbol_lookup);
	int new_symbol = 0;
	index_pair temp[256];
	for (int i = 0; i < 256; i++) {
		temp[i].index = i;
		temp[i].value = count_table[i];
	}
	qsort(temp, 256, sizeof(index_pair), cmp_indexed_desc_func);
	s_adap->saved_indices = malloc(256);
	for (int i = 0; i < 256; i++)
	{
		uint8_t key = temp[i].index;
		s_adap->symbol_coeff_table[key] = new_symbol_lookup[new_symbol];
		s_adap->reverse_table[new_symbol_lookup[new_symbol]] = key;
		s_adap->saved_indices[new_symbol++] = key;
	}
	s_adap->max_index = new_symbol - 1;
}
void s_zero_adap_init_from_saved(more_zero_adaptor* s_adap, uint8_t* saved_indices, uint8_t max_index) {
	uint8_t new_symbol_lookup[256];
	get_lookup(new_symbol_lookup);
	for (int i = 0; i <= max_index; ++i)
	{
		s_adap->symbol_coeff_table[saved_indices[i]] = new_symbol_lookup[i];
		s_adap->reverse_table[new_symbol_lookup[i]] = saved_indices[i];
	}
	s_adap->saved_indices = saved_indices;
	s_adap->max_index = max_index;
}
uint32_t symbol2output(more_zero_adaptor* s_adap, uint32_t symbol) {
	return s_adap->reverse_table[symbol];
}
uint32_t inpput2symbol(more_zero_adaptor* s_adap, uint32_t input) {
	return s_adap->symbol_coeff_table[input];
}