#pragma once
#include <stdint.h>
#include <stdlib.h>

typedef struct {
	uint8_t symbol_coeff_table[256];
	uint8_t reverse_table[256];
	uint8_t* saved_indices;
	uint8_t max_index;
} more_zero_adaptor;

void s_zero_adap_init(more_zero_adaptor* s_adap, uint32_t* count_table);
void s_zero_adap_init_from_saved(more_zero_adaptor* s_adap, uint8_t* saved_indices, uint8_t max_index);
uint32_t symbol2output(more_zero_adaptor* s_adap, uint32_t symbol);
uint32_t inpput2symbol(more_zero_adaptor* s_adap, uint32_t input);
