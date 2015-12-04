#pragma once
#include "stdafx.h"

const uint32_t shiftvalue;

typedef struct {
	uint32_t* symbol_cdf;
	uint32_t maxIndex;
} probability_model;
typedef struct {
	const uint32_t symbol_update_f;
	probability_model p_model;
	uint32_t update_count;
	uint32_t cumulatived_update;
} adaptive_model;

void bin_p_init(probability_model* p_model, uint32_t low_count, uint32_t total_count);
void bin_adapt(adaptive_model* p_model, uint32_t symbol);
uint32_t bin_cdf(probability_model* p_model);
uint32_t get_entropy(uint32_t* count, uint32_t max_symbol);
void get_count(uint32_t* count, uint32_t max_symbol, char* input, uint32_t size);
uint32_t get_low_count(char* input, uint32_t size);