#pragma once
#include "threading.h"
#include "s_adaptor.h"
#include "p_model.h"

extern const uint32_t byte_size;
extern const uint32_t input_size;
extern const uint32_t precision;
extern const uint32_t MSB_pos; //try increased precision for larger CDF
extern const uint32_t max;
extern const uint32_t half;

typedef struct {
	probability_model* p_model;
	more_zero_adaptor* s_adap;
	uint32_t low, range;
	uint32_t processed_count;
	uint8_t current_output;
	consumer_producer_buffer buffer;
} range_coder;
typedef struct {
	range_coder rc;
	bool carry;
	uint32_t underflow;
} encoder;
typedef struct {
	range_coder rc;
	uint32_t file_size; //in bit
	uint32_t code, input_count;
	uint8_t current_input;
} decoder;

void coder_input(range_coder* rc, uint32_t input);
void en_init(encoder* en, probability_model* p, more_zero_adaptor* s, pthread_cond_t* cond_in, void* output_consumer, void(*consume)(void*, uint32_t input));
void en_finalise(encoder* en);
void de_init(decoder* de, probability_model* p, more_zero_adaptor* s, uint32_t file_size, pthread_cond_t* cond_in, void* output_consumer, void(*consume)(void*, uint32_t input));
void de_finalise(decoder* de);