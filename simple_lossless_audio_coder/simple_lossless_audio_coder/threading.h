#pragma once
#include <pthread.h>
#include "stdafx.h"

extern const uint32_t buffersize;

typedef struct {
	pthread_t consumer, producer;
	pthread_mutex_t in_mutex, out_mutex;
	pthread_cond_t condc, condp;
	pthread_cond_t* cond_in;
	void(*consume)(void*, uint32_t input);
	void* output_consumer;
	uint32_t in_write_pos;
	uint32_t in_read_pos;
	uint32_t out_write_pos;
	uint32_t out_read_pos;
	bool recieve_complete;
	bool send_complete;
	bool in_cycle;
	bool out_cycle;
	uint32_t in_buffer[1024], out_buffer[1024];
} consumer_producer_buffer;

void thread_buffer_init(consumer_producer_buffer* b, pthread_cond_t* cond_in, void* output_consumer, void(*consume)(void*, uint32_t input));
bool thread_in_buffer_full(consumer_producer_buffer* b);
bool thread_out_buffer_full(consumer_producer_buffer* b);
bool thread_in_buffer_empty(consumer_producer_buffer* b);
bool thread_out_buffer_empty(consumer_producer_buffer* b);
uint32_t thread_in_buffer_read(consumer_producer_buffer* b);
void thread_in_buffer_write(consumer_producer_buffer* b, uint32_t input);
uint32_t thread_out_buffer_read(consumer_producer_buffer* b);
void thread_out_buffer_write(consumer_producer_buffer* b, uint32_t input);
void produce(consumer_producer_buffer* b);