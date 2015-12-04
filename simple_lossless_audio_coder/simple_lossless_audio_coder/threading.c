#include "stdafx.h"
#include "threading.h"

const uint32_t buffersize = 1024;

bool thread_in_buffer_full(consumer_producer_buffer* b) {
	return (b->in_cycle && b->in_write_pos == b->in_read_pos);
}
bool thread_out_buffer_full(consumer_producer_buffer* b){
	return (b->out_cycle && b->out_write_pos == b->out_read_pos);
}
bool thread_in_buffer_empty(consumer_producer_buffer* b) {
	return (!b->in_cycle && b->in_write_pos == b->in_read_pos);
}
bool thread_out_buffer_empty(consumer_producer_buffer* b) {
	return (!b->out_cycle && b->out_write_pos == b->out_read_pos);
}
uint32_t thread_in_buffer_read(consumer_producer_buffer* b) {
	uint32_t result = b->in_buffer[b->in_read_pos++];
	if (b->in_read_pos == buffersize)
	{
		b->in_read_pos = 0;
		b->in_cycle = false;
	}
	return result;
}
void thread_in_buffer_write(consumer_producer_buffer* b, uint32_t input) {
	b->in_buffer[b->in_write_pos++] = input;
	if (b->in_write_pos == buffersize)
	{
		b->in_write_pos = 0;
		b->in_cycle = true;
	}
}
uint32_t thread_out_buffer_read(consumer_producer_buffer* b) {
	uint32_t result = b->out_buffer[b->out_read_pos++];
	if (b->out_read_pos == buffersize)
	{
		b->out_read_pos = 0;
		b->out_cycle = false;
	}
	return result;
}
void thread_out_buffer_write(consumer_producer_buffer* b, uint32_t input) {
	b->out_buffer[b->out_write_pos++] = input;
	if (b->out_write_pos == buffersize)
	{
		b->out_write_pos = 0;
		b->out_cycle = true;
	}
}
void produce(consumer_producer_buffer* b) {
	while (!b->send_complete && !thread_in_buffer_empty(b))
	{
		pthread_mutex_lock(b->out_mutex);

		while (!b->send_complete && thread_in_buffer_empty(b))
			pthread_cond_wait(b->condp, b->out_mutex);

		if (!thread_in_buffer_empty(b))
			b->consume(b->output_consumer, thread_out_buffer_read(b));

		pthread_cond_signal(b->condp);
		pthread_mutex_unlock(b->out_mutex);
	}
	pthread_exit(0);
}

void thread_buffer_init(consumer_producer_buffer* b, void* output_consumer, void(*consume)(void*, uint32_t input)) {
	pthread_cond_init(&b->condc, NULL);
	pthread_cond_init(&b->condp, NULL);
	pthread_mutex_init(&b->in_mutex, NULL);
	pthread_mutex_init(&b->out_mutex, NULL);
	b->output_consumer = output_consumer;
	b->consume = consume;
	b->in_read_pos = 0; b->out_read_pos = 0;
	b->in_write_pos = 0; b->out_write_pos = 0;
	b->recieve_complete = false;
	b->send_complete = false;
	b->in_cycle = false;
	b->out_cycle = false;
	pthread_create(&b->producer, NULL, produce, b);
}