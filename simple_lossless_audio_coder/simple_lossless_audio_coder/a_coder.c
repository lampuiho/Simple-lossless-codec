#include "stdafx.h"
#include "a_coder.h"

const uint32_t byte_size = 8;
const uint32_t input_size = 16;
const uint32_t precision = 32;
const uint32_t MSB_pos = 31; //try increased precision for larger CDF
const uint32_t max = -1;
const uint32_t half = 0x80000000;

void emit_output(range_coder* rc) {
	thread_out_buffer_write(&rc->buffer, rc->current_output);
	rc->current_output = 0;
}
uint32_t get_mid(range_coder* rc) {
	uint64_t real_range = (uint64_t)rc->range + 1;
	uint32_t result = (uint32_t)((*(rc->p_model->symbol_cdf) * real_range) >> precision);
	return result;
}
void emit_zero(range_coder* rc)
{
	//check overflow
	if (((++rc->processed_count) & 7) == 0)
		emit_output(rc);
	else
		rc->current_output >>= 1;
}
void emit_one(range_coder* rc)
{
	rc->current_output |= 0x80;
	emit_zero(rc);
}
void en_clear_carry_underflow(encoder* en) {
	emit_one(&en->rc);
	for (; en->underflow > 0; en->underflow--)
		emit_zero(&en->rc);
}
void en_clear_underflow(encoder* en) {
	emit_zero(&en->rc);
	for (; en->underflow > 0; en->underflow--)
		emit_one(&en->rc);
}
void en_shift(range_coder* rc) {
	rc->low <<= 1; rc->range = (rc->range << 1) | 1;
}
void en_renorm(encoder* en) {
	while (en->rc.range < half)
	{
		if (en->carry) //previously output bit is 1
		{
			en->carry = (en->rc.low >> MSB_pos) == 1;
			en_clear_carry_underflow(en);
		}
		else if (en->rc.low < half) //high and low fist bit 0
		{
			en_clear_underflow(en);
		}
		else
			en->underflow++; //first bit of low is 1, increment underflow for possible carry

		en_shift(&en->rc);
	}
}
void encode(encoder* en, uint32_t input) {
	uint8_t temp[2];
	for (uint32_t j = 0; j < input_size / byte_size; ++j)
	{
		temp[j] = input >> (byte_size * j);
		temp[j] = inpput2symbol(en->rc.s_adap, temp[j]);
		for (uint32_t i = 0; i < byte_size; ++i)
		{
			int processing_bit = temp[j] & 1;
			temp[j] >>= 1;

			uint8_t range_adj = get_mid(&(en->rc));

			if (processing_bit == 1)
			{
				en->rc.range -= range_adj;

				if (range_adj > (max - en->rc.low))
					en->carry = true;
				en->rc.low += range_adj;
			}
			else
				en->rc.range = range_adj - 1;

			en_renorm(en);
			//p_model.Add(processing_bit);
		}
	}
}
void en_consume(encoder* en) {
	pthread_mutex_lock(&en->rc.buffer.in_mutex);
	while (!en->rc.buffer.recieve_complete && thread_in_buffer_empty(&en->rc.buffer))
		pthread_cond_wait(&en->rc.buffer.condc, &en->rc.buffer.in_mutex);

	if (!thread_in_buffer_empty(&en->rc.buffer))
		encode(en, thread_in_buffer_read(&en->rc.buffer));

	pthread_cond_signal(en->rc.buffer.cond_in);
	pthread_mutex_unlock(&en->rc.buffer.in_mutex);
}
void* en_main(encoder* en) {
	while (!en->rc.buffer.recieve_complete)
	{
		en_consume(en);
	}
	pthread_exit(0);
	return NULL;
}
void de_input(decoder* de) {
	de->current_input = thread_in_buffer_read(&de->rc.buffer);
}
void de_consume(decoder* de) {
	if ((de->input_count++ & 7) == 0) {
		pthread_mutex_lock(&de->rc.buffer.in_mutex);
		while (!de->rc.buffer.recieve_complete && thread_in_buffer_empty(&de->rc.buffer))
			pthread_cond_wait(&de->rc.buffer.condc, &de->rc.buffer.in_mutex);

		if (!thread_in_buffer_empty(&de->rc.buffer))
			de_input(de);

		pthread_cond_signal(de->rc.buffer.cond_in);
		pthread_mutex_unlock(&de->rc.buffer.in_mutex);
	}
	de->code <<= 1;
	de->code |= (de->current_input & 1);
	de->current_input >>= 1;
}
void de_shift(decoder* de) {
	en_shift(&de->rc);
	de_consume(de);
}
void de_renorm(decoder* de) {
	while (de->rc.range < half)
		de_shift(de);
}
void find_symbol(decoder* de) {
	uint32_t dist;
	if (de->code >= de->rc.low)
		dist = de->code - de->rc.low;
	else
		dist = max - de->rc.low + de->code;
	uint32_t mid = get_mid(&de->rc);
	if (dist < mid)
	{
		emit_zero(&de->rc);
		de->rc.range = mid - 1;
		//p_model.Add(0);
	}
	else
	{
		emit_one(&de->rc);
		de->rc.low = de->rc.low + mid;
		de->rc.range -= mid;
		//p_model.Add(1);
	}
	de_renorm(de);
}
void* de_main(decoder* de) {
	while (de->input_count < precision && (!de->rc.buffer.recieve_complete || !thread_in_buffer_empty(&de->rc.buffer)))
		de_consume(de);
	de->code <<= (int)(precision - de->input_count);
	while (de->rc.processed_count < de->file_size)
		find_symbol(de);

	pthread_mutex_lock(&de->rc.buffer.in_mutex);
	de->rc.buffer.recieve_complete = true;
	pthread_cond_signal(de->rc.buffer.cond_in);
	pthread_mutex_unlock(&de->rc.buffer.in_mutex);
	de->rc.buffer.send_complete = true;

	return NULL;
}

void rc_init(range_coder* rc, probability_model* p, more_zero_adaptor* s, pthread_cond_t* cond_in, void* output_consumer, void(*consume)(void*, uint32_t input)) {
	rc->p_model = p;
	rc->s_adap = s;
	rc->low = 0;
	rc->range = 0;
	rc->processed_count = 0;
	rc->current_output = 0;
	thread_buffer_init(&rc->buffer, cond_in, output_consumer, consume);
}
void en_init(encoder* en, probability_model* p, more_zero_adaptor* s, pthread_cond_t* cond_in, void* output_consumer, void(*consume)(void*, uint32_t input)) {
	en->carry = false;
	en->underflow = 0;
	rc_init(&en->rc, p, s, cond_in, output_consumer, consume);
	pthread_create(&en->rc.buffer.consumer, NULL, en_main, en);
}
void de_init(decoder* de, probability_model* p, more_zero_adaptor* s, uint32_t file_size, pthread_cond_t* cond_in, void* output_consumer, void(*consume)(void*, uint32_t input)) {
	de->code = 0; de->current_input = 0; de->file_size = file_size; de->input_count = 0;
	rc_init(&de->rc, p, s, cond_in, output_consumer, consume);
	pthread_create(&de->rc.buffer.consumer, NULL, de_main, de);
}

void en_finalise(encoder* en) {
	pthread_mutex_lock(&en->rc.buffer.in_mutex);
	en->rc.buffer.recieve_complete = true;
	pthread_cond_signal(&en->rc.buffer.condc);
	pthread_mutex_unlock(&en->rc.buffer.in_mutex);

	pthread_join(en->rc.buffer.consumer, NULL);

	if (!thread_in_buffer_empty(&en->rc.buffer))
		encode(en, thread_in_buffer_read(&en->rc.buffer));

	en->underflow++;
	if (!en->carry && (en->rc.low < half))
		en_clear_underflow(en);
	else
		en_clear_carry_underflow(en);

	en->rc.current_output >>= (byte_size - en->rc.processed_count - 1);
	emit_output(&en->rc);
	en->rc.buffer.send_complete = true;
}
void de_finalise(decoder* de) {
	pthread_mutex_lock(&de->rc.buffer.in_mutex);
	de->rc.buffer.recieve_complete = true;
	pthread_cond_signal(&de->rc.buffer.condc);
	pthread_mutex_unlock(&de->rc.buffer.in_mutex);
	pthread_join(de->rc.buffer.consumer, NULL);
}

void coder_input(range_coder* rc, uint32_t input) {
	pthread_mutex_lock(&rc->buffer.in_mutex);
	while (thread_in_buffer_full(&rc->buffer) && !rc->buffer.recieve_complete)
		pthread_cond_wait(rc->buffer.cond_in, &rc->buffer.in_mutex);
	thread_in_buffer_write(&rc->buffer, input);
	pthread_cond_signal(&rc->buffer.condc);
	pthread_mutex_unlock(&rc->buffer.in_mutex);
}