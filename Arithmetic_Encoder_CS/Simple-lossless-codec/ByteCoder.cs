using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Simple_lossless_codec
{
    public abstract class ByteCoder: IDataPacketHandler<byte>
    {
        protected const int byte_size = 8;
        protected const int precision = sizeof(uint) * byte_size;
        protected const int MSB_pos = precision - byte_size; //try increased precision for larger CDF
        protected const uint max_range = 0x1 << MSB_pos;
        protected const uint half = (uint)0x1 << (precision - 1);

        protected ProbabilityAdaptor p_adap = new IncrementProbabilityAdaptor(byte.MaxValue);
        protected uint low = 0, range = uint.MaxValue;
        protected bool complete = false;
        protected Buffer<byte> input_buffer = new Buffer<byte>();
        protected System.Threading.Thread main_loop;

        public List<byte> result;

        public virtual void Finalise()
        {
            lock (input_buffer)
            {
                complete = true;
                System.Threading.Monitor.Pulse(input_buffer);
            }
            main_loop.Join();
        }
        public virtual void Input(byte data)
        {
            lock (input_buffer)
            {
                while (input_buffer.Full())
                    System.Threading.Monitor.Wait(input_buffer);

                input_buffer.Write(data);
                System.Threading.Monitor.Pulse(input_buffer);
            }
        }
        public virtual void Initialise(IDataPacketSender<byte> sender)
        {
            main_loop = new System.Threading.Thread(main);
            main_loop.Start();
            sender.RegisterPacketHandler(this);
        }
        protected virtual void main() //consumer/producer pattern
        {
            result = new List<byte>(1024);
            while (!complete)
            {
                consume();
            }
        }
        protected virtual void consume()
        {
            lock (input_buffer)
            {
                while (!complete && input_buffer.Empty())
                    System.Threading.Monitor.Wait(input_buffer);

                if (!input_buffer.Empty())
                    process_byte(input_buffer.Read());

                System.Threading.Monitor.Pulse(input_buffer);
            }
        }
        protected abstract void process_byte(byte input);
        protected void emit_byte(byte output)
        {
            result.Add(output);
        }
    }

    public class ByteEncoder : ByteCoder
    {
        int processed_count = 0;
        uint underflow = 0;
        byte buffer;
        bool carry = false;

        public override void Finalise()
        {
            base.Finalise();

            while(!input_buffer.Empty())
                process_byte(input_buffer.Read());

            if (carry)
            {
                clear_carry_underflow();
                buffer = (byte)(low >> MSB_pos);
                low <<= byte_size;
                if (low >= half)
                {
                    if (buffer == 0xFF)
                        result[result.Count - 1] += 1;
                    emit_byte((byte)(buffer +1));
                    emit_byte(0);
                }
                else
                {
                    emit_byte(buffer);
                    emit_byte(0xFF);
                }
            }
            else
            {
                byte temp = (byte)(low >> MSB_pos);
                low <<= byte_size;
                if (low >= half)
                {
                    if (buffer == 0xFF)
                        clear_carry_underflow();
                    else
                        clear_underflow();

                    emit_byte((byte)(temp + 1));
                    emit_byte(0);
                }
                else
                {
                    clear_underflow();
                    emit_byte(temp);
                    emit_byte(0xFF);
                }
            }

            result.RemoveAt(0);
        }

        protected override void process_byte(byte input)
        {
            uint range_scaling = range / p_adap.CDF_T;
            uint err = (range % p_adap.CDF_T) + 1;

            if (input == 0)
            {
                uint high_temp = p_adap.CDF(input + 1) * range_scaling - 1;
                uint high_err_temp = (p_adap.CDF(input + 1) * err) / p_adap.CDF_T;
                uint high_range_adj = high_temp + high_err_temp;
                range = high_range_adj;
            }
            else
            {
                uint temp = p_adap.CDF(input) * range_scaling;
                uint err_temp = (p_adap.CDF(input) * err) / p_adap.CDF_T;
                uint range_adj = temp + err_temp;
                if (input == byte.MaxValue)
                    range = range - range_adj;
                else
                {
                    uint high_temp = p_adap.CDF(input + 1) * range_scaling - 1;
                    uint high_err_temp = (p_adap.CDF(input + 1) * err) / p_adap.CDF_T;
                    uint high_range_adj = high_temp + high_err_temp;
                    range = high_range_adj - range_adj;
                }

                try
                {
                    low = checked(low + range_adj);
                }
                catch
                {
                    low = unchecked(low + range_adj);
                    carry = true;
                }
            }
            //range = range_scaling * p_adap.P(input) + (p_adap.P(input) * err) / p_adap.CDF_T - 1;

            rescale();
            p_adap.Add(input);
            processed_count++;
        }
        void rescale()
        {
            if (range < max_range)
            {
                //range is small enough for first 8 bits of low to keep same value or to be addition of 1 depanding on next input range
                //put it into buffer and output last buffered value
                if (carry)
                {
                    carry = false;
                    clear_carry_underflow();
                    buffer = (byte)(low >> MSB_pos);
                }
                else if (low < ((uint)0xFF << MSB_pos)) //ensure buffered output is not 0xFF since a carry means 0x100 (i.e. carried to the next buffered output)
                {
                    clear_underflow();
                    buffer = (byte)(low >> MSB_pos);
                }
                else
                    underflow++; //case 0xFF, wait until next carry occurs

                range <<= byte_size; range |= 0xFF;
                low <<= byte_size;
                rescale();
            }
        }
        void clear_carry_underflow()
        {
            emit_byte((byte)(buffer+1));
            for (; underflow > 0; underflow--)
                emit_byte(0);
        }
        void clear_underflow()
        {
            emit_byte((byte)(buffer));
            for (; underflow > 0; underflow--)
                emit_byte(0xFF);
        }
    }

    public class ByteDecoder : ByteCoder
    {
        //original file size in no. of bytes
        readonly int file_size; readonly List<byte> original;
        int read_pos = 0; uint code = 0;

        public ByteDecoder(int file_size, byte[] original) : base()
        {
            this.file_size = file_size;
            this.original = new List<byte>(original);
        }
        protected override void main()
        {
            result = new List<byte>(file_size);
            while (read_pos < 4 && !complete)
            {
                consume();
            }
            while (result.Count < file_size)
            {
                Find_Symbol();
                rescale();
            }
        }
        protected override void consume()
        {
            code <<= 8; read_pos++;
            base.consume();
        }
        protected override void process_byte(byte input)
        {
            code |= input;
        }
        void Find_Symbol()
        {
            uint last_temp_dist = 0;
            uint range_scaling = range / p_adap.CDF_T;
            uint err = (range % p_adap.CDF_T) + 1;
            uint dist;
            if (code >= low)
                dist = code - low;
            else
                dist = (uint)(uint.MaxValue - low + code);
            for (int i = 1; i <= p_adap.SymbolMax; ++i)
            {
                uint temp = p_adap.CDF(i) * range_scaling;
                uint err_temp = (p_adap.CDF(i) * err) / p_adap.CDF_T;
                uint temp_dist = temp - 1 + err_temp;
                if (dist <= temp_dist)
                {
                    emit_byte((byte)(i - 1));
                    if (original[result.Count - 1] != i - 1)
                        throw new Exception("mismatch");
                    p_adap.Add(i - 1);
                    low = low + last_temp_dist;
                    range = temp_dist - last_temp_dist;
                    return;
                }
                last_temp_dist = temp_dist + 1;
            }
            throw new Exception("Symbol not found");
        }
        void rescale()
        {
            if (range <= (max_range - 1))
            {
                range <<= byte_size; range |= 0xFF;
                low <<= byte_size;
                consume();
                rescale();
            }
        }
    }
}
