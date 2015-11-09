using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Threading.Tasks;

namespace Simple_lossless_codec
{
    public class Data
    {
        public BitArray data { get; }
        public int entropy { get; }

        public Data(BitArray input)
        {
            data = input;
            double low = 0, high = 0;
            foreach (bool bit in data)
            {
                if (bit)
                    low += 1.0;
                else
                    high += 1.0;
            }
            double total = low + high;
            low = low / (total); high = high / (total);
            entropy = (int)Math.Ceiling((-(Math.Log(low, 2) + Math.Log(high, 2))));
        }
    }

    public class EncodedData
    {
        public int o_data_length { get; }
        public BitArray data { get; }

        public EncodedData(int o_data_length, BitArray output)
        {
            this.o_data_length = o_data_length; this.data = output;
        }
    }
	
    public class BitWise_Byte_Buffer
    {
        byte buffer = 0, pos = 0;
        
        public bool push(bool input)
        {
            buffer >>= 1;

            if (input)
                buffer |= 0x80;

            if (pos++ == 7)
                return true;

            return false;
        }
        public bool read()
        {
            bool result = ((buffer >> (8 - pos--)) & 0x01) > 0;
            return result;
        }
    }
	public class Buffer<T>
    {
        const int buf_size = 1048576;
        int pos = -1;
        List<T> data = new List<T>(buf_size);

        public T Read()
        {
            var result = data[++pos];
                if (pos == buf_size - 1)
                {
                    data.Clear();
                    pos = -1;
                }
            return result;
        }
        public void Write(T input)
        {
            data.Add(input);
        }
        public bool Empty()
        {
            return pos == data.Count - 1;
        }
        public bool Full()
        {
            return data.Count == buf_size;
        }

	}
	
	public interface IDataPacketHandler<T>
	{
		void Input(T data); //for conversion to C sake, not turned into actual EventHandler
		void Finalise();
		void Initialise(IDataPacketSender<T> sender);		
	}
    public interface IDataPacketSender<T>
    {
        void RegisterPacketHandler(IDataPacketHandler<T> handler);
    }
	
    public class BinaryCoder
    {
        const uint precision = sizeof(int)*8; const uint half = (uint.MaxValue >> 1) + 1;
        const uint quarter = (uint.MaxValue >> 2) + 1; const uint quarter_3 = (uint)0x3 << (sizeof(int)*8 - 2);

        internal abstract class Binary_Coder_Worker : IDisposable
        {
            internal const int byte_size = 8;

            protected ProbabilityAdaptor adaptation_module = new BufferedBIPA();
            protected List<byte> output;
            protected byte current_byte = 0;
            protected uint position_count = 0, underflow = 0, low = 0, high = uint.MaxValue;

            internal Binary_Coder_Worker() { }
            internal Binary_Coder_Worker(ProbabilityAdaptor adap){
				adaptation_module = adap;
			}

            protected void emit_output()
            {
                //check overflow
                if (++position_count == byte_size)
                {
                    output.Add(current_byte);
                    current_byte = 0;
                    position_count = 0;
                }
                else
                    current_byte >>= 1;
            }
            protected void emit_zeros()
            {
                //lazy implementation
                while (underflow > 0)
                {
                    emit_output();
                    underflow--;
                }
            }
            protected void emit_ones()
            {
                while (underflow > 0)
                {
                    current_byte |= 0x80;
                    emit_output();
                    underflow--;
                }
            }
            protected virtual void shift_range()
            {
                low <<= 1; high = (high << 1) + 1;
            }

            protected virtual uint obtain_mid()
            {
                uint range = high - low;
                uint CDF_Low = adaptation_module.CDF(true);
                //feedback error for multiplication
                uint temp = ((range % adaptation_module.CDF_T + 1) * CDF_Low) / adaptation_module.CDF_T;
                //divide first to avoid overflow
                return low + range / adaptation_module.CDF_T * CDF_Low + temp;
            }
            public virtual void Dispose()
            {
                output = null; adaptation_module = null;
            }
        }

        internal class Binary_Encoder_Worker : Binary_Coder_Worker
        {
            Action[] obtain_range; int processed = 0;
			
			internal Binary_Encoder_Worker() : base()
			{
                obtain_range = new Action[2] { low_range, high_range };
			}

            protected void process_bit(byte input)
            {
                obtain_range[input]();
                check_range();
                processed++;
            }
            internal EncodedData run(Data target)
            {
                output = new List<byte>(target.entropy);
                foreach (bool bit in target.data)
                {
                    process_bit(Convert.ToByte(bit));
                }

                return finalise();
            }
            void check_range()
            {
                //since first bit of low always 0, and that of high always 1 when not equal
                if (high < half)//equal when high bit == 0 or low bit == 1
                {
                    //rescaling operations
                    emit_output();
                    emit_ones();
                    shift_range();
                    check_range();
                }
                else if (low >= half)
                {
                    current_byte |= 0x80; //endianess
                    emit_output();
                    emit_zeros();
                    shift_range();//2(x-half) = 2x - (maxValue + 1) <- automatic overflow
                    check_range();
                }
                else
                    check_underflow();
            }
            void check_underflow()
            {
                if (low >= quarter && high < quarter_3)
                {
                    //pre-scaling to prevent overflow
                    underflow++;
                    low = low - quarter;
                    high = high - quarter;
                    shift_range();
                    check_underflow();
                }
            }
            internal virtual EncodedData finalise()
            {
                underflow++;
                //flush last bits to get within range
                if (low >= quarter)
                {
                    current_byte |= 0x80;
                    emit_output();
                    emit_zeros();
                }
                else
                {
                    emit_output();
                    emit_ones();
                }

                if (position_count != 0) //filler
                {
                    current_byte >>= (byte)(byte_size - position_count - 1);
                    output.Add(current_byte);
                }

                BitArray output_bits = new BitArray(output.ToArray());
                EncodedData result = new EncodedData(processed, output_bits);
                return result;
            }
            void low_range()
            {
                high = obtain_mid() - 1;
                //bit_count[0]++; bit_count[1]++;//update probability for next step for each byte
                //to be replaced by other adaptative algorithms (such aspartial matching)
                adaptation_module.Add(false);
            }
            void high_range()
            {
                low = obtain_mid();
                //bit_count[1]++;
                adaptation_module.Add(true);
            }
            public override void Dispose(){
                base.Dispose();
                obtain_range = null;
            }
        }

        internal class Binary_Decoder_Worker : Binary_Coder_Worker
        {
            EncodedData target;
            uint code = 0;
            int i = 0;

            internal BitArray run(EncodedData target)
            {
                this.target = target;
                output = new List<byte>(target.o_data_length / 8);
                while (i < precision && i < target.data.Length)
                {
                    next_input();
                }
                while (output.Count * byte_size < target.o_data_length)
                {
                    Find_Symbol();
                    check_range();
                    check_underflow();
                }
                return new BitArray(output.ToArray());
            }
            void Find_Symbol()
            {
                //obtain next symbol
                uint temp_bound = obtain_mid();
                if (code < temp_bound)
                {
                    emit_output();
                    high = temp_bound - 1;
                    //add low count
                    adaptation_module.Add(false);
                }
                else
                {
                    current_byte |= 0x80;
                    emit_output();
                    low = temp_bound;
                    //add high count
                    adaptation_module.Add(true);
                }
            }
            void check_range()
            {
                if (high < half)
                {
                    shift_range();
                    check_range();
                }
                else if (low >= half)
                {
                    shift_range();
                    check_range();
                }
                else
                    check_underflow();
            }
            void check_underflow()
            {
                if (low >= quarter && high < quarter_3)
                {
                    //underflow++;
                    code = code - quarter;
                    low = low - quarter;
                    high = high - quarter;
                    shift_range();
                    check_underflow();
                }
            }
            void next_input()
            {
                code <<= 1;
                if (i < target.data.Length)
                    code += Convert.ToByte(target.data[i++]);
            }
            protected override void shift_range()
            {
                base.shift_range();
                next_input();
            }
        }
        
        public EncodedData encode(Data input)
        {
            using (Binary_Encoder_Worker worker = new Binary_Encoder_Worker())
            {
                return worker.run(input);
            }
        }

        public BitArray decode(EncodedData input)
        {
            using (Binary_Decoder_Worker worker = new Binary_Decoder_Worker())
            {
                return worker.run(input);
            }
        }
    }
	
	public class StreamBinaryEncoder<T> : BinaryCoder, IDataPacketHandler<T>
	{
        Streaming_Binary_Encoder_Worker<T> worker;
        public EncodedData result { get; set; }
        internal class Streaming_Binary_Encoder_Worker<X> : Binary_Encoder_Worker
        {
            Buffer<X> buff; bool complete = false;
            System.Threading.Thread main_loop;
            internal Streaming_Binary_Encoder_Worker() : base()
            {
                buff = new Buffer<X>();
                main_loop = new System.Threading.Thread(main);
                main_loop.Start();
            }
            internal void main() //consumer/producer pattern
            {
                output = new List<byte>(1024);
                while (true)
                {
                    System.Threading.Thread.Sleep(10);
                    lock (buff)
                    {
                        while (!complete && buff.Empty())
                            System.Threading.Monitor.Wait(buff);

                        if (complete && buff.Empty())
                            break;

                        process_bits(buff.Read());
                        System.Threading.Monitor.Pulse(buff);
                    }
                }
            }
            internal void input(X data)
            {
                lock (buff)
                {
                    while (buff.Full())
                        System.Threading.Monitor.Wait(buff);

                    buff.Write(data);
                    System.Threading.Monitor.Pulse(buff);
                }
            }
            void process_bits(dynamic input)
            {
                var bytes = BitConverter.GetBytes(input);
                BitArray temp = new BitArray(bytes);
                foreach (bool bit in temp)
                {
                    process_bit(Convert.ToByte(bit));
                }
            }
            internal override EncodedData finalise()
            {
                lock (buff)
                {
                    complete = true;
                    System.Threading.Monitor.Pulse(buff);
                }
                main_loop.Join();
                return base.finalise();
            }
            public override void Dispose()
            {
                base.Dispose();
                main_loop = null;
                buff = null;
            }
        }

        public void Initialise(IDataPacketSender<T> sender)
		{
            worker = new Streaming_Binary_Encoder_Worker<T>();
            sender.RegisterPacketHandler(this);
        }
        public void Input(T data)
		{
            worker.input(data);
        }
        public void Finalise()
		{
            result = worker.finalise();
            worker.Dispose();
            worker = null;
        }

    }
}
