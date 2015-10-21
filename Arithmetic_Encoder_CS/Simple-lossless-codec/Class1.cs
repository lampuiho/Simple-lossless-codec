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

    public class BinaryCoder
    {
        const uint precision = sizeof(int)*8; const uint half = (uint.MaxValue >> 1) + 1;
        const uint quarter = (uint.MaxValue >> 2) + 1; const uint quarter_3 = (uint)0x3 << (sizeof(int)*8 - 2);

        internal class Coder_Worker : IDisposable
        {
            protected List<int> output;
            protected uint[] bit_count = new uint[2] { 3, 4 }; //keep size of bit_count half of sizeof(uint)*8 or overflow can occur
            protected uint current_int = 0, position_count = 0, underflow = 0, low = 0, high = uint.MaxValue;

            protected void emit_output()
            {
                //check overflow
                if (++position_count == sizeof(int)*8)
                {
                    output.Add(unchecked((int)current_int));
                    current_int = 0;
                    position_count = 0;
                }
                else
                    current_int >>= 1;
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
                    current_int |= half;
                    emit_output();
                    underflow--;
                }
            }
            protected virtual void shift_range()
            {
                low <<= 1; high = (high << 1) + 1;
            }
            public virtual void Dispose()
            {
                output = null; bit_count = null;
            }
        }

        internal class Binary_Encoder_Worker : Coder_Worker
        {
            Data target;
            Action[] obtain_range;
            
            internal Binary_Encoder_Worker(Data input)
            {
                this.target = input;
                obtain_range = new Action[2] { low_range, high_range };
            }

            internal EncodedData run()
            {
                output = new List<int>(target.entropy);
                foreach (bool bit in target.data)
                {
                    obtain_range[Convert.ToByte(bit)]();
                    check_range();
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
                    current_int |= half; //endianess
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
            EncodedData finalise()
            {
                underflow++;
                //flush last bits to get within range
                if (low >= quarter)
                {
                    current_int |= half;
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
                    current_int >>= (int)(sizeof(int)*8 - position_count - 1);
                    output.Add(unchecked((int)current_int));
                }

                BitArray output_bits = new BitArray((int[])(object)output.ToArray());
                EncodedData result = new EncodedData(target.data.Length, output_bits);
                return result;
            }
            void low_range()
            {
                long dif = (long)high - (long)low + 1;
                //divide first to avoid overflow
                high = (uint)(low + dif / bit_count[1] * bit_count[0] + (dif % bit_count[1]) * bit_count[0] - 1);
                //feedback error for multiplication
                //bit_count[0]++; bit_count[1]++;//update probability for next step for each byte
                                               //to be replaced by other adaptative algorithms (such aspartial matching)
            }
            void high_range()
            {
                long dif = (long)high - (long)low + 1;
                low = (uint)(low + dif / bit_count[1] * bit_count[0] + (dif % bit_count[1]) * bit_count[0]);
                //bit_count[1]++;
            }
            public override void Dispose(){
                base.Dispose();
                target = null; obtain_range = null;
            }
        }

        internal class Binary_Decoder_Worker : Coder_Worker
        {
            EncodedData target;
            uint code;
            int i = 0;

            internal Binary_Decoder_Worker(EncodedData input)
            {
                target = input;
            }

            internal BitArray run()
            {
                output = new List<int>(target.o_data_length);
                while (i < precision)
                {
                    next_input();
                }
                while (output.Count*sizeof(int)*8 < target.o_data_length)
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
                long dif = (long)high - (long)low + 1;
                uint temp_bound = (uint)(low + dif / bit_count[1] * bit_count[0] + (dif % bit_count[1]) * bit_count[0]);
                if (code < temp_bound)
                {
                    emit_output();
                    high = temp_bound - 1;

                    //add low count
                }
                else
                {
                    current_int |= half;
                    emit_output();
                    low = temp_bound;

                    //add high count
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
                    underflow++;
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
            using (Binary_Encoder_Worker worker = new Binary_Encoder_Worker(input))
            {
                return worker.run();
            }
        }

        public BitArray decode(EncodedData input)
        {
            using (Binary_Decoder_Worker worker = new Binary_Decoder_Worker(input))
            {
                return worker.run();
            }
        }
    }
}
