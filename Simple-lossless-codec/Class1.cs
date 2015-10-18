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
            entropy = (int)input.Length;
        }
    }

    public class EncodedData
    {
        public int o_data_length { get; }
        public uint[] range { get; }
        public BitArray data { get; }

        public EncodedData(int o_data_length, uint[] range, BitArray output)
        {
            this.o_data_length = o_data_length; this.range = range; this.data = output;
        }
    }

    public class BinaryCoder
    {
        const uint precision = sizeof(uint); const uint half = uint.MaxValue >> 1 + 1;
        const uint quarter = uint.MaxValue >> 2 + 1; const uint quarter_3 = 0x3 << (sizeof(uint) - 2);

        internal class Output_Worker:IDisposable
        {
            Data target;
            Action[] obtain_range;
            List<uint> output;
            uint[] range, bit_count; //keep size of bit_count half of sizeof(uint) or overflow can occur
            uint current_int = 0, position_count = 0, underflow = 0;
            
            internal Output_Worker(Data input)
            {
                this.target = input;
                obtain_range = new Action[2] { low_range, high_range };
            }

            internal EncodedData run()
            {
                range = new uint[2] { 0, uint.MaxValue };
                bit_count = new uint[2] { 1, 2 }; //initial condition
                output = new List<uint>(target.entropy);
                foreach (var bit in target.data)
                {
                    obtain_range[(uint)bit]();
                    check_output();
                }

                return finalise();
            }
            void check_output()
            {
                //since first bit of low always 0, and that of high always 1 when not equal
                if (range[1] < half)//equal when high bit == 0 or low bit == 1
                {
                    emit_output();
                    emit_ones();
                    range[0] <<= 1; range[1] = (range[1] << 1) + 1;
                    check_output();
                }
                else if (range[0] >= half)
                {
                    current_int |= half;
                    emit_output();
                    emit_zeros();
                    range[0] <<= 1; range[1] = (range[1] << 1) + 1;
                    check_output();
                }
                else
                    check_underflow();
            }
            void check_underflow()
            {
                if (range[0] >= quarter && range[1] < quarter_3)
                {
                    //second bit equals, discard it
                    underflow++;
                    range[0] = range[0] - quarter;
                    range[1] = range[1] - quarter;
                    range[0] <<= 1; range[1] = (range[1] << 1) + 1;
                    check_underflow();
                }
            }
            void emit_output()
            {
                //check overflow
                if (++position_count == sizeof(int))
                {
                    output.Add(current_int);
                    current_int = 0;
                    position_count = 0;
                }
                else
                    current_int >>= 1;
            }
            void emit_zeros()
            {
                while (underflow > 0)
                {
                    emit_output();
                    underflow--;
                }
            }
            void emit_ones()
            {
                while (underflow > 0)
                {
                    current_int |= half;
                    emit_output();
                    underflow--;
                }
            }
            EncodedData finalise()
            {
                underflow++;
                if (range[0] >= quarter)
                {
                    current_int |= half;//flush bits to know when to stop decoder from looping
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
                    current_int >>= (int)(sizeof(int) - position_count - 1);
                    output.Add(current_int);
                }

                BitArray output_bits = new BitArray((int[])(object)output.ToArray());
                EncodedData result = new EncodedData(target.data.Length, (uint[])range.Clone(), output_bits);
                return result;
            }
            void low_range()
            {
                uint dif = range[1] - range[0] + 1;
                //divide first to avoid overflow
                range[1] = range[0] + dif / bit_count[1] * bit_count[0] + (dif % bit_count[1]) * bit_count[0] - 1;
                //feedback error for multiplication
                bit_count[0]++; bit_count[1]++;//update probability for next step (FIR)
                                               //to be replaced by prediction algorithm (partial matching)
            }
            void high_range()
            {
                uint dif = range[1] - range[0] + 1;
                range[1] = range[0] + dif - 1;
                range[0] = range[0] + dif / bit_count[1] * bit_count[0] + (dif % bit_count[1]) * bit_count[0];
                bit_count[1]++;
            }
            public void Dispose(){
                target = null; obtain_range = null; output = null; range = null; bit_count = null;
            }
        }
        
        EncodedData encode(Data input)
        {
            using (Output_Worker worker = new Output_Worker(input))
            {
                return worker.run();
            }
        }

        BitArray decode(EncodedData input)
        {
                List<bool> output = new List<bool>(input.o_data_length);
            uint code;
            //load precision number of bits to code
                BitArray return_var = new BitArray(output.ToArray());
                return return_var;
        }
    }
}
