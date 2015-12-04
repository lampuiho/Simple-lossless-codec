using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Simple_lossless_codec
{
    public abstract class ProbabilityModel
    {
        //dictionary size can only be byte or short due to limitation of C# array
        protected uint[] symbol_count; //cumulative
        protected int maxIndex;

        public uint CDF_T
        {
            get
            {
                return symbol_count[maxIndex];
            }
        }
        public int SymbolMax
        {
            get
            {
                return maxIndex;
            }
        }

        public ProbabilityModel(int MaxSymbolValue)
        {
            symbol_count = new uint[MaxSymbolValue + 1];
            maxIndex = MaxSymbolValue;
            Init();
        }
        protected abstract void Init();
        public virtual uint CDF(dynamic symbol)
        {
            uint temp = (uint)(symbol);
            if (temp == 0)
                return 0;
            else
                return symbol_count[--temp];
        }
        public virtual uint P(dynamic symbol)
        {
            return CDF(symbol + 1) - CDF(symbol);
        }
    }

    public class ProbabilityAdaptor : ProbabilityModel
    {
        public const int shift_value = 24;
        public const uint max_value = 1 << shift_value;
        public const uint bit_mask = max_value - 1;

        public ProbabilityAdaptor(int MaxSymbolValue) : base(MaxSymbolValue)
        {
        }
        public ProbabilityAdaptor(int MaxSymbolValue, uint[] count) : base(MaxSymbolValue)
        {
            symbol_count = count;
            for (int i = 1; i < symbol_count.Length; ++i)
            {
                symbol_count[i] += symbol_count[i - 1];
            }
            for (int i = 0; i < MaxSymbolValue; ++i)
            {
                symbol_count[i] = (uint)((ulong)max_value * symbol_count[i] / this.CDF_T);
            }
            symbol_count[MaxSymbolValue] = max_value;
        }

        protected override void Init()
        {
            //uint avg_chance = (uint)((max_value + 1) / 2 / symbol_count.Length);
            uint avg_chance = 1;
            symbol_count[0] = avg_chance;
            for (int i = 1; i < symbol_count.Length; ++i)
            {
                symbol_count[i] = avg_chance + symbol_count[i - 1];
            }
        }
        public virtual void Add(dynamic symbol) { }
    }
    public class IncrementProbabilityAdaptor : ProbabilityAdaptor
    {

        public IncrementProbabilityAdaptor(int MaxSymbolValue) : base(MaxSymbolValue)
        {
        }

        public override void Add(dynamic symbol)
        {
            //scale down a bit when approaching overflow
            if (CDF_T == max_value)
            {
                for (int i = 0; i < symbol_count.Length; ++i)
                    symbol_count[i] = (symbol_count[i] >> 1);// + symbol_count[i] % 2);
            }

            for (int i = System.Convert.ToInt32(symbol); i < symbol_count.Length; i++)
                symbol_count[i]++;
        }
    }
    public class BinaryProbabilityModel : ProbabilityAdaptor
    {
        public new const int shift_value = 32;
        public new const uint max_value = uint.MaxValue;


        public BinaryProbabilityModel() : base(0)
        {
            symbol_count[0] = uint.MaxValue >> 1;
        }
        public BinaryProbabilityModel(uint low_count, uint total_count) : base(0)
        {
            symbol_count[0] = (uint)((ulong)low_count * uint.MaxValue / total_count);
        }
        protected override void Init() {
        }
        public new uint CDF
        {
            get
            {
                return symbol_count[0];
            }
        }

        public static uint get_low_count(byte[] array)
        {
            uint low = 0;
            foreach (byte b in array)
            {
                var temp = b;
                for (int i = 0; i < 8; ++i)
                {
                    if ((temp & 1) == 0)
                        low++;
                    temp >>= 1;
                }
            }
            return low;
        }
        public static uint get_entropy(uint low, uint total)
        {
            double high = total - low;
            return (uint)Math.Ceiling(low * Math.Log((double)total / (double)low, 2) + high * Math.Log((double)total / high, 2));
        }
    }
    public class BinaryIncrementProbilityAdpator : IncrementProbabilityAdaptor
    {
        public BinaryIncrementProbilityAdpator() : base(1)
        {
        }
        protected override void Init()
        {
            symbol_count[0] = 3;
            symbol_count[1] = 4;
        }
        public override uint CDF(dynamic symbol)
        {
            if (symbol)
                return symbol_count[0];
            else
                return 0;
        }
    }
    public class BufferedBIPA : BinaryProbabilityModel
    {
        BitWise_Byte_Buffer bit_buffer = new BitWise_Byte_Buffer();

        public override void Add(dynamic symbol)
        {
            bool full = bit_buffer.push(symbol);
            if (full)
                for (int i = 0; i < 8; ++i)
                    base.Add(bit_buffer.read());
        }
    }
    public class BinaryExpIncProbilityAdaptor : BinaryProbabilityModel
    {
        const int power_factor = 4;
        //const int scale_factor = 1 << power_factor;

        public BinaryExpIncProbilityAdaptor() : base() { }
        public BinaryExpIncProbilityAdaptor(uint low_count, uint total_count):base(low_count, total_count) { }
        public override void Add(dynamic symbol)
        {
            if (symbol == 1)
                symbol_count[0] -= symbol_count[0] >> power_factor;
            else
                symbol_count[0] += (max_value - symbol_count[0] + 1) >> power_factor;
        }
    }
    public class ExpoIncProbabilityAdaptor : ProbabilityAdaptor //clamping not fixed
    {
        protected const int power_factor = 4;

        public ExpoIncProbabilityAdaptor(int MaxSymbolValue) : base(MaxSymbolValue)
        {
        }
        
        protected override void Init()
        {
            /*
            uint avg_chance = (uint)((max_value + 1) / symbol_count.Length);
            symbol_count[0] = avg_chance;
            for (int i = 1; i < symbol_count.Length; ++i)
            {
                symbol_count[i] = avg_chance + symbol_count[i - 1];
            }
            */
            maxIndex = 80;
            ProbilityDistributions.g_cdf.CopyTo(symbol_count, 0);
        }

        public override void Add(dynamic symbol)
        {
            for(int i = 0; i < maxIndex; ++i)
                symbol_count[i] -= (symbol_count[i] >> power_factor);

            for(int i = symbol; i < maxIndex; ++i)
                symbol_count[i] += (max_value >> power_factor);
        }
    }
    public class MixedDistExpIncProbabilityAdaptor : ExpoIncProbabilityAdaptor
    {
        uint[,] distributions;
        public MixedDistExpIncProbabilityAdaptor(int MaxSymbolValue) : base(MaxSymbolValue)
        {
        }

        protected override void Init()
        {
            base.Init();
            distributions = new uint[maxIndex + 1, maxIndex + 1];
            for (int i = 0; i <= maxIndex; ++i)
                for (int j = 0; j <= maxIndex; ++j)
                    distributions[i, j] = ProbilityDistributions.u_cdf[i, j];
        }

        public override void Add(dynamic symbol)
        {
           for (int i = 0; i < maxIndex; ++i)
                symbol_count[i] = (uint)(symbol_count[i] + (((int)distributions[symbol,i] - symbol_count[i]) >> power_factor));
        }
    }
}
