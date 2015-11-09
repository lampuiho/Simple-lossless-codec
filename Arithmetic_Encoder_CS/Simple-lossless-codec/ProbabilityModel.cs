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
            int temp = System.Convert.ToInt32(symbol);
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
        protected const uint max_value = 0x1000000;

        public ProbabilityAdaptor(int MaxSymbolValue) : base(MaxSymbolValue)
        {
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
    public class BufferedBIPA : BinaryIncrementProbilityAdpator
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
    public class BinaryExpIncProbilityAdaptor : BinaryIncrementProbilityAdpator
    {
        const int power_factor = 4;
        //const int scale_factor = 1 << power_factor;

        public override void Add(dynamic symbol)
        {
            if ((bool)symbol == true)
                symbol_count[0] -= symbol_count[0] >> power_factor;
            else
                symbol_count[0] += (max_value - symbol_count[0]) >> power_factor;
        }
    }
    public class ExpoIncProbabilityAdaptor : ProbabilityAdaptor
    {
        protected const int power_factor = 4;

        public ExpoIncProbabilityAdaptor(int MaxSymbolValue) : base(MaxSymbolValue)
        {
        }

        public override void Add(dynamic symbol)
        {
            for(int i = 0; i < symbol_count.Count(); ++i)
                symbol_count[i] -= (symbol_count[i] >> power_factor);
            symbol_count[symbol] += (max_value >> power_factor);
        }
    }
    public class MixedDistExpIncProbabilityAdaptor : ExpoIncProbabilityAdaptor
    {
        uint[][] distributions;

        public MixedDistExpIncProbabilityAdaptor(int MaxSymbolValue) : base(MaxSymbolValue)
        {
            init_distribution();
        }

        public override void Add(dynamic symbol)
        {
            var distribution = distributions[symbol];
            for (int i = 0; i < symbol_count.Count(); ++i)
                symbol_count[i] += (distribution[i] - symbol_count[i]) >> power_factor;
        }

        public void init_distribution() { }
    }
}
