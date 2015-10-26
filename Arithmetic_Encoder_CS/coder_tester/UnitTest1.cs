using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Simple_lossless_codec;
using System.Linq;
using System.Collections;

namespace coder_tester
{
    [TestClass]
    public class UnitTest1
    {
        public class Packet_Sender : IDataPacketSender<short>
        {
            public short[] randomData;
            System.Collections.Generic.List<IDataPacketHandler<short>> handlers = new System.Collections.Generic.List<IDataPacketHandler<short>>();

            public Packet_Sender()
            {
                int Min = 0;
                int Max = 20;
                Random randNum = new Random();
                randomData = Enumerable
                    .Repeat<short>(0, 157)
                    .Select(i => (short)randNum.Next(Min, Max))
                    .ToArray();
            }
            public void RegisterPacketHandler(IDataPacketHandler<short> handler)
            {
                handlers.Add(handler);
            }
            public void SendDataPeriodically(int delay_millisecond)
            {
                foreach(var data in randomData)
                {
                    System.Threading.Thread.Sleep(1);
                    foreach(var handler in handlers)
                    {
                        handler.Input(data);
                    }
                }
            }
        }

        [TestMethod]
        public void TestMethod1()
        {
            BitArray test_array = randomBitGenerator();
            Data test_data = new Data(test_array);
            bool test_result = test_array.Length > test_data.entropy;
            Assert.IsTrue(test_result, "Entropy: ({0}) is not smaller than {1}", test_data.entropy, test_array.Length);
        }


        [TestMethod]
        public void TestMethod2()
        {
            BitArray test_array = randomBitGenerator();
            Data test_data = new Data(test_array);
            BinaryCoder coder = new BinaryCoder();
            var encoded = coder.encode(test_data);
            var decoded = coder.decode(encoded);
            bool violation = false;
            Assert.IsTrue(test_array.Length == decoded.Length, "length not equal. decoded length: {0} as opposed to {1}", decoded.Length, test_array.Length);
            for(int i = 0; i < test_array.Length; i++)
            {
                violation = test_array[i] != decoded[i];
                if (violation)
                    break;
            }
            Assert.IsFalse(violation, "Not all bits are equal");
        }

        [TestMethod]
        public void TestMethod3()
        {
            Packet_Sender bitstream = new Packet_Sender();
            StreamBinaryEncoder<short> coder = new StreamBinaryEncoder<short>();
            coder.Initialise(bitstream);
            bitstream.SendDataPeriodically(1);
            coder.Finalise();
            var encoded = coder.result;
            var decoded = coder.decode(encoded);
            bool violation = false;
            var bytes = bitstream.randomData.SelectMany<short,byte>(x=>BitConverter.GetBytes(x)).ToArray();
            BitArray test_array = new BitArray(bytes);
            Assert.IsTrue(test_array.Length == decoded.Length, "length not equal. decoded length: {0} as opposed to {1}", decoded.Length, test_array.Length);
            for (int i = 0; i < bitstream.randomData.Length; i++)
            {
                violation = test_array[i] != decoded[i];
                if (violation)
                    break;
            }
            Assert.IsFalse(violation, "Not all bits are equal");
        }

        BitArray randomBitGenerator()
        {
            int Min = 0;
            int Max = 20;
            Random randNum = new Random();
            int[] test2 = Enumerable
                .Repeat(0, 10000)
                .Select(i => randNum.Next(Min, Max))
                .ToArray();

            BitArray test_array = new BitArray(test2);

            return test_array;
        }
    }
}
