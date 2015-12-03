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
        public class Packet_Sender : IDataPacketSender<int>
        {
            public int[] randomData;
            System.Collections.Generic.List<IDataPacketHandler<int>> handlers = new System.Collections.Generic.List<IDataPacketHandler<int>>();

            public Packet_Sender()
            {
            }
            public Packet_Sender(int[] data)
            {
                randomData = data;
            }
            public Packet_Sender(byte[] data)
            {
                randomData = new int[data.Length / 2];

                for (int i = 0; i < data.Length; i += 2)
                    randomData[i >> 1] = BitConverter.ToInt16(data, i);
            }
            public void GenerateRandom(int amount)
            {
                int Min = -20;
                int Max = 20;
                Random randNum = new Random();
                randomData = Enumerable
                    .Repeat<int>(0, amount)
                    .Select(i => randNum.Next(Min, Max))
                    .ToArray();
            }
            public void RegisterPacketHandler(IDataPacketHandler<int> handler)
            {
                handlers.Add(handler);
            }
            public void SendDataPeriodically(int delay_millisecond)
            {
                foreach(var data in randomData)
                {
                    //System.Threading.Thread.Sleep(1);
                    foreach(var handler in handlers)
                    {
                        handler.Input(data);
                    }
                }
            }
            public static bool SaveBytes(string _FileName, byte[] array)
            {
                try
                {
                    // Open file for reading
                    System.IO.FileStream _FileStream =
                       new System.IO.FileStream(_FileName, System.IO.FileMode.Create,
                                                System.IO.FileAccess.Write);
                    // Writes a block of bytes to this stream using data from
                    // a byte array.
                    _FileStream.Write(array, 0, array.Length);

                    // close file stream
                    _FileStream.Close();

                    return true;
                }
                catch (Exception _Exception)
                {
                    // Error
                    Console.WriteLine("Exception caught in process: {0}",
                                      _Exception.ToString());
                }

                // error occured, return false
                return false;
            }
            public static bool ReadBytes(string _FileName, ref byte[] array)
            {
                try
                {
                    // Open file for reading
                    System.IO.FileStream _FileStream =
                       new System.IO.FileStream(_FileName, System.IO.FileMode.Open,
                                                System.IO.FileAccess.Read);
                    System.IO.FileInfo info = new System.IO.FileInfo(_FileName);

                    // Writes a block of bytes to this stream using data from
                    // a byte array.
                    array = new byte[info.Length];
                    _FileStream.Read(array, 0, (int)info.Length);

                    // close file stream
                    _FileStream.Close();

                    return true;
                }
                catch (Exception _Exception)
                {
                    // Error
                    Console.WriteLine("Exception caught in process: {0}",
                                      _Exception.ToString());
                }

                // error occured, return false
                return false;
            }
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
            System.Diagnostics.Trace.Listeners.Add(new System.Diagnostics.ConsoleTraceListener());
            System.Diagnostics.Trace.WriteLine(string.Format("DataSize: {0}; Entropy: {1}; Encoded Size: {2}", test_data.data.Length, test_data.entropy, encoded.data.Length));
            System.Diagnostics.Trace.Flush();
        }
        /*
        [TestMethod]
        public void TestMethod3()
        {
            Packet_Sender bitstream = new Packet_Sender();
            StreamBinaryEncoder<byte> coder = new StreamBinaryEncoder<byte>();
            coder.Initialise(bitstream);
            bitstream.SendDataPeriodically(1);
            coder.Finalise();
            var encoded = coder.result;
            var decoded = coder.decode(encoded);
            bool violation = false;
            var bytes = bitstream.randomData.SelectMany<byte, byte>(x=>BitConverter.GetBytes(x)).ToArray();
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
        */
        [TestMethod]
        public void ByteCoderTest()
        {
            Packet_Sender bitstream = new Packet_Sender();
            bitstream.GenerateRandom(10000);
            //bitstream.ReadBytes(@"bytes");
            ByteEncoder encoder = new ByteEncoder();
            encoder.Initialise(bitstream);
            bitstream.SendDataPeriodically(1);
            encoder.Finalise();
            ByteDecoder decoder = new ByteDecoder(bitstream.randomData.Count());
            Packet_Sender bitstream2 = new Packet_Sender(encoder.result.Select(x=>(int)x).ToArray());
            decoder.Initialise(bitstream2);
            bitstream2.SendDataPeriodically(1);
            decoder.Finalise();

            bool violation = false;
            Assert.IsTrue(decoder.result.Count == bitstream.randomData.Count(), "length not equal. decoded length: {0} as opposed to {1}", decoder.result.Count, bitstream.randomData.Count());
            for (int i = 0; i < bitstream.randomData.Length; i++)
            {
                violation = bitstream.randomData[i] != decoder.result[i];
                if (violation)
                {
                    //bitstream.SaveBytes(@"bytes");
                    break;
                }
            }
            Assert.IsFalse(violation, "Not all bits are equal");
        }

        [TestMethod]
        public void BinaryCoderTest()
        {
            byte[] original_data = null;
            Packet_Sender.ReadBytes(@"bytes", ref original_data);
            int[] table = new int[byte.MaxValue + 1];
            foreach (byte i in original_data)
                table[i]++;
            BinarySymbolBinaryRangeCoder.MoreZeroAdaptor s_adap = new BinarySymbolBinaryRangeCoder.MoreZeroAdaptor(table);
            var sent_array = original_data.Select(x => (byte) s_adap.get_symbol(x)).ToArray();
            Packet_Sender bitstream = new Packet_Sender(sent_array);
            //bitstream.GenerateRandom(10000);
            uint low_count = BinaryProbabilityModel.get_low_count(sent_array);
            uint total_count = (uint)(sent_array.Length) * 8;
            sent_array = null;
            BinaryProbabilityModel p_model = new BinaryProbabilityModel(low_count, total_count);
            //BinaryProbabilityModel p_model = new BinaryExpIncProbilityAdaptor(low_count, total_count);
            BinaryEncoder encoder = new BinaryEncoder(p_model);
            encoder.Initialise(bitstream);
            bitstream.SendDataPeriodically(1);
            encoder.Finalise();
            uint entropy = BinaryProbabilityModel.get_entropy(low_count, total_count);
            System.Diagnostics.Trace.Listeners.Add(new System.Diagnostics.ConsoleTraceListener());
            System.Diagnostics.Trace.WriteLine(string.Format("DataSize: {0}; Entropy: {1}; Encoded Size: {2}", total_count, entropy, encoder.result.Count * 8));
            System.Diagnostics.Trace.Flush();
            //p_model = new BinaryExpIncProbilityAdaptor(low_count, total_count);
            BinaryDecoder decoder = new BinaryDecoder(p_model, bitstream.randomData.Length * 2);
            Packet_Sender bitstream2 = new Packet_Sender(encoder.result.Select(x => (int)x).ToArray());
            decoder.Initialise(bitstream2);
            bitstream2.SendDataPeriodically(1);
            decoder.Finalise();
            encoder = null;
            var decoded = decoder.result.Select(x => (byte) s_adap.get_output(x)).ToArray();
            Assert.IsTrue(decoded.Length == original_data.Length, "length not equal. decoded length: {0} as opposed to {1}", decoded.Length, original_data.Length);

            bool violation = false;
            for (int i = 0; i < original_data.Length; i++)
            {
                violation = original_data[i] != decoded[i];
                if (violation)
                    break;
            }
            Assert.IsFalse(violation, "Not all bits are equal");
        }

        [TestMethod]
        public void Get16bitEntropy()
        {
            byte[] original_data = null;
            Packet_Sender.ReadBytes(@"bytes", ref original_data);
            int total_count = 0;
            int[] table = new int[ushort.MaxValue + 1];
            foreach(int i in original_data)
            {
                table[(ushort)i]++;
                total_count++;
            }
            for (int i=0; i<= ushort.MaxValue;++i)
                if (table[i] == 0)
                {
                    table[i]++;
                    total_count++;
                }

            double entropy = 0;
            for (int i = 0; i <= ushort.MaxValue; ++i)
            {
                entropy+=table[i] * Math.Log((double)total_count / (double)table[i], 2);
            }

            uint result = (uint)Math.Ceiling(entropy);
        }

        [TestMethod]
        public void Get8bitEntropy()
        {
            byte[] original_data = null;
            Packet_Sender.ReadBytes(@"bytes", ref original_data);

            int total_count = 0;
            int[] table = new int[byte.MaxValue+1];
            foreach (byte i in original_data)
            {
                table[i]++;
                total_count++;
            }
            for (int i = 0; i <= byte.MaxValue; ++i)
                if (table[i] == 0)
                {
                    table[i]++;
                    total_count++;
                }

            double entropy = 0;
            for (int i = 0; i <= byte.MaxValue; ++i)
            {
                entropy += table[i] * Math.Log((double)total_count / (double)table[i], 2);
            }

            uint result = (uint)Math.Ceiling(entropy);
        }
    }
}
