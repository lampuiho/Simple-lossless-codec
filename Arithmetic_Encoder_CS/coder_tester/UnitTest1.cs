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
            public int[] data;
            System.Collections.Generic.List<IDataPacketHandler<int>> handlers = new System.Collections.Generic.List<IDataPacketHandler<int>>();

            public Packet_Sender()
            {
            }
            public Packet_Sender(int[] data)
            {
                this.data = data;
            }
            public Packet_Sender(byte[] data)
            {
                this.data = new int[data.Length];

                for (int i = 0; i < data.Length; i += 1)
                    this.data[i] = data[i];
            }
            public static int[] GenerateRandom(int amount)
            {
                int Min = -1;
                int Max = 1;
                Random randNum = new Random();
                var data = Enumerable
                    .Repeat<int>(0, amount)
                    .Select(i => randNum.Next(Min, Max))
                    .ToArray();
                /*
                var data = Enumerable
                    .Repeat<int>(0, amount).ToArray();
                data[data.Length - 1] = 1;*/
                return data;
            }
            public void RegisterPacketHandler(IDataPacketHandler<int> handler)
            {
                handlers.Add(handler);
            }
            public void SendDataPeriodically(int delay_millisecond)
            {
                foreach(var data in data)
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
        [TestMethod]
        public void TestMethod2()
        {
            BitArray test_array = new BitArray(Packet_Sender.GenerateRandom(1000));
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
            byte[] original_data = null;
            Packet_Sender.ReadBytes(@"bytes", ref original_data);
            //bitstream.GenerateRandom(10000);
            uint[] table = new uint[byte.MaxValue + 1];
            foreach (byte i in original_data)
                table[i]++;
            ProbabilityAdaptor p_model = new ProbabilityAdaptor(byte.MaxValue, table);
            ByteEncoder encoder = new ByteEncoder(p_model);
            Packet_Sender bitstream = new Packet_Sender(original_data);
            encoder.Initialise(bitstream);
            bitstream.SendDataPeriodically(1);
            encoder.Finalise();
            ByteDecoder decoder = new ByteDecoder(bitstream.data.Count(), p_model);
            Packet_Sender bitstream2 = new Packet_Sender(encoder.result.ToArray());
            decoder.Initialise(bitstream2);
            bitstream2.SendDataPeriodically(1);
            decoder.Finalise();

            bool violation = false;
            Assert.IsTrue(decoder.result.Count == bitstream.data.Count(), "length not equal. decoded length: {0} as opposed to {1}", decoder.result.Count, bitstream.data.Count());
            for (int i = 0; i < bitstream.data.Length; i++)
            {
                violation = bitstream.data[i] != decoder.result[i];
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
            //Packet_Sender bitstream = new Packet_Sender(sent_array);
            byte[] original_data = Packet_Sender.GenerateRandom(10000).Select(x=>(byte)x).ToArray();
            int[] table = new int[byte.MaxValue + 1];
            foreach (byte i in original_data)
                table[i]++;
            BinarySymbolBinaryRangeCoder.MoreZeroAdaptor s_adap = new BinarySymbolBinaryRangeCoder.MoreZeroAdaptor(table);
            var sent_array = original_data.Select(x => (byte)s_adap.get_symbol(x)).ToArray();
            //Packet_Sender.ReadBytes(@"bytes", ref original_data);
            uint low_count = BinaryProbabilityModel.get_low_count(sent_array);
            uint total_count = (uint)(sent_array.Length) * 8;
            //sent_array = null;
            BinaryProbabilityModel p_model = new BinaryProbabilityModel(low_count, total_count);
            //BinaryProbabilityModel p_model = new BinaryExpIncProbilityAdaptor(low_count, total_count);
            BinaryEncoder encoder = new BinaryEncoder(p_model);
            Packet_Sender bitstream = new Packet_Sender(sent_array);
            encoder.Initialise(bitstream);
            bitstream.SendDataPeriodically(1);
            encoder.Finalise();
            uint entropy = BinaryProbabilityModel.get_entropy(low_count, total_count);
            System.Diagnostics.Trace.Listeners.Add(new System.Diagnostics.ConsoleTraceListener());
            System.Diagnostics.Trace.WriteLine(string.Format("DataSize: {0}; Entropy: {1}; Encoded Size: {2}", total_count, entropy, encoder.result.Count * 8));
            System.Diagnostics.Trace.Flush();
            //p_model = new BinaryExpIncProbilityAdaptor(low_count, total_count);
            BinaryDecoder decoder = new BinaryDecoder(p_model, bitstream.data.Length);
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
            ushort[] converted = new ushort[original_data.Length / 2];
            for(int i =0;i < original_data.Length; i += 2)
            {
                converted[i / 2] = BitConverter.ToUInt16(original_data, i);
            }
            int[] table = new int[ushort.MaxValue + 1];
            foreach(int i in converted)
            {
                table[i]++;
                total_count++;
            }
            double entropy = 0;
            for (int i = 0; i <= ushort.MaxValue; ++i)
                if (table[i] > 0)
                    entropy+=table[i] * Math.Log((double)total_count / (double)table[i], 2);

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

            double entropy = 0;
            for (int i = 0; i <= byte.MaxValue; ++i)
                if (table[i] > 0)
                    entropy += table[i] * Math.Log((double)total_count / (double)table[i], 2);

            uint result = (uint)Math.Ceiling(entropy);
        }
    }
}
