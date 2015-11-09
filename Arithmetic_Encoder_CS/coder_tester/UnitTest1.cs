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
        public class Packet_Sender : IDataPacketSender<byte>
        {
            public byte[] randomData;
            System.Collections.Generic.List<IDataPacketHandler<byte>> handlers = new System.Collections.Generic.List<IDataPacketHandler<byte>>();

            public Packet_Sender()
            {
                
                int Min = 0;
                int Max = 20;
                Random randNum = new Random();
                randomData = Enumerable
                    .Repeat<byte>(0, 100000)
                    .Select(i => (byte)randNum.Next(Min, Max))
                    .ToArray();
                    
            }
            public void RegisterPacketHandler(IDataPacketHandler<byte> handler)
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
            public bool SaveBytes(string _FileName)
            {
                try
                {
                    // Open file for reading
                    System.IO.FileStream _FileStream =
                       new System.IO.FileStream(_FileName, System.IO.FileMode.Create,
                                                System.IO.FileAccess.Write);
                    // Writes a block of bytes to this stream using data from
                    // a byte array.
                    _FileStream.Write(randomData, 0, randomData.Length);

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
            public bool ReadBytes(string _FileName)
            {
                try
                {
                    // Open file for reading
                    System.IO.FileStream _FileStream =
                       new System.IO.FileStream(_FileName, System.IO.FileMode.Open,
                                                System.IO.FileAccess.Read);
                    // Writes a block of bytes to this stream using data from
                    // a byte array.
                    _FileStream.Read(randomData, 0, randomData.Length);

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

        [TestMethod]
        public void ByteCoderTest()
        {
            Packet_Sender bitstream = new Packet_Sender();
            //bitstream.ReadBytes(@"bytes");
            ByteCoder encoder = new ByteEncoder();
            encoder.Initialise(bitstream);
            bitstream.SendDataPeriodically(1);
            encoder.Finalise();
            ByteCoder decoder = new ByteDecoder(bitstream.randomData.Count(), bitstream.randomData);
            Packet_Sender bitstream2 = new Packet_Sender();
            bitstream2.randomData = encoder.result.ToArray();
            decoder.Initialise(bitstream2);
            bitstream2.SendDataPeriodically(1);
            decoder.Finalise();

            bool violation = false;
            Assert.IsTrue(decoder.result.Count == bitstream.randomData.Count(), "length not equal. decoded length: {0} as opposed to {1}", decoder.result.Count, bitstream.randomData.Count());
            for (int i = 0; i < bitstream.randomData.Length; i++)
            {
                violation = bitstream.randomData[i] != decoder.result[i];
                if (violation)
                    break;
            }
            Assert.IsFalse(violation, "Not all bits are equal");
        }
    }
}
