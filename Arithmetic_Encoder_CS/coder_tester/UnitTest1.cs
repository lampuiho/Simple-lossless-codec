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

        BitArray randomBitGenerator()
        {
            int Min = 0;
            int Max = 20;
            Random randNum = new Random();
            int[] test2 = Enumerable
                .Repeat(0, 10)
                .Select(i => randNum.Next(Min, Max))
                .ToArray();

            BitArray test_array = new BitArray(test2);

            return test_array;
        }
    }
}
