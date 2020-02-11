using System;
using System.Xml.Serialization;

namespace SPOTL
{
    public class LoveNumbers
    {
        [XmlArray("h")]
        public double[] h = new double[3] { Constants.h1Default, Constants.h2Default, Constants.h3Default };
        [XmlArray("k")]
        public double[] k = new double[3] { Constants.k1Default, Constants.k2Default, Constants.k3Default };
        [XmlArray("l")]
        public double[] l = new double[3] { Constants.l1Default, Constants.l2Default, Constants.l3Default };
    }
}
