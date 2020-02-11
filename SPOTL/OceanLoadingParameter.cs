using System;
using System.Collections;
using System.Xml.Serialization;
using System.IO;
using System.Numerics;

namespace SPOTL
{
    public class OceanLoadingParameter
    {
        [XmlElement("Location")]
        public Location Location = null;

        [XmlElement("GrennsFunction")]
        public OceanLoadingProperties.GreenFunctionEntry GreensFunction = null;

        [XmlElement("OceanModel")]
        public SPOTL.OceanLoadingProperties.OceanModelEntry[] OceanModel = null;
      
        [XmlElement("LocalOceanModel")]
        public SPOTL.OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanModel = null;

        public OceanLoadingParameter()
        {
            list = new ArrayList();
        }

        [XmlIgnore()]
        internal ArrayList list = null;

        [XmlArray("OceanLoadWaveItems")]
        public OceanLoadWave[] waveItems
        {
            get
            {
                OceanLoadWave[] items = new OceanLoadWave[list.Count];
                list.CopyTo(items);
                return items;
            }
            set
            {
                if (value == null) return;
                OceanLoadWave[] items = (OceanLoadWave[])value;
                list.Clear();
                foreach (OceanLoadWave item in items)
                    list.Add(item);
            }
        }

        public int AddEntry(OceanLoadWave item)
        {
            return list.Add(item);
        }

        public void RemoveEntryAt(int i)
        {
            list.RemoveAt(i);
        }

        public static void serialisieren(OceanLoadingParameter nout, string path)
        {
            XmlSerializer s = new XmlSerializer(typeof(OceanLoadingParameter));
            TextWriter w = new StreamWriter(path);
            s.Serialize(w, nout);
            w.Close();
        }

        public static OceanLoadingParameter deserialisieren(string path)
        {
            XmlSerializer s = new XmlSerializer(typeof(OceanLoadingParameter));
            TextReader r = new StreamReader(path);
            OceanLoadingParameter nout = (OceanLoadingParameter)s.Deserialize(r);
            r.Close();
            return nout;
        }
    }

    public class OceanLoadWave
    {
        [XmlAttribute("DarwinSymbol")]
        public string dsym = null;

        [XmlArray("CartwrightCode")]
        public int[] icte;

        [XmlElement("OceanLoad")]
        public OceanLoad oceanLoad = null;

        public OceanLoadWave()
        {
            list = new ArrayList();
        }

        [XmlIgnore()]
        internal ArrayList list = null;

        [XmlArray("OceanLoadItems")]
        public OceanLoad[] nloadfItems
        {
            get
            {
                OceanLoad[] items = new OceanLoad[list.Count];
                list.CopyTo(items);
                return items;
            }
            set
            {
                if (value == null) return;
                OceanLoad[] items = (OceanLoad[])value;
                list.Clear();
                foreach (OceanLoad item in items)
                    list.Add(item);
            }
        }

        public int AddEntry(OceanLoad item)
        {
            return list.Add(item);
        }

        public void RemoveEntryAt(int i)
        {
            list.RemoveAt(i);
        }
    }

    public class OceanLoad
    {
        [XmlAttribute("OceanModel")]
        public string OceanModel = null;

        [XmlAttribute("LocalModel")]
        public string LocalModel = null;

        [XmlElement("header")]
        public string[] header = null;

        [XmlElement("g")]
        public Load g = new Load();

        [XmlElement("gdef")]
        public Load gdef = new Load();

        [XmlElement("p")]
        public Load p = new Load();

        [XmlArray("d")]
        public Load[] d = new Load[3] { new Load(), new Load(), new Load() };

        [XmlArray("t")]
        public Load[] t = new Load[2] { new Load(), new Load() };

        [XmlArray("tdef")]
        public Load[] tdef = new Load[2] { new Load(), new Load() };

        [XmlArray("s")]
        public Load[] s = new Load[3] { new Load(), new Load(), new Load() };

        [XmlElement("o")]
        public Load o = new Load();

        public OceanLoad()
        {

        }
    }

    public class Load
    {
        [XmlAttribute("Amplitude")]
        public double Amplitude = double.NaN;

        [XmlAttribute("Phase")]
        public double Phase = double.NaN;

        public Load()
        {
            this.Amplitude = 0.0d;
            this.Phase = 0.0d;
        }

        public Load(double Amplitude, double Phase)
        {
            this.Amplitude = Amplitude;
            this.Phase = Phase;
        }

        public Load(Complex c)
        {
            this.Amplitude = c.Real;
            this.Phase = c.Imaginary;
        }
    }
}
