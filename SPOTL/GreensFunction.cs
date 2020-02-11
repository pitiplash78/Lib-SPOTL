using System;
using System.Collections;
using System.IO;

namespace SPOTL
{
    internal class GreensFunction
    {
        public static GreensFunction readGreensfunction(string path)
        {
            FileStream sourceFileStream = File.OpenRead(path);
            System.IO.Compression.GZipStream decompressingStream = new System.IO.Compression.GZipStream(sourceFileStream, System.IO.Compression.CompressionMode.Decompress);

            GreensFunction tmpgf = new GreensFunction();

            using (var re = new StreamReader(decompressingStream))
            {
                tmpgf.grname = re.ReadLine().Trim(); //Header

                string str = null;
                while ((str = re.ReadLine()) != null)
                {
                    // BlockHeader
                    int ngreen = Convert.ToInt32(str.Substring(0, 1));
                    int num = Convert.ToInt32(str.Substring(1, 3));
                    int ntot = Convert.ToInt32(str.Substring(4, 4));
                    int ngr = Convert.ToInt32(str.Substring(8, 4));
                    double beg = Convert.ToDouble(str.Substring(12, 10), Constants.NumberFormatEN);
                    double end = Convert.ToDouble(str.Substring(22, 10), Constants.NumberFormatEN);
                    double spc = Convert.ToDouble(str.Substring(32, 10), Constants.NumberFormatEN);
                    GreensFunction.FinGrd fingrd = (GreensFunction.FinGrd)Enum.Parse(typeof(GreensFunction.FinGrd), str.Substring(42, 6).Trim()[0].ToString());

                    double[,] grfn = new double[ngr, ngreen];

                    for (int i = 0; i < ngr; i++)
                    {
                        str = re.ReadLine();
                        for (int j = 0; j < ngreen; j++)
                        {
                            double v = double.Parse(str.Substring(j * 13, 13).Trim(), Constants.NumberFormatEN);
                            grfn[i, j] = v;
                        }
                    }
                    tmpgf.AddElementEntry(new GreensFunction.GFelememt(ngreen, num, ntot, ngr, beg, end, spc, fingrd, grfn));
                }
            }
            return tmpgf;
        }

        public GreensFunction()
        {
            gf = new ArrayList();
        }

        /// <summary>
        /// Name of the Greens function.
        /// </summary>
        public string grname = null;

        internal ArrayList gf;

        public GFelememt[] GrennsFunctionItems
        {
            get
            {
                GFelememt[] items = new GFelememt[gf.Count];
                gf.CopyTo(items);
                return items;
            }
            set
            {
                if (value == null) return;
                GFelememt[] items = (GFelememt[])value;
                gf.Clear();
                foreach (GFelememt item in items)
                    gf.Add(item);
            }
        }

        public int AddElementEntry(GFelememt item)
        {
            return gf.Add(item);
        }

        public enum FinGrd
        {
            None,
            C,
            G,
            F,
            L,
        }

        public class GFelememt
        {
            public int ngreen = 0;
            public int num = 0;
            public int ntot = 0;
            public int ngr = 0;
            public double beg = double.NaN;
            public double end = double.NaN;
            public double spc = double.NaN;
            public FinGrd fingrd = FinGrd.None;
            public double[,] grfn = null;

            public GFelememt(int ngreen, int num, int ntot, int ngr,
                             double beg, double end, double spc, FinGrd fingrd, double[,] grfn)
            {
                this.ngreen = ngreen;
                this.num = num;
                this.ntot = ntot;
                this.ngr = ngr;
                this.beg = beg;
                this.end = end;
                this.spc = spc;
                this.fingrd = fingrd;
                this.grfn = grfn;
            }
        }
    }
}
