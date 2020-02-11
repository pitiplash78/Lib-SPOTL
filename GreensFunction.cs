using System;
using System.Collections;
using System.IO;

namespace SPOTL
{
    internal class GreensFunction
    {
        public static GreensFunction readGreensfunction(string model)
        {
            GreensFunction tmpgf = new GreensFunction();

            FileStream fs = new FileStream(model, FileMode.Open);
            StreamReader re = new StreamReader(fs);

            string grname = re.ReadLine().Trim(); //Header

            string str = null;

            while ((str = re.ReadLine()) != null)
            {
                // BlockHeader
                int ngreen = Convert.ToInt32(str.Substring(0, 1));
                int num = Convert.ToInt32(str.Substring(1, 3));
                int ntot = Convert.ToInt32(str.Substring(4, 4));
                int ngr = Convert.ToInt32(str.Substring(8, 4));
                double beg = Convert.ToDouble(str.Substring(12, 10), Utilities.Utility.NumberFormatEN);
                double end = Convert.ToDouble(str.Substring(22, 10), Utilities.Utility.NumberFormatEN);
                double spc = Convert.ToDouble(str.Substring(32, 10), Utilities.Utility.NumberFormatEN);
                string fingrd = str.Substring(42, 6).Trim();

                double[,] grfn = new double[ngr, ngreen];

                for (int i = 0; i < ngr; i++)
                {
                    str = re.ReadLine();
                    for (int j = 0; j < ngreen; j++)
                    {
                        double v = double.Parse(str.Substring(j * 13, 13).Trim(), Utilities.Utility.NumberFormatEN);
                        grfn[i, j] = v;
                    }
                }
                tmpgf.AddElementEntry(new GreensFunction.GFelememt(ngreen, num, ntot, ngr, beg, end, spc, fingrd, grfn));
            }
            re.Close();

            return tmpgf;
        }

        public GreensFunction()
        {
            gf = new ArrayList();
        }
        
        private ArrayList gf;

        public GFelememt[] GrennFunctionItems
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
        
        public class GFelememt
        {
            public int ngreen = 0;
            public int num = 0;
            public int ntot = 0;
            public int ngr = 0;
            public double beg = double.NaN;
            public double end = double.NaN;
            public double spc = double.NaN;
            public string fingrd = null;
            public double[,] grfn = null;

            public GFelememt(int _ngreen, int _num, int _ntot, int _ngr,
                             double _beg, double _end, double _spc, string _fingrd, double[,] _grfn)
            {
                this.ngreen = _ngreen;
                this.num = _num;
                this.ntot = _ntot;
                this.ngr = _ngr;
                this.beg = _beg;
                this.end = _end;
                this.spc = _spc;
                this.fingrd = _fingrd;
                this.grfn = _grfn;
            }
        }
    }
}
