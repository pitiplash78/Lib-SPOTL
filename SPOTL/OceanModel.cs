using System;
using System.Collections;
using System.IO;
using System.IO.Compression;

namespace SPOTL
{
    /// <summary>
    /// Contains all information about a ocean model (reading routine inclusive).
    /// </summary>
    internal class OceanModel
    {
        /// <summary>
        /// Reads a compressed ocean model.
        /// </summary>
        /// <param name="modelPaths">File path of the ocean model to be read.</param>
        /// <returns>The OceanModel class containing all information about a ocean model.</returns>
        public static OceanModel readTidalModel(string modelPaths)
        {
            OceanModel oceanModel = new OceanModel();

            FileStream sourceFileStream = File.OpenRead(modelPaths);
            System.IO.Compression.GZipStream decompressingStream = new System.IO.Compression.GZipStream(sourceFileStream, System.IO.Compression.CompressionMode.Decompress);

            using (var re = new StreamReader(decompressingStream))
            {
                string dsym = re.ReadLine().Trim();
                string tmp = re.ReadLine();
                int[] icte = new int[6];
                for (int i = 0; i < icte.Length; i++)
                    icte[i] = Convert.ToInt32(tmp.Substring(i * 3, 3));
                tmp = re.ReadLine();
                double tlat = double.Parse(tmp.Substring(0, 8), Constants.NumberFormatEN) +
                              double.Parse(tmp.Substring(8, 8), Constants.NumberFormatEN) / 1000d;
                tmp = re.ReadLine();
                double blat = double.Parse(tmp.Substring(0, 8), Constants.NumberFormatEN) +
                              double.Parse(tmp.Substring(8, 8), Constants.NumberFormatEN) / 1000d;
                tmp = re.ReadLine();
                double elong = double.Parse(tmp.Substring(0, 8), Constants.NumberFormatEN) +
                               double.Parse(tmp.Substring(8, 8), Constants.NumberFormatEN) / 1000d;
                tmp = re.ReadLine();
                double wlong = double.Parse(tmp.Substring(0, 8), Constants.NumberFormatEN) +
                               double.Parse(tmp.Substring(8, 8), Constants.NumberFormatEN) / 1000d;
                tmp = re.ReadLine();
                int latc = int.Parse(tmp.Substring(0, 8));
                int longc = int.Parse(tmp.Substring(8, 8));
                string mdnam = re.ReadLine().Trim();
                int[] ir1 = new int[latc * longc];
                int j = 0;
                while ( j < ir1.Length)
                {
                    tmp = re.ReadLine();
                    string[] parts = tmp.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    for (int k = 0; k < parts.Length; k++)
                        ir1[j++] = int.Parse(parts[k]);
                }

                int[] im1 = new int[latc * longc];
                j = 0;
                while (j < im1.Length)
                {
                    tmp = re.ReadLine();
                    string[] parts = tmp.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    for (int k = 0; k < parts.Length; k++)
                        im1[j++] = int.Parse(parts[k]);
                }
                oceanModel = new OceanModel(dsym, icte, tlat, blat, elong, wlong, latc, longc, mdnam, ir1, im1);
            }
            return oceanModel;
        }

        /// <summary>
        /// Darwin symbol. 
        /// </summary>
        public string dsym = null;

        /// <summary>
        /// Cartwright code.
        /// </summary>
        public int[] icte = null;

        /// <summary>
        /// Latitude (N) of top edge of model grid.
        /// </summary>
        public double tlat = double.NaN;

        /// <summary>
        /// Latitude (N) of bottom edge of model grid.
        /// </summary>
        public double blat = double.NaN;

        /// <summary>
        /// Longitude (E) of "right" edge of model grid
        /// </summary>
        public double elong = double.NaN;
        
        /// <summary>
        /// Longitude (E) of "left" edge of model grid
        /// </summary>
        public double wlong = double.NaN;

        /// <summary>
        /// Number of cells EW.
        /// </summary>
        public int latc = 0;

        /// <summary>
        /// Number of cells NS.
        /// </summary>
        public int longc = 0;

        /// <summary>
        /// Modelname.
        /// </summary>
        public string mdnam = null;

        /// <summary>
        /// Array of the real part for each cell.
        /// </summary>
        public int[] ir1 = null;

        /// <summary>
        /// Array of the imaginary part for each cell.
        /// </summary>
        public int[] im1 = null;

        /// <summary>
        /// Constructor for a Ocean Model.
        /// </summary>
        public OceanModel()
        {
        }

        /// <summary>
        /// Constructor for a Ocean Model.
        /// </summary>
        /// <param name="dsym">Darwin symbol.</param>
        /// <param name="icte">Cartwright code.</param>
        /// <param name="tlat">Latitude (N) of top edge of model grid.</param>
        /// <param name="blat">Latitude (N) of bottom edge of model grid.</param>
        /// <param name="elong">Longitude (E) of "right" edge of model grid</param>
        /// <param name="wlong">Longitude (E) of "left" edge of model grid</param>
        /// <param name="latc">Number of cells EW.</param>
        /// <param name="longc">Number of cells NS.</param>
        /// <param name="mdnam">Modelname.</param>
        /// <param name="ir1">Array of the real part for each cell.</param>
        /// <param name="im1">Array of the imaginary part for each cell.</param>
        public OceanModel(string dsym, int[] icte, double tlat, double blat, double elong,
                double wlong, int latc, int longc, string mdnam, int[] ir1, int[] im1)
            {
                this.dsym = dsym;
                this.icte = icte;
                this.tlat = tlat;
                this.blat = blat;
                this.elong = elong;
                this.wlong = wlong;
                this.latc = latc;
                this.longc = longc;
                this.mdnam = mdnam;
                this.ir1 = ir1;
                this.im1 = im1;
            }

        /// <summary>
        /// Returns the cell index (as i and j also ind), and fractional position within a cell, for a given long and lat (in degrees)
        /// </summary>
        /// <param name="rlong">Latitude [degree].</param>
        /// <param name="rlato">Longitude [degree].</param>
        /// <param name="x">Fractional position within a cell x.</param>
        /// <param name="y">Fractional position within a cell y.</param>
        /// <param name="i">Cell index i.</param>
        /// <param name="j">Cell index j.</param>
        /// <returns>Content of the cell.</returns>
        public int celfnd(double rlong, double rlato, ref double x, ref double y, ref int i, ref int j)
        {
            /*     subroutine celfnd(rlong,rlato,x,y,i,j,ind)
             * 
             *     Revision 1.1  2011/09/23 22:39:01  agnew  -  Initial revision */

            /*  returns the cell index (as i and j also ind), and fractional position within a cell, for a given long
             *  and lat (in degrees) */

            int ind = 0;

            /*  check to see that we are within model */
            if (rlato > this.tlat || rlato < this.blat)                     //       if(rlato.gt.tlat.or.rlato.lt.blat) then
            {
                ind = 0;                                                    // 	 ind = 0
                return ind;                                                 // 	 return
            }                                                               //       endif

            double dlong = (rlong - this.wlong) % 360d;                     // amod(rlong-wlong,360.)

            if (dlong < 0.0)                                                //       if(dlong.lt.0) dlong = dlong + 360
                dlong += 360;

            if (dlong > this.elong - this.wlong)                            //       if(dlong.gt.elong-wlong) then
            {
                ind = 0;                                                    // 	 ind = 0
                return ind;                                                 // 	 return
            }                                                               //       endif

            /*  put rlong into the range of longitudes that are specified. */
            x = this.longc * dlong / (this.elong - this.wlong);             //       x = (longc*dlong)/(elong-wlong)
            y = this.latc * (this.tlat - rlato) / (this.tlat - this.blat);  //       y = (latc*(tlat-rlato))/(tlat-blat)
            i = (int)x + 1;                                                 //       i= int(x) + 1
            j = (int)y;                                                     //       j= int(y)
            x = x - i + 0.5d;                                               //       x = x - i + .5
            y = j - y + 0.5d;                                               //       y = j - y + .5
            ind = i + this.longc * j;                                       //       ind = i + longc*j
            return ind;
        }
    }
}
