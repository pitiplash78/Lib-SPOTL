using System;
using System.IO;
using System.Collections;

namespace SPOTL
{
    /// <summary>
    /// Class holds polygon file information for local ocean models.
    /// </summary>
    internal class PolInf
    {
        /// <summary>
        /// Constructor of PolInf class (holding all informations about the polygons).
        /// </summary>
        /// <param name="numberOfLocalModels">Number of the used local ocean models </param>
        public PolInf(int numberOfLocalModels)
        {
            this.npoly = new Data[numberOfLocalModels];
        }

        /// <summary>
        /// Array of the used polygons, containing all information, as the use mode (excluding/including)
        /// </summary>
        public Data[] npoly;

        /// <summary>
        /// Number of the polygons used in culculation, especialy when more than one is used.
        /// </summary>
        public int np = 0;

        /// <summary>
        /// Class containing all informations about one polygon.
        /// </summary>
        public class Data
        {
            /// <summary>
            /// Name of the polygon.
            /// </summary>
            public string polynm;

            /// <summary>
            /// Length of the coordinate array, describing the polygon path.
            /// </summary>
            public int npoly = 0;

            /// <summary>
            /// Array of the X coordinate, describing the polygon path.
            /// </summary>
            public double[] xpoly;

            /// <summary>
            /// Array of the Y coordinate, describing the polygon path.
            /// </summary>
            public double[] ypoly;

            /// <summary>
            /// Flag of using state exclude/ include.
            /// </summary>
            public bool use = true;
        }

        /// <summary>
        /// Internal class for loading the polygon coordinates.
        /// </summary>
        internal class Coordinates
        {
            public double x = double.NaN;
            public double y = double.NaN;

            public Coordinates(double x, double y)
            {
                this.x = x;
                this.y = y;
            }
        }

        /// <summary>
        /// Read one polygon file and returns data containing all informations about the polygon.
        /// </summary>
        /// <param name="path">File path of the polygon file to be read.</param>
        /// <returns>Returns the 'Data' containing the information about the readed polygon file.</returns>
        public static Data readPolgonInformation(string path)
        {
            Data data = new Data();

            FileStream fs = new FileStream(path, FileMode.Open);

            ArrayList polygonCoordinates = new ArrayList();

            using (StreamReader re = new StreamReader(fs))
            {
                data.polynm = re.ReadLine();
                string tmp = null;

                string[] parts = null;
                char[] delimeter = new char[] { ' ' };

                while ((tmp = re.ReadLine()) != null)
                {
                    parts = tmp.Split(delimeter, StringSplitOptions.RemoveEmptyEntries);
                    if(parts.Length == 2)
                        polygonCoordinates.Add(new Coordinates(double.Parse(parts[0], Constants.NumberFormatEN),
                                                               double.Parse(parts[1], Constants.NumberFormatEN)));
                }
            }

            data.npoly = polygonCoordinates.Count;
            data.xpoly = new double[data.npoly];
            data.ypoly = new double[data.npoly];

            for (int i = 0; i < data.npoly; i++)
            {
                Coordinates c = (Coordinates)polygonCoordinates[i];
                data.xpoly[i] = c.x;
                data.ypoly[i] = c.y;
            }

            return data;
        }

        /// <summary>
        /// Checks the existing polygons, according the the inclsion state, of the given point ist inside/outside the polygon.
        /// </summary>
        /// <param name="rlong">Longitude [degree] of the Coordinate to be checked.</param>
        /// <param name="rlat">Latitude [degree] of the Coordinate to be checked.</param>
        /// <returns>
        /// <para>False if</para>
        /// <para>1. there is any polygon (of all those read in) in which the coordinate should not fall, and it does or</para>
        /// <para>2. the coordinate should fall in at least one polygon (of those read in) and it does not</para>
        /// <para>otherwise it returns True.</para>
        /// </returns>
        public bool chkgon(double rlong, double rlat)// subroutine chkgon(rlong,rlat,iok)
        {
            /* returns iok=0 if
             * 1. there is any polygon (of all those read in) in which the coordinate should not fall, and it does or
             * 2. the coordinate should fall in at least one polygon (of those read in) and it does not
             * otherwise returns iok=1 */

            bool iok = false;
           
            /* first convert rlong to -180/180 (or leave in this, or 0/360) */
            if (rlong < -180.0)                                                     //       if(rlong.lt.-180) rlong = rlong + 360
                rlong += 360;

            /* now make rlong2 fall into the other system than rlong (unless rlong is between 0 and 180) */
            double rlong2 = rlong;
            if (rlong < 0.0)                                                        //       if(rlong.lt.0) rlong2 = rlong + 360
                rlong2 = rlong + 360;

            if (rlong > 180.0)                                                      //       if(rlong.gt.180) rlong2 = rlong - 360
                rlong2 = rlong - 360;

            for (int ip = 0; ip < this.np; ip++)                                   //       do 5 ip=1,np
            {
                if (!this.npoly[ip].use)                                        //       if(.not.use(ip)) then
                {
                    /* polygon is one we should not be in; test to see if we are, and if so set
                     * iok to 0 and return - check all such before seeing if we might also be in a polygon we should be */
                    int l = this.npoly[ip].npoly;
                    double[] xu = new double[l];
                    double[] yu = new double[l];
                    for (int i = 0; i < l; i++)                                    //          do 3 i=1,npoly(ip)
                    {
                        xu[i] = this.npoly[ip].xpoly[i];                //             xu(i) = xpoly(i,ip)
                        yu[i] = this.npoly[ip].ypoly[i];                //             yu(i) = ypoly(i,ip)
                    }

                    if (inpoly(ref rlong, ref rlat, ref xu, ref yu, ref l) != 0 ||
                        inpoly(ref rlong2, ref rlat, ref xu, ref yu, ref l) != 0)   // 	 if(inpoly(rlong,rlat,xu,yu,npoly(ip)).ne.0.or.inpoly(rlong2,rlat,xu,yu,npoly(ip)).ne.0)
                    {
                        iok = false;                                                // 	    iok=0
                        return iok;                                                 // 	    return
                    }
                }
            }

            int ianyok = 0;                                                         //       ianyok=0
            for (int ip = 0; ip < this.np; ip++)                                   //       do 15 ip=1,np
            {
                if (this.npoly[ip].use)                                         //       if(use(ip)) then
                {
                    /* polygon is one we should be in; test to see if we are, and if so set iok to 1 and return */
                    ianyok++;                                                       //          ianyok = ianyok+1
                    int l = this.npoly[ip].npoly;
                    double[] xu = new double[l];
                    double[] yu = new double[l];
                    for (int i = 0; i < l; i++)                                    //          do 13 i=1,npoly(ip)
                    {
                        xu[i] = this.npoly[ip].xpoly[i];                //             xu(i) = xpoly(i,ip)
                        yu[i] = this.npoly[ip].ypoly[i];                //             yu(i) = ypoly(i,ip)
                    }
                    if (inpoly(ref rlong, ref rlat, ref xu, ref yu, ref  l) != 0 ||
                        inpoly(ref rlong2, ref rlat, ref xu, ref yu, ref l) != 0)   // 	 if(inpoly(rlong,rlat,xu,yu,npoly(ip)).ne.0.or.inpoly(rlong2,rlat,xu,yu,npoly(ip)).ne.0)
                    {
                        iok = true;                                                 // 	    iok=1
                        return iok;                                                 // 	    return
                    }
                }
            }

            /* not inside any polygons; set iok to 0 if there are any we should have been in */
            iok = true;                                                             //       iok = 1
            if (ianyok > 0)                                                         //       if(ianyok.gt.0) iok = 0
                iok = false;

            return iok;                                                             //       return
        }

        /// <summary>
        /// Checkes vor a given polygon if the given location is inside,on the edge, or outside the polygon.
        /// </summary>
        /// <param name="x">Longitude [degree] of the Coordinate to be checked.</param>
        /// <param name="y">Latitude [degree] of the Coordinate to be checked.</param>
        /// <param name="xv">Array of the X coordinates,describing the polygon.</param>
        /// <param name="yv">Array of the Y coordinates,describing the polygon</param>
        /// <param name="nv">Length of the X and Y coordinte array.</param>
        /// <returns>
        /// <para>Returns 1 if point at (x,y) is inside polygon whose nv vertices are at xv(1),yv(1);...,xv(nv),yv(nv)</para>
        /// <para>Returns 0 if point is outside</para>
        /// <para>Returns 2 if point is on edge or vertex</para>
        /// </returns>
        private int inpoly(ref double x, ref double y, ref double[] xv, ref double[] yv, ref int nv) // function inpoly(x,y,xv,yv,nv)
        {
            /*     Revision 1.1  2011/09/23 22:39:01  agnew
             *     Initial revision
             *     Rewritten by D. Agnew from the version by Godkin and Pulli, in BSSA, Vol 74, pp 1847-1848 (1984) */

            int inpoly = 0;

            int isc = 0;
            for (int i = 0; i < nv - 1; i++)                                        //       do 5 i=1,nv-1
            {
                isc = ncross(xv[i] - x, yv[i] - y, xv[i + 1] - x, yv[i + 1] - y);   //       isc = ncross(xv(i)-x,yv(i)-y,xv(i+1)-x,yv(i+1)-y)
                if (isc == 4)                                                       //       if(isc.eq.4) then
                {
                    /*  on edge - know the answer */
                    inpoly = 2;                                                     //          inpoly = 2
                    return inpoly;                                                  //          return
                }
                inpoly += isc;                                                      //       inpoly = inpoly + isc
            }

            /*  check final segment */
            isc = ncross(xv[nv - 1] - x, yv[nv - 1] - y, xv[0] - x, yv[0] - y);     //       isc = ncross(xv(nv)-x,yv(nv)-y,xv(1)-x,yv(1)-y)

            if (isc == 4)                                                           //       if(isc.eq.4) then
            {
                inpoly = 2;                                                         //          inpoly = 2
                return inpoly;                                                      //          return
            }

            inpoly += isc;                                                          //       inpoly = inpoly + isc
            inpoly /= 2;                                                            //       inpoly = inpoly/2

            /*  convert to all positive (a departure from the original) */
            inpoly = Math.Abs(inpoly);                                              //       inpoly = iabs(inpoly)
            return inpoly;                                                          //       return
        }

        /// <summary>
        /// finds whether the segment from point 1 to point 2 crosses the negative x-axis or goes through the origin (this is the signed crossing number)
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="x2"></param>
        /// <param name="y2"></param>
        /// <returns>
        /// <para>return value       nature of crossing </para>
        /// <para>    4               segment goes through the origin </para>
        /// <para>    2               segment crosses from below </para>
        /// <para>    1               segment ends on -x axis from below or starts on it and goes up </para>
        /// <para>    0               no crossing </para>
        /// <para>   -1               segment ends on -x axis from above or starts on it and goes down </para>
        /// <para>   -2               segment crosses from above </para>
        /// </returns>
        private int ncross(double x1, double y1, double x2, double y2)// function ncross(x1,y1,x2,y2)
        {
            int ncross = 0;
            /*   finds whether the segment from point 1 to point 2 crosses the negative x-axis or goes through the origin (this is the signed crossing number)
             *     return value       nature of crossing
             *         4               segment goes through the origin
             *         2               segment crosses from below
             *         1               segment ends on -x axis from below or starts on it and goes up
             *         0               no crossing
             *        -1               segment ends on -x axis from above or starts on it and goes down
             *        -2               segment crosses from above */

            if (y1 * y2 > 0.0)                                  //       if(y1*y2.gt.0) then
            {
                /* all above (or below) axis */
                ncross = 0;                                     //          ncross = 0
                return ncross;                                  //          return
            }

            double c12 = x1 * y2;                               //       c12 = x1*y2
            double c21 = x2 * y1;                               //       c21 = x2*y1

            if (c12 == c21 && x1 * x2 <= 0.0)                   //       if(c12.eq.c21.and.x1*x2.le.0.) then
            {
                /* through origin */
                ncross = 4;                                     //          ncross = 4
                return ncross;                                  //          return
            }
            if (y1 == 0.0 && x1 > 0.0 || y2 == 0.0 && x2 > 0.0 ||
                y1 < 0.0 && c12 > c21 || y1 > 0.0 && c12 < c21 ||
                y1 == 0.0 && y2 == 0.0 && x1 < 0.0 && x2 < 0.0) // if((y1.eq.0.and.x1.gt.0).or.(y2.eq.0.and.x2.gt.0)
                                                                //     .or.((y1.lt.0).and.(c12.gt.c21)).or.((y1.gt.0).and.(c12.lt.c21))
                                                                //     .or.(y1.eq.0.and.y2.eq.0.and.x1.lt.0.and.x2.lt.0)) then
            {
                /* touches +x axis; crosses +x axis; lies entirely on -x axis */
                ncross = 0;                                     //          ncross = 0
                return ncross;                                  //          return
            }

            if (y1 != 0.0 && y2 != 0.0)                         //       if(y1.ne.0.and.y2.ne.0) then
            {
                /*  cross axis */
                if (y1 < 0.0)                                   //          if(y1.lt.0) ncross = 2
                    ncross = 2;
                if (y1 > 0.0)                                   //          if(y1.gt.0) ncross = -2
                    ncross = -2;
                return ncross;                                  //          return
            }

            /* one end touches -x axis - goes which way? */
            if (y1 == 0.0)                                      //       if(y1.eq.0) then
            {
                if (y2 < 0.0)                                   //          if(y2.lt.0) ncross = -1
                    ncross = -1;
                if (y2 > 0.0)                                   //          if(y2.gt.0) ncross = 1
                    ncross = 1;
            }
            else                                                //       else
            {
                /* y2=0 - ends on x-axis */
                if (y1 < 0.0)                                   //          if(y1.lt.0) ncross = 1
                    ncross = 1;
                if (y1 > 0.0)                                   //          if(y1.gt.0) ncross = -1
                    ncross = -1;
            }
            return ncross;                                      //       return
        }
    }
}
