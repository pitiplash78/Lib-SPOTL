using System;
using System.IO;

namespace SPOTL
{
    /// <summary>
    /// Class containing all informations about the water density (reading routine inclusive).
    /// </summary>
    internal class SeaDensity
    {
        /// <summary>
        /// Contains the data of the water density according to the location.
        /// </summary>
        internal int[,] idens;

        public double rho = 0.0;

        /// <summary>
        /// Reads the water density file.
        /// </summary>
        /// <param name="path">File path for the water density.</param>
        /// <returns>The SeaDensity, containing all informations about the water density.</returns>
        public static SeaDensity readWaterDensity(string path)
        {
            FileStream sourceFileStream = File.OpenRead(path);
            System.IO.Compression.GZipStream decompressingStream = new System.IO.Compression.GZipStream(sourceFileStream, System.IO.Compression.CompressionMode.Decompress);

            SeaDensity sead = new SeaDensity();

            sead.idens = new int[72, 36];

            using (var sr = new StreamReader(decompressingStream))
            {
                for (int j = 0; j < sead.idens.GetLength(1); j++)
                {
                    string[] parts = sr.ReadLine().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                    for (int i = 0; i < sead.idens.GetLength(0); i++)
                        sead.idens[i, j] = int.Parse(parts[i]);
                }
            }
            return sead;
        }

        /// <summary>
        /// Checks the surface density table to returns the value at given location
        /// </summary>
        /// <param name="rlat">Latitude.</param>
        /// <param name="rlong">Longitude.</param>
        /// <returns>Returns the surface density at given location.</returns>
        public double seadens(double rlat, double rlong) // subroutine seadens(rlat,rlong,rho)
        {
            /*     $Log: seadens.f,v $ */
            /*     Revision 1.3  2013/03/11 02:03:16  agnew      force local longitude variable to be 0 to 360 */
            /*     Revision 1.2  2011/12/27 21:26:51  agnew      Fixed Arctic Ocean cells */
            /*     Revision 1.1  2011/11/19 19:10:20  agnew      Initial revision */

            /*    Values are densities averaged over monthly values, computed from T and S grids of statistical means in the 2009 World Ocean Atlas
             *    (annual values of T and S used if no others available) */

            /*  If no grid data, value is found by averaging over adjacent nonzero cells. This means that nonzero values will be returned over many
             *    land areas. */

            double rlo = rlong;
            if (rlong < 0.0)                                        //       if(rlong.lt.0) rlo=rlong+360
                rlo = rlong + 360;
            if (rlong >= 360.0)                                     //       if(rlong.ge.360) rlo=rlong-360
                rlo = rlong - 360;

            double x = rlo / 5.0;                                   //       x=rlo/5.
            double y = (rlat + 90.0) / 5.0;                         //       y=(rlat+90.)/5.
            int nx = (int)x;                                        //       nx=x+1 
            int ny = (int)y;                                        //       ny=y+1

            if (ny >= idens.GetLength(1))  // HACK: ist eingebaut weil am rand vom array gibt error
                ny--;

            int id = (int)idens[nx, ny];                            //       id=idens(nx,ny)
            
            if (id < -800)                                          //       if(id.lt.-800) then
            {
                /* no data; look at adjacent cells and average as many as found; crude but probably adequate */
                id = 0;                                             //         id=0
                int no = 0;                                         //       no=0
                if (ny > 1 && idens[nx, ny - 1] > -800)             //         if(ny.gt.1.and.idens(nx,ny-1).gt.-800) then
                {
                    id += idens[nx, ny - 1];                        //           id=id+idens(nx,ny-1)
                    no++;                                           //           no=no+1
                }
                if (ny < 35 && idens[nx, ny + 1] > -800)            //         if(ny.lt.36.and.idens(nx,ny+1).gt.-800) then
                {
                    id += idens[nx, ny + 1];                        //           id=id+idens(nx,ny+1)
                    no++;                                           //           no=no+1
                }
                if (nx > 1 && idens[nx - 1, ny] > -800)             //         if(nx.gt.1.and.idens(nx-1,ny).gt.-800) then
                {
                    id += idens[nx - 1, ny];                        //           id=id+idens(nx-1,ny)
                    no++;                                           //           no=no+1
                }
                if (nx < 71 && idens[nx + 1, ny] > -800)            //         if(nx.lt.72.and.idens(nx+1,ny).gt.-800) then
                {
                    id += idens[nx + 1, ny];                        //           id=id+idens(nx+1,ny)
                    no++;                                           //           no=no+1
                }

                if (no > 1)                                         //         if(no.gt.1) id=id/no
                    id /= no;

                /*  no nonzero cells found; return a value of zero */
                if (no == 0)                                        //         if(no.eq.0) id=-100000
                    id = -100000;
            }
            this.rho = (double)id / 100.0 + 1e3;                    //       rho=1000.+id/100.
            return rho;                        
        }
    }
}
