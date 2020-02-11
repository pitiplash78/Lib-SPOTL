using System;
using System.Collections.Generic;
using System.Collections;
using System.IO;
using System.IO.Compression;
using System.Text;
using System.Xml.Serialization;
using System.Reflection;
using System.Numerics;

namespace SPOTL
{
    // TODO: mit neuester version vergleichen & updaten | sauber machen | überall sinnvolle summaries hin | inputparameter überprüfen ggf. verbessern
    // TODO: Unit bei jeder berechnung mit angeben oder bei Construtor mit setzen ? -> gleich mit solid earth tides Berechnung 
    /// <summary>
    /// Calculation of the Ocean Tide Loading depending on the (class) OceanLoadingProperties.
    /// </summary>
    public class OceanLoading 
    {
        #region Location
        public Location _Location = null;
        public class Location
        {
            internal double Longitude = double.NaN;
            internal double Latitude = double.NaN;
            internal double Height = double.NaN;

            public Location(double Longitude, double Latitude, double Height)
            {
                this.Longitude = Longitude;
                this.Latitude = Latitude;
                this.Height = Height;
            }
        }
        #endregion

        #region Enum Component
        public enum Component
        {
            None,
            Gravity,
            GravityDeformation,
            GravityAttraction,
            OceanTide,
            Tilt,
            TiltDeformation,
            TiltAttraction,
            HorizontalStrain,
            ShearStrain,
            ArealStrain,
            VolumeStrain,
            VerticalDisplacement,
            HorizontalDisplacement,
        }
        #endregion

        #region Parameter for constructor

        public string PathOceanLoading
        {
            get 
            {
                return _PathOceanLoading; 
            }
            set
            {
                _PathOceanLoading = value;
                if (_PathOceanLoading != null)
                {
                    if (filePath == null)
                        filePath = new FilePath();

                    filePath.EarthModels = _PathOceanLoading + "EarthModels" + Path.DirectorySeparatorChar;
                    filePath.OceanModels = _PathOceanLoading + "OceanModels" + Path.DirectorySeparatorChar;
                    filePath.LocalModels = _PathOceanLoading + "LocalModels" + Path.DirectorySeparatorChar;

                    filePath.WaterDensity = _PathOceanLoading + "WaterDensity" + Path.DirectorySeparatorChar + "seadens.gz";
                    filePath.LndSeaMatrix = _PathOceanLoading + "LandSeaMatrix" + Path.DirectorySeparatorChar + "landsea.matrix.gz";
                }
            }
        }

        internal string _PathOceanLoading = Environment.CurrentDirectory + Path.DirectorySeparatorChar;

        private FilePath filePath = null;
       
        public OceanLoadingProperties.GreenFunctionEntry GreensFunctionModel { get; private set; }

        public OceanLoadingProperties.OceanModelEntry[] OceanTideModel { get; private set; }

        public OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanTideModel { get; private set; }

        /// <summary>
        /// Unit parameter for contructor for calculation the gravity.
        /// </summary>
        public Utilities.Units.UnitNamesEnum OutputUnit_Gravity
        {
            get
            {
                return _OutputUnit_Gravity;
            }
            set
            {
                _OutputUnit_Gravity = value;
                dConversion_Gravity = Utilities.Units.getConversion(BaseUnits.UnitGravity, _OutputUnit_Gravity).Factor;
            }
        }
        internal Utilities.Units.UnitNamesEnum _OutputUnit_Gravity = BaseUnits.UnitGravity;
        internal double dConversion_Gravity = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the potential height.
        /// </summary>
        public Utilities.Units.UnitNamesEnum OutputUnit_PotentialHeight
        {
            get
            {
                return _OutputUnit_PotentialHeight;
            }
            set
            {
                _OutputUnit_PotentialHeight = value;
                dConversion_PotentialHeight = Utilities.Units.getConversion(BaseUnits.UnitGravity, _OutputUnit_PotentialHeight).Factor;
            }
        }
        internal Utilities.Units.UnitNamesEnum _OutputUnit_PotentialHeight = BaseUnits.UnitPotentialHeight;
        internal double dConversion_PotentialHeight = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the tilt.
        /// </summary>
        public Utilities.Units.UnitNamesEnum OutputUnit_Tilt
        {
            get
            {
                return _OutputUnit_Tilt;
            }
            set
            {
                _OutputUnit_Tilt = value;
                dConversion_Tilt = Utilities.Units.getConversion(BaseUnits.UnitTilt, _OutputUnit_Tilt).Factor;
            }
        }
        internal Utilities.Units.UnitNamesEnum _OutputUnit_Tilt = BaseUnits.UnitTilt;
        internal double dConversion_Tilt = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the strain.
        /// </summary>
        public Utilities.Units.UnitNamesEnum OutputUnit_Strain
        {
            get
            {
                return _OutputUnit_Strain;
            }
            set
            {
                _OutputUnit_Strain = value;
                dConversion_Strain = Utilities.Units.getConversion(BaseUnits.UnitStrain, _OutputUnit_Strain).Factor;
            }
        }
        internal Utilities.Units.UnitNamesEnum _OutputUnit_Strain = BaseUnits.UnitStrain;
        internal double dConversion_Strain = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the displacement.
        /// </summary>
        public Utilities.Units.UnitNamesEnum OutputUnit_Displacement
        {
            get
            {
                return _OutputUnit_Displacement;
            }
            set
            {
                _OutputUnit_Displacement = value;
                dConversion_Displacement = Utilities.Units.getConversion(BaseUnits.UnitStrain, _OutputUnit_Displacement).Factor;
            }
        }
        internal Utilities.Units.UnitNamesEnum _OutputUnit_Displacement = BaseUnits.UnitDisplacement;
        internal double dConversion_Displacement = 1d;

        public static class BaseUnits
        {
            public static Utilities.Units.UnitNamesEnum UnitGravity = Utilities.Units.UnitNamesEnum.MicroGal;
            public static Utilities.Units.UnitNamesEnum UnitPotentialHeight = Utilities.Units.UnitNamesEnum.Millimeter;
            public static Utilities.Units.UnitNamesEnum UnitDisplacement = Utilities.Units.UnitNamesEnum.Millimeter;
            public static Utilities.Units.UnitNamesEnum UnitTilt = Utilities.Units.UnitNamesEnum.NanoRadian;
            public static Utilities.Units.UnitNamesEnum UnitStrain = Utilities.Units.UnitNamesEnum.NanoStrain;
        }
        #endregion

        #region Wave parameter list
        private ArrayList waveParamList = null;
        private class waveParam
        {
            public string wave;
            public int[] icte;
            public Complex grav;
            public Complex gravdef;
            public Complex pothi;
            public Complex[] disp;
            public Complex[] tilt;
            public Complex[] tiltdef;
            public Complex[] strain;

            public waveParam(string wave, int[] icte, Complex grav, Complex gravdef,
                                                      Complex pothi,
                                                      Complex[] disp,
                                                      Complex[] tilt, Complex[] tiltdef,
                                                      Complex[] strain)
            {
                this.wave = wave;
                this.icte = icte;
                this.grav = grav;
                this.gravdef = gravdef;
                this.pothi = pothi;
                this.disp = disp;
                this.tilt = tilt;
                this.tiltdef = tiltdef;
                this.strain = strain;
            }
        }

        private class hartid_wPM
        {
            public int[] icte;
            public double amp;
            public double phase;

            public hartid_wPM(int[] icte, double amp, double phase)
            {
                this.icte = icte;
                this.amp = amp;
                this.phase = phase;
            }
        }
        #endregion

        public OceanLoading(Location location, 
                            OceanLoadingProperties.GreenFunctionEntry GreensFunctionModel,  
                            OceanLoadingProperties.OceanModelEntry[] OceanTideModel,
                            OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanTideModel = null)
        {
            if (_Location == null)
                _Location = new Location(location.Longitude,
                                         location.Latitude,
                                         location.Height);

            this.GreensFunctionModel = GreensFunctionModel;
            this.OceanTideModel = OceanTideModel;
            if (LocalOceanTideModel != null)
                this.LocalOceanTideModel = LocalOceanTideModel;
        }
          
        // Public methods for caluculation
        /// <summary>
        /// Calculates value for ocean loading of the according component for given time and Unit and azimuth (if needed for component)
        /// </summary>
        /// <param name="mjd">Time in modified Julian day to calculate the the needed component value.</param>
        /// <param name="component">Tidal component to be calculated.</param>
        /// <param name="azimuth">Optional parameter for components (Horizontal Strain, Tilt having directions  (default 0°== north) </param>
        /// <returns>Calculated value of the according component for given time and Unit and azimuth (if needed for component).</returns>
        public double calculateComponent(double mjd, Component component, double azimuth = 0.0)
        {
            return calculateComponent(Utilities.TimeConversion.ModJulDay2Time(mjd), component, azimuth);
        }

        /// <summary>
        /// Calculates value for ocean loading of the according component for given time and Unit and azimuth (if needed for component)
        /// </summary>
        /// <param name="dt">Time to calculate the the needed component value.</param>
        /// <param name="component">Tidal component to be calculated.</param>
        /// <param name="azimuth">Optional parameter for components (Horizontal Strain, Tilt having directions  (default 0°== north) </param>
        /// <returns>Calculated value of the according component for given time and Unit and azimuth (if needed for component).</returns>
        public double calculateComponent(DateTime dt, Component component, double azimuth = 0.0)
        {
            if(waveParamList == null)
                waveParamList = nfload(GreensFunctionModel, OceanTideModel);

            /* adjust to a new azimuth, rotated azim deg (clockwise in the command line,
             * converted to counterclockwise below to match the usual formulae) from the
             * EW-NS one that the loads are originally in. */
            double caz = Math.Cos(-azimuth * Utilities.Constants.radian2degree);
            double saz = Math.Sin(-azimuth * Utilities.Constants.radian2degree);

            // return value
            double result = double.NaN;

            if (component == Component.Gravity)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, wpTmp.grav.Real, wpTmp.grav.Imaginary);
                }    
                result = (hartid(dt, hartidPM) * dConversion_Gravity);
            }
            else if (component == Component.GravityDeformation)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, wpTmp.gravdef.Real, wpTmp.gravdef.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Gravity);
            }
            else if (component == Component.GravityAttraction)
            {
                hartid_wPM[] hartidPMGravity = new hartid_wPM[waveParamList.Count];
                hartid_wPM[] hartidPMGravityDeformation = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    hartidPMGravity[i] = new hartid_wPM(wpTmp.icte, wpTmp.grav.Real, wpTmp.grav.Imaginary);
                    hartidPMGravityDeformation[i] = new hartid_wPM(wpTmp.icte, wpTmp.gravdef.Real, wpTmp.gravdef.Imaginary);
                }
                result = (hartid(dt, hartidPMGravity) - hartid(dt, hartidPMGravityDeformation)) * dConversion_Gravity;
            }
            else if( component == Component.OceanTide)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, wpTmp.pothi.Real, wpTmp.pothi.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_PotentialHeight);
            }
            else if (component == Component.Tilt)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    Complex tilt = phasor(-saz * argand(wpTmp.tilt[1]) + caz * argand(wpTmp.tilt[0]));
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, tilt.Real, tilt.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Tilt);
            }
            else if (component == Component.TiltDeformation)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    Complex tiltdef = phasor(-saz * argand(wpTmp.tiltdef[1]) + caz * argand(wpTmp.tiltdef[0]));
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, tiltdef.Real, tiltdef.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Tilt);
            }
            else if (component == Component.TiltAttraction)
            {
                hartid_wPM[] hartidPMTilt = new hartid_wPM[waveParamList.Count];
                hartid_wPM[] hartidPMTiltDef = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    Complex tilt = phasor(-saz * argand(wpTmp.tilt[1]) + caz * argand(wpTmp.tilt[0]));
                    hartidPMTilt[i] = new hartid_wPM(wpTmp.icte, tilt.Real, tilt.Imaginary);
                    Complex tiltdef = phasor(-saz * argand(wpTmp.tiltdef[1]) + caz * argand(wpTmp.tiltdef[0]));
                    hartidPMTiltDef[i] = new hartid_wPM(wpTmp.icte, tiltdef.Real, tiltdef.Imaginary);
                }
                result = ( hartid(dt, hartidPMTilt) - hartid(dt, hartidPMTiltDef)) * dConversion_Tilt;
            }

            else if (component == Component.HorizontalStrain)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    Complex strain = phasor(caz * caz * argand(wpTmp.strain[1]) + saz * saz * argand(wpTmp.strain[0]) - 2 * caz * saz * argand(wpTmp.strain[2]));
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, strain.Real, strain.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Strain);
            }
            else if (component == Component.ShearStrain)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    Complex strain = phasor(argand(wpTmp.strain[2]) * (caz * caz - saz * saz) - caz * saz * (argand(wpTmp.strain[0]) - argand(wpTmp.strain[1])));
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, strain.Real, strain.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Strain);
            }
            else if (component == Component.ArealStrain)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    Complex strain = phasor((argand(wpTmp.strain[0]) + argand(wpTmp.strain[1]))); // TODO: checken ob das wirklich richtig ist
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, strain.Real, strain.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Strain);
            }
            else if (component == Component.VolumeStrain)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    Complex strain = phasor(0.6667d * (argand(wpTmp.strain[0]) + argand(wpTmp.strain[1])));
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, strain.Real, strain.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Strain);
            }
            else if (component == Component.VerticalDisplacement)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < hartidPM.Length; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, wpTmp.disp[0].Real, wpTmp.disp[0].Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Displacement);
            }
            else if (component == Component.HorizontalDisplacement)
            {
                hartid_wPM[] hartidPM = new hartid_wPM[waveParamList.Count];
                for (int i = 0; i < waveParamList.Count; i++)
                {
                    waveParam wpTmp = (waveParam)waveParamList[i];
                    Complex displacement = phasor(-saz * argand(wpTmp.disp[2]) + caz * argand(wpTmp.disp[1]));
                    hartidPM[i] = new hartid_wPM(wpTmp.icte, displacement.Real, displacement.Imaginary);
                }
                result = (hartid(dt, hartidPM) * dConversion_Displacement);
            }
            return result;
        }

        // TODO was ist mit den Lokal Modellen 
        #region calculates the wave parameters with respect to the oceanloading
        private ArrayList nfload(OceanLoadingProperties.GreenFunctionEntry gfModel, OceanLoadingProperties.OceanModelEntry[] oclModel)
        {
            // return 
            ArrayList wpl = new ArrayList();

            double rlong = _Location.Longitude;
            double rlato = _Location.Latitude;

            ////////////////////////
            LandseaMatrix lm = LandseaMatrix.readLandSea(filePath.LndSeaMatrix);
            GreensFunction gf = GreensFunction.readGreensfunction(filePath.EarthModels + gfModel.FileName);

            ArrayList oclModelList = new ArrayList();
            for (int i = 0; i < oclModel.Length; i++)
                for ( int j = 0; j< oclModel[i].Waves.Length; j++)
                    oclModelList.Add(filePath.OceanModels + oclModel[i].Waves[j] +oclModel[i].FilePath);
            string[] oclModelPath = (string[])oclModelList.ToArray();
            OceanModel tm = OceanModel.readTidalModel(oclModelPath);

            double ct = Math.Cos(Utilities.Constants.degree2radian * (90d - _Location.Latitude));
            double st = Math.Sin(Utilities.Constants.degree2radian * (90d - _Location.Latitude));

            const int minaz = 150;
            //  const double rho = 1.025e3; - removed in Revision 1.2  2012/03/06 03:32:43  agnew
 

            double clat = double.NaN;
            double clong = double.NaN;

            double close = double.NaN;
            int nclose = 0;

            for (int i = 0; i < tm.TideModelItems.Length; i++)
            {
                OceanModel.TMelement TMe = (OceanModel.TMelement)tm.TideModelItems[i];

                Complex zero = Complex.Zero;
                Complex[] disp = new Complex[3] { zero, zero, zero };
                Complex grav = Complex.Zero;
                Complex pothi = Complex.Zero;
                Complex[] tilt = new Complex[2] { zero, zero };
                Complex[] str = new Complex[3] { zero, zero, zero };

                Complex gravdef = Complex.Zero;
                Complex[] tiltdef = new Complex[2] { zero, zero };

                Complex amp = Complex.Zero;
                Complex cmp = Complex.Zero;

                double dist = 0;
                double numnf = 0;
                double clnf = 0;
                double farnf = 0;

                //  convolutions start here - read in Green functions, which come in
                //  ntot sections; num is the current one, fingrd is whether or not to
                //  use the detailed land-sea grid 
                for (int num = 0; num < gf.GrennFunctionItems.Length; num++)
                {
                    GreensFunction.GFelememt GFe = (GreensFunction.GFelememt)gf.GrennFunctionItems[num];
                    double stp = Utilities.Constants.degree2radian * GFe.spc;
                    //  if this is the first convolution interval see if (1) there is ocean at
                    //  the station location and (2) if it is below sea level. If both are true,
                    //  compute the gravity load (only) for a disk from distance 0 to the start
                    //  of the first integrated Green function 

                    if (num == 0)
                    {
                        amp = ocmodl(rlong, rlato, TMe, GFe, lm, ref dist, ref numnf, ref clnf, ref farnf);

                        if (!Complex.Equals(amp, zero) && _Location.Height < 0)
                        {
                            //amp = rho * amp;
                            double del = Utilities.Constants.degree2radian * GFe.beg / 2;
                            double g = 0d;
                            double t = 0;
                            double pot = 0;
                            ignewt(del, 2 * del, _Location.Longitude, _Location.Latitude, _Location.Height, ref g, ref t, ref pot);
                            grav = grav + g * amp * 2d * Math.PI;
                            close = 0d;
                            clat = _Location.Latitude;
                            clong = _Location.Longitude;
                            nclose = 1;
                        }

                    }
                    //  loop through distance begins here
                    for (int ii = 0; ii < GFe.grfn.GetLength(0); ii++)
                    {
                        dist = GFe.beg + ((double)ii) * GFe.spc;
                        double del = Utilities.Constants.degree2radian * dist;
                        double cd = Math.Cos(del);
                        double sd = Math.Sin(del);

                        double g = 0;
                        double t = 0;
                        double pot = 0;
                        ignewt(del, stp, _Location.Longitude, _Location.Latitude, _Location.Height, ref g, ref t, ref pot);

                        g = g + GFe.grfn[ii, 2];
                        t = t + GFe.grfn[ii, 3];

                        double gdef = GFe.grfn[ii, 2];
                        double tdef = GFe.grfn[ii, 3];

                        pot = pot + GFe.grfn[ii, 6];
                        double u = GFe.grfn[ii, 0];
                        double v = GFe.grfn[ii, 1];
                        double ett = GFe.grfn[ii, 4];
                        double ell = GFe.grfn[ii, 5];

                        // find the number of azimuth steps (scaled with sin(del), but kept above some minimum)
                        int naz = (int)(360d * sd);
                        if (dist < 180d)
                            naz = Math.Max(naz, minaz);
                        double azstp = (2d * Math.PI) / naz;
                        double azstpd = 360d / naz;
                        double azimuth = azstpd / 2d;
                        double caz = Math.Cos(azstp / 2);
                        double saz = Math.Sin(azstp / 2);
                        double saztp = Math.Sin(azstp);
                        double caztp = Math.Cos(azstp);
                        double stpfac = 2 * saz;

                        //  loop through azimuth begins here

                        for (int jj = 0; jj < naz; jj++)
                        {
                            //  do the spherical trigonometry to find the cell lat and long
                            double cb = cd * ct + sd * st * caz;
                            if (Math.Abs(cb) > 1d)
                                cb = cb / Math.Abs(cb);
                            double sb = Math.Sqrt(1d - cb * cb);
                            rlato = 90 - Utilities.Constants.radian2degree * Math.Acos(cb);
                            rlong = 0d;
                            //  if near the north pole leave longitude at zero- otherwise find it
                            if (sb > 1e-3)
                            {
                                sb = 1d / sb;
                                double sg = sd * saz * sb;
                                double cg = (st * cd - sd * ct * caz) * sb;
                                rlong = _Location.Longitude + Utilities.Constants.radian2degree * Math.Atan2(sg, cg);
                            }

                            amp = ocmodl(rlong, rlato, TMe, GFe, lm, ref dist, ref numnf, ref clnf, ref farnf);
                           // amp = rho * amp;
                            if (!Complex.Equals(amp, zero))
                            {
                                // if this is the first nonzero amplitude, save the distance and location
                                if (nclose == 0)
                                {
                                    close = dist;
                                    clat = rlato;
                                    clong = rlong;
                                    nclose = 1;
                                }
                                cmp = amp * stpfac;
                                //  compute the loads.   
                                grav = grav + g * amp * azstp;
                                gravdef = gravdef + gdef * amp * azstp;
                                pothi = pothi + pot * amp * azstp;
                                // disp is vert,ns,ew;
                                disp[0] = disp[0] + u * amp * azstp;
                                disp[1] = disp[1] + v * caz * cmp;
                                disp[2] = disp[2] + v * saz * cmp;
                                // tilt is ns,ew;
                                tilt[0] = tilt[0] + t * caz * cmp;
                                tilt[1] = tilt[1] + t * saz * cmp;
                                tiltdef[0] = tiltdef[0] + tdef * caz * cmp;
                                tiltdef[1] = tiltdef[1] + tdef * saz * cmp;
                                // strain is ns,ew,ne shear
                                double aa = 0.5d * (caz * caz - saz * saz);
                                double bb = 0.5d * azstp + aa * stpfac;
                                double cc = 0.5d * azstp - aa * stpfac;
                                str[0] = str[0] + (ett * bb + ell * cc) * amp;
                                str[1] = str[1] + (ell * bb + ett * cc) * amp;
                                str[2] = str[2] + (ett - ell) * saz * caz * saztp * caztp * amp;
                            }
                            //  use recursion to step sin and cos of the azimuth
                            azimuth = azimuth + azstpd;
                            double xx = saz * caztp + caz * saztp;
                            caz = caz * caztp - saz * saztp;
                            saz = xx;
                        }
                    }
                }
                //c  done with convolution - add to disk file, and convert to amp and phase
                double isp = (double)TMe.icte[0];
                Complex toloc = new Complex(Math.Cos(Utilities.Constants.degree2radian * isp * _Location.Longitude),
                                            Math.Sin(Utilities.Constants.degree2radian * isp * _Location.Longitude));
                toloc = new Complex(1d, 0d);

                //      gravity
                grav = phasor(Complex.Conjugate(1.0e8 * toloc * grav));
                gravdef = phasor(Complex.Conjugate(1.0e8 * gravdef));

                //      potential height
                pothi = phasor(Complex.Conjugate(1.0e3 * toloc * pothi));

                //      displacement (output as Z, N, E) change sign of N and E (see CJ II-143, 5 Jan 82)
                disp[0] = phasor(Complex.Conjugate(1.0e3 * toloc * disp[0]));
                disp[1] = phasor(Complex.Conjugate(-1.0e3 * toloc * disp[1]));
                disp[2] = phasor(Complex.Conjugate(-1.0e3 * toloc * disp[2]));

                //      tilt (output as N, E)
                tilt[0] = phasor(Complex.Conjugate(1.0e9 * toloc * tilt[0]));
                tilt[1] = phasor(Complex.Conjugate(1.0e9 * toloc * tilt[1]));
                tiltdef[0] = phasor(Complex.Conjugate(1.0e9 * toloc * tiltdef[0]));
                tiltdef[1] = phasor(Complex.Conjugate(1.0e9 * toloc * tiltdef[1]));

                //      strain (output as N, E, EN shear)
                str[0] = phasor(Complex.Conjugate(1.0e9 * toloc * str[0]));
                str[1] = phasor(Complex.Conjugate(1.0e9 * toloc * str[1]));
                str[2] = phasor(Complex.Conjugate(1.0e9 * toloc * str[2]));

                wpl.Add(new waveParam(TMe.dsym, TMe.icte, grav, gravdef, pothi, disp, tilt, tiltdef, str));
            }
            Complex cz = Complex.Zero;
            wpl.Add(new waveParam("0", new int[6] { -1, 0, 0, 0, 0, 0 },
                cz, cz, cz, new Complex[3] { cz, cz, cz }, new Complex[] { cz, cz }, new Complex[] { cz, cz }, new Complex[3] { cz, cz, cz }));

            return wpl;
        }

        #region functions for Loadings
        private Complex ocmodl(double rlong, double rlato,
                          OceanModel.TMelement TMe, GreensFunction.GFelememt GFe, LandseaMatrix lsM,
                                ref double dist, ref double numnf, ref double clnf, ref double farnf)
        {
            //  For given latitude and east longitude (in degrees and >= 0) returns
            //  the complex amplitude of the tide (in meters, with greenwich phase).
            //  this amplitude is of course zero on land.

            double x = 0;
            double y = 0;
            int i = 0;
            int j = 0;
            int ind = 0;

            Complex amp = Complex.Zero;

            celfnd(rlong, rlato, TMe, ref x, ref y, ref i, ref j, ref ind);

            if (ind > (TMe.latc * TMe.longc) || ind <= 0)
                return amp;

            if (GFe.fingrd != "F")// fingrd.ne.'F')
            {
                //  large distance - return the model value (0 if on land)
                amp = 0.001d * (new Complex((double)TMe.ir1[ind - 1], (double)TMe.im1[ind - 1]));
                //         amp = .001*cmplx(float(ir1(ind)),float(im1(ind)))
                return amp;
            }
            else
            {
                //  close distance - use the
                //  5-minute map file to decide if there is land or not. If not, bilinearly
                //  interpolate from the ocean model, filling out empty cells to make this
                //  possible (though if all the cells needed are zero, return this).

                if (lsM.getLandSea(rlato, rlong))
                {
                    //  detailed grid says land - return 0 
                    return amp = Complex.Zero;
                }
                else
                {
                    //  detailed map says ocean -load an interpolating array of four elements
                    //  with values from adjacent cells, depending on which quadrant of the cell
                    //  we are in. In the process note how many nonzero values we have; save one
                    int ix = 1;
                    if (x <= 0) ix = -ix;
                    int iy = -TMe.longc;
                    if (y <= 0) iy = -iy;
                    x = Math.Abs(x);
                    y = Math.Abs(y);

                    int nnz = 0;
                    int ibind = 0;
                    Complex[] bg = new Complex[4];
                    Complex bgs = new Complex();

                    for (int ii = 0; ii < 4; ii++)
                    {
                        if (ii == 0) ibind = ind;
                        if (ii == 1) ibind = ind + ix;
                        if (ii == 2) ibind = ind + ix + iy;
                        if (ii == 3) ibind = ind + iy;
                        if (ibind <= 0 || ibind > TMe.latc * TMe.longc)
                            bg[ii] = Complex.Zero;
                        if (ibind > 0 && ibind <= TMe.latc * TMe.longc)
                        {
                            bg[ii] = 0.001d * (new Complex((double)TMe.ir1[ibind - 1],
                                                (double)TMe.im1[ibind - 1]));
                            //  bg(i) = 0.001*cmplx(float(ir1(ibind)),float(im1(ibind)))
                        }
                        if (!Complex.Equals(bg[ii], Complex.Zero))
                        {
                            nnz++;
                            bgs = bg[ii];
                        }
                    }
                    //                      now we have a grid for interpolation: the values are arranged as
                    //
                    //           4      3
                    //           1      2
                    //
                    //  where (1) is the value we are closest to.
                    //
                    //  we now fill out the zero values (if any) from the nonzero ones
                    //c  (Except that if there are 0 or 4 nonzero values we needn't bother)

                    if (nnz == 1)
                    {
                        //  the one saved value is the nonzero one: fill out the whole array with this
                        for (int ii = 0; ii < bg.Length; ii++)
                            bg[ii] = bgs;
                    }
                    else if (nnz == 2)
                    {
                        //  two values: either in a row, a column, or diagonal:
                        if (!Complex.Equals(bg[0], Complex.Zero) &&
                            !Complex.Equals(bg[1], Complex.Zero))
                        {
                            bg[3] = bg[0];
                            bg[2] = bg[1];
                        }
                        if (!Complex.Equals(bg[1], Complex.Zero) &&
                            !Complex.Equals(bg[2], Complex.Zero))
                        {
                            bg[3] = bg[2];
                            bg[0] = bg[1];
                        }
                        if (!Complex.Equals(bg[2], Complex.Zero) &&
                            !Complex.Equals(bg[3], Complex.Zero))
                        {
                            bg[0] = bg[3];
                            bg[1] = bg[2];
                        }
                        if (!Complex.Equals(bg[3], Complex.Zero) &&
                            !Complex.Equals(bg[0], Complex.Zero))
                        {
                            bg[2] = bg[3];
                            bg[1] = bg[0];
                        }
                        //     diagonals
                        if (!Complex.Equals(bg[3], Complex.Zero) &&
                            !Complex.Equals(bg[1], Complex.Zero))
                        {
                            bg[2] = 0.5d * (bg[3] + bg[1]);
                            bg[0] = bg[2];
                        }
                        if (!Complex.Equals(bg[0], Complex.Zero) &&
                            !Complex.Equals(bg[2], Complex.Zero))
                        {
                            bg[3] = 0.5d * (bg[0] + bg[2]);
                            bg[1] = bg[2];
                        }
                    }
                    else if (nnz == 3)
                    {
                        // one missing value
                        if (Complex.Equals(bg[0], Complex.Zero)) bg[0] = bg[3] + bg[1] - bg[2];
                        if (Complex.Equals(bg[1], Complex.Zero)) bg[1] = bg[0] + bg[2] - bg[3];
                        if (Complex.Equals(bg[2], Complex.Zero)) bg[2] = bg[3] + bg[1] - bg[0];
                        if (Complex.Equals(bg[3], Complex.Zero)) bg[3] = bg[0] + bg[2] - bg[1];
                    }
                    //  bilinear interpolation
                    amp = bg[0] * (1 - x) * (1 - y) +
                          bg[1] * x * (1 - y) +
                          bg[3] * (1 - x) * y +
                          bg[2] * x * y;
                    //        amp = bg(1)*(1-x)*(1-y) + bg(2)*x*(1-y) +  bg(4)*(1-x)*y     + bg(3)*x*y

                    //  if all the values are zero, we will return zero - but we pass some
                    //  information about this through common block coninf; the first time this
                    //  happens, set clnf and farnf to the current distance from the station,
                    //  and afterwards only farnf (since we always increase the distance).
                    if (nnz == 0)
                    {
                        if (numnf == 0)
                        {
                            clnf = dist;
                            farnf = dist;
                        }
                        else
                            farnf = dist;
                        numnf = numnf + 1;
                    }

                }
                return amp;
            }
        }

        #region
        /*  common block holds polygon file information */
        internal class PolInf
        {
            // public bool ispoly = false;
            public Data[] npoly; // length & model
            public class Data
            {
                public int npoly; // length of xoply & ypoly
                public double[] xpoly;
                public double[] ypoly;
                public string polyf;
                public string polynm;
                public bool use = true;
            }
        }

        internal PolInf polinf = new PolInf();

        /* Common Block Declarations */

        //struct {
        //    char polnam[80], alpoly[1];
        //} polmor_;

        //#define polmor_1 polmor_


        /*     Revision 1.3  2012/03/13 21:18:10  agnew Modified polygon common block to include polygon names 
         *     Revision 1.2  2011/11/27 04:32:37  agnew larger dimension for polygon information
         *     Revision 1.1  2011/09/23 22:39:02  agnew Initial revision */


        private void readLocalModels(string path) // int plstart_()
        {
            FileStream sourceFileStream = File.OpenRead(path);
            GZipStream decompressingStream = new GZipStream(sourceFileStream, CompressionMode.Decompress);

            /*K1  
                  1  1  0  0  0  0
                     -30    -983
                     -45     -17
                     155      17
                     129     983
                     421     751
                OSU Regional Models (2010): Tasmania and SE Austra
                      0      0      0      0      0      0      0      0      0      0 */

            using (var sr = new StreamReader(decompressingStream))
            {
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
                Console.WriteLine(sr.ReadLine());
            }


            /* Format strings */
            string fmt_102 = "(a)";
            string fmt_104 = "(a1)";

            /* System generated locals */
            int i__1, i__2;

            /* Local variables */
            int i__, j, llu;
            char iuse;

            /*    Reads in, if requested, information about */
            /*  polygonal areas (up to 10) within which the point either must or */
            /*  must not fall. */

            //       logical*1 ispoly,use
            ////       character*1 iuse,alpoly
            ////       character*80 polyf,polnam,polynm
            ///*  common block holds polygon file information */
            ////       
            ////       common/polmor/polnam,alpoly
            ////       np=0
            //    polinf_1.np = 0;
            ////       llu = 4
            //    llu = 4;
            ////       open(llu,file=polyf,status='old')
            //    o__1.oerr = 0;
            //    o__1.ounit = llu;
            //    o__1.ofnmlen = 80;
            //    o__1.ofnm = polinf_1.polyf;
            //    o__1.orl = 0;
            //    o__1.osta = "old";
            //    o__1.oacc = 0;
            //    o__1.ofm = 0;
            //    o__1.oblnk = 0;
            //    f_open(&o__1);
            ///* first read name for whole file, and number of polygons; for each one, */
            ///* then read name (not used), number of points, + for point must be in, */
            ///* - for may not be, and then the points (long and lat) themselves. */
            ///*   The polygon coordinates may not cross a 360-degree discontinuity */
            ///* and must be in the range -180 to 180 or 0-360 (E +) */
            ////       read(llu,102) polnam
            //    io___2.ciunit = llu;
            //    s_rsfe(&io___2);
            //    do_fio(&c__1, polmor_1.polnam, (ftnlen)80);
            //    e_rsfe();
            ////  102  format(a)
            ////       read(llu,*) np
            //    io___3.ciunit = llu;
            //    s_rsle(&io___3);
            //    do_lio(&c__3, &c__1, (char *)&polinf_1.np, (ftnlen)sizeof(int));
            //    e_rsle();
            ////       do i=1,np
            //    i__1 = polinf_1.np;
            //    for (i__ = 1; i__ <= i__1; ++i__) {
            ////         read(llu,102) polynm(i)
            //    io___5.ciunit = llu;
            //    s_rsfe(&io___5);
            //    do_fio(&c__1, polinf_1.polynm + (i__ - 1) * 80, (ftnlen)80);
            //    e_rsfe();
            ////         read(llu,*) npoly(i)
            //    io___6.ciunit = llu;
            //    s_rsle(&io___6);
            //    do_lio(&c__3, &c__1, (char *)&polinf_1.npoly[i__ - 1], (ftnlen)sizeof(
            //        int));
            //    e_rsle();
            ////         read(llu,104) iuse
            //    io___7.ciunit = llu;
            //    s_rsfe(&io___7);
            //    do_fio(&c__1, iuse, (ftnlen)1);
            //    e_rsfe();
            ////  104    format(a1)
            ////         if(iuse.eq.'+') use(i) = .true.
            //    if (*(unsigned char *)iuse == '+') {
            //        polinf_1.use[i__ - 1] = TRUE_;
            //    }
            ////         if(iuse.eq.'-') use(i) = .false.
            //    if (*(unsigned char *)iuse == '-') {
            //        polinf_1.use[i__ - 1] = FALSE_;
            //    }
            ///* global override from variable alpoly, read from command line */
            ////         if(alpoly.eq.'+') use(i) = .true.
            //    if (*(unsigned char *)polmor_1.alpoly == '+') {
            //        polinf_1.use[i__ - 1] = TRUE_;
            //    }
            ////         if(alpoly.eq.'-') use(i) = .false.
            //    if (*(unsigned char *)polmor_1.alpoly == '-') {
            //        polinf_1.use[i__ - 1] = FALSE_;
            //    }
            ////         do j=1,npoly(i)
            //    i__2 = polinf_1.npoly[i__ - 1];
            //    for (j = 1; j <= i__2; ++j) {
            ////           read(llu,*) xpoly(j,i),ypoly(j,i)
            //        io___10.ciunit = llu;
            //        s_rsle(&io___10);
            //        do_lio(&c__4, &c__1, (char *)&polinf_1.xpoly[j + i__ * 500 - 501],
            //             (ftnlen)sizeof(double));
            //        do_lio(&c__4, &c__1, (char *)&polinf_1.ypoly[j + i__ * 500 - 501],
            //             (ftnlen)sizeof(double));
            //        e_rsle();
            ////         enddo
            //    }
            ////       enddo
            //    }
            ////       close(llu)
            //    cl__1.cerr = 0;
            //    cl__1.cunit = llu;
            //    cl__1.csta = 0;
            //    f_clos(&cl__1);
            //       return
            //return 0;
            //       end
        } /* plstart_ */

        /// <summary>
        /// 
        /// </summary>
        /// <param name="rlong"></param>
        /// <param name="rlat"></param>
        /// <returns></returns>
        private bool chkgon(double rlong, double rlat)//, int iok)     //       subroutine chkgon(rlong,rlat,iok)
        {
            /* returns iok=0 if
             *  1. there is any polygon (of all those read in) in which the
             *     coordinate should not fall, and it does
             *     or
             * 2. the coordinate should fall in at least one polygon
             *    (of those read in) and it does not otherwise returns iok=1 */
            bool iok = false;

            /* Local variables */
            double[] xu = new double[500];
            double[] yu = new double[500];

            //       
            /* first convert rlong to -180/180 (or leave in this, or 0/360) */
            if (rlong < -180.0)                 //       if(rlong.lt.-180) rlong = rlong + 360
                rlong += 360;

            /* now make rlong2 fall into the other system than rlong (unless rlong is between 0 and 180) */
            double rlong2 = rlong;
            if (rlong < 0.0)    //       if(rlong.lt.0) rlong2 = rlong + 360
                rlong2 = rlong + 360;

            if (rlong > 180.0) //       if(rlong.gt.180) rlong2 = rlong - 360
                rlong2 = rlong - 360;

            for (int ip = 1; ip <= polinf.npoly.Length; ip++) //       do 5 ip=1,np
            {
                if (!polinf.npoly[ip - 1].use) //       if(.not.use(ip)) then
                {
                    /* polygon is one we should not be in; test to see if we are, and if so set
                     * iok to 0 and return - check all such before seeing if we might also be in a polygon we should be */
                    for (int i = 1; i <= polinf.npoly[ip - 1].npoly; i++)   //          do 3 i=1,npoly(ip)
                    {
                        xu[i - 1] = polinf.npoly[ip - 1].xpoly[i - 1]; //             xu(i) = xpoly(i,ip)
                        yu[i - 1] = polinf.npoly[ip - 1].ypoly[i - 1]; //             yu(i) = ypoly(i,ip)
                    }
                    // 	 if(i
                    if (inpoly(ref rlong, ref rlat, ref xu, ref yu, ref polinf.npoly[ip - 1].npoly) != 0 ||
                        inpoly(ref rlong2, ref rlat, ref xu, ref yu, ref polinf.npoly[ip - 1].npoly) != 0)
                    {
                        iok = true; // 	    iok=0
                        return iok;   // 	    return
                    }
                }
            }

            int ianyok = 0;     //       ianyok=0

            for (int ip = 1; ip <= polinf.npoly.Length; ++ip)   //       do 15 ip=1,np
            {
                if (polinf.npoly[ip - 1].use)       //       if(use(ip)) then
                {
                    /* polygon is one we should be in; test to see if we are, and if so set iok to 1 and return */
                    ianyok++;       //          ianyok = ianyok+1

                    for (int i = 1; i <= polinf.npoly[ip - 1].npoly; i++) //          do 13 i=1,npoly(ip)
                    {
                        xu[i - 1] = polinf.npoly[ip - 1].xpoly[i - 1];     //             xu(i) = xpoly(i,ip)
                        yu[i - 1] = polinf.npoly[ip - 1].ypoly[i - 1];     //             yu(i) = ypoly(i,ip)
                    }
                    //       TODO hier fhelt was
                    if (inpoly(ref rlong, ref rlat, ref xu, ref yu, ref  polinf.npoly[ip - 1].npoly) != 0 ||
                        inpoly(ref rlong2, ref rlat, ref xu, ref yu, ref  polinf.npoly[ip - 1].npoly) != 0)
                    {
                        iok = false;    // 	    iok=1
                        return iok;     // 	    return
                    }
                }
            }
            /* not inside any polygons; set iok to 0 if there are any we should have been in */
            iok = false;    //       iok = 1

            if (ianyok > 0) //       if(ianyok.gt.0) iok = 0
                iok = true;

            return iok; //       return
        }

        /*  Original date of RCS archived version is  May 13  1989 */


        /*     $Id: inpoly.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $ */

        /*     $Log: inpoly.f,v $ */
        /*     Revision 1.1  2011/09/23 22:39:01  agnew */
        /*     Initial revision */


        /* ------------------------------------- */

        private int inpoly(ref double x, ref double y, ref double[] xv, ref double[] yv, ref int nv) // function inpoly(x,y,xv,yv,nv)
        {
            int inpoly = 0;

            /*   Returns 1 if point at (x,y) is inside polygon whose nv vertices are at xv(1),yv(1);...,xv(nv),yv(nv)
             *   Returns 0 if point is outside
             *   Returns 2 if point is on edge or vertex */

            /*     Rewritten by D. Agnew from the version by Godkin and Pulli, in BSSA, Vol 74, pp 1847-1848 (1984) */

            int isc = 0;
            for (int i = 0; i < nv; i++)   //       do 5 i=1,nv-1
            {
                isc = ncross(xv[i] - x, yv[i] - y, xv[i + 1] - x, yv[i + 1] - y);   //       isc = ncross(xv(i)-x,yv(i)-y,xv(i+1)-x,yv(i+1)-y)
                //       if(isc.eq.4) then
                if (isc == 4)
                {
                    /*  on edge - know the answer */
                    inpoly = 2;                                                     //          inpoly = 2
                    return inpoly;                                                  //          return
                }
                inpoly += isc;                                                      //       inpoly = inpoly + isc
            }

            /*  check final segment */
            isc = ncross(xv[nv - 1] - x, yv[nv - 1] - y, xv[0] - x, yv[0] - y);    //       isc = ncross(xv(nv)-x,yv(nv)-y,xv(1)-x,yv(1)-y)

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

        private int ncross(double x1, double y1, double x2, double y2)     // function ncross(x1,y1,x2,y2)
        {
            int ncross = 0;
            /*   finds whether the segment from point 1 to point 2 crosses the negative x-axis or goes through the origin (this is
             *   the signed crossing number)
             *     return value       nature of crossing
             *         4               segment goes through the origin
             *         2               segment crosses from below
             *         1               segment ends on -x axis from below or starts on it and goes up
             *         0               no crossing
             *        -1               segment ends on -x axis from above or starts on it and goes down
             *        -2               segment crosses from above */

            if (y1 * y2 > 0.0)                  //       if(y1*y2.gt.0) then
            {
                /* all above (or below) axis */
                ncross = 0;                     //          ncross = 0
                return ncross;                  //          return
            }

            double c12 = x1 * y2;               //       c12 = x1*y2
            double c21 = x2 * y1;               //       c21 = x2*y1

            if (c12 == c21 && x1 * x2 <= 0.0)   //       if(c12.eq.c21.and.x1*x2.le.0.) then
            {
                /* through origin */
                ncross = 4;                     //          ncross = 4
                return ncross;                  //          return
            }
            if (y1 == 0.0 && x1 > 0.0 ||
                y2 == 0.0 && x2 > 0.0 ||
                y1 < 0.0 && c12 > c21 ||
                y1 > 0.0 && c12 < c21 ||
                y1 == 0.0 && y2 == 0.0 && x1 < 0.0 && x2 < 0.0)
            // if((y1.eq.0.and.x1.gt.0).or.(y2.eq.0.and.x2.gt.0)
            //     .or.((y1.lt.0).and.(c12.gt.c21)).or.((y1.gt.0).and.(c12.lt.c21))
            //     .or.(y1.eq.0.and.y2.eq.0.and.x1.lt.0.and.x2.lt.0)) then
            {
                /* touches +x axis; crosses +x axis; lies entirely on -x axis */
                ncross = 0;                     //          ncross = 0
                return ncross;                  //          return
            }

            if (y1 != 0.0 && y2 != 0.0)         //       if(y1.ne.0.and.y2.ne.0) then
            {
                /*  cross axis */
                if (y1 < 0.0)                   //          if(y1.lt.0) ncross = 2
                    ncross = 2;
                if (y1 > 0.0)                   //          if(y1.gt.0) ncross = -2
                    ncross = -2;
                return ncross;                  //          return
            }

            /* one end touches -x axis - goes which way? */
            if (y1 == 0.0)                      //       if(y1.eq.0) then
            {
                if (y2 < 0.0)                   //          if(y2.lt.0) ncross = -1
                    ncross = -1;
                if (y2 > 0.0)                   //          if(y2.gt.0) ncross = 1
                    ncross = 1;
            }
            else                                //       else
            {
                /* y2=0 - ends on x-axis */
                if (y1 < 0.0)                   //          if(y1.lt.0) ncross = 1
                    ncross = 1;
                if (y1 > 0.0)                   //          if(y1.gt.0) ncross = -1
                    ncross = -1;
            }
            return ncross;                      //       return
        }

        #endregion

        private Complex ocmodl1(double rlong, double rlato,
                          OceanModel.TMelement TMe, GreensFunction.GFelememt GFe, LandseaMatrix lsM,
                                ref double dist, ref double numnf, ref double clnf, ref double farnf)
        {
            Complex amp = Complex.Zero; // return value


            //struct {
            //    real xpoly[5000]	/* was [500][10] */, ypoly[5000]	/* was [500][
            //        10] */;
            //    integer np, npoly[10];
            //    char polyf[80], polynm[800];
            //    logical1 use[10], ispoly;
            //} polinf_;

            //#define polinf_1 polinf_

            //struct {
            //    real dist, close, clat, clong;
            //    integer numnf;
            //    real clnf, farnf;
            //} coninf_;

            //#define coninf_1 coninf_

            //struct {
            //    real rhol;
            //} ldens_;

            //#define ldens_1 ldens_


            /*     Revision 1.6  2012/06/25 01:29:40  agnew added density common block, corrected comments
             *     Revision 1.5  2012/03/13 21:17:32  agnew Modified polygon common block
             *     Revision 1.4  2012/03/06 03:38:09  agnew added new grid options, density, and call to seadens
             *     Revision 1.3  2011/11/27 04:31:57  agnew larger dimension for polygon information
             *     Revision 1.2  2011/11/18 16:52:49  agnew corrected error in interpolation (incorrect indexing of surrounding cells in one case) 
             *     Revision 1.1  2011/09/23 22:39:02  agnew
             *     Initial revision */

            //       subroutine ocmodl(rlato,rlong,fingrd,amp)
            /* Subroutine: int ocmodl_(real *rlato, real *rlong, char *fingrd, complex *
                amp, ftnlen fingrd_len)
            {
                        /* Initialized data */
            int newf = 1;

            /* Local variables */

            SeaDensity seaDensity = new SeaDensity(); // TODO muss dann aber wo anders hin

            /*  For given latitude and east longitude (in degrees and ge 0) returns the complex amplitude of the tide (in kg/m**2, with greenwich phase).
             *  this amplitude is of course zero on land. If fingrd is "F", then the land-sea grid used is a detailed one;
             *  otherwise it is the one used by the ocean model.
             *  
             *  On the first call, this calls a separate routine (ocstart) which gets auxiliary information (passed to other routines in common) and 
             *  returns the two arrays of real and imaginary amplitude.  Later calls simply go to the appropriate point of the arrays.
             *
             *  The first call also reads in, if requested, information about polygonal areas (up to 10) within which the point either must or
             *  must not fall. */

            //       save ir1,im1,newf
            //       integer*2 ir1,im1
            //       logical*1 ispoly,use
            //       character*80 polyf,polynm

            /*  dimensioning of arrays below must be large enough to handle the largest */
            /*  models */

            //       dimension ir1(5000000),im1(5000000)
            //       common /modlim/tlat,blat,wlong,elong,latc,longc
            /*  common block holds polygon file information */
            //       
            //       common/coninf/dist,close,clat,clong,numnf,clnf,farnf
            //       common/ldens/rhol

            double x = 0;
            double y = 0;
            int i = 0;
            int j = 0;
            int ind = 0;

            if (newf == 1)          //       if(newf.eq.1) then
            {
                /* new files - load arrays (polygon information loaded into common) */
                //          newf = 0
                newf = 0;
                // 	 call ocstart(ir1,im1)
               // ocstart(ir1, im1);
                // 	 if(ispoly) call plstart
                if (polinf != null)
                {
                    //TODO plstart_();
                }
                //       endif
            }
            /* ind is cell number; x and y are positions within cell, relative to center, in E and N, in fractions of a cell dimension */
            celfnd(rlong, rlato, TMe,
                ref x, ref y, ref  i, ref j, ref ind);              //       call celfnd(rlong,rlato,x,y,i,j,ind)

            if (ind > TMe.latc * TMe.longc || ind <= 0)             //       if(ind.gt.latc*longc.or.ind.le.0) then
            {
                amp = Complex.Zero;                                 //          amp = cz
                return amp;                                         //          return
            }
            /* if cell is inside a polygon we don't want, or outside one we do, reject */
            /* it without going further */

            if (polinf != null)                                    //       if(ispoly) then
            {
                bool iok = chkgon(rlong, rlato);                         // 	 call chkgon(rlong,rlato,iok)
                if (iok)// == 0)                                       // 	 if(iok.eq.0) then
                {
                    amp = Complex.Zero;                             //             amp = cz
                    return amp;                                     //          return
                }
            }

            double rho = 0; // TODO ??
            if (GFe.fingrd == "C" || GFe.fingrd == "G")             //       if(fingrd.eq.'C'.or.fingrd.eq.'G') then
            {
                /*  large distance - return the model value (0 if on land for tidal models) */
                amp = 0.001d * (new Complex(
                                (double)TMe.ir1[ind - 1],
                                (double)TMe.im1[ind - 1]));         //            amp = .001*cmplx(float(ir1(ind)),float(im1(ind)))

                if (GFe.fingrd == "C")                              //          if(fingrd.eq.'C') call seadens(rlato,rlong,rho)
                    rho = seaDensity.seadens(rlato, rlong);

                if (GFe.fingrd == "G")                              //          if(fingrd.eq.'G') rho = 1.e3
                    rho = 1e3;

                amp *= rho;                                         //          amp=rho*amp
                seaDensity.rhol = rho;                              //          rhol=rho
                
                return amp;                                         //          return
            }

            if (GFe.fingrd == "F" || GFe.fingrd == "L")             //       if(fingrd.eq.'F'.or.fingrd.eq.'L') then
            {
                /* close distance, or land load - use the 5-minute map file to decide if there is land or not. If there is
                 * not (for F) or is (for L), bilinearly interpolate from the model, filling out empty cells to make this possible (though if all
                 * the cells needed are zero, return this). */
                int lnd = 0;
                if (lsM.getLandSea(rlato, rlong))                   //          call lndsea(rlato,rlong,lnd)
                    lnd = 1; // Land
                else
                    lnd = 2; // Water

                if (lnd == 1 && GFe.fingrd == "F" || 
                    lnd == 2 && GFe.fingrd == "L")                  // if((lnd.eq.1.and.fingrd.eq.'F').or.(lnd.eq.2.and.fingrd.eq.'L')) then
                {
                    /*  detailed grid says land (for F) or sea (for L) - return zero */
                    amp = Complex.Zero;                             //             amp = cz
                    return amp;                                     //             return
                }

                /* for L option, do not try to interpolate cells, density is 1.e3 */
                if (GFe.fingrd == "L")                              //          if(fingrd.eq.'L')then
                {
                    amp = 0.001d * (new Complex(
                                (double)TMe.ir1[ind - 1],
                                (double)TMe.im1[ind - 1]));         //            amp = .001*cmplx(float(ir1(ind)),float(im1(ind)))
                    amp *= 1e3;                                     //            amp = 1.e3*amp
                    return amp;                                     // return
                }

                /*  now left with option F, and detailed map says ocean load an interpolating array of four elements with values from adjacent
                 *  cells, depending on which quadrant of the cell we are in. In the process note how many nonzero values we have; save one */
                int ix = 1;                                         //   ix = 1
                if (x <= 0) ix = -ix;                               // 	 if(x.le.0) ix = -ix
                int iy = -TMe.longc;                                //   iy = -longc
                if (y <= 0) iy = -iy;                               // 	 if(y.le.0) iy = -iy
                x = Math.Abs(x);                                    // 	 x=abs(x)
                y = Math.Abs(y);                                    // 	 y=abs(y)

                int nnz = 0;
                int ibind = 0;
                Complex[] bg = new Complex[4];
                Complex bgs = new Complex();

                for (int ii = 0; ii < 4; ii++)                      // 	 do i=1,4
                {
                    if (ii == 0) ibind = ind;                       // 	   if(i.eq.1) ibind=ind
                    if (ii == 1) ibind = ind + ix;                  // 	   if(i.eq.2) ibind=ind+ix
                    if (ii == 2) ibind = ind + ix + iy;             // 	   if(i.eq.3) ibind=ind+ix+iy
                    if (ii == 3) ibind = ind + iy;                  // 	   if(i.eq.4) ibind=ind+iy

                    if (ibind <= 0 || 
                        ibind > TMe.latc * TMe.longc)               //     if(ibind.le.0.or.ibind.gt.latc*longc) 
                        bg[ii] = Complex.Zero;                      //         bg(i) = cz
                    if (ibind > 0 && 
                        ibind <= TMe.latc * TMe.longc)              //      if(ibind.gt.0.and.ibind.le.latc*longc)
                    {
                        bg[ii] = 0.001d * (new Complex(
                            (double)TMe.ir1[ibind - 1],
                            (double)TMe.im1[ibind - 1]));           //         bg(i) = 0.001*cmplx(float(ir1(ibind)),float(im1(ibind)))
                    }
                    if (!Complex.Equals(bg[ii], Complex.Zero))
                    {
                        nnz++;                                      //  if(bg(i).ne.cz) nnz = nnz+1
                        bgs = bg[ii];                               //  if(bg(i).ne.cz) bgs = bg(i)
                    }
                }

                /*  now we have a grid for interpolation: the values are arranged as 
                 *           4      3
                 *           1      2
                 *  where (1) is the value we are closest to.
                 *  
                 * we now fill out the zero values (if any) from the nonzero ones (Except that if there are 0 or 4 nonzero values we needn't bother) */
                if (nnz == 1)                                       //          if(nnz.eq.1) then
                {
                    /*  the one saved value is the nonzero one: fill out the whole array with this */
                    for (int ii = 0; ii < 4; ii++)                  //             do i=1,4
                        bg[ii] = bgs;                               //               bg(i) = bgs
                }

                if (nnz == 2)                                       //          if(nnz.eq.2) then
                {
                    /*  two values: either in a row, a column, or diagonal: */

                    if (!Complex.Equals(bg[0], Complex.Zero) &&
                        !Complex.Equals(bg[1], Complex.Zero))       //             if(bg(1).ne.cz.and.bg(2).ne.cz) then
                    {
                        bg[3] = bg[0];                              //                bg(4) = bg(1)
                        bg[2] = bg[1];                              //                bg(3) = bg(2)
                    }
                    else if (!Complex.Equals(bg[1], Complex.Zero) &&
                             !Complex.Equals(bg[2], Complex.Zero))  //             elseif(bg(2).ne.cz.and.bg(3).ne.cz) then
                    {
                        bg[3] = bg[2];                              //                bg(4) = bg(3)
                        bg[0] = bg[1];                              //                bg(1) = bg(2)
                    }
                    else if (!Complex.Equals(bg[2], Complex.Zero) &&
                             !Complex.Equals(bg[3], Complex.Zero))  //             elseif(bg(3).ne.cz.and.bg(4).ne.cz) then
                    {
                        bg[0] = bg[3];                              //                bg(1) = bg(4)
                        bg[1] = bg[2];                              //                bg(2) = bg(3)
                    }
                    else if (!Complex.Equals(bg[3], Complex.Zero) &&
                             !Complex.Equals(bg[0], Complex.Zero))  //             elseif(bg(4).ne.cz.and.bg(1).ne.cz) then
                    {
                        bg[2] = bg[3];                              //                bg(3) = bg(4)
                        bg[1] = bg[0];                              //                bg(2) = bg(1)
                    }
                    /*     diagonals */
                    else if (!Complex.Equals(bg[3], Complex.Zero) &&
                             !Complex.Equals(bg[1], Complex.Zero))  //             elseif(bg(4).ne.cz.and.bg(2).ne.cz) then
                    {
                        bg[2] = 0.5d * (bg[3] + bg[1]);             //                bg(3) = 0.5*(bg(4)+bg(2))
                        bg[0] = bg[2];                              //                bg(1) = bg(3)
                    }
                    else if (!Complex.Equals(bg[0], Complex.Zero) &&
                             !Complex.Equals(bg[2], Complex.Zero))  //             elseif(bg(1).ne.cz.and.bg(3).ne.cz) then
                    {
                        bg[3] = 0.5d * (bg[0] + bg[2]);             //                bg(4) = 0.5*(bg(1)+bg(3))
                        bg[1] = bg[2];                              //                bg(2) = bg(3)                                 
                    }
                }

                if (nnz == 3)                                       // 	 if(nnz.eq.3) then
                {
                    /* one missing value */
                    if (Complex.Equals(bg[0], Complex.Zero)) 
                        bg[0] = bg[3] + bg[1] - bg[2];              // 	    if(bg(1).eq.cz) bg(1) = bg(4)+bg(2)-bg(3)
                    if (Complex.Equals(bg[1], Complex.Zero)) 
                        bg[1] = bg[0] + bg[2] - bg[3];              // 	    if(bg(2).eq.cz) bg(2) = bg(1)+bg(3)-bg(4)
                    if (Complex.Equals(bg[2], Complex.Zero))
                        bg[2] = bg[3] + bg[1] - bg[0];              // 	    if(bg(3).eq.cz) bg(3) = bg(4)+bg(2)-bg(1)
                    if (Complex.Equals(bg[3], Complex.Zero)) 
                        bg[3] = bg[0] + bg[2] - bg[1];              // 	    if(bg(4).eq.cz) bg(4) = bg(1)+bg(3)-bg(2)
                }
                /*  bilinear interpolation */
                amp = bg[0] * (1 - x) * (1 - y) + bg[1] * x * (1 - y) +
                      bg[3] * (1 - x) * y + bg[2] * x * y;          //  amp = bg(1)*(1-x)*(1-y) + bg(2)*x*(1-y) + bg(4)*(1-x)*y + bg(3)*x*y

                /*  if all the values are zero, we will return zero - but we pass some information about this through common block coninf; the first time this
                 *  happens, set clnf and farnf to the current distance from the station, and afterwards only farnf (since we always increase the distance). */
                if (nnz == 0) // 	 if(nnz.eq.0) then
                {
                    if (numnf == 0)                                 // 	    if(numnf.eq.0) then
                    {
                        clnf = dist;                                // 	       clnf = dist
                        farnf = dist;                               // 	       farnf = dist
                    }
                    else                                            // 	    else
                    {
                        farnf = dist;                               // 	       farnf = dist
                    }
                    numnf++;                                        // 	    numnf=numnf+1
                }

                /*  scale the amplitude by the sea-water density */
                rho = seaDensity.seadens(rlato, rlong);             //       call seadens(rlato,rlong,rho)
                amp *= rho;                                         //       amp=rho*amp
                seaDensity.rhol = rho;                              //       rhol=rho
                return amp;                                         //       return
            }
            return amp;                                             //       return
        }

        private void celfnd(double rlong, double rlato, OceanModel.TMelement TMe,
                    ref double x, ref double y, ref int i, ref int j, ref int ind)
        {
            //              subroutine celfnd(rlong,rlato,x,y,i,j,ind)
            //c  returns the cell index (as i and j also ind), and
            //c  fractional position within a cell, for a given long
            //c  and lat (in degrees)
            //c
            //      common/modlim/tlat,blat,wlong,elong,latc,longc
            //c  check to see that we are within model
            if (rlato > TMe.tlat || rlato < TMe.blat)
            {
                ind = 0;
                return;
            }

            double dlong = (rlong - TMe.wlong) % 360d; // amod(rlong-wlong,360.)
            if (dlong < 0d) dlong += 360; //      if(dlong.lt.0) dlong = dlong +60
            if (dlong > (TMe.elong - TMe.wlong))
            {
                ind = 0;
                return;
            }
            //c  put rlong into the range of longitudes that are specified.

            x = (TMe.longc * dlong) / (TMe.elong - TMe.wlong);
            y = (TMe.latc * (TMe.tlat - rlato)) / (TMe.tlat - TMe.blat);
            i = (int)x + 1;  //verändert wegen 0 basiertem array         nee
            j = (int)y;// -1;     // verändert wegen 0 basiertem array             neee
            x = x - i + 0.5d;
            y = j - y + 0.5d;
            ind = i + TMe.longc * j;
        }

        private void ignewt(double del, double stp,
                            double stationLongitude, double stationLatitude, double stationHeight,
                            ref double grav, ref double tilt, ref double pot)
        {
            //  returns the integrated green function for newtonian gravity, tilt, and potential.
            //  the function is the integral over the interval centered at del with 
            //     width stp (ie, the interval [del-stp/2,del+stp/2], del and stp both being in radians
            //  the height correction is included in the green functions,
            //    the station height in meters being passed as ht in the common block stloc

            const double gn = 6.67e-11;
            const double an = 6.371e+6;

            //  height correction, and find constant parts
            double eps = stationHeight / an;
            double eps1 = 1d + eps;
            double eps2 = eps * eps;
            double eps3 = 2 * eps1 * eps1;
            double g2 = 2 * gn / (1d + 1.5d * eps);
            double ct = Math.Cos(Utilities.Constants.degree2radian * (90d - stationLatitude));
            double g = 9.7803327d * (1d + 0.005279d * ct * ct) - 3.08e-6 * stationHeight;
            double em = gn / g;
            double plc = 4 * an * em;

            double c1 = Math.Sin(0.25d * stp);
            double sthl = 0.5d * stp;

            // sign for gravity makes + acceleration up
            // use different approximations if < or > 0.05 rad (2.9 deg)
            if (del >= 0.05d)
            {
                grav = -g2 * Math.Cos(del / 2d) * c1 * (1d + eps / (Math.Cos(sthl) - Math.Cos(del)));
                tilt = em * (-2d * Math.Sin(del / 2d) * c1 +
                    Math.Log(Math.Tan((del + sthl) / 4d) / Math.Tan((del - sthl) / 4d)));
            }
            else
            { //      if(del.lt.0.05) then
                double d1 = del + stp / 2d;
                double d2 = del - stp / 2d;
                grav = -gn * (eps1 * d1 * d1 - 2 * eps) / (eps3 * Math.Sqrt(eps1 * d1 * d1 + eps2));
                grav += gn * (eps1 * d2 * d2 - 2 * eps) / (eps3 * Math.Sqrt(eps1 * d2 * d2 + eps2));
                tilt = em * (d2 / Math.Sqrt(eps2 + d2 * d2) - d1 / Math.Sqrt(eps2 + d1 * d1) +
                  Math.Log((d1 + Math.Sqrt(eps2 + d1 * d1)) / (d2 + Math.Sqrt(eps2 + d2 * d2))));
            }
            pot = plc * Math.Cos(del / 2d) * c1;
        }

        private Complex phasor(Complex c)
        {
            Complex tmp = Complex.Zero;
            double absolut = (double)Math.Sqrt(c.Real * c.Real + c.Imaginary * c.Imaginary);
            if (!Complex.Equals(c, Complex.Zero))
                tmp = new Complex(absolut, Utilities.Constants.radian2degree * Math.Atan2(c.Imaginary, c.Real));
            return tmp;
        }
        private Complex argand(Complex c)
        {
            double p = c.Imaginary / Utilities.Constants.radian2degree;
            double a = c.Real;
            return new Complex(a * Math.Cos(p), a * Math.Sin(p));
        }
        #endregion
        #endregion

        #region calculates the value for oceanloading of the respective component 
        private double hartid(DateTime dt, hartid_wPM[] _hartidPM)
        {
            // convert times to Julian days (UT) then to Julian centuries from J2000.00 (ET)
            double jd = juldat(dt);
            double dayfr = dt.Hour / 24d + dt.Minute / 1440d + dt.Second / 84600d;
            double year = dt.Year + (dt.DayOfYear + dayfr) / (365 + leap(dt.Year));
            double delta = etutc(year);// d);
            double djd = jd - 0.5d + dayfr;
            double time = (djd - 2451545d + delta / 86400d) / 36525d;

            //
            double[] a = new double[141];
            double[] f = new double[141];
            double[] p = new double[141];
            double[] hc = new double[282];
            double[] scr = new double[423];
            double[,] idt = new double[6, _hartidPM.Length];
            double[] amp = new double[_hartidPM.Length];
            double[] phase = new double[_hartidPM.Length];
            int nt = 141;

            //   get harmonics information
            for (int i = 0; i < _hartidPM.Length; i++)
            {
                hartid_wPM wptt = _hartidPM[i];

                for (int j = 0; j < 6; j++)
                    idt[j, i] = wptt.icte[j];

                amp[i] = wptt.amp;
                phase[i] = wptt.phase;
            }

            //  interpolate to larger set of harmonics
            admint(amp, idt, phase, a, ref f, ref p, _hartidPM.Length, ref nt, time, dayfr);

            for (int i = 1; i <= nt; i++)
            {
                //  set up for first recursion, and normalize frequencies
                p[i - 1] = Utilities.Constants.degree2radian * p[i - 1];
                f[i - 1] = Math.PI * f[i - 1] / 43200.0d;
                //c set up harmonic coefficients, compute tide, and write out
                hc[2 * i - 2] = a[i - 1] * Math.Cos(p[i - 1]);
                hc[2 * i - 1] = -a[i - 1] * Math.Sin(p[i - 1]);
            }

            return recurs(hc, nt, f, scr);
        }

        private double hartid(DateTime dt, ArrayList _hartidPM)
        {
            // convert times to Julian days (UT) then to Julian centuries from J2000.00 (ET)
            double jd = juldat(dt);
            double dayfr = dt.Hour / 24d + dt.Minute / 1440d + dt.Second / 84600d;
            double year = dt.Year + (dt.DayOfYear + dayfr) / (365 + leap(dt.Year));
            double delta = etutc(year);// d);
            double djd = jd - 0.5d + dayfr;
            double time = (djd - 2451545d + delta / 86400d) / 36525d;

            //
            double[] a = new double[141];
            double[] f = new double[141];
            double[] p = new double[141];
            double[] hc = new double[282];
            double[] scr = new double[423];
            double[,] idt = new double[6, _hartidPM.Count];
            double[] amp = new double[_hartidPM.Count];
            double[] phase = new double[_hartidPM.Count];
            int nt = 141;

            //   get harmonics information
            for (int i = 0; i < _hartidPM.Count; i++)
            {
                hartid_wPM wptt = (hartid_wPM)_hartidPM[i];

                for (int j = 0; j < 6; j++)
                    idt[j, i] = wptt.icte[j];

                amp[i] = wptt.amp;
                phase[i] = wptt.phase;
            }

            //  interpolate to larger set of harmonics
            admint(amp, idt, phase, a, ref f, ref p, _hartidPM.Count, ref nt, time, dayfr);

            for (int i = 1; i <= nt; i++)
            {
                //  set up for first recursion, and normalize frequencies
                p[i - 1] = Utilities.Constants.degree2radian * p[i - 1];
                f[i - 1] = Math.PI * f[i - 1] / 43200.0d;
                //c set up harmonic coefficients, compute tide, and write out
                hc[2 * i - 2] = a[i - 1] * Math.Cos(p[i - 1]);
                hc[2 * i - 1] = -a[i - 1] * Math.Sin(p[i - 1]);
            }

            return recurs(hc, nt, f, scr);
        }

        #region function für hartid
        private int leap(int iy)
        {//  function leap(iy)
            // returns 1 if year is a leap year, by Clavian (Gregorian) rule

            int leap = 1 - ((iy % 4) + 3) / 4;
            if ((iy % 100) == 0 && (iy % 400) != 0)
                leap = 0;
            return leap;
        }
        private int juldat(DateTime dt)
        {//     function juldat(it)
            //  Julian Date from Gregorian date, Algorithm from p. 604, Explanatory
            //    Supplement Amer Ephemeris & Nautical Almanac (cf Comm CACM, 11, 657 (1968)
            //    and 15, 918 (1972)) Valid for all positive values of Julian Date
            return ((1461 * (dt.Year + 4800 + (dt.Month - 14) / 12)) / 4
               + (367 * (dt.Month - 2 - 12 * ((dt.Month - 14) / 12))) / 12
               - (3 * ((dt.Year + 4900 + (dt.Month - 14) / 12) / 100)) / 4 + dt.Day - 32075);
        }
        private double etutc(double year)
        { //subroutine etutc(year,delta)
            double delta = 0d;
            //   input year (and fraction thereof), 1700-2006
            //   output delta = et - utc in seconds
            //   tables from p 90, Explanatory Supplement to the A.E., Dr. R A Broucke,
            //   JPL, and the leap.sec table in GAMIT
            //   utc (and ut) is time usually used (eg in time signals)
            double[] d = new double[142]
                {
                    5.15d,4.64d,5.36d,3.49d,3.27d,2.45d,4.03d,1.76d,3.3d,1.0d,2.42d,.94d,2.31d,
                    2.27d,-.22d,.03d,-.05d,-.06d,-.57d,.03d,-.47d,.98d,-.86d,2.45d,.22d,.37d,2.79d,
                    1.2d,3.52d,1.17d,2.67d,3.06d,2.66d,2.97d,3.28d,3.31d,3.33d,3.23d,3.6d,3.52d,
                    4.27d,2.68d,2.75d,2.67d,1.94d,1.39d,1.66d,.88d,.33d,-.17d,-1.88d,-3.43d,-4.05d,
                    -5.77d,-7.06d,-7.36d,-7.67d,-7.64d,-7.93d,-7.82d,-8.35d,-7.91d,-8.03d,-9.14d,
                    -8.18d,-7.88d,-7.62d,-7.17d,-8.14d,-7.59d,-7.17d,-7.94d,-8.23d,-7.88d,-7.68d,
                    -6.94d,-6.89d,-7.11d,-5.87d,-5.04d,-3.9d,-2.87d,-.58d,.71d,1.8d,3.08d,4.63d,
                    5.86d,7.21d,8.58d,10.5d,12.1d,12.49d,14.41d,15.59d,15.81d,17.52d,19.01d,18.39d,
                    19.55d,20.36d,21.01d,21.81d,21.76d,22.35d,22.68d,22.94d,22.93d,22.69d,22.94d,
                    23.2d,23.31d,23.63d,23.47d,23.68d,23.62d,23.53d,23.59d,23.99d,23.8d,24.2d,
                    24.99d,24.97d,25.72d,26.21d,26.37d,26.89d,27.68d,28.13d,28.94d,29.42d,29.66d,
                    30.29d,30.96d,31.09d,31.59d,32.06d,31.82d,32.69d,33.05d,33.16d,33.59
                };
            double[] tx = new double[39]
                {
                    61.5d   ,62d    , 62.5d   ,63d     ,63.5d   ,64d     ,64.5d   ,65d ,
                    65.5d   ,66d    , 66.5d   ,67d     ,67.5d   ,68d     ,68.25d  ,68.5d,
                    68.75d  ,69d    , 69.25d  ,69.5d   ,69.75d  ,70d     ,70.25d  ,70.5d,
                    70.75d  ,71d    , 71.085d ,71.162d ,71.247d ,71.329d ,71.414d ,71.496d,
                    71.581d ,71.666d, 71.748d ,71.833d ,71.915d ,71.999d ,72.0d
                };
            double[] ty = new double[39]
                {
                    33.59d  ,34.032d ,34.235d ,34.441d ,34.644d ,34.95d  ,35.286d ,35.725d,
                    36.16d  ,36.498d ,36.968d ,37.444d ,37.913d ,38.39d  ,38.526d ,38.76d ,
                    39d     ,39.238d ,39.472d ,39.707d ,39.946d ,40.185d ,40.42d  ,40.654d,
                    40.892d ,41.131d ,41.211d ,41.284d ,41.364d ,41.442d ,41.522d ,41.600d,
                    41.680d  ,41.761d,41.838d ,41.919d ,41.996d ,42.184d ,42.184d
                };
            //      dimension d(142),tx(39),ty(39),st(23),si(23)

            //c   update st,si, and nstep to allow for added leap seconds
            double[] st = new double[23]
                {
                    1972.5d,1973.0d,1974.0d,1975.0d,1976.0d,1977.0d,1978.0d,1979.0d,
                    1980.0d,1981.5d,1982.5d,1983.5d,1985.5d,1988.0d,1990.0d,1991.0d,
                    1992.5d,1993.5d,1994.5d,1996.0d,1997.5d,1999.0d,2006.0
                };
            double[] si = new double[23]
                {
                    1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d, 1d,
                    1d, 1d, 1d, 1d, 1d, 1d, 1d 
                };
            int nstep = 23;

            //   for 1820.5 to 1961.5, data is spaced at yearly intervals
            if (year > 1972.0d)// cas = 6;
            {
                delta = 42.184;
                for (int i = 0; i < nstep; i++)
                {
                    if (year >= st[i])
                        delta = delta + si[i];
                    if (year < st[i])
                        return delta;
                }
            }
            else if (year > 1961.5d)// cas = 2;
            {
                for (int i = 0; i < 38; i++)
                {  // 2    do 3 i = 1,38
                    if (year - 1900.0d == tx[i])
                    {
                        delta = ty[i];
                        return delta;
                    }
                    if (year - 1900.0d < tx[i])
                    {
                        delta = ty[i - 1] + (ty[i] - ty[i - 1]) * ((year - 1900d - tx[i - 1]) / (tx[i] - tx[i - 1]));
                        return delta;
                    }
                }
            }
            else if (year > 1820.5d)// cas = 1;
            {
                double n = year - 1819.5d - 1;
                double frac = year - (1819.5d + n);
                delta = (d[(int)(n + 1)] - d[(int)n]) * frac + d[(int)n];
                return delta;
            }
            else
            {
                //  for oldest epochs, use approximation
                if (year >= 1785.0d)
                    return delta = 6;
                if (year < 1785.0d)
                    return delta = (year - 1750.0d) / 5.0d;
                if (year < 1700.0d)
                    return delta = 0d;
            }
            return delta;
        }

        private void admint(double[] ampin, double[,] idtin, double[] phin,
                            double[] amp, ref double[] f, ref double[] p, int nin, ref int nout, double time, double dayfr)
        {
            //   subroutine admint(ampin,idtin,phin,amp,f,p,nin,nout)
            //c  Returns the amplitude, frequency, and phase of a set of tidal
            //c  constituents. n is input as the number wanted, and returned as the number
            //c  provided.  The constituents used are stored in the arrays idd (doodson
            //c  number) and tamp (Cartwright-Edden amplitude).  The actual amp and
            //c  phase of each of these are determined by spline interpolation of the
            //c  real and imaginary part of the admittance, as specified at a subset
            //c  of the constituents.
            //c  The phase is determined for a time set in common block /date/ (see
            //c  subroutine tdfrph), outside of this subroutine.
            //c   Note that the arrays f and p must be specified as double precision.

            double fr = double.NaN;
            double pr = double.NaN;

            double[] rl = new double[nin];
            double[] aim = new double[nin];
            double[] rf = new double[nin];
            int[] key = new int[nin];
            double[] scr = new double[nin];
            double[] zdi = new double[nin];
            double[] zdr = new double[nin];
            double[] di = new double[nin];
            double[] dr = new double[nin];
            double[] sdi = new double[nin];
            double[] sdr = new double[nin];
            double[] tdi = new double[nin];
            double[] tdr = new double[nin];

            const int nt = 141;

            double[] tamp = new double[141]
                {
                  0.63221d,  0.29411d,  0.12105d,  0.07991d,  0.02382d, -0.02359d,  0.02299d,  0.01933d,
                 -0.01787d,  0.01719d,  0.01602d,  0.00467d, -0.00466d, -0.00452d,  0.00447d,  0.00447d,
                  0.00259d, -0.00246d, -0.00217d,  0.00197d,  0.00195d,  0.00191d, -0.00190d,  0.00180d,
                  0.00130d,  0.00117d,  0.00113d,  0.00106d, -0.00102d, -0.00102d,  0.00101d,  0.00090d,
                 -0.00086d,  0.00085d,  0.00085d,  0.00077d,  0.00074d,  0.00074d, -0.00072d,  0.00070d,
                  0.00066d,  0.00065d, -0.00065d,  0.00063d,  0.00063d, -0.00060d,  0.00059d,  0.00054d,
                  0.00048d, -0.00046d,  0.00041d,  0.36864d, -0.26223d, -0.12200d, -0.05021d,  0.05003d,
                 -0.04947d,  0.02062d,  0.02061d,  0.01128d, -0.00953d, -0.00947d, -0.00801d,  0.00741d,
                 -0.00730d,  0.00723d, -0.00713d, -0.00664d,  0.00525d,  0.00414d,  0.00409d,  0.00394d,
                  0.00394d,  0.00342d,  0.00342d,  0.00289d,  0.00288d,  0.00216d, -0.00194d,  0.00193d,
                 -0.00180d,  0.00169d,  0.00169d,  0.00152d,  0.00151d, -0.00151d,  0.00138d,  0.00137d,
                  0.00137d, -0.00125d, -0.00107d,  0.00102d,  0.00090d,  0.00087d, -0.00079d,  0.00079d,
                  0.00078d, -0.00075d, -0.00075d,  0.00067d, -0.00060d, -0.00060d,  0.00054d,  0.00054d,
                 -0.00054d, -0.00047d, -0.00044d,  0.00044d,  0.00042d,  0.00041d, -0.06661d, -0.03518d,
                 -0.03099d,  0.02793d, -0.02762d, -0.01275d, -0.00673d, -0.00584d, -0.00529d, -0.00492d,
                 -0.00288d, -0.00258d, -0.00242d,  0.00231d,  0.00228d, -0.00204d,  0.00188d, -0.00181d,
                 -0.00169d, -0.00100d, -0.00093d, -0.00084d,  0.00077d,  0.00077d, -0.00070d, -0.00050d,
                 -0.00049d,  0.00049d,  0.00048d,  0.00044d, -0.00042 
                };

            int[,] idd = new int[141, 6] 
                {
                    { 2, 0, 0, 0, 0, 0}, { 2, 2,-2, 0, 0, 0}, { 2,-1, 0, 1, 0, 0}, { 2, 2, 0, 0, 0, 0}, 
                    { 2, 2, 0, 0, 1, 0}, { 2, 0, 0, 0,-1, 0}, { 2,-1, 2,-1, 0, 0}, { 2,-2, 2, 0, 0, 0},
                    { 2, 1, 0,-1, 0, 0}, { 2, 2,-3, 0, 0, 1}, { 2,-2, 0, 2, 0, 0}, { 2,-3, 2, 1, 0, 0},  
                    { 2, 1,-2, 1, 0, 0}, { 2,-1, 0, 1,-1, 0}, { 2, 3, 0,-1, 0, 0}, { 2, 1, 0, 1, 0, 0}, 
                    { 2, 2, 0, 0, 2, 0}, { 2, 2,-1, 0, 0,-1}, { 2, 0,-1, 0, 0, 1}, { 2, 1, 0, 1, 1, 0},
                    { 2, 3, 0,-1, 1, 0}, { 2, 0, 1, 0, 0,-1}, { 2, 0,-2, 2, 0, 0}, { 2,-3, 0, 3, 0, 0},  
                    { 2,-2, 3, 0, 0,-1}, { 2, 4, 0, 0, 0, 0}, { 2,-1, 1, 1, 0,-1}, { 2,-1, 3,-1, 0,-1},
                    { 2, 2, 0, 0,-1, 0}, { 2,-1,-1, 1, 0, 1}, { 2, 4, 0, 0, 1, 0}, { 2,-3, 4,-1, 0, 0},
                    { 2,-1, 2,-1,-1, 0}, { 2, 3,-2, 1, 0, 0}, { 2, 1, 2,-1, 0, 0}, { 2,-4, 2, 2, 0, 0},  
                    { 2, 4,-2, 0, 0, 0}, { 2, 0, 2, 0, 0, 0}, { 2,-2, 2, 0,-1, 0}, { 2, 2,-4, 0, 0, 2},
                    { 2, 2,-2, 0,-1, 0}, { 2, 1, 0,-1,-1, 0}, { 2,-1, 1, 0, 0, 0}, { 2, 2,-1, 0, 0, 1},
                    { 2, 2, 1, 0, 0,-1}, { 2,-2, 0, 2,-1, 0}, { 2,-2, 4,-2, 0, 0}, { 2, 2, 2, 0, 0, 0},  
                    { 2,-4, 4, 0, 0, 0}, { 2,-1, 0,-1,-2, 0}, { 2, 1, 2,-1, 1, 0}, { 1, 1, 0, 0, 0, 0},
                    { 1,-1, 0, 0, 0, 0}, { 1, 1,-2, 0, 0, 0}, { 1,-2, 0, 1, 0, 0}, { 1, 1, 0, 0, 1, 0},
                    { 1,-1, 0, 0,-1, 0}, { 1, 2, 0,-1, 0, 0}, { 1, 0, 0, 1, 0, 0}, { 1, 3, 0, 0, 0, 0},  
                    { 1,-2, 2,-1, 0, 0}, { 1,-2, 0, 1,-1, 0}, { 1,-3, 2, 0, 0, 0}, { 1, 0, 0,-1, 0, 0},
                    { 1, 1, 0, 0,-1, 0}, { 1, 3, 0, 0, 1, 0}, { 1, 1,-3, 0, 0, 1}, { 1,-3, 0, 2, 0, 0},
                    { 1, 1, 2, 0, 0, 0}, { 1, 0, 0, 1, 1, 0}, { 1, 2, 0,-1, 1, 0}, { 1, 0, 2,-1, 0, 0},  
                    { 1, 2,-2, 1, 0, 0}, { 1, 3,-2, 0, 0, 0}, { 1,-1, 2, 0, 0, 0}, { 1, 1, 1, 0, 0,-1},
                    { 1, 1,-1, 0, 0, 1}, { 1, 4, 0,-1, 0, 0}, { 1,-4, 2, 1, 0, 0}, { 1, 0,-2, 1, 0, 0},
                    { 1,-2, 2,-1,-1, 0}, { 1, 3, 0,-2, 0, 0}, { 1,-1, 0, 2, 0, 0}, { 1,-1, 0, 0,-2, 0},  
                    { 1, 3, 0, 0, 2, 0}, { 1,-3, 2, 0,-1, 0}, { 1, 4, 0,-1, 1, 0}, { 1, 0, 0,-1,-1, 0},
                    { 1, 1,-2, 0,-1, 0}, { 1,-3, 0, 2,-1, 0}, { 1, 1, 0, 0, 2, 0}, { 1, 1,-1, 0, 0,-1},
                    { 1,-1,-1, 0, 0, 1}, { 1, 0, 2,-1, 1, 0}, { 1,-1, 1, 0, 0,-1}, { 1,-1,-2, 2, 0, 0},  
                    { 1, 2,-2, 1, 1, 0}, { 1,-4, 0, 3, 0, 0}, { 1,-1, 2, 0, 1, 0}, { 1, 3,-2, 0, 1, 0},
                    { 1, 2, 0,-1,-1, 0}, { 1, 0, 0, 1,-1, 0}, { 1,-2, 2, 1, 0, 0}, { 1, 4,-2,-1, 0, 0},
                    { 1,-3, 3, 0, 0,-1}, { 1,-2, 1, 1, 0,-1}, { 1,-2, 3,-1, 0,-1}, { 1, 0,-2, 1,-1, 0},  
                    { 1,-2,-1, 1, 0, 1}, { 1, 4,-2, 1, 0, 0}, { 0, 2, 0, 0, 0, 0}, { 0, 1, 0,-1, 0, 0},
                    { 0, 0, 2, 0, 0, 0}, { 0, 0, 0, 0, 1, 0}, { 0, 2, 0, 0, 1, 0}, { 0, 3, 0,-1, 0, 0},
                    { 0, 1,-2, 1, 0, 0}, { 0, 2,-2, 0, 0, 0}, { 0, 3, 0,-1, 1, 0}, { 0, 0, 1, 0, 0,-1},  
                    { 0, 2, 0,-2, 0, 0}, { 0, 2, 0, 0, 2, 0}, { 0, 3,-2, 1, 0, 0}, { 0, 1, 0,-1,-1, 0},
                    { 0, 1, 0,-1, 1, 0}, { 0, 4,-2, 0, 0, 0}, { 0, 1, 0, 1, 0, 0}, { 0, 0, 3, 0, 0,-1},
                    { 0, 4, 0,-2, 0, 0}, { 0, 3,-2, 1, 1, 0}, { 0, 3,-2,-1, 0, 0}, { 0, 4,-2, 0, 1, 0},  
                    { 0, 0, 2, 0, 1, 0}, { 0, 1, 0, 1, 1, 0}, { 0, 4, 0,-2, 1, 0}, { 0, 3, 0,-1, 2, 0},
                    { 0, 5,-2,-1, 0, 0}, { 0, 1, 2,-1, 0, 0}, { 0, 1,-2, 1,-1, 0}, { 0, 1,-2, 1, 1, 0},
                    { 0, 2,-2, 0,-1, 0}
                };

            int k = 0;
            int nlp = 0;
            int ndi = 0;
            int nsd = 0;
            int ntd = 0;
            double[] d = new double[6];
            double[] dd = new double[6];
            int[] itmsave = new int[6];

            for (int ll = 0; ll < nin; ll++)
            {
                //  see if Doodson numbers match
                for (int kk = 0; kk < nt; kk++)
                {
                    int ii = 0;
                    for (int i = 0; i < 6; i++)
                        ii = ii + (int)Math.Abs(idd[kk, i] - idtin[i, ll]);

                    if (ii == 0 && k < nin /*ncon*/)
                    {
                        rl[k] = ampin[ll] * Math.Cos(Utilities.Constants.degree2radian * phin[ll]) / Math.Abs(tamp[kk]);
                        aim[k] = ampin[ll] * Math.Sin(Utilities.Constants.degree2radian * phin[ll]) / Math.Abs(tamp[kk]);
                        // Now have real and imaginary parts of admittance, scaled by C-E
                        // amplitude. Admittance phase is whatever was used in the original
                        // expression. (Usually phase is given relative to some reference
                        // but amplitude is in absolute units). Next get frequency.
                        double[] tt = new double[6] {idd[kk,0],idd[kk,1],idd[kk,2],
                                                         idd[kk,3],idd[kk,4],idd[kk,5]};
                        tdfrph(tt, ref fr, ref pr, time, dayfr, ref d, ref dd);
                        rf[k] = fr;
                        k++;
                    }
                }
            }

            //  done going through constituents--there are k of them
            //  have specified admittance at a number of points. sort these by frequency
            //  and separate diurnal and semidiurnal, recopying admittances to get them in order
            shells(ref rf, ref key, k);

            for (int i = 0; i < k; i++)
            {
                if (rf[i] < 0.5d) nlp++;
                if (rf[i] < 1.5d && rf[i] > 0.5d) ndi++;
                if (rf[i] < 2.5d && rf[i] > 1.5d) nsd++;
                if (rf[i] < 3.5d && rf[i] > 2.5d) ntd++;
                scr[i] = rl[(int)key[i]];
            }
            for (int i = 0; i < k; i++)
            {
                rl[i] = scr[i];
                scr[i] = aim[key[i]];
            }
            for (int i = 0; i < k; i++)
                aim[i] = scr[i];

            //  now set up splines (8 cases - four species, each real and imag)
            //  we have to allow for the case when there are no constituent amplitudes
            //  for the long-period or terdiurnal.

            if (nlp != 0)
            {
                zdr = spline(nlp, rf, rl, scr);
                zdi = spline(nlp, rf, aim, scr);
            }

            dr = spline(ndi, cpfrom(rf, nlp), cpfrom(rl, nlp), scr);
            di = spline(ndi, cpfrom(rf, nlp), cpfrom(aim, nlp), scr);
            sdr = spline(nsd, cpfrom(rf, nlp + ndi), cpfrom(rl, nlp + ndi), scr);
            sdi = spline(nsd, cpfrom(rf, nlp + ndi), cpfrom(aim, nlp + ndi), scr);
            if (ntd != 0)
            {
                tdr = spline(ntd, cpfrom(rf, nlp + ndi + nsd), cpfrom(rl, nlp + ndi + nsd), scr);
                tdi = spline(ntd, cpfrom(rf, nlp + ndi + nsd), cpfrom(aim, nlp + ndi + nsd), scr);
            }

            //  evaluate all harmonics using the interpolated admittance
            int j = 0;
            int istart = 0;
            for (int i = 0; i < nt; i++)
            {
                if (!((idd[i, 0] == 0 && nlp == 0) || (idd[i, 0] == 3 && ntd == 0)))
                {
                    double[] tt = new double[6] {idd[i,0],idd[i,1],idd[i,2],
                                                         idd[i,3],idd[i,4],idd[i,5]};
                    tdfrph(tt, ref f[j], ref p[j], time, dayfr, ref d, ref dd);
                    //  equilibrium phase corrections (including to local phase) (LP tides corrected 1 Aug 2009)
                    if (idd[i, 0] == 0) p[j] = p[j] + 180d;
                    if (idd[i, 0] == 1) p[j] = p[j] + 90d;
                    double sf = f[j];
                    double re = 0d;
                    double am = 0d;

                    if (idd[i, 0] == 0)
                    {
                        re = eval(sf, nlp, rf, rl, zdr, ref istart);
                        am = eval(sf, nlp, rf, aim, zdi, ref istart);
                    }
                    if (idd[i, 0] == 1)
                    {
                        re = eval(sf, ndi, cpfrom(rf, nlp), cpfrom(rl, nlp), dr, ref istart);
                        am = eval(sf, ndi, cpfrom(rf, nlp), cpfrom(aim, nlp), di, ref istart);
                    }
                    if (idd[i, 0] == 2)
                    {
                        re = eval(sf, nsd, cpfrom(rf, nlp + ndi), cpfrom(rl, nlp + ndi), sdr, ref istart);
                        am = eval(sf, nsd, cpfrom(rf, nlp + ndi), cpfrom(aim, nlp + ndi), sdi, ref istart);
                    }
                    if (idd[i, 0] == 3)
                    {
                        re = eval(sf, ntd, cpfrom(rf, nlp + ndi + nsd), cpfrom(rl, nlp + ndi + nsd), tdr, ref istart);
                        am = eval(sf, ntd, cpfrom(rf, nlp + ndi + nsd), cpfrom(aim, nlp + ndi + nsd), tdi, ref istart);
                    }
                    amp[j] = tamp[i] * Math.Sqrt(re * re + am * am);
                    p[j] = p[j] + Math.Atan2(am, re) / Utilities.Constants.degree2radian;
                    j++;
                }
            }
            nout = j - 1;
        }

        private void tdfrph(double[] idood, ref double freq, ref double phase, double time, double dayfr,
            ref double[] d, ref double[] dd)//, ref int[] itmsave)
        {
            //   Given the Doodson number of a tidal constituent (in idood), returns
            //  the frequency and phase.  Phase is returned in degrees and frequency
            //  in cycles/day.

            //   Note that phases must be decreased by 90 degrees if the sum of the order
            //  and the species number is odd (as for the 2nd degree diurnals, and 3rd
            //  degree low-frequency and semidiurnals).
            //   These phases may need further adjustment to allow for the spherical
            //  harmonic normalization used; e.g. for that used for the potential by
            //  Cartwright and Tayler, 180 degrees must be added for (species,order)
            //  = (1,2), (1,3), or (3,3). 

            // IERS expressions for the Delauney arguments, in degrees
            double f1 = 134.9634025100d + time * (477198.8675605000d + time * (0.0088553333d +
                                          time * (0.0000143431d + time * (-0.0000000680d))));
            double f2 = 357.5291091806d + time * (35999.0502911389d + time * (-0.0001536667d +
                                          time * (0.0000000378d + time * (-0.0000000032d))));
            double f3 = 93.2720906200d + time * (483202.0174577222d + time * (-0.0035420000d +
                                          time * (-0.0000002881d + time * (0.0000000012d))));
            double f4 = 297.8501954694d + time * (445267.1114469445d + time * (-0.0017696111d +
                                          time * (0.0000018314d + time * (-0.0000000088d))));
            double f5 = 125.0445550100d + time * (-1934.1362619722d + time * (0.0020756111d +
                                          time * (0.0000021394d + time * (-0.0000000165d))));

            //  convert to Doodson (Darwin) variables
            d[0] = 360.0d * dayfr - f4;
            d[1] = f3 + f5;
            d[2] = d[1] - f4;
            d[3] = d[1] - f1;
            d[4] = -f5;
            d[5] = d[2] - f2;

            //find frequencies of Delauney variables (in cycles/day), and from these
            //the same for the Doodson arguments

            double fd1 = 0.0362916471d + 0.0000000013d * time;
            double fd2 = 0.0027377786d;
            double fd3 = 0.0367481951d - 0.0000000005d * time;
            double fd4 = 0.0338631920d - 0.0000000003d * time;
            double fd5 = -0.0001470938d + 0.0000000003d * time;
            dd[0] = 1.0d - fd4;
            dd[1] = fd3 + fd5;
            dd[2] = dd[1] - fd4;
            dd[3] = dd[1] - fd1;
            dd[4] = -fd5;
            dd[5] = dd[2] - fd2;

            //compute phase and frequency of the given tidal constituent

            freq = 0d;
            phase = 0d;
            for (int i = 0; i < 6; i++)
            {//  do i = 1,6
                freq = freq + idood[i] * dd[i];
                phase = phase + idood[i] * d[i];
            }
            //adjust phases to fall in the positive range 0 to 360
            phase = phase % 360d;
            if (phase < 0)
                phase = phase + 360d;
        }

        private void shells(ref double[] x, ref int[] k, int n)
        {//                  subroutine shells(x,k,n)
            //  sorts an array x, of length n, sorting upward, and returns
            //  an array k which may be used to key another array to the
            //  sorted pattern (i.e., if we had an array f to which x 
            //  corresponded before sorting, then after calling shells,
            //  f(k(1)) will be the element of f corresponding to the
            //  smallest x, f(k(2)) the next smallest, and so on.
            //   revised 29-dec-82 so k is sorted in turn, the values of
            //  k that point to identical values of x being put in increasing order

            double[] tmp = new double[n];// x.Clone();
            for (int i = 0; i < n; i++)
                tmp[i] = x[i];

            //sortieren
            Array.Sort(tmp);

            //indexe suchen und in key eintragen
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    if (x[i] == tmp[j])
                        k[j] = i;
                }
            }
            // zurückschreiben
            for (int i = 0; i < n; i++)
                x[i] = tmp[i];
        }

        private double[] spline(int nn, double[] x, double[] u, double[] a)
        {
            //                      subroutine spline(nn,x,u,s,a)
            //c  finds array  s  for spline interpolator  eval.
            //c  nn  number of data points supplied (may be negative, see below)
            //c  x  array containing x-coordinates where function is sampled.  xx(1),xx(2),... must be a strictly increasing sequence.
            //c  u  array containing sample values that are to be interpolated.
            //c  s  output array of 2nd derivative at sample points.
            //c  a  working space array of dimension at least  nn.
            //c  if the user wishes to force the derivatives at the ends of the series to     
            //c  assume specified values, he should put du(1)/dx and du(n)/dx in s1,s2
            //c  and call the routine with nn=-number of terms in series.  normally a parabola
            //c  is fitted through the 1st and last 3 points to find the slopes.
            //c  if less than 4 points are supplied, straight lines are fitted.
            //      dimension x(*),u(*),s(*),a(*)
            //c
            //      q(u1,x1,u2,x2)=(u1/x1**2-u2/x2**2)/(1.0/x1-1.0/x2)
            //c

            double[] s = new double[Math.Abs(nn) + 1];

            int n = Math.Abs(nn) - 1;

            //  series too short for cubic spline - use straight lines.
            if (n < 3)
            {
                for (int i = 0; i < n; i++)
                    s[i] = 0.0d;
                return s;
            }

            double q1 = qFormula(u[1] - u[0], x[1] - x[0], u[2] - u[0], x[2] - x[0]);
            double qn = qFormula(u[n - 1] - u[n], x[n - 1] - x[n], u[n - 2] - u[n], x[n - 2] - x[n]);

            if (nn <= 0)
            {
                q1 = s[0];
                qn = s[1];
            }

            s[0] = 6.0d * ((u[1] - u[0]) / (x[1] - x[0]) - q1);
            int n1 = n - 1;
            for (int i = 1; i <= n1; i++)
                s[i] = (u[i - 1] / (x[i] - x[i - 1]) - u[i] *
                        (1.0d / (x[i] - x[i - 1]) + 1.0d / (x[i + 1] - x[i])) +
                        u[i + 1] / (x[i + 1] - x[i]))
                       * 6.0d;

            s[n] = 6.0 * (qn + (u[n1] - u[n]) / (x[n] - x[n1]));
            a[0] = 2.0 * (x[1] - x[0]);
            a[1] = 1.5d * (x[1] - x[0]) + 2.0d * (x[2] - x[1]);
            s[1] = s[1] - 0.5d * s[0];
            double c = double.NaN;
            for (int i = 2; i <= n1; i++)
            {
                c = (x[i] - x[i - 1]) / a[i - 1];
                a[i] = 2.0 * (x[i + 1] - x[i - 1]) - c * (x[i] - x[i - 1]);
                s[i] = s[i] - c * s[i - 1];
            }
            c = (x[n] - x[n1]) / a[n1];
            a[n] = (2.0d - c) * (x[n] - x[n1]);
            s[n] = s[n] - c * s[n1];

            s[n] = s[n] / a[n];
            for (int j = 1; j <= n1 + 1; j++)
            {
                int i = n - j;
                s[i] = (s[i] - (x[i + 1] - x[i]) * s[i + 1]) / a[i];
            }
            return s;
        }
        private double qFormula(double u1, double x1, double u2, double x2)
        {
            return (u1 / (x1 * x1) - u2 / (x2 * x2)) / (1.0d / x1 - 1.0d / x2);
        }
        private double[] cpfrom(double[] src, int ind)
        {
            double[] tmp = new double[src.Length];

            for (int i = ind; i < src.Length; i++)
                tmp[i - ind] = src[i];

            return tmp;
        }

        private double eval(double y, int nn, double[] x, double[] u, double[] s, ref int istart)
        {   //      function eval(y,nn,x,u,s)
            //  performs cubic spline interpolation of a function sampled at unequally
            //  spaced intervals.  the routine spline  should be called to set up the array s
            //  y  the coordinate at which function value is desired.
            //  nn  number of samples of original function.
            //  x  array containing sample coordinates. the sequence x(1),x(2).....x(nn)     
            //     must be strictly increasing.
            //  u  array containing samples of function at the coordinates x(1),x(2)...      
            //  s  array containing the 2nd derivatives at the sample points.  found by the  
            //     routine  spline, which must be called once before beginning interpolation.
            //  if  y  falls outside range(x(1),x(nn))  the value at the nearest endpoint    
            //  of the series is used.

            double eval = double.NaN;
            int k = 0;
            int k1 = 0;

            nn = Math.Abs(nn) - 1;

            if (y < x[0]) // out of range. substitue end values
            {
                eval = u[0];
                return eval;
            }
            else if (y > x[nn]) // out of range. substitue end values
            {
                eval = u[nn];
                return eval;
            }
            else if (y - x[istart] < 0) // locate interval (x(k1),x(k))  containing y
            {
                //  scan downwards in x array
                for (k = 1; k <= istart; k++)
                {// 1200 do 1300 k=1,istart
                    k1 = istart - k;
                    if (x[k1] <= y)
                        goto gt1350;
                }// 1300 continue
            gt1350:
                k = k1 + 1;
                goto gt1500;
            }
            else if (y - x[istart] >= 0)  // locate interval (x(k1),x(k))  containing y
            {
                //  scan up the x array
                for (k = istart; k <= nn; k++)
                {// 1000 do 1100 k=istart,nn
                    if (x[k] > y)
                        goto gt1150;
                }// 1100 continue
            gt1150:
                k1 = k - 1;
                goto gt1500;
            }
        gt1500:
            istart = k1;
            //c  evaluate interpolate
            double dy = x[k] - y;
            double dy1 = y - x[k1];
            double dk = x[k] - x[k1];
            eval = (((s[k1] * dy * dy * dy) + (s[k] * dy1 * dy1 * dy1)) * (1d / (6d * dk))) +
                (dy1 * ((u[k] / dk) - (s[k] * dk) / 6.0d)) + (dy * ((u[k1] / dk) - (s[k1] * dk) / 6.0d));

            return eval;
        }

        private double recurs(double[] hc, int nf, double[] om, double[] scr)
        {   //              subroutine recurs(x,n,hc,nf,om,scr)

            //c  does sin and cos recursion to fill in data x, of length n, for
            //c  nf sines and cosines with frequenciies om (normalized so the
            //c  nyquist is pi). hc contains alternating cosine and sine coefficients
            //c  scr is a scratch array of length 3*nf (n.b. - it is double precision)
            //c
            //c  set up for start of recursion by computing harmonic values
            //c  at start point and just before it
            for (int i = 1; i <= nf; i++)
            {
                scr[3 * i - 2 - 1] = hc[2 * i - 1 - 1];
                scr[3 * i - 1 - 1] = hc[2 * i - 1 - 1] * Math.Cos(om[i - 1]) - hc[2 * i - 1] * Math.Sin(om[i - 1]);
                scr[3 * i - 1] = 2d * Math.Cos(om[i - 1]);
            }

            double x = 0d;
            //  do recursive computation for each harmonic
            for (int j = 1; j <= nf; j++)
            {
                x += +scr[3 * j - 2 - 1];
                double sc = scr[3 * j - 2 - 1];
                scr[3 * j - 2 - 1] = scr[3 * j - 1] * sc - scr[3 * j - 1 - 1];
                scr[3 * j - 1 - 1] = sc;
            }
            return x;
        }
        #endregion
        #endregion
    }
}
