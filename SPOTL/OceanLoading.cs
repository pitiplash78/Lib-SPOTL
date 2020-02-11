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
    /* Original Code from D.C. Agnew down loaded "Current Version 3.3.0.2 (September 2013)"
     * at june 2014 from http://igppweb.ucsd.edu/~agnew/Spotl/spotlmain.html
     *
     * The whole functionatilty is not implemented. But everything for calculation of each component!! */

    /// <summary>
    /// Calculation of the Ocean Tide Loading depending on the (class) OceanLoadingProperties.
    /// </summary>
    public class OceanLoading
    {
        #region Components to be calculable
        public enum Component
        {
            None,
            Gravity,
            GravityDeformation,
            PotentialHeight,
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
        /// <summary>
        /// Path of the parameter and used files for ocean loading calculation after SPOTL.
        /// </summary>
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

                    filePath.EarthModels = _PathOceanLoading + FilePath._EarthModels + Path.DirectorySeparatorChar;
                    filePath.OceanModels = _PathOceanLoading + FilePath._OceanModels + Path.DirectorySeparatorChar;
                    filePath.LocalModels = _PathOceanLoading + FilePath._LocalModels + Path.DirectorySeparatorChar;

                    filePath.WaterDensity = _PathOceanLoading + FilePath._WaterDensityDirectory + Path.DirectorySeparatorChar + FilePath._WaterDensityFileName;
                    filePath.LndSeaMatrix = _PathOceanLoading + FilePath._LndSeaMatrixDirectory + Path.DirectorySeparatorChar + FilePath._LndSeaMatrixFileName;
                }
            }
        }
        /// <summary>
        /// Internal directory base path for ocean loading calculation
        /// </summary>
        private string _PathOceanLoading = null;

        /// <summary>
        /// Path for the XML-Leap Second Table
        /// </summary>
        private string PathLeapSecondTable = null;

        /// <summary>
        /// Object getting the difference DDT = ET-UTC or TDT - UTC
        /// </summary>
        private TimeConversion.TAI_UTC TAI_UTC = null;

        /// <summary>
        /// Object for the internal class of directory and file paths used in library
        /// </summary>
        private FilePath filePath = null;

        /// <summary>
        /// 
        /// </summary>
        public OceanLoadingProperties.GreenFunctionEntry GreensFunctionModel { get; private set; }

        /// <summary>
        /// 
        /// </summary>
        public OceanLoadingProperties.OceanModelEntry[] OceanTideModel { get; private set; }

        public OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanTideModel { get; private set; }

        /// <summary>
        /// Holds the information about the calculated parameters.
        /// </summary>
        public OceanLoadingParameter oceanLoadingParameter { get; private set; }

        /// <summary>
        /// Unit parameter for contructor for calculation the gravity.
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_Gravity
        {
            get
            {
                return _OutputUnit_Gravity;
            }
            set
            {
                _OutputUnit_Gravity = value;
                dConversion_Gravity = Units.Units.getConversion(BaseUnits.UnitGravity, _OutputUnit_Gravity).Factor;
            }
        }
        private Units.Units.UnitNamesEnum _OutputUnit_Gravity = BaseUnits.UnitGravity;
        private double dConversion_Gravity = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the potential height.
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_PotentialHeight
        {
            get
            {
                return _OutputUnit_PotentialHeight;
            }
            set
            {
                _OutputUnit_PotentialHeight = value;
                dConversion_PotentialHeight = Units.Units.getConversion(BaseUnits.UnitGravity, _OutputUnit_PotentialHeight).Factor;
            }
        }
        private Units.Units.UnitNamesEnum _OutputUnit_PotentialHeight = BaseUnits.UnitPotentialHeight;
        private double dConversion_PotentialHeight = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the tilt.
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_Tilt
        {
            get
            {
                return _OutputUnit_Tilt;
            }
            set
            {
                _OutputUnit_Tilt = value;
                dConversion_Tilt = Units.Units.getConversion(BaseUnits.UnitTilt, _OutputUnit_Tilt).Factor;
            }
        }
        private Units.Units.UnitNamesEnum _OutputUnit_Tilt = BaseUnits.UnitTilt;
        private double dConversion_Tilt = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the strain.
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_Strain
        {
            get
            {
                return _OutputUnit_Strain;
            }
            set
            {
                _OutputUnit_Strain = value;
                dConversion_Strain = Units.Units.getConversion(BaseUnits.UnitStrain, _OutputUnit_Strain).Factor;
            }
        }
        private Units.Units.UnitNamesEnum _OutputUnit_Strain = BaseUnits.UnitStrain;
        private double dConversion_Strain = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the displacement.
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_Displacement
        {
            get
            {
                return _OutputUnit_Displacement;
            }
            set
            {
                _OutputUnit_Displacement = value;
                dConversion_Displacement = Units.Units.getConversion(BaseUnits.UnitDisplacement, _OutputUnit_Displacement).Factor;
            }
        }
        private Units.Units.UnitNamesEnum _OutputUnit_Displacement = BaseUnits.UnitDisplacement;
        private double dConversion_Displacement = 1d;

        /// <summary>
        /// Unit parameter for contructor for calculation the ocean tide.
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_OceanTide
        {
            get
            {
                return _OutputUnit_OceanTide;
            }
            set
            {
                _OutputUnit_OceanTide = value;
                dConversion_OceanTide = Units.Units.getConversion(BaseUnits.UnitOceanTide, _OutputUnit_OceanTide).Factor;
            }
        }
        private Units.Units.UnitNamesEnum _OutputUnit_OceanTide = BaseUnits.UnitOceanTide;
        private double dConversion_OceanTide = 1d;

        public static class BaseUnits
        {
            public static Units.Units.UnitNamesEnum UnitGravity = Units.Units.UnitNamesEnum.MicroGal;
            public static Units.Units.UnitNamesEnum UnitPotentialHeight = Units.Units.UnitNamesEnum.Millimeter;
            public static Units.Units.UnitNamesEnum UnitDisplacement = Units.Units.UnitNamesEnum.Millimeter;
            public static Units.Units.UnitNamesEnum UnitTilt = Units.Units.UnitNamesEnum.NanoRadian;
            public static Units.Units.UnitNamesEnum UnitStrain = Units.Units.UnitNamesEnum.NanoStrain;
            public static Units.Units.UnitNamesEnum UnitOceanTide = Units.Units.UnitNamesEnum.Meter;
        }

        public void setUnit(Component component, Units.Units.UnitNamesEnum Unit)
        {
            switch (component)
            {
                case SPOTL.OceanLoading.Component.Gravity:
                case SPOTL.OceanLoading.Component.GravityDeformation:
                    OutputUnit_Gravity = Unit;
                    break;
                case SPOTL.OceanLoading.Component.HorizontalDisplacement:
                case SPOTL.OceanLoading.Component.VerticalDisplacement:
                    OutputUnit_Displacement = Unit;
                    break;
                case SPOTL.OceanLoading.Component.HorizontalStrain:
                case SPOTL.OceanLoading.Component.ShearStrain:
                case SPOTL.OceanLoading.Component.ArealStrain:
                case SPOTL.OceanLoading.Component.VolumeStrain:
                    OutputUnit_Strain = Unit;
                    break;
                case SPOTL.OceanLoading.Component.Tilt:
                case SPOTL.OceanLoading.Component.TiltAttraction:
                case SPOTL.OceanLoading.Component.TiltDeformation:
                    OutputUnit_Tilt = Unit;
                    break;
                case SPOTL.OceanLoading.Component.PotentialHeight:
                    OutputUnit_PotentialHeight = Unit;
                    break;
                case SPOTL.OceanLoading.Component.OceanTide:
                    OutputUnit_OceanTide = Unit;
                    break;
            }
        }
        #endregion

        /// <summary>
        /// Constructor for oceanloading prediction, calculating the ocean loading parameter for given settings.
        /// </summary>
        /// <param name="pathOceanLoading">Path of the parameter and used files for ocean loading calculation after SPOTL.</param>
        /// <param name="oceanLoadingProperties">Ocean loading properties.</param>
        /// <param name="Location">Location at where the ocean loading has to be predicted</param>
        /// <param name="GreensFunctionModel">Greensfunction to be used for oceanloading prediction.</param>
        /// <param name="OceanTideModel">Array of ocean tidal models for ocean loading prediction.</param>
        /// <param name="LocalOceanTideModel">Array of local ocean tidal models for ocean loading prediction.</param>
        public OceanLoading(string pathOceanLoading,
                            string pathLeapSecondTable,
                            OceanLoadingProperties oceanLoadingProperties,
                            OceanLoadingProperties.Location Location,
                            OceanLoadingProperties.GreenFunctionEntry GreensFunctionModel,
                            OceanLoadingProperties.OceanModelEntry[] OceanTideModel,
                            OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanTideModel = null)
        {
            this.PathOceanLoading = pathOceanLoading;
            this.PathLeapSecondTable = pathLeapSecondTable;

            TAI_UTC = TimeConversion.TAI_UTC.deserialisieren(PathLeapSecondTable);

            this.stloc = new StationLocation(Location.Name, Location.Longitude, Location.Latitude, Location.Height);

            this.GreensFunctionModel = GreensFunctionModel;

            this.OceanTideModel = OceanTideModel;

            if (LocalOceanTideModel != null)
                this.LocalOceanTideModel = LocalOceanTideModel;

            greensFunction = GreensFunction.readGreensfunction(filePath.EarthModels + GreensFunctionModel.FileName);
            landSeaMatrix = LandSeaMatrix.readLandSea(filePath.LndSeaMatrix);
            seaDensity = SeaDensity.readWaterDensity(filePath.WaterDensity);

            // Preparation of polygon information if local model is selected
            if (this.LocalOceanTideModel != null)
            {
                this.polInf = new PolInf(this.LocalOceanTideModel.Length);
                for (int i = 0; i < this.LocalOceanTideModel.Length; i++)
                    polInf.npoly[i] = PolInf.readPolgonInformation(filePath.LocalModels + this.LocalOceanTideModel[i].FilePath + Path.DirectorySeparatorChar + this.LocalOceanTideModel[i].FilePath + ".polygon");
            }

            for (int j = 0; j < oceanLoadingProperties.WaveSequenceOceanModel.Length; j++)
            {
                string dsym = oceanLoadingProperties.WaveSequenceOceanModel[j].ToLower();
                OceanLoadWave oceanLoadWave = null;
                OceanModel oceanModel = null;

                if (polInf != null)
                {
                    polInf.np = 0;
                    for (int i = 0; i < polInf.npoly.Length; i++)
                        polInf.npoly[i].use = true;
                }

                // Local Ocean Models
                if (LocalOceanTideModel != null)
                {
                    for (int i = 0; i < LocalOceanTideModel.Length; i++)
                    {
                        OceanLoad oceanLoadLocal = null;
                        for (int k = 0; k < LocalOceanTideModel[i].DarwinSymbol.Length; k++)
                        {
                            if (LocalOceanTideModel[i].DarwinSymbol[k].ToLower() == dsym)
                            {
                                oceanModel = OceanModel.readTidalModel(filePath.LocalModels + LocalOceanTideModel[i].FilePath + Path.DirectorySeparatorChar +
                                                                       LocalOceanTideModel[i].DarwinSymbol[k] + "." + LocalOceanTideModel[i].FilePath + ".asc.gz");
                                polInf.np = i + 1;
                                oceanLoadLocal = nloadf(oceanModel);
                                oceanLoadLocal.LocalModel = true.ToString();
                                oceanLoadLocal.OceanModel = LocalOceanTideModel[i].FilePath;
                                polInf.npoly[i].use = false;
                            }
                        }
                        if (oceanLoadLocal != null)
                        {
                            if (oceanLoadWave == null)
                                oceanLoadWave = new OceanLoadWave();
                            oceanLoadWave.AddEntry(oceanLoadLocal);
                        }
                    }
                }
                // Glogal Ocean Model
                if (OceanTideModel != null)
                {
                    for (int i = 0; i < OceanTideModel.Length; i++)
                    {
                        OceanLoad oceanLoadGlobal = null;
                        for (int k = 0; k < OceanTideModel[i].DarwinSymbol.Length; k++)
                        {
                            if (OceanTideModel[i].DarwinSymbol[k].ToLower() == dsym)
                            {
                                oceanModel = OceanModel.readTidalModel(filePath.OceanModels + OceanTideModel[i].FilePath + Path.DirectorySeparatorChar +
                                                                       OceanTideModel[i].DarwinSymbol[k] + "." + OceanTideModel[i].FilePath + ".asc.gz");
                                oceanLoadGlobal = nloadf(oceanModel);
                                oceanLoadGlobal.LocalModel = false.ToString();
                                oceanLoadGlobal.OceanModel = OceanTideModel[i].FilePath;
                            }
                        }
                        if (oceanLoadGlobal != null)
                        {
                            if (oceanLoadWave == null)
                                oceanLoadWave = new OceanLoadWave();
                            oceanLoadWave.AddEntry(oceanLoadGlobal);
                        }
                    }
                }

                if (oceanLoadWave != null)
                {
                    oceanLoadWave.dsym = oceanModel.dsym;
                    oceanLoadWave.icte = oceanModel.icte;
                    loadcomb(ref oceanLoadWave);

                    if (oceanLoadingParameter == null)
                        oceanLoadingParameter = new OceanLoadingParameter();

                    oceanLoadingParameter.AddEntry(oceanLoadWave);
                }
            }

            if (oceanLoadingParameter != null)
            {
                oceanLoadingParameter.Location = new Location(stloc.Name, stloc.rlam, stloc.rlat, stloc.ht);
                oceanLoadingParameter.GreensFunction = new OceanLoadingProperties.GreenFunctionEntry(GreensFunctionModel);
                oceanLoadingParameter.OceanModel = OceanTideModel;
                oceanLoadingParameter.LocalOceanModel = LocalOceanTideModel;
            }
        }

        #region private used objects containing the information for the calculation
        /// <summary>
        /// Object contains the intenal needed information about the station location
        /// </summary>
        private StationLocation stloc = null;

        /// <summary>
        /// Class contains the intenal needed information about the station location
        /// </summary>
        private class StationLocation
        {
            /// <summary>
            /// Name of the station.
            /// </summary>
            public string Name;

            /// <summary>
            /// Latitude in degree [-90..90].
            /// </summary>
            public double rlat;

            /// <summary>
            /// Longitude in degree [-180..180].
            /// </summary>
            public double rlam;

            /// <summary>
            /// Height im metern.
            /// </summary>
            public double ht;

            /// <summary>
            /// Cosinus of the co-latitude in rad.
            /// </summary>
            public double ct;

            /// <summary>
            /// Sinus of the co-latitude in rad.
            /// </summary>
            public double st;

            /// <summary>
            /// Constructor for station information
            /// </summary>
            /// <param name="StationName">Name of the station.</param>
            /// <param name="Longitude">Longitude [degree] of the station.</param>
            /// <param name="Latitude">Latitude [degree] of the station.</param>
            /// <param name="Height">Heigth [meter] of the station.</param>
            public StationLocation(string StationName, double Longitude, double Latitude, double Height)
            {
                this.Name = StationName;
                this.rlat = Latitude;
                this.rlam = Longitude;
                this.ht = Height;

                this.ct = Math.Cos(Constants.degree2radian * (90.0 - Latitude));
                this.st = Math.Sin(Constants.degree2radian * (90.0 - Latitude));
            }
        }

        /// <summary>
        /// Object of the GreensFunction class containing all information (meta-data, data).
        /// </summary>
        private GreensFunction greensFunction = null;

        /// <summary>
        /// Object of the LandSeaMatrix class containing all information (meta-data, data).
        /// </summary>
        private LandSeaMatrix landSeaMatrix = null;

        /// <summary>
        /// Object of the SeaDensity class containing all information (meta-data, data).
        /// </summary>
        private SeaDensity seaDensity = null;

        /// <summary>
        /// Object containing additional informations about the load calculation.
        /// </summary>
        private ConInf conInf = new ConInf();

        /// <summary>
        /// The struct contains additional informations about the load calculation.
        /// </summary>
        private struct ConInf
        {
            /// <summary>
            /// Actual distance for the loop over the distance.
            /// </summary>
            public double dist;

            /// <summary>
            /// Closest nonzero load distance.
            /// </summary>
            public double close;

            /// <summary>
            /// Flag if closest nonzero load distance is found.
            /// </summary>
            public bool nclose;

            /// <summary>
            /// Latitude of the closest nonzero load distance.
            /// </summary>
            public double clat;

            /// <summary>
            /// Longitude of the closest nonzero load distance.
            /// </summary>
            public double clong;

            /// <summary>
            /// Count of zero loads.
            /// </summary>
            public int numnf;

            /// <summary>
            /// Range, where zero load was occured (nearest distance).
            /// </summary>
            public double clnf;

            /// <summary>
            /// Range, where zero load was occured (farest distance).
            /// </summary>
            public double farnf;
        };

        /// <summary>
        /// Object of the PolInf class containing all information about the polygons for local ocean models included (meta-data, data).
        /// </summary>
        private PolInf polInf = null;

        /// <summary>
        /// Object holds the varibles calculated in function ignewt.
        /// </summary>
        private IgnewtVariables ignewtVariables;

        /// <summary>
        /// Class holds the varibles calculated in function ignewt.
        /// </summary>
        private struct IgnewtVariables
        {
            public bool iflg;
            public double eps;
            public double eps1;
            public double eps2;
            public double g2;
            public double em;
            public double plc;
        }
        #endregion

        /// <summary>
        /// Constructor for oceanloading prediction, on give ocean loading parameter.
        /// </summary>
        /// <param name="nout"></param>

        public OceanLoading(string pathLeapSecondTable, OceanLoadingParameter oceanLoadingParameter)
        {
            TAI_UTC = TimeConversion.TAI_UTC.deserialisieren(pathLeapSecondTable);

            this.oceanLoadingParameter = oceanLoadingParameter;
        }

        #region calculates the wave parameters with respect to the oceanloading
        private enum ModeOfOperation
        {
            None,
            /// <summary>
            /// l - means that the the phases are in local phase at the site, again with lags negative 
            /// </summary>
            l,  
            /// <summary>
            /// g - means that the phases are Greenwich phase, BUT with lags negative (the reverse of the usual phase convention in ocean tides) 
            /// </summary>
            g, 
            /// <summary>
            /// m - causes the program to write to standard output the coor dinates of the corners  
            /// <para>of the cells used for the convolution, in a form (lat and long separated by</para>
            /// <para>impossible values) that is often suitable for plotting. This output is quite voluminous and should</para>
            /// <para>always be redirected to a file.</para>
            /// </summary>
            m,  
        }


        /// <summary>
        /// Computes the load tide caused by the oceans, using a station-centered grid and integrated green functions. */
        /// </summary>
        private OceanLoad nloadf(OceanModel oceanModel)
        {
            /*     Revision 1.2  2012/03/06 03:32:43  agnew - removed density (now done in ocmodl)
             *     Revision 1.1  2011/09/23 22:39:02  agnew - Initial revision */

            /*   computes the load tide caused by the oceans, using a station-centered grid and integrated green functions.
             *    calls:  ignewt (for the newtonian green functions)
             *            ocmodl (for the ocean load - calls ocstart to read model in)
             *            getcl  (to get command-line parameters) - obesolete
             *            getgrf (to get Green functions) - obsolete
             *            lodout (to write results to standard output) */

            /* mode of operation: 'l' - means that the the phases are in local phase at the site, again with lags negative.
             *                    'g' - means that the phases are Greenwich phase, BUT with lags negative (the reverse of
             *                          the usual phase convention in ocean tides).
             *                    'm' - m causes the program to write to standard output the coor dinates of the corners 
             *                          of the cells used for the convolution, in a form (lat and long separated by
             *                          impossible values) that is often suitable for plotting. This output is quite voluminous and should
             *                          always be redirected to a file.
             *                          
             * modo is set to 'l' caused by the generall wanted operation mode. */
            ModeOfOperation modo = ModeOfOperation.l;
            // distance for mode 'm'
            double distm = 0;

            Complex amp = Complex.Zero;
            Complex grav = Complex.Zero;
            Complex gravdef = Complex.Zero;
            Complex pothi = Complex.Zero;
            Complex[] tilt = new Complex[2] { Complex.Zero, Complex.Zero };
            Complex[] tiltdef = new Complex[2] { Complex.Zero, Complex.Zero };
            Complex[] str = new Complex[3] { Complex.Zero, Complex.Zero, Complex.Zero };
            Complex[] disp = new Complex[3] { Complex.Zero, Complex.Zero, Complex.Zero };
            conInf.nclose = false;
            /* convolutions start here - read in Green functions, which come in ntot sections; num is the current one, fingrd is whether or not to
             * use the detailed land-sea grid */
            for (int num = 0; num < greensFunction.GrennsFunctionItems.Length; num++)// while (num < ntot) //L3:
            {
                // call getgrf(grfile,num,ntot,nf,fingrd)                                   - obsolete due to feeding in at function call
                double stp = Constants.degree2radian *
                   greensFunction.GrennsFunctionItems[num].spc;                             //       stp = dr*spc

                /* if this is the first convolution interval see if (1) there is ocean at the station location and (2) if it is below sea level. If both are true,
                 * compute the gravity load (only) for a disk from distance 0 to the start of the first integrated Green function */
                double g = 0.0d;
                double t = 0.0d;
                double pot = 0.0d;
                if (num == 0)                                                               //       if(num.eq.1) then
                {
                    amp = ocmodl(stloc.rlat, stloc.rlam,
                                 ref oceanModel,
                                 ref greensFunction.GrennsFunctionItems[num].fingrd);              //          call ocmodl(rlat,rlam,fingrd,amp)
                    if ((amp != Complex.Zero) && stloc.ht < 0.0)                            // 	 if(amp.ne.cz.and.ht.lt.0) then
                    {
                        double del = Constants.degree2radian *
                            greensFunction.GrennsFunctionItems[num].beg / 2;                // 	    del=dr*beg/2
                        ignewt(del, del * 2, ref g, ref t, ref pot);                        //             call ignewt(del,2*del,g,t,pot)
                        grav += g * amp * Math.PI * 2;                                      //             grav = grav + g*amp*tpi
                        conInf.close = 0.0;                                                 // 	    close = 0.
                        conInf.clat = stloc.rlat;                                           // 	    clat = rlat
                        conInf.clong = stloc.rlam;                                          // 	    clong = rlam
                        conInf.nclose = true;                                               // 	    nclose = 1
                    }                                                                       // 	 endif
                }                                                                           //       endif

                /*  loop through distance begins here */
                for (int ii = 0; ii < greensFunction.GrennsFunctionItems[num].ngr; ii++)    //       do 13 ii=1,nf
                {
                    conInf.dist = greensFunction.GrennsFunctionItems[num].beg + ii *
                        greensFunction.GrennsFunctionItems[num].spc;                        //       dist = beg + (ii-1.)*spc
                    double del = Constants.degree2radian * conInf.dist;                     //       del = dr*dist
                    double cd = Math.Cos(del);                                              //       cd = cos(del)
                    double sd = Math.Sin(del);                                              //       sd = sin(del)
                    ignewt(del, stp, ref g, ref t, ref pot);                                //       call ignewt(del,stp,g,t,pot)
                    g += greensFunction.GrennsFunctionItems[num].grfn[ii, 2];               //       g = g + grfn(ii,3)
                    t += greensFunction.GrennsFunctionItems[num].grfn[ii, 3];               //       t = t + grfn(ii,4)
                    double gdef = greensFunction.GrennsFunctionItems[num].grfn[ii, 2];
                    double tdef = greensFunction.GrennsFunctionItems[num].grfn[ii, 3];
                    pot += greensFunction.GrennsFunctionItems[num].grfn[ii, 6];             //       pot = pot + grfn(ii,7)
                    double u = greensFunction.GrennsFunctionItems[num].grfn[ii, 0];         //       u = grfn(ii,1)
                    double v = greensFunction.GrennsFunctionItems[num].grfn[ii, 1];         //       v = grfn(ii,2)
                    double ett = greensFunction.GrennsFunctionItems[num].grfn[ii, 4];       //       ett = grfn(ii,5)
                    double ell = greensFunction.GrennsFunctionItems[num].grfn[ii, 5];       //       ell = grfn(ii,6)

                    /* find the number of azimuth steps (scaled with sin(del), but kept above some minimum) */
                    int naz = (int)(sd * 360.0);                                            //       naz = 360.*sd

                    if (conInf.dist < 180.0)                                                //       if(dist.lt.180.) naz = max0(naz,minaz)
                        naz = Math.Max(naz, Constants.minaz);

                    double azstp = (Math.PI * 2.0d) / naz;                                  //       azstp = tpi/naz
                    double azstpd = 360.0 / naz;                                            //       azstpd = 360./naz
                    double azimuth = azstpd / 2.0d;                                         //       azimuth = azstpd/2
                    double caz = Math.Cos(azstp / 2.0d);                                    //       caz = cos(azstp/2)
                    double saz = Math.Sin(azstp / 2.0d);                                    //       saz = sin(azstp/2)
                    double saztp = Math.Sin(azstp);                                         //       saztp = sin(azstp)
                    double caztp = Math.Cos(azstp);                                         //       caztp = cos(azstp)
                    double stpfac = saz * 2;                                                //       stpfac = 2*saz

                    /*  loop through azimuth begins here */
                    for (int jj = 1; jj <= naz; jj++)                                       //       do 11 jj=1,naz
                    {
                        /* do the spherical trigonometry to find the cell lat and long */
                        double cb = cd * stloc.ct + sd * stloc.st * caz;                    //       cb = cd*ct + sd*st*caz
                        if (Math.Abs(cb) > 1.0d)                                            //       if(abs(cb).gt.1) cb = cb/abs(cb)
                            cb /= Math.Abs(cb);
                        double sb = Math.Sqrt(1.0d - cb * cb);                              //       sb = sqrt(1.-cb*cb)
                        double rlato = 90.0d - Constants.radian2degree * Math.Acos(cb);     //       rlato = 90 - rd*acos(cb) 
                        double rlong = 0.0d;                                                //       rlong = 0.
                        /*  if near the north pole leave longitude at zero- otherwise find it */
                        if (sb > 0.001d)                                                    //       if(sb.gt.1.e-3) then
                        {
                            sb = 1.0 / sb;                                                  //          sb = 1./sb
                            double sg = sd * saz * sb;                                      //          sg = sd*saz*sb
                            double cg = (stloc.st * cd - sd * stloc.ct * caz) * sb;         //          cg = (st*cd - sd*ct*caz)*sb
                            rlong = stloc.rlam +
                                Constants.radian2degree * Math.Atan2(sg, cg);               //          rlong = rlam + rd*atan2(sg,cg)
                        }                                                                   //       endif
                        amp = ocmodl(rlato, rlong, ref oceanModel,
                                     ref greensFunction.GrennsFunctionItems[num].fingrd);   //       call ocmodl(rlato,rlong,fingrd,amp)

                        if (amp != Complex.Zero && modo != ModeOfOperation.m)               //       if(amp.ne.(0.,0.).and.modo.ne.'m') then
                        {
                            /* if this is the first nonzero amplitude, save the distance and location */
                            if (!conInf.nclose)                                             //     if(nclose.eq.0) then
                            {
                                conInf.close = conInf.dist;                                 // 	    close = dist
                                conInf.clat = rlato;                                        // 	    clat = rlato
                                conInf.clong = rlong;                                       // 	    clong = rlong
                                conInf.nclose = true;                                       // 	    nclose = 1
                            }                                                               //          endif
                            Complex cmp = amp * stpfac;                                     //          cmp = amp*stpfac
                            /*  compute the loads.  disp is vert,ns,ew; tilt is ns,ew; strain is ns,ew,ne shear */
                            disp[0] += u * amp * azstp;                                     //          disp(1) = disp(1) + u*amp*azstp
                            disp[1] += v * caz * cmp;                                       //          disp(2) = disp(2) + v*caz*cmp
                            disp[2] += v * saz * cmp;                                       //          disp(3) = disp(3) + v*saz*cmp

                            grav += g * amp * azstp;                                        //          grav = grav + g*amp*azstp
                            gravdef += gdef * amp * azstp;

                            pothi += pot * amp * azstp;                                     // 	        pothi = pothi + pot*amp*azstp

                            tilt[0] += t * caz * cmp;                                       //          tilt(1) = tilt(1) + t*caz*cmp
                            tilt[1] += t * saz * cmp;                                       //          tilt(2) = tilt(2) + t*saz*cmp
                            tiltdef[0] += tdef * caz * cmp;
                            tiltdef[1] += tdef * saz * cmp;

                            double aa = (caz * caz - saz * saz) * 0.5;                      //          aa = .5*(caz*caz - saz*saz)
                            double bb = azstp * 0.5 + aa * stpfac;                          //          bb = .5*azstp + aa*stpfac
                            double cc = azstp * 0.5 - aa * stpfac;                          //          cc = .5*azstp - aa*stpfac
                            str[0] += (ett * bb + ell * cc) * amp;                          //          str(1) = str(1) + (ett*bb + ell*cc)*amp
                            str[1] += (ell * bb + ett * cc) * amp;                          //          str(2) = str(2) + (ell*bb + ett*cc)*amp
                            str[2] += (ett - ell) * saz * caz * saztp * caztp * amp;        //          str(3) = str(3) + (ett-ell)*saz*caz*saztp*caztp*amp
                        }
                        else if ((amp != Complex.Zero) &&
                                 modo ==  ModeOfOperation.m && conInf.dist <= distm)        //       elseif(amp.ne.(0.,0.).and.modo.eq.'m'.and.dist.le.distm) then
                        {
                            /*  otherwise, output locations of cell corners, for later plotting */
                            ldbxdr(azimuth, conInf.dist, azstpd,
                                greensFunction.GrennsFunctionItems[num].spc);               //          call ldbxdr(azimuth,dist,azstpd,spc)
                        }                                                                   //       endif

                        /*  use recursion to step sin and cos of the azimuth */
                        azimuth += azstpd;                                                  //       azimuth = azimuth + azstpd
                        double xx = saz * caztp + caz * saztp;                              //       xx = saz*caztp + caz*saztp
                        caz = caz * caztp - saz * saztp;                                    //       caz = caz*caztp - saz*saztp
                        saz = xx;                                                           //       saz = xx
                    }                                                                       //  11   continue
                }                                                                           //  13   continue
            }
            //       if(num.lt.ntot) go to 3                                                - obsolete due to used while loop over the Greens function elements 

            // Calculate ocean tide waves taken from oclook.f
            Complex oamp = Complex.Zero;
            {
                GreensFunction.FinGrd fingrd = GreensFunction.FinGrd.F;                     //         fingrd = 'F'
                oamp = ocmodl(stloc.rlat, stloc.rlam, ref oceanModel, ref fingrd);                 //          call ocmodl(rlato,rlong,fingrd,amp)
                if (oamp != Complex.Zero)                                                   //          if(amp.ne.cz) amp=amp/rho
                    oamp /= seaDensity.rho;

                if (polInf != null)                                                         // 	 if(ispoly) then
                {
                    bool iok = polInf.chkgon(stloc.rlam, stloc.rlat);                       // 	    call chkgon(rlong,rlato,iok)
                    if (!iok) // == 0)                                                      // 	    if(iok.eq.0) amp = cz
                        oamp = Complex.Zero;
                }                                                                           //          endif

                oamp = Complex.Conjugate(oamp);                                             //             amp = conjg(amp)
                oamp = phasor(oamp);                                                        //             amp = phasor(amp)
            }


            if (modo == ModeOfOperation.m)                                                  //       if(modo.eq.'m') stop
                return new OceanLoad();

            /*  done with convolution - add to disk file, and convert to amp and phase */
            Complex[] cload = new Complex[14] { Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero,
                                                Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero };
            cload[0] = grav;                                                                //       cload(1) = grav
            cload[10] = gravdef;
            cload[9] = pothi;                                                               //       cload(10) = pothi
            cload[4] = tilt[0];                                                             //       cload(5) = tilt(1)
            cload[5] = tilt[1];                                                             //       cload(6) = tilt(2)
            cload[11] = tiltdef[0];
            cload[12] = tiltdef[1];
            for (int i = 0; i < 3; i++)                                                     //       do 15 i=1,3
            {
                cload[i + 1] = disp[i];                                                     //       cload(i+1) = disp(i)
                cload[i + 6] = str[i];                                                      //       cload(i+6) = str(i)
            }                                                                               //  15   continue
            cload[13] = oamp;

            return lodout(cload, modo, ref greensFunction, ref oceanModel, ref polInf);     //       call lodout(cload,stnam,rlat,rlam,ht,modo)
        }

        private void oclook(OceanModel oceanModel, ref OceanLoad oceanload)
        {
            /* program oclook - looks at model values: 
             *   either prints out values for a particular place outputs grid within an area (for plotting) 
             *   outputs values within an area */

            /*     Revision 1.7  2013/03/11 02:01:28  agnew  - Fixed divide-by-zero if amp (so possibly rho) is zero
             *     Revision 1.6  2012/11/20 16:29:32  agnew  - fixed bug with coarse flag 
             *     Revision 1.5  2012/06/25 01:30:11  agnew  - corrected use of density output from ocmodl by adding density common block
             *     Revision 1.4  2012/03/13 21:17:50  agnew  - Modified polygon common block 
             *     Revision 1.3  2011/12/23 00:54:45  agnew  - fixed index error for g a p r i option, modernized loops, cleaned up code for e option
             *     Revision 1.2  2011/11/27 04:31:33  agnew  - larger dimension for polygon information 
             *     Revision 1.1  2011/09/23 22:39:02  agnew  - Initial revision */

            Complex amp = Complex.Zero;

            char modo = 'o';

            if (modo == 'o')                                                    //       elseif(modo.eq.'o') then
            {
                GreensFunction.FinGrd fingrd = GreensFunction.FinGrd.F;         //         fingrd = 'F'
                ocmodl(stloc.rlat, stloc.rlam,
                    ref oceanModel,
                    ref fingrd);               //          call ocmodl(rlato,rlong,fingrd,amp)
                if (amp != Complex.Zero)                                        //          if(amp.ne.cz) amp=amp/rho
                    amp /= seaDensity.rho;

                if (polInf != null)                                             // 	 if(ispoly) then
                {
                    bool iok = polInf.chkgon(stloc.rlam, stloc.rlat);           // 	    call chkgon(rlong,rlato,iok)
                    if (!iok) // == 0)                                          // 	    if(iok.eq.0) amp = cz
                        amp = Complex.Zero;
                }                                                               //          endif

                //if (amp == Complex.Zero)         // 	 if(amp.eq.cz) then
                {
                    //    //             write(6,119)
                    //    //  119        format('Tide is zero or excluded by polygon')

                }       // 	 endif

                //if (amp != Complex.Zero)                                        // 	 if(amp.ne.cz) then
                {
                    //             write(6,121) 'S',stnam, rlato, rlong, ht
                    //  121        format(a1,2x,a40,3x,f10.4,f10.4,f10.0)
                    //             write(6,123) 'O',dsym,(icte(j),j=1,6),mdnam
                    //  123        format(a1,1x,a4,5x,6i2,5x,a50)
                    //             write(6,125)
                    //  125        format('L g          Phases are Greenwich, lags negative')
                    //             write(6,127)
                    //  127        format('X')

                    amp = Complex.Conjugate(amp);                               //             amp = conjg(amp)
                    amp = phasor(amp);                                          //             amp = phasor(amp)

                    oceanload.o = new Load(amp);                                //             write(6,129) amp //  129        format('o',9x,2f10.4)
                }       // 	 endif
            }
            return;               //       end
        }
        
        // private static class modlim
        // {
        //    public static double tlat, blat, wlong, elong;
        //    public static int latc, longc;
        //}
        //private void oclook_orig(OceanModel oceanModel)
        //{
        //    /* program oclook - looks at model values: 
        //     *   either prints out values for a particular place outputs grid within an area (for plotting) 
        //     *   outputs values within an area */

        //    /*     Revision 1.7  2013/03/11 02:01:28  agnew  - Fixed divide-by-zero if amp (so possibly rho) is zero
        //     *     Revision 1.6  2012/11/20 16:29:32  agnew  - fixed bug with coarse flag 
        //     *     Revision 1.5  2012/06/25 01:30:11  agnew  - corrected use of density output from ocmodl by adding density common block
        //     *     Revision 1.4  2012/03/13 21:17:50  agnew  - Modified polygon common block 
        //     *     Revision 1.3  2011/12/23 00:54:45  agnew  - fixed index error for g a p r i option, modernized loops, cleaned up code for e option
        //     *     Revision 1.2  2011/11/27 04:31:33  agnew  - larger dimension for polygon information 
        //     *     Revision 1.1  2011/09/23 22:39:02  agnew  - Initial revision */

        //    Complex cz = Complex.Zero;
        //    Complex oamp = Complex.Zero;
        //    char fingrd = 'C';
        //    string stnam = "  Interpolated ocean-tide value         ";

        //    /* Format strings */
        //    string fmt_100 = "(\002usage: oclook modelfile e\002,/,\002 or   "
        //        + " oclook modelfile lat long d\002,/,\002 or    oclook modelfile l"
        //        + "at long o\002,/,\002 or    oclook modelfile slat nlat wlong elon"
        //        + "g g\002,/,\002 or    oclook modelfile slat nlat wlong elong a"
        //        + " \002,/,\002 or    oclook modelfile slat nlat wlong elong p \002"
        //        + ",/,\002 or    oclook modelfile slat nlat wlong elong r \002,/"
        //       + ",\002 or    oclook modelfile slat nlat wlong elong i \002,/,\002"
        //        + " (in all cases but the first, a polyfile may be put \002,/,\002 "
        //        + "     on the end)\002)";
        //    string fmt_102 = "(f10.0)";
        //    string fmt_104 = "(\002Cell indices are \002,i5,i5,\002 (\002,i7"
        //         + ",\002)\002,\002 with fractional locations \002,2f6.3)";
        //    string fmt_106 = "(\002Cell amps (R&I,[amp&ph]) are \002,2f10.4"
        //         + ",/,\002 [\002,2f10.4,\002 ]\002)";
        //    string fmt_107 = "(\002No ocean in this cell\002)";
        //    string fmt_108 = "(\002Fine-grid (interpolated) amps (R&I,[amp&ph"
        //         + "]) are \002,2f10.4,/,\002 [\002,2f10.4,\002 ]\002)";
        //    string fmt_109 = "(\002Interpolated tide is zero\002)";
        //    string fmt_110 = "(\002Point would be excluded by polygon\002)";
        //    string fmt_119 = "(\002Tide is zero or excluded by polygon\002)";
        //    string fmt_121 = "(a1,2x,a40,3x,f10.4,f10.4,f10.0)";
        //    string fmt_123 = "(a1,1x,a4,5x,6i2,5x,a50)";
        //    string fmt_125 = "(\002L g          Phases are Greenwich, lags ne"
        //         + "gative\002)";
        //    string fmt_127 = "(\002X\002)";
        //    string fmt_129 = "(\002o\002,9x,2f10.4)";
        //    string fmt_132 = "(\002Area is outside model\002)";
        //    string fmt_142 = "(2f10.4)";
        //    string fmt_140 = "(2f10.4,g14.6)";

        //    /* Local variables */
        //    double x, y;
        //    int[] i1 = new int[2];
        //    int[] i2 = new int[2];
        //    int ih, jh, il, jl;
        //    //double ht;
        //    Complex amp;
        //    char modo;//, dumm;
        //    double slat;
        //    int ngrs;
        //    Complex toloc;
        //    double rlato, rlong, rnlat;
        //    double ftlong, rilong;

        //    //       if(iargc().lt.2.or.iargc().gt.7) then
        //    //if (iargc_() < 2 || iargc_() > 7)
        //    //    {
        //    // 	 write(6,100)
        //    //  100  
        //    // 	 stop
        //    //       endif
        //    //    }
        //    //       call getarg(1,mdfile)
        //    //getarg_(&c__1, modpar_1.mdfile, (ftnlen)80);
        //    //       if(iargc().eq.2) then
        //    //   if (iargc_() == 2)
        //    //   {
        //    //          call getarg(2,modo)
        //    //getarg_(&c__2, modo, (ftnlen)1);

        //    // }   //       endif
        //    //       if(iargc().eq.4.or.iargc().eq.5) then
        //    //  if (iargc_() == 4 || iargc_() == 5) 
        //    //  {
        //    //          call getarg(2,dumm)
        //    //	getarg_(&c__2, dumm, (ftnlen)80);
        //    //          read(dumm,102) rlato
        //    //  102     format(f10.0)
        //    //          call getarg(3,dumm)
        //    //	getarg_(&c__3, dumm, (ftnlen)80);
        //    //          read(dumm,102) rlong
        //    //          call getarg(4,modo)
        //    //	getarg_(&c__4, modo, (ftnlen)1);
        //    //          if(iargc().eq.5) then
        //    //	if (iargc_() == 5)
        //    //{
        //    //             ispoly = .true.
        //    //  polinf_1.ispoly = true;
        //    //             call getarg(5,polyf)
        //    //	    getarg_(&c__5, polinf_1.polyf, (ftnlen)80);
        //    //          endif
        //    //	}
        //    //       endif
        //    //    }
        //    //       if(iargc().ge.6) then
        //    //    if (iargc_() >= 6) {
        //    //          call getarg(2,dumm)
        //    //	getarg_(&c__2, dumm, (ftnlen)80);
        //    //          read(dumm,102) slat
        //    //          call getarg(3,dumm)
        //    //	getarg_(&c__3, dumm, (ftnlen)80);
        //    //          read(dumm,102) rnlat
        //    //          call getarg(4,dumm)
        //    //	getarg_(&c__4, dumm, (ftnlen)80);
        //    //          read(dumm,102) ftlong
        //    //          call getarg(5,dumm)
        //    //	getarg_(&c__5, dumm, (ftnlen)80);
        //    //          read(dumm,102) rilong
        //    //          call getarg(6,modo)
        //    //	getarg_(&c__6, modo, (ftnlen)1);
        //    //          if(iargc().eq.7) then
        //    //	if (iargc_() == 7) {
        //    //             ispoly = .true.
        //    //	    polinf_1.ispoly = TRUE_;
        //    //             call getarg(7,polyf)
        //    //	    getarg_(&c__7, polinf_1.polyf, (ftnlen)80);
        //    //          endif
        //    //	}
        //    //       endif
        //    //    }

        //    //  done with getting arguments

        //    for (int num = 0; num < greensFunction.GrennsFunctionItems.Length; num++)// while (num < ntot) //L3:
        //    {
        //        if (modo == 'd')        //       if(modo.eq.'d') then
        //        {
        //            fingrd = 'C'; // 	 fingrd = 'C'
        //            ocmodl(stloc.rlat, stloc.rlam,
        //                ref oceanModel,
        //                ref greensFunction.GrennsFunctionItems[num]);
        //            //    rlato, rlong, fingrd, amp);       //          call ocmodl(rlato,rlong,fingrd,amp)


        //            if (amp != Complex.Zero) //          if(amp.ne.cz) amp=amp/rho
        //                amp /= seaDensity.rho;


        //            int ind = oceanModel.celfnd(rlong, rlato,
        //                        ref x, ref y, ref  i, ref j);      //          call celfnd(rlong,rlato,x,y,i,j,ind)
        //            // 	 write(6,104) i,j,ind,x,y
        //            //  104  

        //            if (amp != Complex.Zero) // 	 if(amp.ne.cz) write(6,106) amp,phasor(amp)
        //            {
        //            }

        //            if (amp == Complex.Zero) // 	 if(amp.eq.cz) write(6,107)
        //            {
        //            }
        //            //  106  
        //            //  107     format('No ocean in this cell')
        //            fingrd = 'F';       // 	 fingrd = 'F'
        //            ocmodl(stloc.rlat, stloc.rlam,
        //                ref oceanModel,
        //                ref greensFunction.GrennsFunctionItems[num]);
        //            // rlato, rlong, fingrd, amp);       //          call ocmodl(rlato,rlong,fingrd,amp)
        //            if (amp != Complex.Zero) //          if(amp.ne.cz) amp=amp/rho
        //                amp /= seaDensity.rho;

        //            if (amp != Complex.Zero)     // 	 if(amp.ne.cz) write(6,108) amp,phasor(amp)
        //            {
        //            }

        //            if (amp == Complex.Zero)        // 	 if(amp.eq.cz) write(6,109)
        //            {
        //            }
        //            //  108  
        //            //  109     format('Interpolated tide is zero')

        //            if (polInf != null)    // 	 if(ispoly) then
        //            {
        //                bool iok = polInf.chkgon(rlong, rlato);           // 	    call chkgon(rlong,rlato,iok)
        //                if (!iok)// == 0)   // 	   if(iok.eq.0) write(6,110)
        //                {
        //                }
        //                //  110       format('Point would be excluded by polygon')
        //                //          endif
        //            }

        //        }
        //        else if (modo == 'o')         //       elseif(modo.eq.'o') then
        //        {

        //            fingrd = 'F';       //          fingrd = 'F'
        //            ocmodl(stloc.rlat, stloc.rlam,
        //                ref oceanModel,
        //                ref greensFunction.GrennsFunctionItems[num]);
        //            // rlato, rlong, fingrd, amp);   //          call ocmodl(rlato,rlong,fingrd,amp)
        //            if (amp != Complex.Zero) //          if(amp.ne.cz) amp=amp/rho
        //                amp /= seaDensity.rho;

        //            if (polInf != null)        // 	 if(ispoly) then
        //            {
        //                bool iok = polInf.chkgon(rlong, rlato);           // 	    call chkgon(rlong,rlato,iok)
        //                if (!iok)// == 0)       // 	    if(iok.eq.0) amp = cz
        //                    amp = Complex.Zero;
        //            }   //          endif

        //            if (amp == Complex.Zero)         // 	 if(amp.eq.cz) then
        //            {
        //                //             write(6,119)
        //                //  119        format('Tide is zero or excluded by polygon')

        //            }       // 	 endif

        //            if (amp != Complex.Zero) // 	 if(amp.ne.cz) then
        //            {
        //                //             write(6,121) 'S',stnam, rlato, rlong, ht
        //                //  121        format(a1,2x,a40,3x,f10.4,f10.4,f10.0)
        //                //             write(6,123) 'O',dsym,(icte(j),j=1,6),mdnam
        //                //  123        format(a1,1x,a4,5x,6i2,5x,a50)
        //                //             write(6,125)
        //                //  125        format('L g          Phases are Greenwich, lags negative')
        //                //             write(6,127)
        //                //  127        format('X')

        //                amp = Complex.Conjugate(amp);       //             amp = conjg(amp)
        //                amp = phasor(amp); //             amp = phasor(amp)
        //                //             write(6,129) amp
        //                //  129        format('o',9x,2f10.4)

        //            }       // 	 endif

        //            /*  for options g a p r i, this section outputs either the grid cells */
        //            /*    or the tidal values in the cells, within the region specified on */
        //            /*    the command line */
        //            /*  for option e, this section outputs the edges of the boundary between */
        //            /*    zero and nonzero cells */


        //        }
        //        else //       else
        //        {
        //            ocmodl(0.0, 0.0,
        //                ref oceanModel,
        //                ref greensFunction.GrennsFunctionItems[num]);
        //            //     0.0, 0.0, fingrd, amp);       //          call ocmodl(0.,0.,fingrd,amp)

        //            if (amp != Complex.Zero)     //          if(amp.ne.cz) amp=amp/rho
        //                amp /= seaDensity.rho;

        //            if (modo == 'e')        //          if(modo.eq.'e') then
        //            {
        //                slat = modlim.blat;       //             slat=blat
        //                rnlat = modlim.tlat;      //             rnlat=tlat
        //                ftlong = modlim.wlong;        //             ftlong=wlong
        //                rilong = modlim.elong;        //             rilong=elong  
        //            }       //          endif

        //            if (rilong < ftlong)        // 	 if(rilong.lt.ftlong) rilong=rilong+360
        //                rilong += 360;


        //            fndgde(ftlong, rilong, slat, rnlat, ref il, ref ih, ref jl, ref jh);       //          call fndgde(ftlong,rilong,slat,rnlat,il,ih,jl,jh)

        //            if (il == 0 && ih == 0) //          if(il.eq.0.and.ih.eq.0) then
        //            {
        //                //             write(6,132)
        //                //  132        format('Area is outside model')

        //                //  	    stop

        //            }       //          endif
        //            // following code takes care of case where longitudes specified span edge of model

        //            if (ih < il)    // 	 if(ih.lt.il) then
        //            {
        //                ngrs = 2;   // 	    ngrs = 2
        //                i1[0] = 1;  // 	    i1(1) = 1
        //                i2[0] = ih; // 	    i2(1) = ih
        //                i1[1] = il; // 	    i1(2) = il
        //                i2[1] = modlim.longc; // 	    i2(2) = longc
        //            }
        //            else if (ih == il) // 	 elseif(ih.eq.il) then
        //            {
        //                ngrs = 1;       // 	    ngrs = 1
        //                i1[0] = 1;      // 	    i1(1) = 1
        //                i2[0] = modlim.longc;       // 	    i2(1) = longc
        //            }
        //            else                // 	 else
        //            {
        //                ngrs = 1;       // 	    ngrs = 1
        //                i1[0] = il;     // 	    i1(1) = il
        //                i2[0] = ih;     // 	    i2(1) = ih
        //            }   // 	 endif

        //            double dlat = (modlim.tlat - modlim.blat) / modlim.latc; // 	 dlat = (tlat-blat)/latc
        //            double dlong = (modlim.elong - modlim.wlong) / modlim.longc; // 	 dlong = (elong-wlong)/longc

        //            /* now scan through cells by columns, then rows; for option e output
        //             * any EW-running edges between zero and nonzero, or at top and bottom */
        //            for (int k = 0; k < ngrs; k++)     // 	 do k=1,ngrs
        //            {
        //                for (int i = i1[k]; i < i2[k]; i++)  // 	   do i=i1(k),i2(k)
        //                {

        //                    for (int j = jl; j <= jh; ++j)      // 	     do j=jl,jh
        //                    {
        //                        double clat = modlim.tlat - (j - 1) * dlat - dlat / 2;      // 	       clat = tlat - (j-1)*dlat - dlat/2
        //                        double clong = modlim.wlong + (i - 1) * dlong + dlong / 2;      // 	       clong = wlong + (i-1)*dlong + dlong/2

        //                        if (clong > 180.0)      // 	       if(clong.gt.180) clong = clong - 360
        //                            clong += -360.0;

        //                        ocmodl(clat, clong,
        //                                ref oceanModel,
        //                                ref greensFunction.GrennsFunctionItems[num]);
        //                        //   clat, clong, fingrd, amp);        //                call ocmodl(clat,clong,fingrd,amp)
        //                        if (amp != Complex.Zero)        //                if(amp.ne.cz) amp=amp/rho
        //                            amp = amp / seaDensity.rho;

        //                        if (polInf != null)        // 	       if(ispoly) then
        //                        {
        //                            bool iok = polInf.chkgon(clong, clat);           //                  call chkgon(clong,clat,iok)
        //                            if (!iok)// == 0)               // 	         if(iok.eq.0) amp=cz
        //                                amp = Complex.Zero;

        //                        }   //                endif

        //                        /* write out EW edges between zero and nonzero cells, including the case
        //                         * where the zero cell is one outside the grid */
        //                        if (modo == 'e')        // 	       if(modo.eq.'e') then
        //                        {
        //                            if (amp == Complex.Zero && oamp != Complex.Zero) //                   if(amp.eq.cz.and.oamp.ne.cz) then
        //                            {
        //                                //                      write(6,142) clong-dlong/2,clat+dlat/2
        //                                //                      write(6,142) clong+dlong/2,clat+dlat/2
        //                                //                      write(6,*) '370 100'
        //                            }       //                   endif
        //                            if (amp != Complex.Zero && (j == jl || oamp == Complex.Zero)) //                   if(amp.ne.cz.and.(j.eq.jl.or.oamp.eq.cz)) then
        //                            {
        //                                //                      write(6,142) clong-dlong/2,clat+dlat/2
        //                                //                      write(6,142) clong+dlong/2,clat+dlat/2
        //                                //                      write(6,*) '370 100'
        //                            }   //                   endif

        //                            if (amp != Complex.Zero && j == jh) //                   if(amp.ne.cz.and.j.eq.jh) then
        //                            {
        //                                //                      write(6,142) clong-dlong/2,clat+dlat/2
        //                                //                      write(6,142) clong+dlong/2,clat+dlat/2
        //                                //                      write(6,*) '370 100'
        //                            }   //                   endif
        //                        }   //                endif

        //                        // all other options, for writing out result in, or coordinates of, individual grid cells
        //                        if (amp != Complex.Zero)        // 	       if(amp.ne.cz) then
        //                        {
        //                            toloc = phasor(amp);        //                   toloc = phasor(amp)
        //                            if (modo == 'r')  // 	          if(modo.eq.'r') write(6,140) clong,clat,real(amp)
        //                            {
        //                            }
        //                            else if (modo == 'i') // 	          if(modo.eq.'i') write(6,140) clong,clat,aimag(amp)
        //                            {
        //                            }
        //                            else if (modo == 'a')      // 	          if(modo.eq.'a') write(6,140) clong,clat,real(toloc)
        //                            {
        //                            }
        //                            else if (modo == 'p')        // 	          if(modo.eq.'p') write(6,140) clong,clat,aimag(toloc)
        //                            {
        //                            }
        //                            //  140              format(2f10.4,g14.6)
        //                            else if (modo == 'g')    // 	          if(modo.eq.'g') then
        //                            {
        //                                //                      write(6,142) clong-dlong/2,clat-dlat/2
        //                                //                      write(6,142) clong-dlong/2,clat+dlat/2
        //                                //                      write(6,142) clong+dlong/2,clat+dlat/2
        //                                //                      write(6,142) clong+dlong/2,clat-dlat/2
        //                                //                      write(6,142) clong-dlong/2,clat-dlat/2
        //                                //                      write(6,*) '370 100'
        //                                //  142                 format(2f10.4)

        //                            }   //                   endif
        //                        }                   // 	       endif
        //                        oamp = amp;         // 	       oamp = amp
        //                    }                       //              enddo
        //                }                           //            enddo
        //            }                               //          enddo   

        //            // for the e option, do a second scan, by rows, to find all the NS-running edges
        //            if (modo == 'e')                        // 	 if(modo.eq.'e') then
        //            {
        //                oamp = Complex.Zero;                //             oamp=cz
        //                for (int k = 1; k < ngrs; k++)         // 	    do k=1,ngrs
        //                {
        //                    for (int j = jl; j < jh; j++)                // 	      do j=jl,jh
        //                    {
        //                        for (int i = i1[k]; i < i2[k]; i++) // 	        do i=i1(k),i2(k)
        //                        {
        //                            double clat = modlim.tlat - (j - 1) * dlat - dlat / 2;           // 	           clat = tlat - (j-1)*dlat - dlat/2
        //                            double clong = modlim.wlong + (i - 1) * dlong + dlong / 2;     // 	           clong = wlong + (i-1)*dlong + dlong/2

        //                            if (clong > 180.0)          // 	           if(clong.gt.180) clong = clong - 360
        //                                clong += -360.0;

        //                            ocmodl(clat, clong,
        //                                ref oceanModel,
        //                                ref greensFunction.GrennsFunctionItems[num]);
        //                            // clat, clong, fingrd, amp);     //                    call ocmodl(clat,clong,fingrd,amp)

        //                            if (amp != Complex.Zero && oamp != Complex.Zero)     //                    if(amp.eq.cz.and.oamp.ne.cz) then
        //                            {
        //                                //                       write(6,142) clong-dlong/2,clat+dlat/2
        //                                //                       write(6,142) clong-dlong/2,clat-dlat/2
        //                                //                       write(6,*) '370 100'

        //                            }               //                    endif
        //                            if (amp != Complex.Zero && (i2[k] == 1 || oamp == Complex.Zero)) // TODO check == 1 oder == 0   //                    if(amp.ne.cz.and.(i.eq.i1(k).or.oamp.eq.cz)) then
        //                            {
        //                                //                       write(6,142) clong-dlong/2,clat+dlat/2
        //                                //                       write(6,142) clong-dlong/2,clat-dlat/2
        //                                //                       write(6,*) '370 100'

        //                            }               //                    endif
        //                            if (amp != Complex.Zero && i2[k] == 1)     // TODO check == 1 oder == 0         //                    if(amp.ne.cz.and.i.eq.i2(k)) then
        //                            {
        //                                //s_wsfe(&io___79);                               //                       write(6,142) clong+dlong/2,clat+dlat/2
        //                                //s_wsfe(&io___80);                                   //                       write(6,142) clong+dlong/2,clat-dlat/2
        //                                //s_wsle(&io___81);                                   //                       write(6,*) '370 100'

        //                            }                   //                    endif
        //                            oamp = amp;         // 	           oamp = amp
        //                        }                       //                 enddo
        //                    }                           //               enddo
        //                }                               //             enddo
        //            }                                   // 	 endif
        //        }                                       //       endif

        //        // s_stop("", (ftnlen)0);      //       stop
        //    }
        //    return;               //       end
        //}


        //private void fndgde(double ftlon, double rilon, double slat, double rnlat, ref int il, ref int ih, ref int jl, ref int jh)   //       subroutine fndgde(ftlon,rilon,slat,rnlat,il,ih,jl,jh)
        //{
        //    /*     Revision 1.1  2011/09/23 22:39:01  agnew  - Initial revision */

        //    /*  returns the cell indices for the corners of a grid
        //     * ftlon,rilon are the left and right longitudes wanted,
        //     * slat and rnlat the south and north limits
        //     * returns the cell indices that cover these, allowing for
        //     * the case when the limits given are outside the model; will not
        //     * give the full limits of indices when the longitudes cross the
        //     * edge of the model. */

        //    /*  put longitudes into model range */
        //    double dlongl = (ftlon - modlim.wlong) % 360.0;     //       dlongl = amod(ftlon-wlong,360.)
        //    if (dlongl < 0.0)                                   //       if(dlongl.lt.0) dlongl = dlongl + 360
        //        dlongl += 360;
            
        //    double dlongr = (rilon - modlim.wlong) % 360.0;     //       dlongr = amod(rilon-wlong,360.)
        //    if (dlongr < 0.0)                                   //       if(dlongr.lt.0) dlongr = dlongr + 360
        //        dlongr += 360;
           
        //    /* allow for limits completely outside the model */
        //    if (slat > modlim.tlat || 
        //        rnlat < modlim.blat || 
        //        dlongl > modlim.elong - modlim.wlong ||
        //        dlongl > modlim.elong -
        //        modlim.wlong)                                   //      if(slat.gt.tlat.or.rnlat.lt.blat.or.dlongl.gt.elong-wlong.or.dlongl.gt.elong-wlong) then       
        //    {
        //        il = 0;                                         //         il=0
        //        ih = 0;                                         // 	ih=0
        //        jl = 0;                                         // 	jl=0
        //        jh = 0;                                         // 	jh=0
        //    }                                                   //       endif

        //    double x = modlim.longc * dlongl / 
        //        (modlim.elong - modlim.wlong);                  //       x = (longc*dlongl)/(elong-wlong)
            
        //    double y = modlim.latc * (modlim.tlat - rnlat) / 
        //        (modlim.tlat - modlim.blat);                    //       y = (latc*(tlat-rnlat))/(tlat-blat)

        //    il = (int)x + 1;                                    //       il= int(x) + 1
        //    jl = (int)y + 1;                                    //       jl= int(y) + 1 

        //    if (rnlat > modlim.tlat)                            //       if(rnlat.gt.tlat) jl = 1
        //        jl = 1;

        //    x = modlim.longc * dlongr / 
        //        (modlim.elong - modlim.wlong);                  //       x = (longc*dlongr)/(elong-wlong)
            
        //    y = modlim.latc * (modlim.tlat - slat) / 
        //        (modlim.tlat - modlim.blat);                    //       y = (latc*(tlat-slat))/(tlat-blat)

        //    ih = (int)x + 1;                                    //       ih= int(x) + 1
        //    jh = (int)y + 1;                                    //       jh= int(y) + 1

        //    if (slat < modlim.blat)                             //       if(slat.lt.blat) jh = latc
        //        jh = modlim.latc;
            
        //    /* if NS edge of model is in grid, ih will be le il (though both are > 0) */
        //}                                                       //       end

        private void loadcomb(ref OceanLoadWave nloadfWave)
        {
            if (nloadfWave.oceanLoad == null)
                nloadfWave.oceanLoad = new OceanLoad();

            Complex g = argand(new Complex(nloadfWave.oceanLoad.g.Amplitude, nloadfWave.oceanLoad.g.Phase));
            Complex gdef = argand(new Complex(nloadfWave.oceanLoad.gdef.Amplitude, nloadfWave.oceanLoad.gdef.Phase));
            Complex p = argand(new Complex(nloadfWave.oceanLoad.p.Amplitude, nloadfWave.oceanLoad.p.Phase));
            Complex d0 = argand(new Complex(nloadfWave.oceanLoad.d[0].Amplitude, nloadfWave.oceanLoad.d[0].Phase));
            Complex d1 = argand(new Complex(nloadfWave.oceanLoad.d[1].Amplitude, nloadfWave.oceanLoad.d[1].Phase));
            Complex d2 = argand(new Complex(nloadfWave.oceanLoad.d[2].Amplitude, nloadfWave.oceanLoad.d[2].Phase));
            Complex t0 = argand(new Complex(nloadfWave.oceanLoad.t[0].Amplitude, nloadfWave.oceanLoad.t[0].Phase));
            Complex t1 = argand(new Complex(nloadfWave.oceanLoad.t[1].Amplitude, nloadfWave.oceanLoad.t[1].Phase));
            Complex tdef0 = argand(new Complex(nloadfWave.oceanLoad.tdef[0].Amplitude, nloadfWave.oceanLoad.tdef[0].Phase));
            Complex tdef1 = argand(new Complex(nloadfWave.oceanLoad.tdef[1].Amplitude, nloadfWave.oceanLoad.tdef[1].Phase));
            Complex s0 = argand(new Complex(nloadfWave.oceanLoad.s[0].Amplitude, nloadfWave.oceanLoad.s[0].Phase));
            Complex s1 = argand(new Complex(nloadfWave.oceanLoad.s[1].Amplitude, nloadfWave.oceanLoad.s[1].Phase));
            Complex s2 = argand(new Complex(nloadfWave.oceanLoad.s[2].Amplitude, nloadfWave.oceanLoad.s[2].Phase));
            Complex o = argand(new Complex(nloadfWave.oceanLoad.o.Amplitude, nloadfWave.oceanLoad.o.Phase));

            for (int i = 0; i < nloadfWave.nloadfItems.Length; i++)
            {
                g += argand(new Complex(nloadfWave.nloadfItems[i].g.Amplitude, nloadfWave.nloadfItems[i].g.Phase));
                gdef += argand(new Complex(nloadfWave.nloadfItems[i].gdef.Amplitude, nloadfWave.nloadfItems[i].gdef.Phase));
                p += argand(new Complex(nloadfWave.nloadfItems[i].p.Amplitude, nloadfWave.nloadfItems[i].p.Phase));
                d0 += argand(new Complex(nloadfWave.nloadfItems[i].d[0].Amplitude, nloadfWave.nloadfItems[i].d[0].Phase));
                d1 += argand(new Complex(nloadfWave.nloadfItems[i].d[1].Amplitude, nloadfWave.nloadfItems[i].d[1].Phase));
                d2 += argand(new Complex(nloadfWave.nloadfItems[i].d[2].Amplitude, nloadfWave.nloadfItems[i].d[2].Phase));
                t0 += argand(new Complex(nloadfWave.nloadfItems[i].t[0].Amplitude, nloadfWave.nloadfItems[i].t[0].Phase));
                t1 += argand(new Complex(nloadfWave.nloadfItems[i].t[1].Amplitude, nloadfWave.nloadfItems[i].t[1].Phase));
                tdef0 += argand(new Complex(nloadfWave.nloadfItems[i].tdef[0].Amplitude, nloadfWave.nloadfItems[i].tdef[0].Phase));
                tdef1 += argand(new Complex(nloadfWave.nloadfItems[i].tdef[1].Amplitude, nloadfWave.nloadfItems[i].tdef[1].Phase));
                s0 += argand(new Complex(nloadfWave.nloadfItems[i].s[0].Amplitude, nloadfWave.nloadfItems[i].s[0].Phase));
                s1 += argand(new Complex(nloadfWave.nloadfItems[i].s[1].Amplitude, nloadfWave.nloadfItems[i].s[1].Phase));
                s2 += argand(new Complex(nloadfWave.nloadfItems[i].s[2].Amplitude, nloadfWave.nloadfItems[i].s[2].Phase));
                o += argand(new Complex(nloadfWave.nloadfItems[i].o.Amplitude, nloadfWave.nloadfItems[i].o.Phase));
            }

            nloadfWave.oceanLoad.g = new Load(phasor(g));
            nloadfWave.oceanLoad.gdef = new Load(phasor(gdef));
            nloadfWave.oceanLoad.p = new Load(phasor(p));
            nloadfWave.oceanLoad.d[0] = new Load(phasor(d0));
            nloadfWave.oceanLoad.d[1] = new Load(phasor(d1));
            nloadfWave.oceanLoad.d[2] = new Load(phasor(d2));
            nloadfWave.oceanLoad.t[0] = new Load(phasor(t0));
            nloadfWave.oceanLoad.t[1] = new Load(phasor(t1));
            nloadfWave.oceanLoad.tdef[0] = new Load(phasor(tdef0));
            nloadfWave.oceanLoad.tdef[1] = new Load(phasor(tdef1));
            nloadfWave.oceanLoad.s[0] = new Load(phasor(s0));
            nloadfWave.oceanLoad.s[1] = new Load(phasor(s1));
            nloadfWave.oceanLoad.s[2] = new Load(phasor(s2));
            nloadfWave.oceanLoad.o = new Load(phasor(o));
        }

        public void writenfloadOut(string path)
        {
            FileStream fs = null;

            if (!File.Exists(path))
                fs = new FileStream(path, FileMode.Create);
            else
                fs = new FileStream(path, FileMode.Append);

            StreamWriter wr = new StreamWriter(fs);

            for (int i = 0; i < oceanLoadingParameter.waveItems.Length; i++)
            {
                for (int j = 0; j < oceanLoadingParameter.waveItems[i].nloadfItems.Length; j++)
                {
                    for (int k = 0; k < oceanLoadingParameter.waveItems[i].nloadfItems[j].header.Length; k++)
                        wr.WriteLine(oceanLoadingParameter.waveItems[i].nloadfItems[j].header[k]);
                    wr.WriteLine(String.Format(Constants.NumberFormatEN, "{0}{1,19:0.0000}{2,10:0.0000}", "g",
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].g.Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].g.Phase));
                    wr.WriteLine(String.Format(Constants.NumberFormatEN, "{0}{1,19:0.0000}{2,10:0.0000}", "p",
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].p.Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].p.Phase));
                    wr.WriteLine(String.Format(Constants.NumberFormatEN, "{0}{1,19:0.0000}{2,10:0.0000}{3,10:0.0000}{4,10:0.0000}{5,10:0.0000}{6,10:0.0000}", "d",
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].d[0].Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].d[0].Phase,
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].d[1].Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].d[1].Phase,
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].d[2].Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].d[2].Phase));
                    wr.WriteLine(String.Format(Constants.NumberFormatEN, "{0}{1,19:0.0000}{2,10:0.0000}{3,10:0.0000}{4,10:0.0000}", "t",
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].t[0].Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].t[0].Phase,
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].t[1].Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].t[1].Phase));
                    wr.WriteLine(String.Format(Constants.NumberFormatEN, "{0}{1,19:0.0000}{2,10:0.0000}{3,10:0.0000}{4,10:0.0000}{5,10:0.0000}{6,10:0.0000}", "s",
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].s[0].Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].s[0].Phase,
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].s[1].Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].s[1].Phase,
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].s[2].Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].s[2].Phase));
                    wr.WriteLine(String.Format(Constants.NumberFormatEN, "{0}{1,19:0.0000}{2,10:0.0000}", "o",
                        oceanLoadingParameter.waveItems[i].nloadfItems[j].o.Amplitude, oceanLoadingParameter.waveItems[i].nloadfItems[j].o.Phase));
                }
            }

            wr.Flush();
            wr.Close();
        }

        private void ldbxdr(double azimuth, double distance, double azstp, double disstp)
        {
            /*     Revision 1.1  2011/09/23 22:39:01  agnew  -  Initial revision */

            /* outputs the coordinates (longitude-latitude) of a load cell centered 
             * at given azimuth and distance, with dimensions azstp and disstp. The 
             * station location is passed to invspt in a common block */

            double[,] xx = new double[2, 2];
            double[,] yy = new double[2, 2];
            double b = 0;
            double rlong = 0;

            for (int i = 1; i <= 2; i++)                                //       do 7 i=1,2
            {
                for (int j = 1; j <= 2; j++)                            //       do 5 j=1,2
                {
                    int idelsg = 2 * i - 3;                             //       idelsg = 2*i - 3
                    int ialpsg = 2 * j - 3;                             //       ialpsg = 2*j - 3
                    double cornaz = azimuth + ialpsg * (azstp / 2.0);   //       cornaz = azimuth+ialpsg*(azstp/2.)
                    double rad = distance + idelsg * (disstp / 2.0);    //       rad = distance+idelsg*(disstp/2.)
                    invspt(cornaz, rad, ref b, ref rlong);              //       call invspt(cornaz,rad,b,rlong)
                    yy[i, j] = 90 - b;                                  //       yy(i,j) = 90-b
                    xx[i, j] = rlong;                                   //       xx(i,j) = rlong
                }                                                       //   5   continue
            }                                                           //   7   continue
        }

        /// <summary>
        /// Solves the inverse spherical triangle problem. 'stloc' holds the needed station parameter.
        /// </summary>
        /// <param name="alp"></param>
        /// <param name="del"></param>
        /// <param name="b">Returns the colatitude in degree.</param>
        /// <param name="rlong">Returns the east longitude of the place which is at a distance of delta degrees and at an azimuth of alp (clockwise from north).</param>
        private void invspt(double alp, double del, ref double b, ref double rlong)
        {
            /*     Revision 1.1  2011/09/23 22:39:01  agnew - Initial revision */

            /*  solves the inverse spherical triangle problem - given a station whose colatitude is t, and east longitude rlam (in degrees), returns the
             *  colatitude b and east longitude rlong of the place which is at a distance of delta degrees and at an azimuth of alp (clockwise from north). 
             *  the common block stloc holds the cosine and sine of the station colatitude, and its east longitude. */

            double ca = Math.Cos(alp * Constants.degree2radian);        //       ca = cos(alp*dtr)
            double sa = Math.Sin(alp * Constants.degree2radian);        //       sa = sin(alp*dtr)
            double cd = Math.Cos(del * Constants.degree2radian);        //       cd = cos(del*dtr)
            double sd = Math.Sin(del * Constants.degree2radian);        //       sd = sin(del*dtr)
            double cb = cd * stloc.ct + sd * stloc.st * ca;             //       cb = cd*ct + sd*st*ca
            double sb = Math.Sqrt(1.0 - cb * cb);                       //       sb = sqrt(1.-cb*cb)
            b = Math.Acos(cb) / Constants.degree2radian;                                    //       b = acos(cb)/dtr
            if (sb <= .001)                                             //       if(sb.le.1.e-3) then
            {
                /*  special case - the point is at the poles */
                rlong = 0.0;                                            // 	        rlong = 0
                return;                                                 // 	        return
            }                                                           //       endif
            double sg = sd * sa / sb;                                   //       sg = sd*sa/sb
            double cg = (stloc.st * cd - sd * stloc.ct * ca) / sb;      //       cg = (st*cd-sd*ct*ca)/sb
            double g = Math.Atan2(sg, cg) / Constants.degree2radian;                        //       g = atan2(sg,cg)/dtr
            rlong = stloc.rlam + g;                                     //       rlong = rlam + g
            return;                                                     //       return
        }

        private OceanLoad lodout(Complex[] cload, ModeOfOperation modo, ref GreensFunction greensFunction, ref OceanModel oceanModel, ref PolInf polInf) //       subroutine lodout(cload,stnam,rlat,rlam,ht,modo)
        {
            /*     Revision 1.7  2013/03/11 02:21:04  agnew  -  updated version number to 3.3.0.2
             *     Revision 1.6  2012/06/25 17:14:29  agnew  -  updated version number to 3.3.0.1
             *     Revision 1.5  2012/03/13 21:19:15  agnew  -  Modified polygon common block and output polygon names
             *     Revision 1.4  2012/03/06 03:36:28  agnew  -  added lines for new grid options
             *     Revision 1.3  2011/12/30 23:29:19  agnew  -  changed version number given in output, to 3.3.0 
             *     Revision 1.2  2011/12/30 23:03:36  agnew  -  removed unneeded statement numbers 
             *     Revision 1.1  2011/09/23 22:39:01  agnew  -  Initial revision */

            /* puts the data into a single ASCII rescaled file for all constituents, in order of frequency. Loads are written out in
             * amp and local or Greenwich phase, lags neg */

            OceanLoad nout = new OceanLoad();
            ArrayList header = new ArrayList();

            header.Add("S  " +
                stloc.Name.PadRight(43) +
                stloc.rlat.ToString("0.0000", Constants.NumberFormatEN).PadLeft(10) +
                stloc.rlam.ToString("0.0000", Constants.NumberFormatEN).PadLeft(10) +
                stloc.ht.ToString("0", Constants.NumberFormatEN).PadLeft(9) +
                ".");            //       write(luo,201) 'S',stnam, rlat, rlam, ht | //  201  format(a1,2x,a40,3x,f10.4,f10.4,f10.0)

            header.Add(string.Format(Constants.NumberFormatEN, "O {0}{1,2:0}{2,2:0}{3,2:0}{4,2:0}{5,2:0}{6,2:0}     {7}",
                oceanModel.dsym.PadRight(9),
                oceanModel.icte[0], oceanModel.icte[1], oceanModel.icte[2], oceanModel.icte[3], oceanModel.icte[4], oceanModel.icte[5],
                oceanModel.mdnam));                                                          //       write(luo,203) 'O',dsym,(icte(j),j=1,6),mdnam | //  203  format(a1,1x,a4,5x,6i2,5x,a50)
            header.Add("G  " + greensFunction.grname);                                      //       write(luo,205) 'G',grname | //  205  format(a1,2x,a78)

            for (int i = 0; i < greensFunction.GrennsFunctionItems.Length; i++)            //       do i=1,nring
            {
                if (greensFunction.GrennsFunctionItems[i].fingrd == GreensFunction.FinGrd.F)                    //         if(statgr(i).eq.'F') write(luo,211) rin(i),rout(i),rsz(i)
                    header.Add(string.Format("G    Rings from {0,6} to {1,6} spacing {2,6} {3}",
                                             greensFunction.GrennsFunctionItems[i].beg.ToString("0.00", Constants.NumberFormatEN),
                                             greensFunction.GrennsFunctionItems[i].end.ToString("0.00", Constants.NumberFormatEN),
                                             greensFunction.GrennsFunctionItems[i].spc.ToString("0.0000", Constants.NumberFormatEN),
                                             " - detailed grid (ocean), seawater"));//211    format('G    Rings from ',f6.2,' to ',f6.2,' spacing ', f4.2,' - detailed grid (ocean), seawater')

                if (greensFunction.GrennsFunctionItems[i].fingrd == GreensFunction.FinGrd.L)                    //         if(statgr(i).eq.'L') write(luo,213) rin(i),rout(i),rsz(i)
                    header.Add(string.Format("G    Rings from {0,6} to {1,6} spacing {2,6} {3}",
                                             greensFunction.GrennsFunctionItems[i].beg.ToString("0.00", Constants.NumberFormatEN),
                                             greensFunction.GrennsFunctionItems[i].end.ToString("0.00", Constants.NumberFormatEN),
                                             greensFunction.GrennsFunctionItems[i].spc.ToString("0.0000", Constants.NumberFormatEN),
                                             " - detailed grid (land), freshwater")); //213    format('G    Rings from ',f6.2,' to ',f6.2,' spacing ', f4.2,' - detailed grid (land), freshwater')

                if (greensFunction.GrennsFunctionItems[i].fingrd == GreensFunction.FinGrd.C)                    //         if(statgr(i).eq.'C') write(luo,215) rin(i),rout(i),rsz(i)
                    header.Add(string.Format("G    Rings from {0,6} to {1,6} spacing {2,6} {3}",
                                             greensFunction.GrennsFunctionItems[i].beg.ToString("0.00", Constants.NumberFormatEN),
                                             greensFunction.GrennsFunctionItems[i].end.ToString("0.00", Constants.NumberFormatEN),
                                             greensFunction.GrennsFunctionItems[i].spc.ToString("0.0000", Constants.NumberFormatEN),
                                             " - model grid, seawater"));             //215    format('G    Rings from ',f6.2,' to ',f6.2,' spacing ', f4.2,' - model grid, seawater')

                if (greensFunction.GrennsFunctionItems[i].fingrd == GreensFunction.FinGrd.G)                    //         if(statgr(i).eq.'G') write(luo,217) rin(i),rout(i),rsz(i)
                    header.Add(string.Format("G    Rings from {0,6} to {1,6} spacing {2,6} {3}",
                                             greensFunction.GrennsFunctionItems[i].beg.ToString("0.00", Constants.NumberFormatEN),
                                             greensFunction.GrennsFunctionItems[i].end.ToString("0.00", Constants.NumberFormatEN),
                                             greensFunction.GrennsFunctionItems[i].spc.ToString("0.0000", Constants.NumberFormatEN),
                                             " - model grid, freshwater"));           //217    format('G    Rings from ',f6.2,' to ',f6.2,' spacing ', f4.2,' - model grid, freshwater')
            }                                                                               //       enddo

            if (polInf != null)                                                             //       if(alpoly.ne.'N') then
            {
                //          write(luo,221) polnam  | //  221     format('P',1x,a78)
                for (int i = 0; i < polInf.np; i++)                               //          do i=1,np
                {
                    if (polInf.npoly[i].use)                                                //             if(use(i).eqv..true.) write(6,223) polynm(i) | //  223        format('P   Included polygon: ',a50)
                        header.Add("P   Included polygon: " + polInf.npoly[i].polynm);

                    if (!polInf.npoly[i].use)                                               //             if(use(i).eqv..false.) write(6,225) polynm(i) | //  225        format('P   Excluded polygon: ',a50)
                        header.Add("P   Excluded polygon: " + polInf.npoly[i].polynm);
                }                                                                           //          enddo
            }                                                                               //       endif

            /* following three lines depend on a specific routine, which returns time information in string   dates */

            //string dates = DateTime.Now.ToString();//       call fdate(dates)
            //       write(luo,231) dates
            //  231  format('C   Version 3.3.0.2 of load program, run at ',a24)

            header.Add("C   closest nonzero load was " +
                conInf.close.ToString("0.00", Constants.NumberFormatEN).PadLeft(6) +
                " degrees away, at" +
                conInf.clat.ToString("0.00", Constants.NumberFormatEN).PadLeft(8) +
                conInf.clong.ToString("0.00", Constants.NumberFormatEN).PadLeft(8)); //       write(luo,233) close,clat,clong | //  233  format('C   closest nonzero load was ',f6.2,' degrees away, at', 2f8.2)
            if (conInf.numnf != 0)                                                        //       if(numnf.ne.0) write(luo,235) numnf,clnf,farnf | //  235  format('C   ',i6,' zero loads found where ocean present, range', f7.2,'-',f7.2,' deg')
                header.Add("C   " +
                    conInf.numnf.ToString("0").PadLeft(6) +
                    " zero loads found where ocean present, range " +
                    conInf.clnf.ToString("0.00", Constants.NumberFormatEN) + " - " +
                    conInf.farnf.ToString("0.00", Constants.NumberFormatEN) + " deg");

            if (modo == ModeOfOperation.g)                                                                //       if(modo.eq.'g') write(luo,241) | //  241  format('L g          Phases are Greenwich, lags negative')
                header.Add("L g          Phases are Greenwich, lags negative");

            if (modo == ModeOfOperation.l)                                                                //       if(modo.eq.'l') write(luo,243) |  //  243  format('L l          Phases are local, lags negative')
                header.Add("L l          Phases are local, lags negative");

            header.Add("X");                                                                //       write(luo,249) | //  249  format('X')

            nout.header = (string[])header.ToArray(typeof(string));

            int isp = oceanModel.icte[0];     //       isp = icte(1)
            Complex toloc = new Complex(Math.Cos(Constants.degree2radian * isp * stloc.rlam),
                                Math.Sin(Constants.degree2radian * isp * stloc.rlam));  //       toloc = cmplx(cos(dr*isp*rlam),sin(dr*isp*rlam))

            if (modo == ModeOfOperation.g)                                                                //       if(modo.eq.'g') toloc = cmplx(1.,0.)
                toloc = new Complex(1.0, 0.0);

            Complex[] cmp = new Complex[3];

            //      gravity
            cmp[0] = 1.0e8 * toloc * cload[0];                                              //       cmp(1) = 1.e8*toloc*cload(1)
            cmp[0] = phasor(Complex.Conjugate(cmp[0]));                                     //       dumm = conjg(cmp(1)) | //       cmp(1) = phasor(dumm) 
            nout.g = new Load(cmp[0]);                                                                //       write(luo,251) cmp(1) | //  251  format('g',9x,2f10.4)
            nout.gdef = new Load(phasor(Complex.Conjugate(1.0e8 * cload[10])));

            //      potential height
            cmp[0] = 1.0e3 * toloc * cload[9];                                              //       cmp(1) = 1.e3*toloc*cload(10)
            cmp[0] = phasor(Complex.Conjugate(cmp[0]));                                     //       dumm = conjg(cmp(1)) | //      cmp(1) = phasor(dumm)
            nout.p = new Load(cmp[0]);                                                                //       write(luo,253) cmp(1) | //  253  format('p',9x,2f10.4)

            //      displacement (output as E N Z) - change sign of N and E (see CJ II-143, 5 Jan 82)
            cmp[0] = -1.0e3 * toloc * cload[3];                                             //       cmp(1) = -1.e3*toloc*cload(4)
            cmp[1] = -1.0e3 * toloc * cload[2];                                             //       cmp(2) = -1.e3*toloc*cload(3)
            cmp[2] = 1.0e3 * toloc * cload[1];                                              //       cmp(3) =  1.e3*toloc*cload(2)
            for (int i = 0; i < 3; i++)                                                     //       do i = 1,3
                cmp[i] = phasor(Complex.Conjugate(cmp[i]));                                 //         dumm = conjg(cmp(i)) | //         cmp(i) = phasor(dumm)
            nout.d = new Load[] { new Load(cmp[0]), new Load(cmp[1]), new Load(cmp[2]) };                              //       write(luo,255) cmp | //  255  format('d',9x,6f10.4)

            //      tilt (output as E N)
            cmp[0] = 1.0e9 * toloc * cload[5];                                              //       cmp(1) = 1.e9*toloc*cload(6)
            cmp[1] = 1.0e9 * toloc * cload[4];                                              //       cmp(2) = 1.e9*toloc*cload(5)
            for (int i = 0; i < 2; i++)                                                     //       do i = 1,2
                cmp[i] = phasor(Complex.Conjugate(cmp[i]));                                 //         dumm = conjg(cmp(i)) | //         cmp(i) = phasor(dumm)

            nout.t = new Load[] { new Load(cmp[0]), new Load(cmp[1]) };                                       //       write(luo,257) cmp(1),cmp(2) | //  257  format('t',9x,4f10.4)
            nout.tdef = new Load[] { new Load(phasor(Complex.Conjugate(1.0e9 * toloc * cload[11]))),
                                     new Load(phasor(Complex.Conjugate(1.0e9 * toloc * cload[12]))) };
            //      strain (output as E, N, EN shear)
            cmp[0] = 1.0e9 * toloc * cload[7];                                              //       cmp(1) = 1.e9*toloc*cload(8)
            cmp[1] = 1.0e9 * toloc * cload[6];                                              //       cmp(2) = 1.e9*toloc*cload(7)
            cmp[2] = 1.0e9 * toloc * cload[8];                                              //       cmp(3) = 1.e9*toloc*cload(9)
            for (int i = 0; i < 3; i++) //       do i = 1,3
                cmp[i] = phasor(Complex.Conjugate(cmp[i]));                                 //         dumm = conjg(cmp(i)) | //         cmp(i) = phasor(dumm)
            nout.s = new Load[] { new Load(cmp[0]), new Load(cmp[1]), new Load(cmp[2]) };                              //       write(luo,259) cmp | //  259  format('s',9x,6f10.4)

            // ocen tide [m]
            nout.o = new Load(cload[13]);

            return nout;                                                                    //       return
        }

        /// <summary>
        /// <para>For given latitude and east longitude (in degrees and ge 0) returns the complex amplitude of the tide (in kg/m**2, with greenwich phase).</para>
        /// <para>this amplitude is of course zero on land. If fingrd is "F", then the land-sea grid used is a detailed one;</para>
        /// <para>otherwise it is the one used by the ocean model.</para>
        /// </summary>
        /// <param name="rlong">East longitude</param>
        /// <param name="rlato">Latitude</param>
        /// <param name="TMe"></param>
        /// <param name="GFe"></param>
        /// <param name="lsM"></param>
        /// <param name="dist"></param>
        /// <param name="numnf"></param>
        /// <param name="clnf"></param>
        /// <param name="farnf"></param>
        /// <returns></returns>
        private Complex ocmodl(double rlato, double rlong, ref OceanModel oceanModel, ref GreensFunction.FinGrd fingrd)
        {
            /*     Revision 1.6  2012/06/25 01:29:40 : agnew : added density common block, corrected comments
             *     Revision 1.5  2012/03/13 21:17:32 : agnew : Modified polygon common block
             *     Revision 1.4  2012/03/06 03:38:09 : agnew : added new grid options, density, and call to seadens
             *     Revision 1.3  2011/11/27 04:31:57 : agnew : larger dimension for polygon information
             *     Revision 1.2  2011/11/18 16:52:49 : agnew : corrected error in interpolation (incorrect indexing of surrounding cells in one case) 
             *     Revision 1.1  2011/09/23 22:39:02 : agnew
             *     Initial revision */

            /*  For given latitude and east longitude (in degrees and ge 0) returns the complex amplitude of the tide (in kg/m**2, with greenwich phase).
            *  this amplitude is of course zero on land. If fingrd is "F", then the land-sea grid used is a detailed one;
            *  otherwise it is the one used by the ocean model.
            *  
            *  On the first call, this calls a separate routine (ocstart) which gets auxiliary information (passed to other routines in common) and 
            *  returns the two arrays of real and imaginary amplitude.  Later calls simply go to the appropriate point of the arrays.
            *
            *  The first call also reads in, if requested, information about polygonal areas (up to 10) within which the point either must or
            *  must not fall. */

            Complex amp = Complex.Zero; // return value

            /* ind is cell number; x and y are positions within cell, relative to center, in E and N, in fractions of a cell dimension */
            double x = 0;
            double y = 0;
            int i = 0;
            int j = 0;
            int ind = oceanModel.celfnd(rlong, rlato,
                                        ref x, ref y, ref  i, ref j);   //       call celfnd(rlong,rlato,x,y,i,j,ind)

            if (ind > oceanModel.latc * oceanModel.longc || ind <= 0)   //       if(ind.gt.latc*longc.or.ind.le.0) then
            {
                amp = Complex.Zero;                                     //          amp = cz
                return amp;                                             //          return
            }
            /* if cell is inside a polygon we don't want, or outside one we do, reject it without going further */

            if (polInf != null)                                         //       if(ispoly) then
            {
                bool iok = polInf.chkgon(rlong, rlato);                 // 	 call chkgon(rlong,rlato,iok)
                if (!iok)// == 0)                                       // 	 if(iok.eq.0) then
                {
                    amp = Complex.Zero;                                 //             amp = cz
                    return amp;                                         //          return
                }
            }

            double rho = 0.0;
            if (fingrd == GreensFunction.FinGrd.C || fingrd == GreensFunction.FinGrd.G)                 //       if(fingrd.eq.'C'.or.fingrd.eq.'G') then
            {
                /*  large distance - return the model value (0 if on land for tidal models) */
                amp = 0.001d * (new Complex(
                                oceanModel.ir1[ind - 1],
                                oceanModel.im1[ind - 1]));              //            amp = .001*cmplx(float(ir1(ind)),float(im1(ind)))

                if (fingrd == GreensFunction.FinGrd.C)                                  //          if(fingrd.eq.'C') call seadens(rlato,rlong,rho)
                    rho = seaDensity.seadens(rlato, rlong);

                if (fingrd == GreensFunction.FinGrd.G)                                  //          if(fingrd.eq.'G') rho = 1.e3
                    rho = 1e3;

                amp *= rho;                                             //          amp=rho*amp

                return amp;                                             //          return
            }

            if (fingrd == GreensFunction.FinGrd.F || fingrd == GreensFunction.FinGrd.L)                 //       if(fingrd.eq.'F'.or.fingrd.eq.'L') then
            {
                /* close distance, or land load - use the 5-minute map file to decide if there is land or not. If there is
                 * not (for F) or is (for L), bilinearly interpolate from the model, filling out empty cells to make this possible (though if all
                 * the cells needed are zero, return this). */
                bool lnd = landSeaMatrix.getLandSea(rlato, rlong);      //          call lndsea(rlato,rlong,lnd) | true (1) = Land; false (2) = Water

                if (lnd && fingrd == GreensFunction.FinGrd.F ||
                    !lnd && fingrd == GreensFunction.FinGrd.L)                          // if((lnd.eq.1.and.fingrd.eq.'F').or.(lnd.eq.2.and.fingrd.eq.'L')) then
                {
                    /*  detailed grid says land (for F) or sea (for L) - return zero */
                    amp = Complex.Zero;                                 //             amp = cz
                    return amp;                                         //             return
                }

                /* for L option, do not try to interpolate cells, density is 1.e3 */
                if (fingrd == GreensFunction.FinGrd.L)                                  //          if(fingrd.eq.'L')then
                {
                    amp = 0.001d * (new Complex(
                                oceanModel.ir1[ind - 1],
                                oceanModel.im1[ind - 1]));              //            amp = .001*cmplx(float(ir1(ind)),float(im1(ind)))
                    amp *= 1e3;                                         //            amp = 1.e3*amp
                    return amp;                                         // return
                }

                /*  now left with option F, and detailed map says ocean load an interpolating array of four elements with values from adjacent
                 *  cells, depending on which quadrant of the cell we are in. In the process note how many nonzero values we have; save one */
                int ix = 1;                                             //   ix = 1
                if (x <= 0) ix = -ix;                                   // 	 if(x.le.0) ix = -ix
                int iy = -oceanModel.longc;                                //   iy = -longc
                if (y <= 0) iy = -iy;                                   // 	 if(y.le.0) iy = -iy
                x = Math.Abs(x);                                        // 	 x=abs(x)
                y = Math.Abs(y);                                        // 	 y=abs(y)

                int nnz = 0;
                int ibind = 0;
                Complex[] bg = new Complex[4];
                Complex bgs = new Complex();

                for (int ii = 0; ii < 4; ii++)                          // 	 do i=1,4
                {
                    if (ii == 0) ibind = ind;                           // 	   if(i.eq.1) ibind=ind
                    if (ii == 1) ibind = ind + ix;                      // 	   if(i.eq.2) ibind=ind+ix
                    if (ii == 2) ibind = ind + ix + iy;                 // 	   if(i.eq.3) ibind=ind+ix+iy
                    if (ii == 3) ibind = ind + iy;                      // 	   if(i.eq.4) ibind=ind+iy

                    if (ibind <= 0 ||
                        ibind > oceanModel.latc * oceanModel.longc)     //     if(ibind.le.0.or.ibind.gt.latc*longc) 
                        bg[ii] = Complex.Zero;                          //         bg(i) = cz
                    if (ibind > 0 &&
                        ibind <= oceanModel.latc * oceanModel.longc)    //      if(ibind.gt.0.and.ibind.le.latc*longc)
                    {
                        bg[ii] = 0.001d * (new Complex(
                            oceanModel.ir1[ibind - 1],
                            oceanModel.im1[ibind - 1]));                //         bg(i) = 0.001*cmplx(float(ir1(ibind)),float(im1(ibind)))
                    }
                    if (!Complex.Equals(bg[ii], Complex.Zero))
                    {
                        nnz++;                                          //  if(bg(i).ne.cz) nnz = nnz+1
                        bgs = bg[ii];                                   //  if(bg(i).ne.cz) bgs = bg(i)
                    }
                }

                /*  now we have a grid for interpolation: the values are arranged as 
                 *           4      3
                 *           1      2
                 *  where (1) is the value we are closest to.
                 * we now fill out the zero values (if any) from the nonzero ones (Except that if there are 0 or 4 nonzero values we needn't bother) */
                if (nnz == 1)                                           //          if(nnz.eq.1) then
                {
                    /*  the one saved value is the nonzero one: fill out the whole array with this */
                    for (int ii = 0; ii < 4; ii++)                      //             do i=1,4
                        bg[ii] = bgs;                                   //               bg(i) = bgs
                }

                if (nnz == 2)                                           //          if(nnz.eq.2) then
                {
                    /*  two values: either in a row, a column, or diagonal: */

                    if (!Complex.Equals(bg[0], Complex.Zero) &&
                        !Complex.Equals(bg[1], Complex.Zero))           //             if(bg(1).ne.cz.and.bg(2).ne.cz) then
                    {
                        bg[3] = bg[0];                                  //                bg(4) = bg(1)
                        bg[2] = bg[1];                                  //                bg(3) = bg(2)
                    }
                    else if (!Complex.Equals(bg[1], Complex.Zero) &&
                             !Complex.Equals(bg[2], Complex.Zero))      //             elseif(bg(2).ne.cz.and.bg(3).ne.cz) then
                    {
                        bg[3] = bg[2];                                  //                bg(4) = bg(3)
                        bg[0] = bg[1];                                  //                bg(1) = bg(2)
                    }
                    else if (!Complex.Equals(bg[2], Complex.Zero) &&
                             !Complex.Equals(bg[3], Complex.Zero))      //             elseif(bg(3).ne.cz.and.bg(4).ne.cz) then
                    {
                        bg[0] = bg[3];                                  //                bg(1) = bg(4)
                        bg[1] = bg[2];                                  //                bg(2) = bg(3)
                    }
                    else if (!Complex.Equals(bg[3], Complex.Zero) &&
                             !Complex.Equals(bg[0], Complex.Zero))      //             elseif(bg(4).ne.cz.and.bg(1).ne.cz) then
                    {
                        bg[2] = bg[3];                                  //                bg(3) = bg(4)
                        bg[1] = bg[0];                                  //                bg(2) = bg(1)
                    }
                    /*     diagonals */
                    else if (!Complex.Equals(bg[3], Complex.Zero) &&
                             !Complex.Equals(bg[1], Complex.Zero))      //             elseif(bg(4).ne.cz.and.bg(2).ne.cz) then
                    {
                        bg[2] = 0.5d * (bg[3] + bg[1]);                 //                bg(3) = 0.5*(bg(4)+bg(2))
                        bg[0] = bg[2];                                  //                bg(1) = bg(3)
                    }
                    else if (!Complex.Equals(bg[0], Complex.Zero) &&
                             !Complex.Equals(bg[2], Complex.Zero))      //             elseif(bg(1).ne.cz.and.bg(3).ne.cz) then
                    {
                        bg[3] = 0.5d * (bg[0] + bg[2]);                 //                bg(4) = 0.5*(bg(1)+bg(3))
                        bg[1] = bg[2];                                  //                bg(2) = bg(3)                                 
                    }
                }

                if (nnz == 3)                                           // 	 if(nnz.eq.3) then
                {
                    /* one missing value */
                    if (Complex.Equals(bg[0], Complex.Zero))
                        bg[0] = bg[3] + bg[1] - bg[2];                  // 	    if(bg(1).eq.cz) bg(1) = bg(4)+bg(2)-bg(3)
                    if (Complex.Equals(bg[1], Complex.Zero))
                        bg[1] = bg[0] + bg[2] - bg[3];                  // 	    if(bg(2).eq.cz) bg(2) = bg(1)+bg(3)-bg(4)
                    if (Complex.Equals(bg[2], Complex.Zero))
                        bg[2] = bg[3] + bg[1] - bg[0];                  // 	    if(bg(3).eq.cz) bg(3) = bg(4)+bg(2)-bg(1)
                    if (Complex.Equals(bg[3], Complex.Zero))
                        bg[3] = bg[0] + bg[2] - bg[1];                  // 	    if(bg(4).eq.cz) bg(4) = bg(1)+bg(3)-bg(2)
                }
                /*  bilinear interpolation */
                amp = bg[0] * (1 - x) * (1 - y) + bg[1] * x * (1 - y) +
                      bg[3] * (1 - x) * y + bg[2] * x * y;              //  amp = bg(1)*(1-x)*(1-y) + bg(2)*x*(1-y) + bg(4)*(1-x)*y + bg(3)*x*y

                /*  if all the values are zero, we will return zero - but we pass some information about this through common block coninf; the first time this
                 *  happens, set clnf and farnf to the current distance from the station, and afterwards only farnf (since we always increase the distance). */
                if (nnz == 0)                                           // 	 if(nnz.eq.0) then
                {
                    if (conInf.numnf == 0)                             // 	    if(numnf.eq.0) then
                    {
                        conInf.clnf = conInf.dist;                    // 	       clnf = dist
                        conInf.farnf = conInf.dist;                   // 	       farnf = dist
                    }
                    else                                                // 	    else
                    {
                        conInf.farnf = conInf.dist;                   // 	       farnf = dist
                    }
                    conInf.numnf++;                                    // 	    numnf=numnf+1
                }

                /*  scale the amplitude by the sea-water density */
                rho = seaDensity.seadens(rlato, rlong);                 //       call seadens(rlato,rlong,rho)
                amp *= rho;                                             //       amp=rho*amp
                return amp;                                             //       return
            }
            return amp;                                                 //       return
        }

        private void ignewt(double del, double stp, ref double grav, ref double tilt, ref double pot)
        {

            /*     Revision 1.3  2012/03/01 21:36:12  agnew - for gravity, added case when height is zero
             *     Revision 1.2  2012/03/01 21:32:23  agnew - gravity calculation modified to use a single formula
             *                                                both gravity and tilt calculations converted to double precision
             *     Revision 1.1  2011/09/23 22:39:01  agnew - Initial revision */

            /* returns the integrated green function for newtonian gravity, tilt, and potential.
             * the function is the integral over the interval cntered at del with width stp (ie, the interval [del-stp/2,del+stp/2],
             * del and stp both being in radians the height correction is included in the green functions, 
             * the station height in meters being passed as ht in the common block stloc */




            /*  on first call do part of height correction, and find constant parts */
            if (!ignewtVariables.iflg)                                      //       if(iflg.eq.0) then
            {
                ignewtVariables.iflg = true;                                //          iflg = 1
                ignewtVariables.eps = stloc.ht / Constants.a;               //          eps = ht/a
                ignewtVariables.eps1 = ignewtVariables.eps + 1.0;           //          eps1=1.+eps
                ignewtVariables.eps2 = ignewtVariables.eps *
                    ignewtVariables.eps;                                    //          eps2=eps*eps
                ignewtVariables.g2 = Constants.gn /
                    (ignewtVariables.eps1 * ignewtVariables.eps1);          //          g2 = gn/(eps1*eps1)
                double g = (stloc.ct * 0.005279 * stloc.ct + 1) *
                    9.7803327 - stloc.ht * 3.08e-6;                         // 	 g = 9.7803327*(1+.005279*ct*ct) - 3.08e-6*ht
                ignewtVariables.em = Constants.gn / g;                      // 	 em = gn/g
                ignewtVariables.plc = Constants.a * 4 * ignewtVariables.em; // 	 plc = 4*a*em
            }                                                               //       endif
            /* sign for gravity makes + acceleration up */
            if (ignewtVariables.eps != 0.0)                                 //       if(eps.ne.0) then
            {
                double s = Math.Sin((del + stp / 2) / 2.0);                 //         s=dsin((del+stp/2)/2.d0)
                double gt = (ignewtVariables.eps1 * 2.0 * (s * s) -
                    ignewtVariables.eps) /
                    Math.Sqrt(ignewtVariables.eps1 * 4 * (s * s) +
                    ignewtVariables.eps2);                                  //         gt=(2.d0*eps1*s**2-eps)/dsqrt(4*eps1*s**2+eps2)
                s = Math.Sin((del - stp / 2) / 2.0);                        //         s=dsin((del-stp/2)/2.d0)
                grav = gt - (ignewtVariables.eps1 * 2.0 * (s * s) -
                    ignewtVariables.eps) /
                    Math.Sqrt(ignewtVariables.eps1 * 4 * (s * s) +
                    ignewtVariables.eps2);                                  //         grav=gt-(2.d0*eps1*s**2-eps)/dsqrt(4*eps1*s**2+eps2)
                grav = -ignewtVariables.g2 * grav;                          //         grav=-g2*grav
            }                                                               //       endif

            if (ignewtVariables.eps == 0.0)                                 //       if(eps.eq.0) then
            {
                grav = -ignewtVariables.g2 * (Math.Sin((del + stp / 2) / 2.0) -
                    Math.Sin((del - stp / 2) / 2.0));                       //         grav=-g2*(dsin((del+stp/2)/2.d0)-dsin((del-stp/2)/2.d0))
            }                                                               //       endif

            /* for tilt, use different approximations if < or > 0.05 rad (2.9 deg) */
            double c1 = Math.Sin(stp * 0.25);                               //       c1 = dsin(0.25d0*stp)
            if (del >= 0.05)                                                //       if(del.ge.0.05) then
            {
                tilt = ignewtVariables.em * (Math.Sin(del / 2.0) * (-2.0) * c1 +
                    Math.Log(Math.Tan((del + stp * 0.5) / 4.0) /
                    Math.Tan((del - stp * 0.5) / 4.0)));                    //        tilt = em*(-2.*dsin(del/2.d0)*c1 + dlog(dtan((del+0.5*stp)/4.d0)/dtan((del-0.5*stp)/4.d0)))
            }                                                               //       endif

            if (del < 0.05)                                                 //       if(del.lt.0.05) then
            {
                double d1 = del + stp / 2;                                  //          d1 = del + stp/2
                double d2 = del - stp / 2;                                  //          d2 = del - stp/2
                tilt = ignewtVariables.em * (d2 / Math.Sqrt(ignewtVariables.eps2 + d2 * d2) - d1 /
                    Math.Sqrt(ignewtVariables.eps2 + d1 * d1) +
                    Math.Log((d1 + Math.Sqrt(ignewtVariables.eps2 + d1 * d1)) /
                    (d2 + Math.Sqrt(ignewtVariables.eps2 + d2 * d2))));     //         tilt = em*(d2/dsqrt(eps2+d2*d2) - d1/dsqrt(eps2+d1*d1) + dlog( (d1+dsqrt(eps2+d1*d1))/(d2+dsqrt(eps2+d2*d2)) ) )
            }                                                               //       endif
            pot = ignewtVariables.plc * Math.Cos(del / 2.0) * c1;           //       pot=plc*cos(del/2.)*c1
        }
        #endregion

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
            return calculateComponent(TimeConversion.TimeConversion.ModJulDay2Time(mjd), component, azimuth);
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
            /* adjust to a new azimuth, rotated azim deg (clockwise in the command line,
             * converted to counterclockwise below to match the usual formulae) from the
             * EW-NS one that the loads are originally in. */
            double caz = Math.Cos(-azimuth * Constants.degree2radian);
            double saz = Math.Sin(-azimuth * Constants.degree2radian);

            // return value
            double result = double.NaN;
            int length = oceanLoadingParameter.waveItems.Length;
            TidalWave[] hartidPM = new TidalWave[length + 1];

            if (component == Component.Gravity)
            {
                for (int i = 0; i < length; i++)
                {
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.g.Amplitude,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.g.Phase);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Gravity;
            }
            else if (component == Component.GravityDeformation)
            {
                for (int i = 0; i < length; i++)
                {
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.gdef.Amplitude,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.gdef.Phase);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Gravity;
            }
            else if (component == Component.PotentialHeight)
            {
                TidalWave[] hartidPMGravityDeformation = new TidalWave[length + 1];
                for (int i = 0; i < length; i++)
                {
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.g.Amplitude,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.g.Phase);
                    hartidPMGravityDeformation[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte,
                                                                  oceanLoadingParameter.waveItems[i].oceanLoad.gdef.Amplitude,
                                                                  oceanLoadingParameter.waveItems[i].oceanLoad.gdef.Phase);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);
                hartidPMGravityDeformation[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = (hartid(dt, hartidPM, ModeOfOperation.l) - hartid(dt, hartidPMGravityDeformation, ModeOfOperation.l)) * dConversion_Gravity;
            }
            else if (component == Component.OceanTide)
            {
                for (int i = 0; i < length; i++)
                {
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.o.Amplitude,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.o.Phase);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);
                
                result = hartid(dt, hartidPM, ModeOfOperation.g) * dConversion_OceanTide;
            }
            else if (component == Component.Tilt)
            {
                for (int i = 0; i < length; i++)
                {
                    Complex tEW = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.t[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.t[0].Phase));
                    Complex tNS = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.t[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.t[1].Phase));

                    Complex tilt = phasor(-saz * tEW + caz * tNS);

                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, tilt.Real, tilt.Imaginary);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Tilt;
            }
            else if (component == Component.TiltDeformation)
            {
                for (int i = 0; i < length; i++)
                {
                    Complex tEW = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.tdef[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.tdef[0].Phase));
                    Complex tNS = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.tdef[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.tdef[1].Phase));

                    Complex tiltdef = phasor(-saz * tEW + caz * tNS);
                    
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, tiltdef.Real, tiltdef.Imaginary);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Tilt;
            }
            else if (component == Component.TiltAttraction)
            {
                TidalWave[] hartidPMTiltDef = new TidalWave[length + 1];
                for (int i = 0; i < length; i++)
                {
                    Complex tEW = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.t[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.t[0].Phase));
                    Complex tNS = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.t[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.t[1].Phase));

                    Complex tilt = phasor(-saz * tEW + caz * tNS);
                    
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, tilt.Real, tilt.Imaginary);
                    //
                    Complex tEWd = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.tdef[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.tdef[0].Phase));
                    Complex tNSd = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.tdef[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.tdef[1].Phase));

                    Complex tiltdef = phasor(-saz * tEWd + caz * tNSd);
                    
                    hartidPMTiltDef[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, tiltdef.Real, tiltdef.Imaginary);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);
                hartidPMTiltDef[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = (hartid(dt, hartidPM, ModeOfOperation.l) - hartid(dt, hartidPMTiltDef, ModeOfOperation.l)) * dConversion_Tilt;
            }

            else if (component == Component.HorizontalStrain)
            {
                for (int i = 0; i < length; i++)
                {
                    Complex sEW = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[0].Phase));
                    Complex sNS = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[1].Phase));
                    Complex sSh = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[2].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[2].Phase));

                    Complex strain = phasor(caz * caz * sNS + saz * saz * sEW - 2 * caz * saz * sSh);

                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, strain.Real, strain.Imaginary);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Strain;
            }
            else if (component == Component.ShearStrain)
            {
                for (int i = 0; i < length; i++)
                {
                    Complex sEW = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[0].Phase));
                    Complex sNS = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[1].Phase));
                    Complex sSh = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[2].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[2].Phase));

                    Complex strain = phasor(sSh * (caz * caz - saz * saz) - caz * saz * (sEW - sNS));

                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, strain.Real, strain.Imaginary);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Strain;
            }
            else if (component == Component.ArealStrain)
            {
                for (int i = 0; i < length; i++)
                {
                    Complex sEW = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[0].Phase));
                    Complex sNS = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[1].Phase));

                    Complex strain = phasor(sEW + sNS);
                    
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, strain.Real, strain.Imaginary);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Strain;
            }
            else if (component == Component.VolumeStrain)
            {
                for (int i = 0; i < length; i++)
                {
                    Complex sEW = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[0].Phase));
                    Complex sNS = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.s[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.s[1].Phase));

                    Complex strain = phasor(0.6667d * (sEW + sNS));

                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, strain.Real, strain.Imaginary);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Strain;
            }
            else if (component == Component.VerticalDisplacement)
            {
                for (int i = 0; i < length; i++)
                {
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.d[2].Amplitude,
                                                oceanLoadingParameter.waveItems[i].oceanLoad.d[2].Phase);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Displacement;
            }
            else if (component == Component.HorizontalDisplacement)
            {
                for (int i = 0; i < length; i++)
                {
                    Complex dEW = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.d[0].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.d[0].Phase));
                    Complex dNS = argand(new Complex(oceanLoadingParameter.waveItems[i].oceanLoad.d[1].Amplitude, oceanLoadingParameter.waveItems[i].oceanLoad.d[1].Phase));

                    Complex displacement = phasor(-saz * dEW + caz * dNS);
                    hartidPM[i] = new TidalWave(oceanLoadingParameter.waveItems[i].icte, displacement.Real, displacement.Imaginary);
                }
                hartidPM[hartidPM.Length - 1] = new TidalWave(new int[6] { -1, 0, 0, 0, 0, 0 }, 0.0, 0.0);

                result = hartid(dt, hartidPM, ModeOfOperation.l) * dConversion_Displacement;
            }
            return result;
        }

        #region Private functions for ocean loading prediction.

        /// <summary>
        /// Parameter of a tidal wave to be predicted.
        /// </summary>
        private class TidalWave
        {
            /// <summary>
            /// Doodsen number of the tidal wave to be predicted.
            /// </summary>
            public int[] DoodsenNumber;

            /// <summary>
            /// Amplitude of the tidal wave to be predicted.
            /// </summary>
            public double Amplitude;

            /// <summary>
            /// Phase of the tidal wave to be predicted.
            /// </summary>
            public double Phase;

            public TidalWave(int[] DoodsenNumber, double Amplitude, double Phase)
            {
                this.DoodsenNumber = DoodsenNumber;
                this.Amplitude = Amplitude;
                this.Phase = Phase;
            }
        }

        /// <summary>
        /// Finds tides by summing harmonics; the amp and phase of these are specified in subroutine admint.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="_hartidPM"></param>
        /// <returns></returns>
        private double hartid(DateTime dt, TidalWave[] _hartidPM, ModeOfOperation modeOfOperation)
        {
            /* Revision 1.2  2011/11/18 20:51:13  agnew - modified for larger number of harmonics (342)
             * Revision 1.1  2011/09/23 22:39:01  agnew - Initial revision */

            // convert times to Julian days (UT) then to Julian centuries from J2000.00 (ET)
            double jd = juldat(dt);
            double dayfr = dt.Hour / 24d + dt.Minute / 1440d + dt.Second / 84600d;
            double year = dt.Year + ((dt.DayOfYear + dayfr) / (365 + leap(dt.Year)));
            //double delta = etutc(year);                                                   // Original version
            double delta = TAI_UTC.GetDifferenceFrom_TDT_UTC(year);
            double djd = jd - 0.5d + dayfr;
            double time = (djd - 2451545d + delta / 86400d) / 36525d;

            /*   get harmonics information */
            TidalWave[] tidalWavesLocalPhase = new TidalWave[_hartidPM.Length];
            for (int i = 0; i < _hartidPM.Length; i++)
            {
                TidalWave wptt = _hartidPM[i];
                if (modeOfOperation == ModeOfOperation.l)
                    tidalWavesLocalPhase[i] =
                        new TidalWave(
                            wptt.DoodsenNumber, wptt.Amplitude, wptt.Phase + wptt.DoodsenNumber[0] * oceanLoadingParameter.Location.Longitude);
                else if (modeOfOperation == ModeOfOperation.g)
                    tidalWavesLocalPhase[i] =
                        new TidalWave(
                            wptt.DoodsenNumber, wptt.Amplitude, wptt.Phase);
            }

            int ntt = Harmonics.nt;                                                         //       ntt=nt
            double[] a = new double[Harmonics.nt];
            double[] f = new double[Harmonics.nt];
            double[] p = new double[Harmonics.nt];
            admint(tidalWavesLocalPhase, ref a, ref f, ref p, ref ntt, time, dayfr);        //       call admint(amp,idt,phase,a,f,p,nin,ntt)

            /*  set up for first recursion, and normalize frequencies */
            double[] wf = new double[Harmonics.nt];
            for (int i = 0; i < ntt; i++)                                                   //       do i=1,ntt
            {
                p[i] = Constants.degree2radian * p[i];                                      //         p(i) = dr*p(i)
                f[i] = Constants.pi * f[i] / 43200.0;                                       //         f(i) = samp*pi*f(i)/43200.d0
                wf[i] = f[i];                                                               //         f(i) = samp*pi*f(i)/43200.d0
            }                                                                               //       enddo

            /* set up harmonic coefficients, compute tide, and write out */
            double[] hc = new double[Harmonics.nt * 2];
            for (int i = 1; i < ntt; i++)                                                   //       do i=1,ntt
            {
                hc[2 * i - 2] = a[i - 1] * Math.Cos(p[i - 1]);                              //         hc(2*i-1) = a(i)*dcos(p(i))
                hc[2 * i - 1] = -a[i - 1] * Math.Sin(p[i - 1]);                             //         hc(2*i)  = -a(i)*dsin(p(i))
            }                                                                               //       enddo
            double[] scr = new double[3 * Harmonics.nt];
            return recurs(hc, ntt, wf, ref scr);                                            //       call recurs(x,np,hc,ntt,wf,scr)
        }

        /// <summary>
        /// Returns the amplitude, frequency, and phase of a set of tidal constituents.
        /// </summary>
        /// <param name="tidalWaveIn">Tidal wave (doodsen number,amplitude, phase) for prediction.</param>
        /// <param name="amp">Calculated mplitude.</param>
        /// <param name="f">Calculated frequency.</param>
        /// <param name="p">Calculated phase.</param>
        /// <param name="nout">Number of used harmonics.</param>
        /// <param name="time">Time in J2000.</param>
        /// <param name="dayfr">Time of day in fractional days.</param>
        private void admint(TidalWave[] tidalWaveIn,
                             ref double[] amp, ref double[] f, ref double[] p,
                             ref int nout, double time, double dayfr)                       //       subroutine admint(ampin,idtin,phin,amp,f,p,nin,nout)
        {
            /* Revision 1.3  2012/03/29 20:34:55 : agnew : moved parameter statements to satisfy gfortran
             * Revision 1.2  2011/11/18 20:51:38 : agnew : modified to 342 harmonics, following hardisp (IERS standard)
             * Revision 1.1  2011/09/23 22:39:00 : agnew : Initial revision */

            /* Returns the amplitude, frequency, and phase of a set of tidal
             * constituents. n is input as the number wanted, and returned as the number
             * provided.  The constituents used are stored in the arrays idd (doodson
             * number) and tamp (Cartwright-Edden amplitude).  The actual amp and
             * phase of each of these are determined by spline interpolation of the
             * real and imaginary part of the admittance, as specified at a subset
             * of the constituents. 
             * The phase is determined for a time set in common block /date/ (see
             * subroutine tdfrph), outside of this subroutine. */

            double[] sdr = new double[tidalWaveIn.Length];
            double[] rl = new double[tidalWaveIn.Length];
            double[] aim = new double[tidalWaveIn.Length];
            double[] rf = new double[tidalWaveIn.Length + 1];
            double[] zdi = new double[tidalWaveIn.Length];
            double[] zdr = new double[tidalWaveIn.Length];
            double[] di = new double[tidalWaveIn.Length];
            double[] dr = new double[tidalWaveIn.Length];
            double[] sdi = new double[tidalWaveIn.Length];
            double[] scr = new double[tidalWaveIn.Length];
            int[] key = new int[tidalWaveIn.Length];

            int k = 0;                                                                      //       k = 0
            for (int ll = 0; ll < tidalWaveIn.Length; ll++)                                 //       do ll=1,nin
            {
                /*  see if Doodson numbers match */
                int ii = 0, kk = 0;
                for (kk = 0; kk < Harmonics.nt; kk++)                                       //          do kk=1,nt
                {
                    ii = 0;                                                                 //             ii = 0
                    for (int i = 0; i < 6; i++)                                             //             do i=1,6
                        ii += Math.Abs(Harmonics.idd[kk, i] -
                              tidalWaveIn[ll].DoodsenNumber[i]);                            //                ii = ii + iabs(idd(i,kk)-idtin(i,ll))

                    if (ii == 0)                                                            //             if(ii.eq.0) go to 5
                        break;
                }                                                                           //          enddo

                /* have a match - put line into array */
                if (ii == 0 && k < tidalWaveIn.Length)                                      //  5       if(ii.eq.0.and.k.lt.ncon) then
                {
                    rl[k] = tidalWaveIn[ll].Amplitude *
                        Math.Cos(Constants.degree2radian * tidalWaveIn[ll].Phase) /
                        Math.Abs(Harmonics.tamp[kk]);                                       //             rl(k) = ampin(ll)*cos(dtr*phin(ll))/abs(tamp(kk))
                    aim[k] = tidalWaveIn[ll].Amplitude *
                        Math.Sin(Constants.degree2radian * tidalWaveIn[ll].Phase) /
                        Math.Abs(Harmonics.tamp[kk]);                                       //             aim(k)= ampin(ll)*sin(dtr*phin(ll))/abs(tamp(kk))

                    /* Now have real and imaginary parts of admittance, scaled by C-E
                     * amplitude. Admittance phase is whatever was used in the original
                     * expression. (Usually phase is given relative to some reference
                     * but amplitude is in absolute units). Next get frequency. */

                    double fr = 0.0;
                    double pr = 0.0;
                    tdfrph(Harmonics.get_idd(kk), ref fr, ref pr, time, dayfr);             //             call tdfrph(idd(1,kk),fr,pr)

                    rf[k] = fr;                                                             //             rf(k) = fr
                    k++;                                                                    //             k = k + 1 - replaced from beginning of the if-condition zo here
                }                                                                           //          endif
            }                                                                               //       enddo

            /* Done going through constituents-there are k of them have specified 
             * admittance at a number of points. Sort these by frequency and 
             * separate diurnal and semidiurnal, recopying admittances to get them in order.*/

            shells(ref rf, ref key, k);                                                     //       call shells(rf,key,k)

            int nlp = 0;                                                                    //       nlp = 0
            int ndi = 0;                                                                    //       ndi = 0
            int nsd = 0;                                                                    //       nsd = 0
            for (int i = 0; i < k; i++)                                                     //       do i=1,k
            {
                if (rf[i] < 0.5)                                                            //          if(rf(i).lt.0.5) nlp = nlp + 1
                    nlp++;
                if (rf[i] < 1.5 && rf[i] > 0.5)                                             //          if(rf(i).lt.1.5.and.rf(i).gt.0.5) ndi = ndi + 1
                    ndi++;
                if (rf[i] < 2.5 && rf[i] > 1.5)                                             //          if(rf(i).lt.2.5.and.rf(i).gt.1.5) nsd = nsd + 1
                    nsd++;
                scr[i] = rl[key[i] - 1];                                                    //          scr(i) = rl(key(i))
            }                                                                               //       enddo

            for (int i = 0; i < k; i++)                                                     //       do i=1,k
            {
                rl[i] = scr[i];                                                             //          rl(i) = scr(i)
                scr[i] = aim[key[i] - 1];                                                   //          scr(i) = aim(key(i))
            }                                                                               //       enddo

            for (int i = 0; i < k; i++)                                                     //       do i=1,k
                aim[i] = scr[i];                                                            //          aim(i) = scr(i)

            /* Now set up splines (8 cases - four species, each real and imag)
             * we have to allow for the case when there are no constituent amplitudes
             * for the long-period. */

            if (nlp != 0)                                                                   //       if(nlp.ne.0) call spline(nlp,rf,rl,zdr,scr)
                spline(nlp, rf, rl, ref zdr, scr);

            if (nlp != 0)                                                                   //       if(nlp.ne.0) call spline(nlp,rf,aim,zdi,scr)
                spline(nlp, rf, aim, ref zdi, scr);

            spline(ndi, cpfrom(rf, nlp), cpfrom(rl, nlp), ref dr, scr);                     //       call spline(ndi,rf(nlp+1),rl(nlp+1),dr,scr)
            spline(ndi, cpfrom(rf, nlp), cpfrom(aim, nlp), ref di, scr);                    //       call spline(ndi,rf(nlp+1),aim(nlp+1),di,scr)
            spline(nsd, cpfrom(rf, nlp + ndi), cpfrom(rl, nlp + ndi), ref sdr, scr);        //       call spline(nsd,rf(nlp+ndi+1),rl(nlp+ndi+1),sdr,scr)
            spline(nsd, cpfrom(rf, nlp + ndi), cpfrom(aim, nlp + ndi), ref sdi, scr);       //       call spline(nsd,rf(nlp+ndi+1),aim(nlp+ndi+1),sdi,scr)

            /*  evaluate all harmonics using the interpolated admittance */
            int j = 0;                                                                      //       j = 1
            for (int i = 0; i < Harmonics.nt; i++)                                          //       do i=1,nt
            {
                if (Harmonics.idd[i, 0] == 0 && nlp == 0)                                   //          if(idd(1,i).eq.0.and.nlp.eq.0) go to 11
                    continue;                                                               // goto L11;

                tdfrph(Harmonics.get_idd(i), ref f[j], ref p[j], time, dayfr);              //          call tdfrph(idd(1,i),f(j),p(j))

                /*  phase corrections to equilibrium tide */
                if (Harmonics.idd[i, 0] == 0)                                               //          if(idd(1,i).eq.0) p(j) = p(j) + 180.
                    p[j] += 180.0;

                if (Harmonics.idd[i, 0] == 1)                                               //          if(idd(1,i).eq.1) p(j) = p(j) + 90. 
                    p[j] += 90.0;

                double sf = f[j];                                                           //          sf = f(j)
                double re = 0.0, am = 0.0;

                if (Harmonics.idd[i, 0] == 0)
                {
                    re = eval(sf, nlp, rf, rl, zdr);                                        //          if(idd(1,i).eq.0) re = eval(sf,nlp,rf,rl,zdr)
                    am = eval(sf, nlp, rf, aim, zdi);                                       //          if(idd(1,i).eq.0) am = eval(sf,nlp,rf,aim,zdi)
                }
                if (Harmonics.idd[i, 0] == 1)
                {
                    re = eval(sf, ndi, cpfrom(rf, nlp), cpfrom(rl, nlp), dr);               //          if(idd(1,i).eq.1) re = eval(sf,ndi,rf(nlp+1),rl(nlp+1),dr)
                    am = eval(sf, ndi, cpfrom(rf, nlp), cpfrom(aim, nlp), di);              //          if(idd(1,i).eq.1) am = eval(sf,ndi,rf(nlp+1),aim(nlp+1),di)
                }
                if (Harmonics.idd[i, 0] == 2)
                {
                    re = eval(sf, nsd, cpfrom(rf, nlp + ndi), cpfrom(rl, nlp + ndi), sdr);  //          if(idd(1,i).eq.2) re = eval(sf,nsd,rf(nlp+ndi+1),rl(nlp+ndi+1),sdr)
                    am = eval(sf, nsd, cpfrom(rf, nlp + ndi), cpfrom(aim, nlp + ndi), sdi); //          if(idd(1,i).eq.2) am = eval(sf,nsd,rf(nlp+ndi+1),aim(nlp+ndi+1),sdi)
                }

                amp[j] = Harmonics.tamp[i] * Math.Sqrt(re * re + am * am);                  //          amp(j) = tamp(i)*sqrt(re**2+am**2)
                p[j] += Math.Atan2(am, re) / Constants.degree2radian;                       //          p(j) = p(j) + atan2(am,re)/dtr
                if (p[j] > 180.0)                                                           //          if(p(j).gt.180) p(j)=p(j)-360.
                    p[j] += -360.0;

                j++;                                                                        //          j = j + 1
            }                                                                               //  11      continue
            nout = j - 1;                                                                   //       nout = j - 1
        }

        /// <summary>
        /// Does sin and cos recursion to fill data x.
        /// </summary>
        /// <param name="hc">Contains alternating cosine and sine coefficients.</param>
        /// <param name="nf">Number of sines and cosines frequenciies.</param>
        /// <param name="om">Sines and cosines with frequenciies.</param>
        /// <param name="scr">Scratch array of length 3*nf</param>
        /// <returns>Sin and cos recursion of data.</returns>
        private double recurs(double[] hc, int nf, double[] om, ref double[] scr)           //       subroutine recurs(x,n,hc,nf,om,scr)
        {
            /* does sin and cos recursion to fill in data x, of length n, for
             * nf sines and cosines with frequenciies om (normalized so the
             * nyquist is pi). hc contains alternating cosine and sine coefficients
             * scr is a scratch array of length 3*nf (n.b. - it is double precision) */

            /* set up for start of recursion by computing harmonic values
             * at start point and just before it */

            for (int i = 1; i <= nf; i++)                                                   //       do 3 i = 1,nf
            {
                scr[i * 3 - 3] = hc[2 * i - 2];                                             //       scr(3*i-2) = hc(2*i-1)
                scr[i * 3 - 2] = hc[2 * i - 1] * Math.Cos(om[i - 1]) -
                                 hc[i * 2 - 1] * Math.Sin(om[i - 1]);                       //       scr(3*i-1) = hc(2*i-1)*cos(om(i)) -hc(2*i)*sin(om(i))
                scr[i * 3 - 1] = Math.Cos(om[i - 1]) * 2.0;                                 //  3    scr(3*i) = 2.*dcos(dble(om(i)))
            }
            /*  do recursion over data */
            double x = 0.0;                                                                 //       x(i) = 0.
            /*  do recursive computation for each harmonic */
            for (int j = 1; j <= nf; j++)                                                   //       do 5 j  = 1,nf
            {
                x += scr[j * 3 - 3];                                                        //       x(i) = x(i) + scr(3*j-2)
                double sc = scr[j * 3 - 2];                                                 //       sc = scr(3*j-2)
                scr[j * 3 - 2] = scr[j * 3] * sc - scr[j * 3 - 1];                          //       scr(3*j-2) = scr(3*j)*sc-scr(3*j-1)
                scr[j * 3 - 1] = sc;                                                        //  5    scr(3*j-1) = sc
            }
            return x;
        }

        /// <summary>
        /// Given the Doodson number of a tidal constituent, returns the frequency and phase.
        /// </summary>
        /// <param name="idood">Doodson number of a tidal constituent</param>
        /// <param name="freq">Returns frequency in cycles/day</param>
        /// <param name="phase">Returns phase in degrees.</param>
        /// <param name="time">Time in J2000.</param>
        /// <param name="dayfr">Time of day in fractional days.</param>
        private void tdfrph(int[] idood,
            ref  double freq, ref double phase, double time, double dayfr)                  //       subroutine tdfrph(idood,freq,phase)
        {
            /* Revision 1.2  2011/11/18 18:38:09 : agnew : added 2008 leap second
             * Revision 1.1  2011/09/23 22:39:02 : agnew : Initial revision */

            /* Given the Doodson number of a tidal constituent (in idood), returns
             * the frequency and phase.  Phase is returned in degrees and frequency
             * in cycles/day. 
             * 
             * Note that phases must be decreased by 90 degrees if the sum of the order 
             * and the species number is odd (as for the 2nd degree diurnals, and 3rd
             * degree low-frequency and semidiurnals).
             * These phases may need further adjustment to allow for the spherical
             * harmonic normalization used; e.g. for that used for the potential by 
             * Cartwright and Tayler, 180 degrees must be added for 
             * (species,order) = (1,2), (1,3), or (3,3). */

            /* set the phases and frequencies for each of the Doodson arguments */
            double[] dd = new double[6];
            double[] d = new double[6];

            /* IERS expressions for the Delauney arguments, in degrees */
            double f1 = time * (time * (time * (time * -6.8e-8 + 1.43431e-5) +
                0.0088553333) + 477198.8675605) + 134.96340251;                             //          f1 =  134.9634025100d0 + t*(477198.8675605000d0 + t*(0.0088553333d0 + t*(0.0000143431d0 + t*(-0.0000000680d0))))  
            double f2 = time * (time * (time * (time * -3.2e-9 + 3.78e-8) -
                1.536667e-4) + 35999.0502911389) + 357.5291091806;                          //          f2 = 357.5291091806d0 + t*(35999.0502911389d0 + t*(-0.0001536667d0 + t*(0.0000000378d0 + t*(-0.0000000032d0))))
            double f3 = time * (time * (time * (time * 1.2e-9 - 2.881e-7) -
                0.003542) + 483202.0174577222) + 93.27209062;                               //          f3 = 93.2720906200d0 + t*(483202.0174577222d0 + t*(-0.0035420000d0 + t*(-0.0000002881d0 + t*(0.0000000012d0))))
            double f4 = time * (time * (time * (time * -8.8e-9 + 1.8314e-6) -
                0.0017696111) + 445267.1114469445) + 297.8501954694;                        //          f4 = 297.8501954694d0 + t*(445267.1114469445d0 + t*(-0.0017696111d0 + t*(0.0000018314d0 + t*(-0.0000000088d0))))
            double f5 = time * (time * (time * (time * -1.65e-8 + 2.1394e-6) +
                0.0020756111) - 1934.1362619722) + 125.04455501;                            //          f5 = 125.0445550100d0 + t*(-1934.1362619722d0 + t*(0.0020756111d0 + t*(0.0000021394d0 + t*(-0.0000000165d0))))     

            /*  convert to Doodson (Darwin) variables */
            d[0] = dayfr * 360.0 - f4;                                                      //         d(1) = 360.d0*dayfr - f4
            d[1] = f3 + f5;                                                                 //         d(2) = f3 + f5
            d[2] = d[1] - f4;                                                               // 	       d(3) = d(2) - f4
            d[3] = d[1] - f1;                                                               //         d(4) = d(2) - f1
            d[4] = -f5;                                                                     //         d(5) = -f5
            d[5] = d[2] - f2;                                                               //         d(6) = d(3) - f2 

            /*   find frequencies of Delauney variables (in cycles/day), and from these the same for the Doodson arguments */
            double fd1 = time * 1.3e-9 + 0.0362916471;                                      //         fd1 =  0.0362916471 + 0.0000000013*t
            double fd2 = 0.0027377786;                                                      //         fd2 =  0.0027377786
            double fd3 = 0.0367481951 - time * 5e-10;                                       //         fd3 =  0.0367481951 - 0.0000000005*t
            double fd4 = 0.033863192 - time * 3e-10;                                        //         fd4 =  0.0338631920 - 0.0000000003*t
            double fd5 = time * 3e-10 - 1.470938e-4;                                        //         fd5 = -0.0001470938 + 0.0000000003*t
            dd[0] = 1.0 - fd4;                                                              //         dd(1) = 1.d0 - fd4
            dd[1] = fd3 + fd5;                                                              //         dd(2) = fd3 + fd5
            dd[2] = dd[1] - fd4;                                                            //         dd(3) = dd(2) - fd4
            dd[3] = dd[1] - fd1;                                                            //         dd(4) = dd(2) - fd1
            dd[4] = -fd5;                                                                   //         dd(5) = -fd5
            dd[5] = dd[2] - fd2;                                                            //         dd(6) = dd(3) - fd2 
            /*   end of intialization  */

            /*  compute phase and frequency of the given tidal constituent */
            freq = 0.0;                                                                     //       freq=0.d0
            phase = 0.0;                                                                    //       phase=0.d0

            for (int i = 0; i < 6; i++)                                                     //       do i = 1,6
            {
                freq += idood[i] * dd[i];                                                   //          freq =   freq + idood(i)*dd(i)
                phase += idood[i] * d[i];                                                   //          phase = phase + idood(i)*d(i)
            }                                                                               //       enddo

            /* adjust phases to fall in the positive range 0 to 360 */
            phase = phase % 360.0;                                                          //       phase = dmod(phase,360.d0)
            if (phase < 0.0)                                                                //       if(phase.lt.0.d0) phase = phase + 360.d0
                phase += 360.0;
        }

        /// <summary>
        /// Determines leap years, by Clavian (Gregorian) rule 
        /// </summary>
        /// <param name="iy">Year.</param>
        /// <returns>Returns 1 if year is a leap year.</returns>
        private int leap(int iy)                                                            //       function leap(iy)
        {
            /* returns 1 if year is a leap year, by Clavian (Gregorian) rule */
            int ret_val = 1 - (iy % 4 + 3) / 4;                                             //       leap = 1 - (mod(iy,4)+3)/4
            if (iy % 100 == 0 && iy % 400 != 0)                                             //       if(mod(iy,100).eq.0.and.mod(iy,400).ne.0) leap=0
                ret_val = 0;
            return ret_val;
        }

        /// <summary>
        /// Calculates Julian Date from Gregorian date.
        /// </summary>
        /// <param name="dt">Gregorian date.</param>
        /// <returns>Julian date.</returns>
        private int juldat(DateTime dt)                                                     //       function juldat(it)
        {
            /* Julian Date from Gregorian date, Algorithm from p. 604, Explanatory
             * Supplement Amer Ephemeris & Nautical Almanac (cf Comm CACM, 11, 657 (1968)
             * and 15, 918 (1972)) Valid for all positive values of Julian Date */

            return (dt.Year + 4800 + (dt.Month - 14) / 12) * 1461 / 4 + (dt.Month - 2 - (
                dt.Month - 14) / 12 * 12) * 367 / 12 - (dt.Year + 4900 + (dt.Month - 14) /
                 12) / 100 * 3 / 4 + dt.Day - 32075;        //      juldat=(1461*(it(1)+4800+(it(2)-14)/12))/4 + (367*(it(2)-2-12*((it(2)-14)/12)))/12 - (3*((it(1)+4900+(it(2)-14)/12)/100))/4+it(3)-32075
        }

        // Obsolete and replaced by the Eterna method - Problems on calculation before year 1960
        public double etutc(double year)                                                    //       subroutine etutc(year,delta)
        {
            double[] d__ = new double[142]
            { 
                5.15,4.64,5.36,3.49,
                3.27,2.45,4.03,1.76,3.3,
                1.0,2.42,0.94,2.31,2.27,-.22,
                0.03,-.05,-.06,-.57,0.03,
                -0.47,0.98,-.86,2.45,0.22,0.37,
                2.79,1.2,3.52,1.17,2.67,
                3.06,2.66,2.97,3.28,3.31,3.33,
                3.23,3.6,3.52,4.27,2.68,
                2.75,2.67,1.94,1.39,1.66,0.88,
                0.33,-.17,-1.88,-3.43,-4.05,
                -5.77,-7.06,-7.36,-7.67,-7.64,
                -7.93,-7.82,-8.35,-7.91,-8.03,
                -9.14,-8.18,-7.88,-7.62,-7.17,
                -8.14,-7.59,-7.17,-7.94,-8.23,
                -7.88,-7.68,-6.94,-6.89,-7.11,
                -5.87,-5.04,-3.9,-2.87,-.58,
                .71,1.8,3.08,4.63,5.86,
                7.21,8.58,10.5,12.1,12.49,
                14.41,15.59,15.81,17.52,19.01,
                18.39,19.55,20.36,21.01,21.81,
                21.76,22.35,22.68,22.94,22.93,
                22.69,22.94,23.2,23.31,23.63,
                23.47,23.68,23.62,23.53,23.59,
                23.99,23.8,24.2,24.99,24.97,
                25.72,26.21,26.37,26.89,27.68,
                28.13,28.94,29.42,29.66,30.29,
                30.96,31.09,31.59,32.06,31.82,
                32.69,33.05,33.16,33.59 
            };
            double[] tx = new double[39]
            { 
                61.5,62.0,62.5,63.0,
                63.5,64.0,64.5,65.0,65.5,
                66.0,66.5,67.0,67.5,68.0,68.25,
                68.5,68.75,69.0,69.25,69.5,
                69.75,70.0,70.25,70.5,70.75,
                71.0,71.085,71.162,71.247,71.329,
                71.414,71.496,71.581,71.666,
                71.748,71.833,71.915,71.999,72.0
            };
            double[] ty = new double[39]
            { 
                33.59,34.032,34.235,
                34.441,34.644,34.95,35.286,35.725,
                36.16,36.498,36.968,37.444,
                37.913,38.39,38.526,38.76,39.0,
                39.238,39.472,39.707,39.946,40.185,
                40.42,40.654,40.892,41.131,
                41.211,41.284,41.364,41.442,41.522,
                41.6,41.68,41.761,41.838,41.919,
                41.996,42.184,42.184 
            };
            double[] st = new double[24]
            { 
                1972.5,1973.0,1974.0, 1975.0,1976.0,1977.0,1978.0,1979.0, 
                1980.0,1981.5,1982.5,1983.5,1985.5,
                1988.0,1990.0,1991.0,1992.5,1993.5,
                1994.5,1996.0,1997.5,1999.0,
                2006.0,2008.0 
            };
            double[] si = new double[24]
            { 
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
            };

            int nstep = 24;

            /* Local variables */
            int i;
            double frac;
            double delta = 0.0;

            /*   input year (and fraction thereof), 1700-2006 */
            /*   output delta = et - utc in seconds */
            /*   tables from p 90, Explanatory Supplement to the A.E., Dr. R A Broucke, */
            /*   JPL, and the leap.sec table in GAMIT */
            /*   utc (and ut) is time usually used (eg in time signals) */


            //       dimension d(142),tx(39),ty(39),st(24),si(24)
            //       
            /*   update st,si, and nstep to allow for added leap seconds */
            //       
            //       data si/24*1./
            //       data nstep/24/

            if (year > 1972.0)   //       if(year.gt.1972.) go to 6
                goto L6;

            if (year > 1961.5) //       if(year.gt.1961.5) go to 2
                goto L2;

            if (year > 1820.5) //       if(year.gt.1820.5) go to 1 
                goto L1;

            /*  for oldest epochs, use approximation */
            if (year >= 1785.0) //       if(year.ge.1785.) delta = 6.
                delta = 6.0;

            if (year < 1785.0) //       if(year.lt.1785.) delta = (year-1750.)/5. 
                delta = (year - 1750.0) / 5.0;

            if (year < 1700.0) //       if(year.lt.1700.) delta = 0.
                delta = 0.0;

            return delta;   //       return

            /*   for 1820.5 to 1961.5, data is spaced at yearly intervals */
        L1:
            int n = (int)(year - 1819.5);  //  1    n = year - 1819.5
            frac = year - (n + 1819.5); //       frac = year - (1819.5 + n)
            delta = (d__[n] - d__[n - 1]) * frac + d__[n - 1];  //       delta = (d(n+1) - d(n))*frac + d(n)
            return delta;   //       return

            /*   for 1961.5 to 1972.0, interpolate between unequispaced data */

        L2:
            for (i = 0; i < 38; i++) //  2    do 3 i = 1,38
            {

                if (year - 1900.0 == tx[i]) //       if(year-1900..eq.tx(i)) go to 4
                    goto L4;

                if (year - 1900.0 < tx[i]) //       if(year-1900..lt.tx(i)) go to 5
                    goto L5;


            }   //  3    continue

        L4:
            delta = ty[i];   //  4    delta = ty(i)

            return delta;   //       return
        //  5    
        L5:
            delta = ty[i - 1] + (ty[i] - ty[i - 1]) * ((year -
                1900.0 - tx[i - 1]) / (tx[i] - tx[i - 1])); //delta=ty(i-1)+(ty(i)-ty(i-1))*((year-1900.-tx(i-1))/(tx(i)-tx(i-1)))

            return delta;     //       return
        /*   after 1972 et-utc has only step offsets. st is the array of step times, si the step sizes (an added second is +1.) */

        L6:
            delta = 42.184;    //  6    delta = 42.184

            for (i = 0; i < nstep; i++) //       do 7 i = 1,nstep
            {

                if (year >= st[i]) //       if(year.ge.st(i)) delta = delta + si(i)
                    delta += si[i];

                if (year < st[i])//       if(year.lt.st(i)) return 
                    return delta;

                /* L7: */
            }   //  7    continue

            return delta;   //       return
        }

        /// <summary>
        /// Sorts an array x, of length n,and returns an array k which may be used to key another array to the sorted pattern.
        /// </summary>
        /// <param name="x">Array to be sorted.</param>
        /// <param name="k">Returns the sorted array.</param>
        /// <param name="n">Length of x.</param>
        private void shells(ref double[] x, ref int[] k, int n)                             //       subroutine shells(x,k,n)
        {
            /* Sorts an array x, of length n, sorting upward, and returns 
             * an array k which may be used to key another array to the
             * sorted pattern (i.e., if we had an array f to which x 
             * corresponded before sorting, then after calling shells, 
             * f(k(1)) will be the element of f corresponding to the 
             * smallest x, f(k(2)) the next smallest, and so on.
             * Revised 29-dec-82 so k is sorted in turn, the values of
             * k that point to identical values of x being put in increasing
             * order. */

            int j, l, ik;
            int iex = 0;
            int imax = 0;
            int igap = n;

            for (int i = 0; i < n; ++i)                                                     //       do 1 i = 1,n
                k[i] = i + 1;                                                               //  1    k(i) = i

        L5:
            if (igap <= 1)                                                                  //  5    if(igap.le.1) go to 25
                goto L25;

            igap /= 2;                                                                      //       igap = igap/2
            imax = n - igap;                                                                //       imax = n - igap

        L10:
            iex = 0;                                                                        //  10   iex = 0

            for (int i = 0; i < imax; i++)                                                  //       do 20 i = 1,imax
            {
                int ipl = i + igap;                                                         //       ipl = i + igap

                if (x[i] <= x[ipl])                                                         //       if(x(i).le.x(ipl)) go to 20 
                    continue;                                                               // goto L20;

                double sv = x[i];                                                           //       sv = x(i)
                ik = k[i];                                                                  //       ik = k(i)
                x[i] = x[ipl];                                                              //       x(i) = x(ipl)
                k[i] = k[ipl];                                                              //       k(i) = k(ipl)
                x[ipl] = sv;                                                                //       x(ipl) = sv
                k[ipl] = ik;                                                                //       k(ipl) = ik
                iex++;                                                                      //       iex = iex + 1
            }                                                                               //  20   continue

            if (iex > 0)                                                                    //       if(iex.gt.0) go to 10
                goto L10;

            goto L5;                                                                        //       go to 5

            //  now sort k's (for identical values of x, if any) 
        L25:
            j = 1;                                                                          //  25   j = 1

        L30:
            if (j >= n)                                                                     //  30   if(j.ge.n) return
                return;

            if (x[j] == x[j + 1])                                                           //       if(x(j).eq.x(j+1)) go to 33
                goto L33;

            j++;                                                                            //       j = j + 1
            goto L30;                                                                       //       go to 30

            //  have at least two x's with the same value - see how long this is true
        L33:
            l = j;                                                                          //  33   l = j

        L35:
            if (x[l] != x[l + 1])                                                           //  35   if(x(l).ne.x(l+1)) go to 38
                goto L38;

            l++;                                                                            //       l = l + 1
            if (l < n)                                                                      //       if(l.lt.n) go to 35
                goto L35;

            //  j and l are the indices within which x(i) does not change - sort k
        L38:
            igap = l - j + 1;                                                               //  38   igap = l - j + 1

        L40:
            if (igap <= 1)                                                                  //  40   if(igap.le.1) j = l + 1
                j = l + 1;

            if (igap <= 1)                                                                  //       if(igap.le.1) go to 30
                goto L30;

            igap /= 2;                                                                      //       igap = igap/2

            imax = l - j + 1 - igap;                                                        //       imax = l-j+1 - igap

        L45:
            iex = 0;                                                                        //  45   iex = 0

            for (int i = 0; i < imax; i++)                                                  //       do 50 i=1,imax
            {
                int ipl = i + igap + j - 1;                                                 //       ipl = i + igap + j - 1
                if (k[i + j - 1] <= k[ipl])                                                 //       if(k(i+j-1).le.k(ipl)) go to 50 
                    goto L50;

                ik = k[i + j - 1];                                                          //       ik = k(i+j-1)
                k[i + j - 1] = k[ipl];                                                      //       k(i+j-1) = k(ipl)
                k[ipl] = ik;                                                                //       k(ipl) = ik
                iex++;                                                                      //       iex = iex + 1

        L50:
                ;
            }                                                                               //  50   continue

            if (iex > 0)                                                                    //       if(iex.gt.0) go to 45
                goto L45;

            goto L40;                                                                       //       go to 40
        }

        /// <summary>
        /// Finds array s for spline interpolator eval.
        /// </summary>
        /// <param name="nn">Number of data points supplied (may be negative, see below)</param>
        /// <param name="x">Array containing x-coordinates where function is sampled.  xx(1),xx(2),... must be a strictly increasing sequence.</param>
        /// <param name="u">Array containing sample values that are to be interpolated.</param>
        /// <param name="s">Output array of 2nd derivative at sample points.</param>
        /// <param name="a">Working space array of dimension at least  nn.</param>
        private void spline(int nn, double[] x, double[] u, ref double[] s, double[] a)     //       subroutine spline(nn,x,u,s,a)
        {
            /*  Original date of RCS archived version is  Feb 18  2005
             *  Revision 1.1  2011/09/23 22:39:02 :  agnew
             *  Initial revision                  :  agnew */

            /*  Finds array s for spline interpolator eval.
             *  nn  Number of data points supplied (may be negative, see below).
             *  x   Array containing x-coordinates where function is sampled.  xx(1),xx(2),... 
             *      must be a strictly increasing sequence.
             *  u   Array containing sample values that are to be interpolated.
             *  s   Output array of 2nd derivative at sample points.
             *  a   Working space array of dimension at least  nn.
             *  If the user wishes to force the derivatives at the ends of the series to
             *  assume specified values, he should put du(1)/dx and du(n)/dx in s1,s2
             *  and call the routine with nn=-number of terms in series.  Normally a parabola
             *  is fitted through the 1st and last 3 points to find the slopes.
             *  If less than 4 points are supplied, straight lines are fitted. */

            /* Inline function:
             * q(u1,x1,u2,x2)=(u1/x1**2-u2/x2**2)/(1.0/x1-1.0/x2) */

            int n = Math.Abs(nn);                                                           //       n=iabs(nn)

            if (n <= 3)                                                                     //       if (n.le.3) go to 5000
                goto L5000;

            double r1 = u[2] - u[1];
            double r2 = x[2] - x[1];
            double r3 = u[3] - u[1];
            double r4 = x[3] - x[1];
            double q1 =
                (r1 / (r2 * r2) - r3 / (r4 * r4)) / (1.0 / r2 - 1.0 / r4);                  //       q1=q(u(2)-u(1),x(2)-x(1),u(3)-u(1),x(3)-x(1))

            r1 = u[n - 1] - u[n];
            r2 = x[n - 1] - x[n];
            r3 = u[n - 2] - u[n];
            r4 = x[n - 2] - x[n];
            double qn =
                (r1 / (r2 * r2) - r3 / (r4 * r4)) / (1.0 / r2 - 1.0 / r4);                  //       qn=q(u(n-1)-u(n),x(n-1)-x(n),u(n-2)-u(n),x(n-2)-x(n))

            if (nn > 0)                                                                     //       if (nn.gt.0) go to 1000
                goto L1000;

            q1 = s[0];                                                                      //       q1=s(1)
            qn = s[1];                                                                      //       qn=s(2)

        L1000:
            s[0] = ((u[1] - u[0]) / (x[1] - x[0]) - q1) * 6.0;                              //  1000 s(1)=6.0*((u(2)-u(1))/(x(2)-x(1)) - q1)

            int n1 = n - 1;                                                                 //       n1= n - 1
            for (int i = 1; i < n1; i++)                                                    //       do 2000 i=2,n1
            {
                s[i] = (u[i - 1] / (x[i] - x[i - 1]) - u[i] *
                       (1.0 / (x[i] - x[i - 1]) +
                       1.0 / (x[i + 1] - x[i])) + u[i + 1] /
                       (x[i + 1] - x[i])) * 6.0;                                            // s(i)= (u(i-1)/(x(i)-x(i-1)) - u(i)*(1.0/(x(i)-x(i-1))+1.0/(x(i+1)-x(i))) + u(i+1)/(x(i+1)-x(i)))*6.0      
            }                                                                               //  2000 continue

            s[n] = (qn + (u[n1] - u[n]) / (x[n] - x[n1])) * 6.0;                            //       s(n)=6.0*(qn + (u(n1)-u(n))/(x(n)-x(n1)))
            a[0] = (x[1] - x[0]) * 2.0;                                                     //       a(1)=2.0*(x(2)-x(1))
            a[1] = (x[1] - x[0]) * 1.5 + (x[3] - x[2]) * 2.0;                               //       a(2)=1.5*(x(2)-x(1)) + 2.0*(x(3)-x(2))
            s[2] -= s[0] * 0.5;                                                             //       s(2)=s(2) - 0.5*s(1)

            double c = 0.0;
            for (int i = 2; i < n1; i++)                                                    //       do 3000 i=3,n1
            {
                c = (x[i] - x[i - 1]) / a[i - 1];                                           //       c=(x(i)-x(i-1))/a(i-1)
                a[i] = (x[i + 1] - x[i - 1]) * 2.0 - c * (x[i] - x[i - 1]);                 //       a(i)=2.0*(x(i+1)-x(i-1)) - c*(x(i)-x(i-1))
                s[i] -= c * s[i - 1];                                                       //       s(i)=s(i) - c*s(i-1)
            }                                                                               //  3000 continue

            c = (x[n] - x[n1]) / a[n1];                                                     //       c=(x(n)-x(n1))/a(n1)
            a[n] = (2.0 - c) * (x[n] - x[n1]);                                              //       a(n)=(2.0-c)*(x(n)-x(n1))
            s[n] -= c * s[n1];                                                              //       s(n)=s(n) - c*s(n1)

            /*  back substitiute */
            s[n] /= a[n];                                                                   //       s(n)= s(n)/a(n)
            for (int j = 1; j <= n1; ++j)                                                   //       do 4000 j=1,n1
            {
                int i = n - j;                                                              //       i=n-j
                s[i] = (s[i] - (x[i + 1] - x[i]) * s[i + 1]) / a[i];                        //       s(i) =(s(i) - (x(i+1)-x(i))*s(i+1))/a(i)
            }                                                                               //  4000 continue
            return;                                                                         //       return
        /*  series too short for cubic spline - use straight lines. */
        L5000:
            for (int i = 0; i < n; i++)                                                     //  5000 do 5500 i=1,n
                s[i] = 0.0;                                                                 //  5500 s(i)=0.0
        }

        /// <summary>
        /// Public variable for function eval perform the speed.
        /// </summary>
        public int istart;

        /// <summary>
        /// Performs cubic spline interpolation of a function sampled at unequally spaced intervals.  The routine spline  should be called to set up the array s.
        /// </summary>
        /// <param name="y">The coordinate at which function value is desired.</param>
        /// <param name="nn">Number of samples of original function. </param>
        /// <param name="x">Array containing sample coordinates. the sequence x(1),x(2).....x(nn) must be strictly increasing.</param>
        /// <param name="u">Array containing samples of function at the coordinates x(1),x(2)... </param>
        /// <param name="s">array containing the 2nd derivatives at the sample points. Found by the routine  spline, which must be called once before beginning interpolation.</param>
        /// <returns>The interpolated value.</returns>
        private double eval(double y, int nn, double[] x, double[] u, double[] s)           //       function eval(y,nn,x,u,s)
        {
            /*  Performs cubic spline interpolation of a function sampled at unequally 
             *  spaced intervals.  The routine spline should be called to set up the array s.
             *  y   the coordinate at which function value is desired.
             *  nn  number of samples of original function. 
             *  x   Array containing sample coordinates. the sequence x(1),x(2).....x(nn) must be strictly increasing.
             *  u   Array containing samples of function at the coordinates x(1),x(2)... 
             *  s   Array containing the 2nd derivatives at the sample points. Found by the 
             *      routine  spline, which must be called once before beginning interpolation.
             *  If  y  falls outside range(x(1),x(nn))  the value at the nearest endpoint 
             *  of the series is used. */
            int k;
            int k1 = 0;
            nn--;                                                                           // Additionally, introduced for index adjustment. 
            nn = Math.Abs(nn);                                                              //       nn=iabs(nn)

            if (y < x[0])                                                                   //       if (y.lt.x(1))  go to 3000
                goto L3000;

            if (y > x[nn])                                                                  //       if (y.gt.x(nn)) go to 3100
                goto L3100;

            /*  locate interval (x(k1),x(k))  containing y */
            if (y - x[istart] >= 0.0)                                                       //       if (y-x(istart)) 1200,1000,1000
                goto L1000;
            else
                goto L1200;

        /*  scan up the x array */
        L1000:
            for (k = istart; k <= nn; ++k)                                                  //  1000 do 1100 k=istart,nn
            {
                if (x[k] > y)                                                               //       if (x(k).gt.y) go to 1150
                    goto L1150;
            }                                                                               //  1100 continue

        L1150:
            k1 = k - 1;                                                                     //  1150 k1=k-1
            goto L1500;                                                                     //       go to 1500
        /*  scan downwards in x array */

        L1200:
            for (k = 1; k <= istart; ++k)                                                   //  1200 do 1300 k=1,istart
            {
                k1 = istart - k;                                                            //       k1=istart-k

                if (x[k1] <= y)                                                             //       if (x(k1).le.y) go to 1350
                    goto L1350;
            }                                                                               //  1300 continue

        L1350:
            k = k1 + 1;                                                                     //  1350 k=k1+1

        L1500:
            istart = k1;                                                                    //  1500 istart=k1
            /*  evaluate interpolate */
            double dy = x[k] - y;                                                           //       dy=x(k)-y
            double dy1 = y - x[k1];                                                         //       dy1=y-x(k1)
            double dk = x[k] - x[k1];                                                       //       dk=x(k)-x(k1)
            double deli = 1.0 / (dk * 6.0);                                                 //       deli=1./(6.0*dk)
            double ff1 = s[k1] * dy * dy * dy;                                              //       ff1=s(k1)*dy*dy*dy
            double ff2 = s[k] * dy1 * dy1 * dy1;                                            //       ff2=s(k)*dy1*dy1*dy1
            double f1 = (ff1 + ff2) * deli;                                                 //       f1=(ff1+ff2)*deli
            double f2 = dy1 * (u[k] / dk - s[k] * dk / 6.0);                                //       f2=dy1*((u(k)/dk)-(s(k)*dk)/6.)
            double f3 = dy * (u[k1] / dk - s[k1] * dk / 6.0);                               //       f3=dy*((u(k1)/dk)-(s(k1)*dk)/6.0)
            return f1 + f2 + f3;                                                            //       eval=f1+f2+f3

            /*  out of range.  substitute end values. */
        L3000:
            return u[0];                                                                    //  3000 eval=u(1)

        L3100:
            return u[nn];                                                                   //  3100 eval=u(nn)
        }

        /// <summary>
        /// Copies an array to a returnded starting at given index.
        /// </summary>
        /// <param name="src">Source array to be copied.</param>
        /// <param name="ind">Index of start copieing the array.</param>
        /// <returns>Copied array.</returns>
        private double[] cpfrom(double[] src, int ind)
        {
            double[] tmp = new double[src.Length];

            for (int i = ind; i < src.Length; i++)
                tmp[i - ind] = src[i];

            return tmp;
        }
        #endregion

        #region Private common fucntions
        /// <summary>
        /// Returns, as a complex number, the amplitude and phase of the complex quantity input.
        /// </summary>
        /// <param name="dumm">Complex number.</param>
        /// <returns>Amplitude and phase of the complex quantity input.</returns>
        private Complex phasor(Complex dumm)                                                // function phasor(dumm)
        {
            Complex phasor = new Complex();
            if (dumm == Complex.Zero)                                                       // if(dumm.eq.cz) phasor=cz
                phasor = Complex.Zero;
            if (dumm != Complex.Zero)                                                       // if(dumm.ne.cz) phasor= cmplx(cabs(dumm),57.2958*atan2(aimag(dumm),real(dumm)))
                phasor = new Complex(
                    Complex.Abs(dumm),
                    Constants.radian2degree * Math.Atan2(dumm.Imaginary, dumm.Real));
            return phasor;
        }

        /// <summary>
        /// Returns, as amplitude and phase, as a complex number.
        /// </summary>
        /// <param name="c">>Amplitude and phase of the complex quantity input</param>
        /// <returns>Complex number.</returns>
        private Complex argand(Complex c)                                                   // function phasor(dumm)
        {
            double p = c.Imaginary / Constants.radian2degree;
            double a = c.Real;
            return new Complex(a * Math.Cos(p), a * Math.Sin(p));
        }
        #endregion
    }
}