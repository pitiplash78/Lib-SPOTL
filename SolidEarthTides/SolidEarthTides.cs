using System;
using System.Collections.Generic;
using System.Text;

namespace SPOTL
{
    /// <summary>
    /// Calculates the Solid Earth Tides, depend on the Lovenumbers, Station properties and Time
    /// </summary>
    public class SolidEarthTides
    {
        /* formaly program erthtd.f
         * 
         * BOMM subroutine written by Jon Berger  November 1969 
         * astronomy revised by Judah Levine (after J C Harrison) Sept 1973
         * updated by Karen Young  March 1977, for PDP 11
         * emended by Duncan Agnew June 1978 to speed execution
         * solar astronomy redone and rigid-earth potential added by Duncan Agnew June 1979
         * tide generating part put in separate subroutine and computation of Munk-Cartwright coefficients added by Duncan Agnew Feb 1981,* July 1982
         * 
         * This version rewritten for export, using F77 calls for I/O, by Duncan Agnew Mar 1987
         * 
         * tides are calculated from harmonics 2 through 3 for the lunar terms, 2 for the solar terms.
         * love numbers h(n), k(n), and l(n) are set by data statements.
         * 
         * gravity tide is in microgals, plus for up acceleration
         * rigid-earth potential is in meters (v/g), plus for up
         * tilt tide is in nanoradians (nr)......  1 nr = 2.063e-4 arcsec
         * strain tide is in parts in e9, plus for extension
         * 
         * azimuths for strain and tilt are reckoned
         * positive clockwise from north (east from north) */

        internal Location location = null;
        internal LoveNumbers loveNumber = null;

        /// <summary>
        /// Constuctor for calcuation Solid Earth Tide
        /// </summary>
        /// <param name="stationlongitude"></param>
        /// <param name="stationlatitude"></param>
        /// <param name="stationheight"></param>
        /// <param name="h">Lovenumbers</param>
        /// <param name="k">Lovenumbers</param>
        public SolidEarthTides(Location Location, LoveNumbers LoveNumbers)
        {
            location = new Location()
            {
                Longitude = Location.Longitude,
                Latitude = Location.Latitude,
                Height = Location.Height,
            };

            if (LoveNumbers != null)
            {
                loveNumber = new LoveNumbers()
                {
                    h = LoveNumbers.h,
                    k = LoveNumbers.k,
                    l = LoveNumbers.l,
                };
            }
            else
            {
                loveNumber = new LoveNumbers()
                {
                    h = new double[] { Constants.h1Default, Constants.h2Default, Constants.h3Default },
                    k = new double[] { Constants.k1Default, Constants.k2Default, Constants.k3Default },
                    l = new double[] { Constants.l1Default, Constants.l2Default, Constants.l2Default },
                };
            }
            // finds geocentric position
            sph();
        }

        #region Enum Component
        public enum Component
        {
            None,
            Gravity,
            GravityPotentialHeight,
            GravityAttraction,
            Tilt,
            TiltDeformation,
            TiltAttraction, 
            HorizontalStrain,
            ArealStrain,
            VolumeStrain,
        }
        #endregion

        #region Parameter for constructor
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
        internal Units.Units.UnitNamesEnum _OutputUnit_Gravity = BaseUnits.UnitGravity;
        internal double dConversion_Gravity = 1d;

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
        internal Units.Units.UnitNamesEnum _OutputUnit_Tilt = BaseUnits.UnitTilt;
        internal double dConversion_Tilt = 1d;

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
        internal Units.Units.UnitNamesEnum _OutputUnit_Strain = BaseUnits.UnitStrain;
        internal double dConversion_Strain = 1d;

        public void setUnit(Component component, Units.Units.UnitNamesEnum Unit)
        {
            switch (component)
            {
                case SPOTL.SolidEarthTides.Component.Gravity:
                case SPOTL.SolidEarthTides.Component.GravityAttraction:
                case SPOTL.SolidEarthTides.Component.GravityPotentialHeight:
                    OutputUnit_Gravity = Unit;
                    break;
                case SPOTL.SolidEarthTides.Component.HorizontalStrain:
               // TODO case SPOTL.SolidEarthTides.Component.ShearStrain:
                case SPOTL.SolidEarthTides.Component.ArealStrain:
                case SPOTL.SolidEarthTides.Component.VolumeStrain:
                    OutputUnit_Strain = Unit;
                    break;
                case SPOTL.SolidEarthTides.Component.Tilt:
                case SPOTL.SolidEarthTides.Component.TiltAttraction:
                case SPOTL.SolidEarthTides.Component.TiltDeformation:
                    OutputUnit_Tilt = Unit;
                    break;
            }
        }
        #endregion

        #region Base Units of Components
        /// <summary>
        /// Base Units for according components 
        /// </summary>
        public static class BaseUnits
        {
            public static Units.Units.UnitNamesEnum UnitGravity = Units.Units.UnitNamesEnum.MeterPerSecondSquared;
            public static Units.Units.UnitNamesEnum UnitTilt = Units.Units.UnitNamesEnum.Radian;
            public static Units.Units.UnitNamesEnum UnitStrain = Units.Units.UnitNamesEnum.Strain;
        }
        #endregion

        #region internal common definitions of used variables
        /// <summary>
        /// contains th observation information
        /// </summary>
        internal static class obs
        {
            internal static double cthet = 0d;
            internal static double sthet = 0d;
            internal static double clong = 0d;
            internal static double slong = 0d;
            internal static double dvert = 0d;
            internal static double radn = 0d;
            internal static double g = 0d;
        }

        internal static class bpos
        {
            /// <summary>
            /// sine of colatitude of sublunar point
            /// </summary>
            internal static double dsz;
            
            /// <summary>
            /// cosine of colatitude of sublunar point
            /// </summary>
            internal static double dcz;
            
            /// <summary>
            /// sine of east longitude of sublunar point
            /// </summary>
            internal static double dsl;
            
            /// <summary>
            /// cosine of east longitude of sublunar point
            /// </summary>
            internal static double dcl;
            
            /// <summary>
            /// sine of colatitude of subsolar point
            /// </summary>
            internal static double ssz;
            
            /// <summary>
            /// cosine of colatitude of subsolar point
            /// </summary>
            internal static double scz;
            
            /// <summary>
            ///  sine of east longitude of subsolar point
            /// </summary>
            internal static double ssl;
            
            /// <summary>
            ///  cosine of east longitude of subsolar point
            /// </summary>
            internal static double scl;
             /// <summary>
            /// solar distance in astronomical units
            /// </summary>
            internal static double dpar;
            /// <summary>
            /// solar distance in astronomical units
            /// </summary>
            internal static double sdist;
        }

        internal static class tdiff
        {
            /// <summary>
            /// containing the difference ephemeris minus universal time, in seconds. if this is not known it should
            /// be set to zero, and the argument to the program should be universal rather than ephemeris time.
            /// </summary>
            internal static double etmut = 0d;
        }

        internal static class constants
        {
            internal static double[] rbor = new double[2] { 1.6592496e-2, 4.2635233e-5 };
            internal static double[] amrat = new double[2] { 78451.25d, 2.1235762e12 };
            internal static double a = 6.37814e6;
            internal static double g1 = 9.79828d;
            internal static double g2 = 9.82022d;
        }
        #endregion

        #region private methods
        private void sph()
        {
            /* for a point at geographical north latitude grlat, east longitude elong      
             * (in degrees), and height ht (in meters), finds geocentric position   
             * and local g using formulae for a spheroid. 
             * returns the caluclated values in class obs */

            double gn = 9.798277692d;
            double ae = 6378140.0d;
            double f = 0.00335281d;
            double rm = 0.00344978d;

            obs.clong = Math.Cos(location.Longitude * Constants.degree2radian);
            obs.slong = Math.Sin(location.Longitude * Constants.degree2radian);
            
            //   latitude difference 
            obs.dvert = f * (1.0d + 0.5d * f) * Math.Sin(2.0d * location.Latitude * Constants.degree2radian) -
                        0.5d * f * f * Math.Sin(4.0d * location.Latitude * Constants.degree2radian);
            double gcclat = (Math.PI / 2.0d) - (location.Latitude * Constants.degree2radian - obs.dvert);
            obs.cthet = Math.Cos(gcclat);
            obs.sthet = Math.Sin(gcclat);
            
            // geocentric radius   
            obs.radn = 1.0d - f * Math.Pow(obs.cthet, 2d) * (1.0d + 1.5d * f * (Math.Pow(obs.sthet, 2d)));
            
            // formulae for g are from jeffreys, 4.022 and 4.023      
            obs.g = gn * (1.0d + f - 1.5d * rm + f * (f - (27.0d / 14.0d) * rm) +
                (2.5d * rm - f - f * (f - (39.0d / 14.0d) * rm)) * Math.Pow(obs.cthet, 2d) -
                (f / 2.0d) * (7.0d * f - 15.0d * rm) * Math.Pow((obs.cthet * obs.sthet), 2d));
            
            // free air correction 
            obs.g = obs.g - obs.g * (2.0d * location.Height * (1.0d + f + rm - 2.0d * f * Math.Pow(obs.cthet, 2d)) / ae);
        }
       
        private void ephem(double t)
        {
           /* t is ephemeris time in julian centuries from 12 hr 0 jan 1900
            * (in ordinary civil reckoning this is greenwich noon on 31 december
            * 1899). if the difference between ephemeris and unversal time is
            * not put in (see common tdiff below), t should be in universal time
            * which is (nearly) the time ordinarily kept by clocks.
            * computes positions of the sun and moon at time t, returning results
            * in common block bpos. the solar ephemeris uses the mean sun.
            * Derived from J. Levines revision (after J. C. Harrison)
            * of an earthtide program by J. Berger and W. E. farrell, with small
            * alterations by D. C. Agnew, partly after M. Wimbush. present
            * subroutine version by D. C. Agnew. */

            // compute universal time in hours
            double ts = 876600.0d * t - 12.0d - (tdiff.etmut / 3600.0d);
            double hr = ts % 24.0d;
            if (hr == 0d) hr = 24d;

            // compute obliquity of the ecliptic
            double w = 0.409319747d - 0.0002271107d * t;

            // compute solar constants for given t
            double t2 = t * t;
            double hs = 4.881627982482d + 628.3319508731d * t + 0.523598775578e-5 * t2;
            double pi20 = 62.8318530717958d;
            hs = ((hs % pi20) + pi20) % pi20;
            double ps = 4.908229466993d + 0.03000526416690d * t + 0.790246300201e-5 * t2;
            double es = 0.01675104d - 0.00004180d * t - 0.000000126d * t2;
            double psig = 0.2617993877971d * (hr - 12.0d) + hs;
            double chmp = Math.Cos(hs - ps);
            double shmp = Math.Sin(hs - ps);
            double ls = hs + shmp * es * (2.0d + 2.5d * es * chmp);
            double cz = Math.Sin(w) * Math.Sin(ls);
            double sz = Math.Sqrt(1.0d - cz * cz);
            double psis = Math.Atan2(Math.Cos(w) * Math.Sin(ls), Math.Cos(ls));
            double rbarr = 1.0d + es * (chmp + es * (chmp - shmp) * (chmp + shmp));
            bpos.ssz = sz;
            bpos.scz = cz;
            bpos.ssl = Math.Sin(psis - psig); 
            bpos.scl = Math.Cos(psis - psig);
            bpos.sdist = 1.0d / rbarr;

            // compute lunar constants for given t
            double hm = 4.7199666d + 8399.7091449d * t - 0.0000198d * t2;
            double pm = 5.83515154d + 71.01804120839d * t - 0.180205e-3 * t2;
            double nm = 4.523601515d - 33.75714624d * t + 0.3626406335e-4 * t2;
            // bl bls bf bd are the fundamental arguments of browns theory
            double bl = hm - pm;
            double bls = hs - ps;
            double bf = hm - nm;
            double bd = hm - hs;
            // Lunar lat long and parallax from brown. latter two from improved lunar ephemeris, latitude from ras paper of 1908...
            double tlongm = hm + 0.10976d * Math.Sin(bl) - 0.02224d * Math.Sin(bl - 2d * bd) +
                0.01149d * Math.Sin(2d * bd) + 0.00373d * Math.Sin(2d * bl) - 0.00324d * Math.Sin(bls) -
                0.00200d * Math.Sin(2d * bf) - 0.00103d * Math.Sin(2d * bl - 2d * bd) -
                0.00100d * Math.Sin(bl + bls - 2d * bd) + 0.00093d * Math.Sin(bl + 2d * bd) -
                0.00080d * Math.Sin(bls - 2d * bd) + 0.00072d * Math.Sin(bl - bls) -
                0.00061d * Math.Sin(bd) - 0.00053d * Math.Sin(bl + bls);
            double tlatm = 0.08950d * Math.Sin(bf) + 0.00490d * Math.Sin(bl + bf) -
                0.00485d * Math.Sin(bf - bl) - 0.00303d * Math.Sin(bf - 2d * bd) +
                0.00097d * Math.Sin(2d * bd + bf - bl) - 0.00081d * Math.Sin(bl + bf - 2d * bd) +
                0.00057d * Math.Sin(bf + 2d * bd);
            double plx = (3422.45d + 186.54d * Math.Cos(bl) + 34.31d * Math.Cos(bl - 2d * bd) +
                28.23d * Math.Cos(2d * bd) + 10.17d * Math.Cos(2d * bl) + 3.09d * Math.Cos(bl + 2d * bd) +
                1.92d * Math.Cos(bls - 2d * bd) + 1.44d * Math.Cos(bl + bls - 2d * bd) +
                1.15d * Math.Cos(bl - bls) - 0.98d * Math.Cos(bd) - 0.95d * Math.Cos(bl + bls) -
                0.71d * Math.Cos(bl - 2d * bf) + 0.62d * Math.Cos(3d * bl) + 0.60d * Math.Cos(bl - 4d * bd));
            double sinmla = Math.Sin(tlatm);
            double cosmla = Math.Cos(tlatm);
            double sinmln = Math.Sin(tlongm);
            double cosmln = Math.Cos(tlongm);

            // convert from celestial lat and long according to explaiMOn suppl of na and le page 26
            bpos.dcz = cosmla * sinmln * Math.Sin(w) + sinmla * Math.Cos(w);
            bpos.dsz = Math.Sqrt(1d - bpos.dcz * bpos.dcz);
            double ll = Math.Atan2(cosmla * sinmln * Math.Cos(w) - sinmla * Math.Sin(w), cosmla * cosmln) - psig;
            bpos.dsl = Math.Sin(ll); 
            bpos.dcl = Math.Cos(ll); 
            bpos.dpar = plx;
        }

        // calculation single component with given component
        private double getSolidEarthTides(double time, Component component, double azimuth)
        {
            //   computes positions of the sun and moon at time t:
            ephem(time);

            double[] az = null;

            //  computes the earth tides on an elastic earth, given the solar and lunar positions
            double value = double.NaN;
            if (component == Component.Gravity)
            {
                value = elastdGravity();
            }
            else if (component == Component.GravityAttraction)
            {
                value = elastdGravity() / (1 - 3 / 2 * loveNumber.k[0] + loveNumber.h[0]);
            }
            else if (component == Component.GravityPotentialHeight)
            {
                value = elastdGravity() / (1 - 3 / 2 * loveNumber.k[0] + loveNumber.h[0]) * (3 / 2 * loveNumber.k[0] + loveNumber.h[0]);
            }
            else if (component == Component.Tilt)
            {
                az = new double[] { azimuth, Math.Cos(azimuth * Constants.degree2radian), Math.Sin(azimuth * Constants.degree2radian) };
                value = elastdTiltSingleComponent(az);
            }
            else if (component == Component.TiltDeformation)
            {
                az = new double[] { azimuth, Math.Cos(azimuth * Constants.degree2radian), Math.Sin(azimuth * Constants.degree2radian) };
                double t = elastdTiltSingleComponent(az);
                value = t * (-loveNumber.h[0] + loveNumber.l[0]) / (1 + loveNumber.k[0] - loveNumber.h[0]);
            }
            else if (component == Component.TiltAttraction)
            {
                az = new double[] { azimuth, Math.Cos(azimuth * Constants.degree2radian), Math.Sin(azimuth * Constants.degree2radian) };
                double t = elastdTiltSingleComponent(az);
                value = t * (1 + loveNumber.k[0] - loveNumber.l[0]) / (1 + loveNumber.k[0] - loveNumber.h[0]);
            }
            else if (component == Component.HorizontalStrain)
            {
                az = new double[] { azimuth, Math.Cos(azimuth * Constants.degree2radian), Math.Sin(azimuth * Constants.degree2radian) };
                value = elastdStrainSingleComponent(az);
            }
            else if (component == Component.ArealStrain)
            {
                az = new double[] { 0.0, Math.Cos(0.0 * Constants.degree2radian), Math.Sin(0.0 * Constants.degree2radian) };
                double value1 = elastdStrainSingleComponent(az);
                az = new double[] { 90.0, Math.Cos(90.0 * Constants.degree2radian), Math.Sin(90.0 * Constants.degree2radian) };
                double value2 = elastdStrainSingleComponent(az);
                value = value1 + value2;
            }
            else if (component == Component.VolumeStrain)
            {
                az = new double[] { 0.0, Math.Cos(0.0 * Constants.degree2radian), Math.Sin(0.0 * Constants.degree2radian) };
                double value1 = elastdStrainSingleComponent(az);
                az = new double[] { 90.0, Math.Cos(90.0 * Constants.degree2radian), Math.Sin(90.0 * Constants.degree2radian) };
                double value2 = elastdStrainSingleComponent(az);
                value = 0.6667 * (value1 + value2);
            }

            return value;
        }

        private double elastdGravity()
        {
            /* dc gravity tide is also known as the Honkasalo correction
             * **note that the love numbers for an elastic earth are used in computing the dc tide as well*/
            double gdc = -3.0481e-7 * (3d * obs.cthet * obs.cthet - 1d) * (1d + (2d / 2d) * loveNumber.h[0] - (3d / 2d) * loveNumber.k[0]) * obs.radn;
            double tnsdc = -9.1445e-7 * obs.cthet * obs.sthet * (1d + loveNumber.k[0] - loveNumber.h[0]) * obs.radn / obs.g;
            double etdc = -1.555e-8 * (loveNumber.h[0] * (3d * obs.cthet * obs.cthet - 1d) - 6d * loveNumber.l[0] * (2d * obs.cthet * obs.cthet - 1d));
            double eldc = -1.555e-8 * (loveNumber.h[0] * (3d * obs.cthet * obs.cthet - 1d) - 6d * loveNumber.l[0] * obs.cthet * obs.cthet);
            double potdc = 0.0992064d * (1d - 3d * obs.cthet * obs.cthet);
            double re = 1d / (obs.radn * constants.a);

            double grav = 0d;
            double gnth = 0d;

            double[] p = new double[3];
            double[] pp = new double[3];

            double[] coor = new double[] { bpos.dsz, bpos.dcz, bpos.dsl, bpos.dcl, bpos.ssz, bpos.scz, bpos.ssl, bpos.scl };

            //  compute normalized parallax
            double[] pa = new double[2] { bpos.dpar / 3422.45d, 1d / bpos.sdist };

            for (int ii = 0; ii < 2; ii++)
            {
                int id = 3;
                if (ii == 1) id = 1;
                int ir = (4 * ii) - 1;
                /*  find cosine of zenith angle, potential constants, legendre polynomials
                 * and their derivatives, and derivatives of the cosine of the zenith angle. */

                double cll = obs.clong * coor[ir + 4] + obs.slong * coor[ir + 3];
                double sll = obs.slong * coor[ir + 4] - obs.clong * coor[ir + 3];
                double cz = coor[ir + 2];
                double sz = coor[ir + 1];
                double cu = obs.cthet * cz + obs.sthet * sz * cll;
                double xi = constants.rbor[ii] * pa[ii] * obs.radn;
                double cc = constants.amrat[ii] * constants.rbor[ii] * pa[ii];
                double[] rkr = new double[3] { cc * Math.Pow(xi, 2), cc * Math.Pow(xi, 3), cc * Math.Pow(xi, 4) };
                p[0] = 0.5d * (3d * cu * cu - 1d);
                pp[0] = 3d * cu;

                if (ii != 1)
                {
                    p[1] = 0.5d * cu * (5d * cu * cu - 3d);
                    p[2] = 0.25d * (7d * cu * p[1] - 3d * p[0]);
                    pp[1] = 1.5d * (5d * cu * cu - 1d);
                    pp[2] = 0.25d * (7d * p[1] + cu * pp[1]) - 3d * pp[0];
                }
                double cut = -obs.sthet * cz + obs.cthet * sz * cll;
                double cutt = -cu;
                double cul = -obs.sthet * sz * sll;
                double cull = -obs.sthet * sz * cll;
                double cutl = -obs.cthet * sz * sll;

                // gravity
                for (int j = 0; j < id; j++)
                    grav += (1d + (2d / ((double)j + 2d)) * loveNumber.h[j] - (((double)j + 3d) / ((double)j + 2d)) * loveNumber.k[j]) *
                            (j + 2) * rkr[j] * p[j] * constants.g1 * re;
                gnth -= (1d + loveNumber.k[0] - loveNumber.h[0]) * rkr[0] * pp[0] * constants.g1 * cut * re;
            }
            // gravity
            grav += (gnth * obs.dvert - gdc);

            return -grav * dConversion_Gravity;
        }

        private double elastdTiltSingleComponent(double[] az)
        {
            /* dc gravity tide is also known as the Honkasalo correction
             * **note that the love numbers for an elastic earth are used in computing the dc tide as well */
            double tnsdc = -9.1445e-7 * obs.cthet * obs.sthet * (1d + loveNumber.k[0] - loveNumber.h[0]) * obs.radn / obs.g;
            double re = 1d / (obs.radn * constants.a);

            double tilt = 0d;
            double tltcor = 0d;
            double[] p = new double[3];
            double[] pp = new double[3];

            double[] coor = new double[] { bpos.dsz, bpos.dcz, bpos.dsl, bpos.dcl, bpos.ssz, bpos.scz, bpos.ssl, bpos.scl };

            //  compute normalized parallax
            double[] pa = new double[2] { bpos.dpar / 3422.45d, 1d / bpos.sdist };

            for (int ii = 0; ii < 2; ii++)
            {
                int id = 3;
                if (ii == 1) id = 1;
                int ir = (4 * ii) - 1;
                /* find cosine of zenith angle, potential constants, legendre polynomials
                 * and their derivatives, and derivatives of the cosine of the zenith angle. */

                double cll = obs.clong * coor[ir + 4] + obs.slong * coor[ir + 3];
                double sll = obs.slong * coor[ir + 4] - obs.clong * coor[ir + 3];
                double cz = coor[ir + 2];
                double sz = coor[ir + 1];
                double cu = obs.cthet * cz + obs.sthet * sz * cll;
                double xi = constants.rbor[ii] * pa[ii] * obs.radn;
                double cc = constants.amrat[ii] * constants.rbor[ii] * pa[ii];
                double[] rkr = new double[3] { cc * Math.Pow(xi, 2), cc * Math.Pow(xi, 3), cc * Math.Pow(xi, 4) };
                p[0] = 0.5d * (3d * cu * cu - 1d);
                pp[0] = 3d * cu;

                if (ii != 1)
                {
                    p[1] = 0.5d * cu * (5d * cu * cu - 3d);
                    p[2] = 0.25d * (7d * cu * p[1] - 3d * p[0]);
                    pp[1] = 1.5d * (5d * cu * cu - 1d);
                    pp[2] = 0.25d * (7d * p[1] + cu * pp[1]) - 3d * pp[0];
                }
                double cut = -obs.sthet * cz + obs.cthet * sz * cll;
                double cul = -obs.sthet * sz * sll;

                for (int j = 0; j < id; j++)
                    tilt = tilt - ((1d + loveNumber.k[j] - loveNumber.h[j]) * rkr[j] * pp[j] * re) * (cut * az[1] - cul * az[2] / obs.sthet);
                tltcor = tltcor + (1d + (2d / (1d + 1d)) * loveNumber.h[0] - ((1 + 2d) / (1 + 1d)) * loveNumber.k[0]) * 2d * rkr[0] * p[0] * re;
            }

            tilt = tilt * (constants.g1 / obs.g) - tltcor * obs.dvert * az[2] - tnsdc * az[2];

            return tilt * dConversion_Tilt;
        }

        private double elastdStrainSingleComponent(double[] az)
        {
            /* dc gravity tide is also known as the Honkasalo correction
             * **note that the love numbers for an elastic earth are used in computing the dc tide as well*/
            double gdc = -3.0481e-7 * (3d * obs.cthet * obs.cthet - 1d) * (1d + (2d / 2d) * loveNumber.h[0] - (3d / 2d) * loveNumber.k[0]) * obs.radn;
            double tnsdc = -9.1445e-7 * obs.cthet * obs.sthet * (1d + loveNumber.k[0] - loveNumber.h[0]) * obs.radn / obs.g;
            double etdc = -1.555e-8 * (loveNumber.h[0] * (3d * obs.cthet * obs.cthet - 1d) - 6d * loveNumber.l[0] * (2d * obs.cthet * obs.cthet - 1d));
            double eldc = -1.555e-8 * (loveNumber.h[0] * (3d * obs.cthet * obs.cthet - 1d) - 6d * loveNumber.l[0] * obs.cthet * obs.cthet);
            double potdc = 0.0992064d * (1d - 3d * obs.cthet * obs.cthet);
            double re = 1d / (obs.radn * constants.a);

            double strain = 0d;
            double[] e = new double[3];

            double[] p = new double[3];
            double[] pp = new double[3];
            double[] ppp = new double[3] { 3d, 0d, 0d };

            double[] coor = new double[] { bpos.dsz, bpos.dcz, bpos.dsl, bpos.dcl, bpos.ssz, bpos.scz, bpos.ssl, bpos.scl };

            //  compute normalized parallax
            double[] pa = new double[2] { bpos.dpar / 3422.45d, 1d / bpos.sdist };

            for (int ii = 0; ii < 2; ii++)
            {
                int id = 3;
                if (ii == 1) id = 1;
                int ir = (4 * ii) - 1;
                /*  find cosine of zenith angle, potential constants, legendre polynomials
                 * and their derivatives, and derivatives of the cosine of the zenith angle. */

                double cll = obs.clong * coor[ir + 4] + obs.slong * coor[ir + 3];
                double sll = obs.slong * coor[ir + 4] - obs.clong * coor[ir + 3];
                double cz = coor[ir + 2];
                double sz = coor[ir + 1];
                double cu = obs.cthet * cz + obs.sthet * sz * cll;
                double xi = constants.rbor[ii] * pa[ii] * obs.radn;
                double cc = constants.amrat[ii] * constants.rbor[ii] * pa[ii];
                double[] rkr = new double[3] { cc * Math.Pow(xi, 2), cc * Math.Pow(xi, 3), cc * Math.Pow(xi, 4) };
                p[0] = 0.5d * (3d * cu * cu - 1d);
                pp[0] = 3d * cu;

                if (ii != 1)
                {
                    p[1] = 0.5d * cu * (5d * cu * cu - 3d);
                    p[2] = 0.25d * (7d * cu * p[1] - 3d * p[0]);
                    pp[1] = 1.5d * (5d * cu * cu - 1d);
                    pp[2] = 0.25d * (7d * p[1] + cu * pp[1]) - 3d * pp[0];
                    ppp[1] = 15d * cu;
                    ppp[2] = 7.5d * (7d * cu * cu - 1d);
                }
                double cut = -obs.sthet * cz + obs.cthet * sz * cll;
                double cutt = -cu;
                double cul = -obs.sthet * sz * sll;
                double cull = -obs.sthet * sz * cll;
                double cutl = -obs.cthet * sz * sll;

                // strain
                for (int j = 0; j < id; j++)
                {
                    e[0] += (rkr[j] * (loveNumber.l[j] * (pp[j] * cut * obs.cthet / obs.sthet + (ppp[j] * cul * cul +
                                pp[j] * cull) / (obs.sthet * obs.sthet)) + loveNumber.h[j] * p[j]));
                    e[1] += (rkr[j] * (loveNumber.l[j] * (ppp[j] * cut * cut + pp[j] * cutt) +
                                loveNumber.h[j] * p[j]));
                }
            }
            // strain
                strain = (e[0] * az[2] * az[2] - e[2] * az[2] * az[1] + e[1] * az[1] * az[1]) *
                             re * (constants.g1 / constants.g2) - etdc * az[1] * az[1] - eldc * az[2] * az[2];

            return strain * dConversion_Strain;
        }

        // calculation for two tilt channels, according to given direction
        private double[] elastdTilt(double[,] azt)
        {
            /* dc gravity tide is also known as the Honkasalo correction
             * **note that the love numbers for an elastic earth are used in computing the dc tide as well */
            double tnsdc = -9.1445e-7 * obs.cthet * obs.sthet * (1d + loveNumber.k[0] - loveNumber.h[0]) * obs.radn / obs.g;
            double re = 1d / (obs.radn * constants.a);

            double[] tilt = new double[2];
            double[] tltcor = new double[2];
            double[] p = new double[3];
            double[] pp = new double[3];

            double[] coor = new double[] { bpos.dsz, bpos.dcz, bpos.dsl, bpos.dcl, bpos.ssz, bpos.scz, bpos.ssl, bpos.scl };

            //  compute normalized parallax
            double[] pa = new double[2] { bpos.dpar / 3422.45d, 1d / bpos.sdist };

            for (int ii = 0; ii < 2; ii++)
            {
                int id = 3;
                if (ii == 1) id = 1;
                int ir = (4 * ii) - 1;
                /*  find cosine of zenith angle, potential constants, legendre polynomials
                 * and their derivatives, and derivatives of the cosine of the zenith angle. */
                double cll = obs.clong * coor[ir + 4] + obs.slong * coor[ir + 3];
                double sll = obs.slong * coor[ir + 4] - obs.clong * coor[ir + 3];
                double cz = coor[ir + 2];
                double sz = coor[ir + 1];
                double cu = obs.cthet * cz + obs.sthet * sz * cll;
                double xi = constants.rbor[ii] * pa[ii] * obs.radn;
                double cc = constants.amrat[ii] * constants.rbor[ii] * pa[ii];
                double[] rkr = new double[3] { cc * Math.Pow(xi, 2), cc * Math.Pow(xi, 3), cc * Math.Pow(xi, 4) };
                p[0] = 0.5d * (3d * cu * cu - 1d);
                pp[0] = 3d * cu;

                if (ii != 1)
                {
                    p[1] = 0.5d * cu * (5d * cu * cu - 3d);
                    p[2] = 0.25d * (7d * cu * p[1] - 3d * p[0]);
                    pp[1] = 1.5d * (5d * cu * cu - 1d);
                    pp[2] = 0.25d * (7d * p[1] + cu * pp[1]) - 3d * pp[0];
                }

                double cut = -obs.sthet * cz + obs.cthet * sz * cll;
                double cul = -obs.sthet * sz * sll;

                for (int kk = 0; kk < 2; kk++)
                {
                    for (int j = 0; j < id; j++)
                    {
                        tilt[kk] = tilt[kk] - ((1d + loveNumber.k[j] - loveNumber.h[j]) * rkr[j] * pp[j] * re) *
                            (cut * azt[1, kk] - cul * azt[2, kk] / obs.sthet);
                    }
                    tltcor[kk] = tltcor[kk] +
                            (1d + (2d / (1d + 1d)) * loveNumber.h[0] - ((1 + 2d) / (1 + 1d)) * loveNumber.k[0]) * 2d * rkr[0] * p[0] * re;
                }
            }

            for (int kk = 0; kk < 2; kk++)
                tilt[kk] = tilt[kk] * (constants.g1 / obs.g) - tltcor[kk] * obs.dvert * azt[2, kk] - tnsdc * azt[2, kk];

            return new double[] {tilt[0] * dConversion_Tilt, tilt[1] * dConversion_Tilt};
        }
        #endregion

        // public methods--------------------------------------------------------------------------------
        /// <summary>
        /// Calculates Component For Solid Earth Tides, depending on the Station location and Lovenumbers 
        /// (given by the constructor) and time.
        /// </summary>
        /// <param name="mjd">time in modified Julian Day</param>
        /// <param name="component">Enumerator 'component'</param>
        /// <param name="azimuth">Azimuth for Tilt and horizontal Strain (neglected for Gravity)</param>
        /// <returns>The computed tidal value in Units given by the 'BaseUnit' class.</returns>
        public double calculateComponent(double mjd, Component component, double azimuth = 0.0)
        {
            // finds geocentric position
            sph();

            //   t is ephemeris time in julian centuries from 12 hr 0 jan 1900
            double t = (mjd - 15019.5d) / 36525.0d;

            return getSolidEarthTides(t, component, azimuth);
        }

        /// <summary>
        /// Calculates Component For Solid Earth Tides, depending on the Station location and Lovenumbers 
        /// (given by the constructor) and time.
        /// </summary>
        /// <param name="mjd">time as DateTime object.</param>
        /// <param name="component">Enumerator 'component'</param>
        /// <param name="azimuth">Azimuth for Tilt and horizontal Strain (neglected for Gravity)</param>
        /// <returns>The computed tidal value in Units given by the 'BaseUnit' class.</returns>
        public double calculateComponent(DateTime dt, Component component, double azimuth = 0.0)
        {
            // finds geocentric position
            sph();

            //   t is ephemeris time in julian centuries from 12 hr 0 jan 1900
            double zhr = (double)dt.Hour + (double)dt.Minute / 60d + (double)dt.Second / 3600d;
            int iyr = dt.Year - 1900;
            double tt = zhr +
                        24d * (double)(dt.DayOfYear - 1) +
                        8760d * iyr +
                        24d * Math.Truncate((iyr - 1d) / 4d);
            tdiff.etmut = (double)(41.184f + iyr - 70f);
            double t = (tt + 12d + (tdiff.etmut / 3600d)) / 876600d;

            return getSolidEarthTides(t, component, azimuth);
        }


        public struct Tilt
        {
            public double NS;
            public double EW;
        }

        /// <summary>
        ///  Calculates Tilt Components For Solid Earth Tides, depending on the Station location and Lovenumbers.
        /// </summary>
        /// <param name="mjd">Time in modified Julian Day</param>
        /// <returns>double[2] (0) NS , (1) EW in [rad] </returns>
        public Tilt SolidEarthTidesNS_EW_Tilt(double mjd)
        {
            sph();

            //   t is ephemeris time in julian centuries from 12 hr 0 jan 1900
            double t = (mjd - 15019.5d) / 36525.0d;

            //   computes positions of the sun and moon at time t:
            ephem(t);

            // --- Tilt in directions:
            double[,] azt = new double[3, 2];
            azt[0, 0] = 0d; // N
            azt[0, 1] = 90d;// E
            // --- Sine and cosine in rad:
            azt[1, 0] = Math.Cos(azt[0, 0] * Constants.degree2radian);
            azt[1, 1] = Math.Cos(azt[0, 1] * Constants.degree2radian);
            azt[2, 0] = Math.Sin(azt[0, 0] * Constants.degree2radian);
            azt[2, 1] = Math.Sin(azt[0, 1] * Constants.degree2radian);

            //  computes the earth tides on an elastic earth, given the solar and lunar positions
            double[] tilt = elastdTilt(azt);

            return new Tilt { NS = tilt[0], EW = tilt[1] };
        }
    }
}