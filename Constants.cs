﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SPOTL
{
    /// <summary>
    /// Class conaining all needed constant values.
    /// </summary>
    public static class Constants
    {
        /// <summary>
        /// Numberformat for english notification (. as decimal separator).
        /// </summary>
        public static System.Globalization.NumberFormatInfo NumberFormatEN = new System.Globalization.CultureInfo("en-US", false).NumberFormat;

        /// <summary>
        /// Numberformat for german notification (, as decimal separator).
        /// </summary>
        public static System.Globalization.NumberFormatInfo NumberFormatDE = new System.Globalization.CultureInfo("de-DE", false).NumberFormat;


        /// <summary>
        /// Convserion of degree to radians. 
        /// </summary>
        public const double degree2radian = Math.PI / 180d;// 0.0174532935d; 
                                 // public const double dtr = 0.01745329252;
                                               // double dr = 0.01745329252;

        public const double pi = Math.PI; // 3.1415926535898;

        /// <summary>
        /// Conversion of radians to degree.
        /// </summary>
        public const double radian2degree = 180d / Math.PI;// 57.29577951d; 

        public static double h1Default = 0.61140d;
        public static double h2Default = 0.28910d;
        public static double h3Default = 0.17500d;

        public static double k1Default = 0.30400d;
        public static double k2Default = 0.09421d;
        public static double k3Default = 0.04300d;

        public static double l1Default = 0.08320d;
        public static double l2Default = 0.01450d;
        public static double l3Default = 0.01030d;

        /// <summary>
        /// Gravitational constant.
        /// </summary>
        public const double gn = 6.67e-11;

        /// <summary>
        /// Mean Earth raduis in meter.
        /// </summary>
        public const double a = 6.371e6;

        /// <summary>
        /// Minumum number of azimuth steps.
        /// </summary>
        public const int minaz = 150;
    }

    /// <summary>
    /// Contains all information about the aharmonics.
    /// </summary>
    public static class Harmonics
    {
        /// <summary>
        /// Number of harmonics.
        /// </summary>
        public const int nt = 342;

        /// <summary>
        /// AMplitude of harmonics.
        /// </summary>
        public static double[] tamp = new double[nt] 
            {
                 0.632208, 0.294107, 0.121046, 0.079915, 0.023818,-0.023589, 0.022994,  0.019333, -0.017871, 0.017192, 0.016018, 0.004671, -0.004662, -0.004519,   0.00447,
                 0.004467, 0.002589,-0.002455,-0.002172, 0.001972, 0.001947, 0.001914, -0.001898,  0.001802, 0.001304,  0.00117,  0.00113,  0.001061, -0.001022, -0.001017,
                 0.001014,  9.01e-4, -8.57e-4,   8.55e-4,  8.55e-4,  7.72e-4,  7.41e-4,   7.41e-4,  -7.21e-4,  6.98e-4,  6.58e-4,  6.54e-4,  -6.53e-4,   6.33e-4,   6.26e-4,
                 -5.98e-4,   5.9e-4,  5.44e-4,   4.79e-4, -4.64e-4,  4.13e-4,  -3.9e-4,   3.73e-4,   3.66e-4,  3.66e-4,  -3.6e-4, -3.55e-4,   3.54e-4,   3.29e-4,   3.28e-4,
                  3.19e-4,  3.02e-4,  2.79e-4,  -2.74e-4, -2.72e-4,  2.48e-4, -2.25e-4,   2.24e-4,  -2.23e-4, -2.16e-4,  2.11e-4,  2.09e-4,   1.94e-4,   1.85e-4,  -1.74e-4,
                 -1.71e-4,  1.59e-4,  1.31e-4,   1.27e-4,   1.2e-4,  1.18e-4,  1.17e-4,   1.08e-4,   1.07e-4,  1.05e-4, -1.02e-4,  1.02e-4,    9.9e-5,   -9.6e-5,    9.5e-5,
                  -8.9e-5,  -8.5e-5,  -8.4e-5,   -8.1e-5,  -7.7e-5,  -7.2e-5,  -6.7e-5,    6.6e-5,    6.4e-5,   6.3e-5,   6.3e-5,   6.3e-5,    6.2e-5,    6.2e-5,     -6e-5,
                   5.6e-5,   5.3e-5,   5.1e-5,      5e-5, 0.368645,-0.262232,-0.121995, -0.050208,  0.050031, -0.04947,  0.02062, 0.020613,  0.011279,  -0.00953, -0.009469, 
                -0.008012, 0.007414,  -0.0073,  0.007227,-0.007131,-0.006644, 0.005249,  0.004137,  0.004087, 0.003944, 0.003943,  0.00342,  0.003418,  0.002885,  0.002884,
                  0.00216,-0.001936, 0.001934, -0.001798,  0.00169, 0.001689, 0.001516,  0.001514, -0.001511, 0.001383, 0.001372, 0.001371, -0.001253, -0.001075,   0.00102, 
                  9.01e-4,  8.65e-4, -7.94e-4,   7.88e-4,  7.82e-4, -7.47e-4, -7.45e-4,    6.7e-4,  -6.03e-4, -5.97e-4,  5.42e-4,  5.42e-4,  -5.41e-4,  -4.69e-4,   -4.4e-4,
                  4.38e-4,  4.22e-4,   4.1e-4,  -3.74e-4, -3.65e-4,  3.45e-4,  3.35e-4,  -3.21e-4,  -3.19e-4,  3.07e-4,  2.91e-4,   2.9e-4,  -2.89e-4,   2.86e-4,   2.75e-4,
                  2.71e-4,  2.63e-4, -2.45e-4,   2.25e-4,  2.25e-4,  2.21e-4, -2.02e-4,     -2e-4,  -1.99e-4,  1.92e-4,  1.83e-4,  1.83e-4,   1.83e-4,   -1.7e-4,   1.69e-4, 
                  1.68e-4,  1.62e-4,  1.49e-4,  -1.47e-4, -1.41e-4,  1.38e-4,  1.36e-4,   1.36e-4,   1.27e-4,  1.27e-4, -1.26e-4, -1.21e-4,  -1.21e-4,   1.17e-4,  -1.16e-4,
                 -1.14e-4, -1.14e-4, -1.14e-4,   1.14e-4,  1.13e-4,  1.09e-4,  1.08e-4,   1.06e-4,  -1.06e-4, -1.06e-4,  1.05e-4,  1.04e-4,  -1.03e-4,     -1e-4,     -1e-4,
                    -1e-4,   9.9e-5,  -9.8e-5,    9.3e-5,   9.3e-5,     9e-5,  -8.8e-5,    8.3e-5,   -8.3e-5,  -8.2e-5,  -8.1e-5,  -7.9e-5,   -7.7e-5,   -7.5e-5,   -7.5e-5,
                  -7.5e-5,   7.1e-5,   7.1e-5,   -7.1e-5,   6.8e-5,   6.8e-5,   6.5e-5,    6.5e-5,    6.4e-5,   6.4e-5,   6.4e-5,  -6.4e-5,     -6e-5,    5.6e-5,    5.6e-5,
                   5.3e-5,   5.3e-5,   5.3e-5,   -5.3e-5,   5.3e-5,   5.3e-5,   5.2e-5,      5e-5, -0.066607,-0.035184,-0.030988, 0.027929, -0.027616, -0.012753, -0.006728,
                -0.005837,-0.005286,-0.004921, -0.002884,-0.002583,-0.002422,  0.00231,  0.002283, -0.002037, 0.001883,-0.001811,-0.001687, -0.001004,  -9.25e-4,  -8.44e-4,
                  7.66e-4,  7.66e-4,    -7e-4,  -4.95e-4, -4.92e-4,  4.91e-4,  4.83e-4,   4.37e-4,  -4.16e-4, -3.84e-4,  3.74e-4, -3.12e-4,  -2.88e-4,  -2.73e-4,   2.59e-4,
                  2.45e-4, -2.32e-4,  2.29e-4,  -2.16e-4,  2.06e-4, -2.04e-4, -2.02e-4,      2e-4,   1.95e-4,  -1.9e-4,  1.87e-4,   1.8e-4,  -1.79e-4,    1.7e-4,   1.53e-4,
                 -1.37e-4, -1.19e-4, -1.19e-4,  -1.12e-4,  -1.1e-4,  -1.1e-4,  1.07e-4,   -9.5e-5,   -9.5e-5,  -9.1e-5,    -9e-5,  -8.1e-5,   -7.9e-5,   -7.9e-5,    7.7e-5,
                  -7.3e-5,   6.9e-5,  -6.7e-5,   -6.6e-5,   6.5e-5,   6.4e-5,  -6.2e-5,      6e-5,    5.9e-5,  -5.6e-5,   5.5e-5,  -5.1e-5
            };

        /// <summary>
        /// Get the doodsen number of the harmics according to the index.
        /// </summary>
        /// <param name="index">Index of the harmonics.</param>
        /// <returns>Array of Doodsen number.</returns>
        public static int[] get_idd(int index)
        {
            return new int[6]
                    {
                        idd[index, 0],
                        idd[index, 1],
                        idd[index, 2], 
                        idd[index, 3], 
                        idd[index, 4],
                        idd[index, 5] 
                    };
        }

        /// <summary>
        /// Doodsen numbers for harmonics.
        /// </summary>
        public static int[,] idd = new int[nt, 6] 
            {
                 { 2, 0, 0, 0, 0, 0},
                 { 2, 2,-2, 0, 0, 0}, 
                 { 2,-1, 0, 1, 0, 0}, 
                 { 2, 2, 0, 0, 0, 0}, 
                 { 2, 2, 0, 0, 1, 0}, 
                 { 2, 0, 0, 0,-1, 0},
                 { 2,-1, 2,-1, 0, 0},
                 { 2,-2, 2, 0, 0, 0}, 
                 { 2, 1, 0,-1, 0, 0}, 
                 { 2, 2,-3, 0, 0, 1},
                 { 2,-2, 0, 2, 0, 0},
                 { 2,-3, 2, 1, 0, 0},
                 { 2, 1,-2, 1, 0, 0},
                 { 2,-1, 0, 1,-1, 0}, 
                 { 2, 3, 0,-1, 0, 0}, 
                 { 2, 1, 0, 1, 0, 0}, 
                 { 2, 2, 0, 0, 2, 0},
                 { 2, 2,-1, 0, 0,-1},
                 { 2, 0,-1, 0, 0, 1}, 
                 { 2, 1, 0, 1, 1, 0},
                 { 2, 3, 0,-1, 1, 0},
                 { 2, 0, 1, 0, 0,-1},
                 { 2, 0,-2, 2, 0, 0},
                 { 2,-3, 0, 3, 0, 0},
                 { 2,-2, 3, 0, 0,-1},
                 { 2, 4, 0, 0, 0, 0},
                 { 2,-1, 1, 1, 0,-1}, 
                 { 2,-1, 3,-1, 0,-1}, 
                 { 2, 2, 0, 0,-1, 0}, 
                 { 2,-1,-1, 1, 0, 1},
                 { 2, 4, 0, 0, 1, 0},
                 { 2,-3, 4,-1, 0, 0}, 
                 { 2,-1, 2,-1,-1, 0}, 
                 { 2, 3,-2, 1, 0, 0}, 
                 { 2, 1, 2,-1, 0, 0}, 
                 { 2,-4, 2, 2, 0, 0},
                 { 2, 4,-2, 0, 0, 0}, 
                 { 2, 0, 2, 0, 0, 0},
                 { 2,-2, 2, 0,-1, 0},
                 { 2, 2,-4, 0, 0, 2},
                 { 2, 2,-2, 0,-1, 0},
                 { 2, 1, 0,-1,-1, 0},
                 { 2,-1, 1, 0, 0, 0},
                 { 2, 2,-1, 0, 0, 1},
                 { 2, 2, 1, 0, 0,-1},
                 { 2,-2, 0, 2,-1, 0},
                 { 2,-2, 4,-2, 0, 0},
                 { 2, 2, 2, 0, 0, 0},
                 { 2,-4, 4, 0, 0, 0}, 
                 { 2,-1, 0,-1,-2, 0}, 
                 { 2, 1, 2,-1, 1, 0},
                 { 2,-1,-2, 3, 0, 0},
                 { 2, 3,-2, 1, 1, 0},
                 { 2, 4, 0,-2, 0, 0},
                 { 2, 0, 0, 2, 0, 0},
                 { 2, 0, 2,-2, 0, 0},
                 { 2, 0, 2, 0, 1, 0},
                 { 2,-3, 3, 1, 0,-1},
                 { 2, 0, 0, 0,-2, 0},
                 { 2, 4, 0, 0, 2, 0},
                 { 2, 4,-2, 0, 1, 0},
                 { 2, 0, 0, 0, 0, 2},
                 { 2, 1, 0, 1, 2, 0}, 
                 { 2, 0,-2, 0,-2, 0},
                 { 2,-2, 1, 0, 0, 1},
                 { 2,-2, 1, 2, 0,-1},
                 { 2,-1, 1,-1, 0, 1},
                 { 2, 5, 0,-1, 0, 0},
                 { 2, 1,-3, 1, 0, 1},
                 { 2,-2,-1, 2, 0, 1},
                 { 2, 3, 0,-1, 2, 0},
                 { 2, 1,-2, 1,-1, 0},
                 { 2, 5, 0,-1, 1, 0},
                 { 2,-4, 0, 4, 0, 0},
                 { 2,-3, 2, 1,-1, 0},
                 { 2,-2, 1, 1, 0, 0},
                 { 2, 4, 0,-2, 1, 0},
                 { 2, 0, 0, 2, 1, 0},
                 { 2,-5, 4, 1, 0, 0},
                 { 2, 0, 2, 0, 2, 0},
                 { 2,-1, 2, 1, 0, 0},
                 { 2, 5,-2,-1, 0, 0},
                 { 2, 1,-1, 0, 0, 0},
                 { 2, 2,-2, 0, 0, 2},
                 { 2,-5, 2, 3, 0, 0},
                 { 2,-1,-2, 1,-2, 0},
                 { 2,-3, 5,-1, 0,-1},
                 { 2,-1, 0, 0, 0, 1},
                 { 2,-2, 0, 0,-2, 0},
                 { 2, 0,-1, 1, 0, 0},
                 { 2,-3, 1, 1, 0, 1},
                 { 2, 3, 0,-1,-1, 0},
                 { 2, 1, 0, 1,-1, 0},
                 { 2,-1, 2, 1, 1, 0},
                 { 2, 0,-3, 2, 0, 1},
                 { 2, 1,-1,-1, 0, 1},
                 { 2,-3, 0, 3,-1, 0},
                 { 2, 0,-2, 2,-1, 0},
                 { 2,-4, 3, 2, 0,-1},
                 { 2,-1, 0, 1,-2, 0},
                 { 2, 5, 0,-1, 2, 0},
                 { 2,-4, 5, 0, 0,-1},
                 { 2,-2, 4, 0, 0,-2},
                 { 2,-1, 0, 1, 0, 2},
                 { 2,-2,-2, 4, 0, 0},
                 { 2, 3,-2,-1,-1, 0},
                 { 2,-2, 5,-2, 0,-1},
                 { 2, 0,-1, 0,-1, 1},
                 { 2, 5,-2,-1, 1, 0},
                 { 1, 1, 0, 0, 0, 0},
                 { 1,-1, 0, 0, 0, 0},
                 { 1, 1,-2, 0, 0, 0},
                 { 1,-2, 0, 1, 0, 0},
                 { 1, 1, 0, 0, 1, 0},
                 { 1,-1, 0, 0,-1, 0},
                 { 1, 2, 0,-1, 0, 0},
                 { 1, 0, 0, 1, 0, 0},
                 { 1, 3, 0, 0, 0, 0}, 
                 { 1,-2, 2,-1, 0, 0}, 
                 { 1,-2, 0, 1,-1, 0},
                 { 1,-3, 2, 0, 0, 0},
                 { 1, 0, 0,-1, 0, 0}, 
                 { 1, 1, 0, 0,-1, 0},
                 { 1, 3, 0, 0, 1, 0},
                 { 1, 1,-3, 0, 0, 1},
                 { 1,-3, 0, 2, 0, 0},
                 { 1, 1, 2, 0, 0, 0},
                 { 1, 0, 0, 1, 1, 0},
                 { 1, 2, 0,-1, 1, 0},
                 { 1, 0, 2,-1, 0, 0},
                 { 1, 2,-2, 1, 0, 0},
                 { 1, 3,-2, 0, 0, 0},
                 { 1,-1, 2, 0, 0, 0},
                 { 1, 1, 1, 0, 0,-1},
                 { 1, 1,-1, 0, 0, 1}, 
                 { 1, 4, 0,-1, 0, 0}, 
                 { 1,-4, 2, 1, 0, 0}, 
                 { 1, 0,-2, 1, 0, 0},
                 { 1,-2, 2,-1,-1, 0},
                 { 1, 3, 0,-2, 0, 0},
                 { 1,-1, 0, 2, 0, 0},
                 { 1,-1, 0, 0,-2, 0},
                 { 1, 3, 0, 0, 2, 0}, 
                 { 1,-3, 2, 0,-1, 0},
                 { 1, 4, 0,-1, 1, 0},
                 { 1, 0, 0,-1,-1, 0}, 
                 { 1, 1,-2, 0,-1, 0},
                 { 1,-3, 0, 2,-1, 0},
                 { 1, 1, 0, 0, 2, 0},
                 { 1, 1,-1, 0, 0,-1},
                 { 1,-1,-1, 0, 0, 1}, 
                 { 1, 0, 2,-1, 1, 0},
                 { 1,-1, 1, 0, 0,-1},
                 { 1,-1,-2, 2, 0, 0}, 
                 { 1, 2,-2, 1, 1, 0},
                 { 1,-4, 0, 3, 0, 0},
                 { 1,-1, 2, 0, 1, 0},
                 { 1, 3,-2, 0, 1, 0},
                 { 1, 2, 0,-1,-1, 0},
                 { 1, 0, 0, 1,-1, 0},
                 { 1,-2, 2, 1, 0, 0},
                 { 1, 4,-2,-1, 0, 0},
                 { 1,-3, 3, 0, 0,-1},
                 { 1,-2, 1, 1, 0,-1},
                 { 1,-2, 3,-1, 0,-1}, 
                 { 1, 0,-2, 1,-1, 0},
                 { 1,-2,-1, 1, 0, 1},
                 { 1, 4,-2, 1, 0, 0},
                 { 1,-4, 4,-1, 0, 0},
                 { 1,-4, 2, 1,-1, 0},
                 { 1, 5,-2, 0, 0, 0},
                 { 1, 3, 0,-2, 1, 0},
                 { 1,-5, 2, 2, 0, 0},
                 { 1, 2, 0, 1, 0, 0},
                 { 1, 1, 3, 0, 0,-1},
                 { 1,-2, 0, 1,-2, 0},
                 { 1, 4, 0,-1, 2, 0},
                 { 1, 1,-4, 0, 0, 2},
                 { 1, 5, 0,-2, 0, 0},
                 { 1,-1, 0, 2, 1, 0},
                 { 1,-2, 1, 0, 0, 0},
                 { 1, 4,-2, 1, 1, 0},
                 { 1,-3, 4,-2, 0, 0}, 
                 { 1,-1, 3, 0, 0,-1},
                 { 1, 3,-3, 0, 0, 1},
                 { 1, 5,-2, 0, 1, 0},
                 { 1, 1, 2, 0, 1, 0},
                 { 1, 2, 0, 1, 1, 0},
                 { 1,-5, 4, 0, 0, 0}, 
                 { 1,-2, 0,-1,-2, 0},
                 { 1, 5, 0,-2, 1, 0},
                 { 1, 1, 2,-2, 0, 0},
                 { 1, 1,-2, 2, 0, 0},
                 { 1,-2, 2, 1, 1, 0},
                 { 1, 0, 3,-1, 0,-1},
                 { 1, 2,-3, 1, 0, 1},
                 { 1,-2,-2, 3, 0, 0},
                 { 1,-1, 2,-2, 0, 0},
                 { 1,-4, 3, 1, 0,-1}, 
                 { 1,-4, 0, 3,-1, 0}, 
                 { 1,-1,-2, 2,-1, 0},
                 { 1,-2, 0, 3, 0, 0},
                 { 1, 4, 0,-3, 0, 0}, 
                 { 1, 0, 1, 1, 0,-1},
                 { 1, 2,-1,-1, 0, 1},
                 { 1, 2,-2, 1,-1, 0},
                 { 1, 0, 0,-1,-2, 0},
                 { 1, 2, 0, 1, 2, 0},
                 { 1, 2,-2,-1,-1, 0},
                 { 1, 0, 0, 1, 2, 0},
                 { 1, 0, 1, 0, 0, 0},
                 { 1, 2,-1, 0, 0, 0}, 
                 { 1, 0, 2,-1,-1, 0},
                 { 1,-1,-2, 0,-2, 0},
                 { 1,-3, 1, 0, 0, 1},
                 { 1, 3,-2, 0,-1, 0},
                 { 1,-1,-1, 0,-1, 1},
                 { 1, 4,-2,-1, 1, 0},
                 { 1, 2, 1,-1, 0,-1},
                 { 1, 0,-1, 1, 0, 1},
                 { 1,-2, 4,-1, 0, 0},
                 { 1, 4,-4, 1, 0, 0},
                 { 1,-3, 1, 2, 0,-1},
                 { 1,-3, 3, 0,-1,-1},
                 { 1, 1, 2, 0, 2, 0},
                 { 1, 1,-2, 0,-2, 0},
                 { 1, 3, 0, 0, 3, 0},
                 { 1,-1, 2, 0,-1, 0},
                 { 1,-2, 1,-1, 0, 1},
                 { 1, 0,-3, 1, 0, 1},
                 { 1,-3,-1, 2, 0, 1},
                 { 1, 2, 0,-1, 2, 0},
                 { 1, 6,-2,-1, 0, 0},
                 { 1, 2, 2,-1, 0, 0},
                 { 1,-1, 1, 0,-1,-1},
                 { 1,-2, 3,-1,-1,-1},
                 { 1,-1, 0, 0, 0, 2},
                 { 1,-5, 0, 4, 0, 0},
                 { 1, 1, 0, 0, 0,-2},
                 { 1,-2, 1, 1,-1,-1},
                 { 1, 1,-1, 0, 1, 1},
                 { 1, 1, 2, 0, 0,-2},
                 { 1,-3, 1, 1, 0, 0},
                 { 1,-4, 4,-1,-1, 0},
                 { 1, 1, 0,-2,-1, 0},
                 { 1,-2,-1, 1,-1, 1},
                 { 1,-3, 2, 2, 0, 0},
                 { 1, 5,-2,-2, 0, 0},
                 { 1, 3,-4, 2, 0, 0},
                 { 1, 1,-2, 0, 0, 2},
                 { 1,-1, 4,-2, 0, 0},
                 { 1, 2, 2,-1, 1, 0},
                 { 1,-5, 2, 2,-1, 0},
                 { 1, 1,-3, 0,-1, 1},
                 { 1, 1, 1, 0, 1,-1},
                 { 1, 6,-2,-1, 1, 0}, 
                 { 1,-2, 2,-1,-2, 0},
                 { 1, 4,-2, 1, 2, 0},
                 { 1,-6, 4, 1, 0, 0},
                 { 1, 5,-4, 0, 0, 0}, 
                 { 1,-3, 4, 0, 0, 0},
                 { 1, 1, 2,-2, 1, 0}, 
                 { 1,-2, 1, 0,-1, 0},
                 { 0, 2, 0, 0, 0, 0},
                 { 0, 1, 0,-1, 0, 0},
                 { 0, 0, 2, 0, 0, 0},
                 { 0, 0, 0, 0, 1, 0},
                 { 0, 2, 0, 0, 1, 0},
                 { 0, 3, 0,-1, 0, 0},
                 { 0, 1,-2, 1, 0, 0},
                 { 0, 2,-2, 0, 0, 0},
                 { 0, 3, 0,-1, 1, 0},
                 { 0, 0, 1, 0, 0,-1},
                 { 0, 2, 0,-2, 0, 0},
                 { 0, 2, 0, 0, 2, 0},
                 { 0, 3,-2, 1, 0, 0},
                 { 0, 1, 0,-1,-1, 0},
                 { 0, 1, 0,-1, 1, 0},
                 { 0, 4,-2, 0, 0, 0},
                 { 0, 1, 0, 1, 0, 0}, 
                 { 0, 0, 3, 0, 0,-1},
                 { 0, 4, 0,-2, 0, 0},
                 { 0, 3,-2, 1, 1, 0},
                 { 0, 3,-2,-1, 0, 0},
                 { 0, 4,-2, 0, 1, 0},
                 { 0, 0, 2, 0, 1, 0},
                 { 0, 1, 0, 1, 1, 0},
                 { 0, 4, 0,-2, 1, 0},
                 { 0, 3, 0,-1, 2, 0},
                 { 0, 5,-2,-1, 0, 0},
                 { 0, 1, 2,-1, 0, 0},
                 { 0, 1,-2, 1,-1, 0},
                 { 0, 1,-2, 1, 1, 0},
                 { 0, 2,-2, 0,-1, 0},
                 { 0, 2,-3, 0, 0, 1},
                 { 0, 2,-2, 0, 1, 0},
                 { 0, 0, 2,-2, 0, 0},
                 { 0, 1,-3, 1, 0, 1},
                 { 0, 0, 0, 0, 2, 0},
                 { 0, 0, 1, 0, 0, 1},
                 { 0, 1, 2,-1, 1, 0},
                 { 0, 3, 0,-3, 0, 0},
                 { 0, 2, 1, 0, 0,-1},
                 { 0, 1,-1,-1, 0, 1}, 
                 { 0, 1, 0, 1, 2, 0},
                 { 0, 5,-2,-1, 1, 0},
                 { 0, 2,-1, 0, 0, 1},
                 { 0, 2, 2,-2, 0, 0},
                 { 0, 1,-1, 0, 0, 0}, 
                 { 0, 5, 0,-3, 0, 0},
                 { 0, 2, 0,-2, 1, 0},
                 { 0, 1, 1,-1, 0,-1},
                 { 0, 3,-4, 1, 0, 0},
                 { 0, 0, 2, 0, 2, 0}, 
                 { 0, 2, 0,-2,-1, 0},
                 { 0, 4,-3, 0, 0, 1}, 
                 { 0, 3,-1,-1, 0, 1},
                 { 0, 0, 2, 0, 0,-2},
                 { 0, 3,-3, 1, 0, 1}, 
                 { 0, 2,-4, 2, 0, 0},
                 { 0, 4,-2,-2, 0, 0},
                 { 0, 3, 1,-1, 0,-1},
                 { 0, 5,-4, 1, 0, 0},
                 { 0, 3,-2,-1,-1, 0},
                 { 0, 3,-2, 1, 2, 0},
                 { 0, 4,-4, 0, 0, 0},
                 { 0, 6,-2,-2, 0, 0},
                 { 0, 5, 0,-3, 1, 0},
                 { 0, 4,-2, 0, 2, 0},
                 { 0, 2, 2,-2, 1, 0},
                 { 0, 0, 4, 0, 0,-2},
                 { 0, 3,-1, 0, 0, 0},
                 { 0, 3,-3,-1, 0, 1},
                 { 0, 4, 0,-2, 2, 0},
                 { 0, 1,-2,-1,-1, 0},
                 { 0, 2,-1, 0, 0,-1},
                 { 0, 4,-4, 2, 0, 0},
                 { 0, 2, 1, 0, 1,-1},
                 { 0, 3,-2,-1, 1, 0},
                 { 0, 4,-3, 0, 1, 1},
                 { 0, 2, 0, 0, 3, 0},
                 { 0, 6,-4, 0, 0, 0}
            };
    }
}
