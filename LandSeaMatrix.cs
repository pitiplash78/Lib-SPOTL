using System;
//using System.Collections.Generic;
using System.IO;
//using System.IO.Compression;
using System.Text;

namespace SPOTL
{
    internal class LandseaMatrix
    {
        private int horCells = 0;
        private int verCells = 0;
        private int mixedCellsLength = 0;

        private int[,] matrix = null;
        const int water = -2; 
        const int land = -1;
        private BitMatrix[] submatrix = null;
        private LandseaMatrix()
        { }

        // TODO Sinnvoll auslagern -> Utilities
        public sealed class BitMatrix
        {
            public int width;
            public int height;
            public int rowSize;
            public int[] bits;

            // A helper to construct a square matrix.
            public BitMatrix(int dimension)
                : this(dimension, dimension)
            {
            }

            public BitMatrix(int width, int height)
            {
                if (width < 1 || height < 1)
                {
                    throw new ArgumentException("Both dimensions must be greater than 0");
                }
                this.width = width;
                this.height = height;
                int rowSize = width >> 5;
                if ((width & 0x1f) != 0)
                {
                    rowSize++;
                }
                this.rowSize = rowSize;
                bits = new int[rowSize * height];
            }

            /**
             * <p>Gets the requested bit, where true means black.</p>
             *
             * @param x The horizontal component (i.e. which column)
             * @param y The vertical component (i.e. which row)
             * @return value of given bit in matrix
             */
            public bool Get(int x, int y)
            {
                int offset = y * rowSize + (x >> 5);
                return ((bits[offset] >> (x & 0x1f)) & 1) != 0;
            }

            /**
             * <p>Sets the given bit to true.</p>
             *
             * @param x The horizontal component (i.e. which column)
             * @param y The vertical component (i.e. which row)
             */
            public void Set(int x, int y)
            {
                int offset = y * rowSize + (x >> 5);
                bits[offset] |= 1 << (x & 0x1f);
            }

            /**
             * <p>Flips the given bit.</p>
             *
             * @param x The horizontal component (i.e. which column)
             * @param y The vertical component (i.e. which row)
             */
            public void Flip(int x, int y)
            {
                int offset = y * rowSize + (x >> 5);
                bits[offset] ^= 1 << (x & 0x1f);
            }

            /**
             * Clears all bits (sets to false).
             */
            public void Clear()
            {
                int max = bits.Length;
                for (int i = 0; i < max; i++)
                {
                    bits[i] = 0;
                }
            }

            /**
             * <p>Sets a square region of the bit matrix to true.</p>
             *
             * @param left The horizontal position to begin at (inclusive)
             * @param top The vertical position to begin at (inclusive)
             * @param width The width of the region
             * @param height The height of the region
             */
            public void SetRegion(int left, int top, int width, int height)
            {
                if (top < 0 || left < 0)
                {
                    throw new ArgumentException("Left and top must be nonnegative");
                }
                if (height < 1 || width < 1)
                {
                    throw new ArgumentException("Height and width must be at least 1");
                }
                int right = left + width;
                int bottom = top + height;
                if (bottom > this.height || right > this.width)
                {
                    throw new ArgumentException("The region must fit inside the matrix");
                }
                for (int y = top; y < bottom; y++)
                {
                    int offset = y * rowSize;
                    for (int x = left; x < right; x++)
                    {
                        bits[offset + (x >> 5)] |= 1 << (x & 0x1f);
                    }
                }
            }


            /**
             * @return The width of the matrix
             */
            public int GetWidth()
            {
                return width;
            }

            /**
             * @return The height of the matrix
             */
            public int GetHeight()
            {
                return height;
            }

            public override String ToString()
            {
                StringBuilder result = new StringBuilder(height * (width + 1));
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        result.Append(Get(x, y) ? "X " : "  ");
                    }
                    result.Append('\n');
                }
                return result.ToString();
            }

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static LandseaMatrix readLandSea(string path)
        {
            FileStream fs = new FileStream(path, FileMode.Open);
            StreamReader re = new StreamReader(fs);

            LandseaMatrix lm = new LandseaMatrix();

            re.ReadLine();
            lm.horCells = int.Parse(re.ReadLine().Substring(16, 10));
            lm.verCells = int.Parse(re.ReadLine().Substring(16, 10));
            lm.water = int.Parse(re.ReadLine().Substring(16, 10));
            lm.land = int.Parse(re.ReadLine().Substring(16, 10));
            re.ReadLine();
            lm.mixedCellsLength = int.Parse(re.ReadLine().Substring(16, 10));
            re.ReadLine();
            re.ReadLine();

            lm.matrix = new int[lm.horCells, lm.verCells];
            string str = null;

            for (int j = 0; j < lm.verCells; j++)
            {
                str = re.ReadLine();
                for (int i = 0; i < lm.horCells; i++)
                {
                    string tmp = str.Substring(i * 6, 6);
                    lm.matrix[i, j] = int.Parse(tmp);
                }
            }

            lm.submatrix = new LandseaMatrix.BitMatrix[lm.mixedCellsLength];

            for (int k = 0; k < lm.mixedCellsLength; k++)
            {
                lm.submatrix[k] = new LandseaMatrix.BitMatrix(32);
                str = re.ReadLine();
                for (int j = 0; j < 32; j++)
                {
                    str = re.ReadLine();
                    for (int i = 0; i < 32; i++)
                        if (str.Substring(i, 1) == "1")
                            lm.submatrix[k].Set(i, j);
                }
            }
            re.Close();
            return lm;
        }

        // überprüfen anhand von coordinate ob land/wassser
        /// <summary>
        /// Determine the land / sea behavior vor given coordinate.
        /// </summary>
        /// <param name="latitude"></param>
        /// <param name="longitude"></param>
        /// <returns> land (true) / water (false) </returns>
        public bool getLandSea(double latitude, double longitude)
        {

            // special branches for the polar regions, to avoid possible indexing trouble
            if (latitude >= 85d) return false;
            if (latitude < -87d) return true;

            //länge := i
            double rlong = longitude;
            if (longitude < 0d) rlong += 360;
            if (longitude > 360d) rlong -= 360;
            int ilong = (int)(2 * rlong);
            // breite := j
            int ilat = (int)(2 * (double)(90d - latitude));

            bool tmp = false;
            int lnd = matrix[ilong, ilat];
            if (lnd == land)
                tmp = true;
            else if (lnd == water)
                tmp = false;
            else
            {
                int i = (int)(64d * (rlong - (0.5d * ((double)ilong))));
                int j = 31 - (int)(64d * (latitude - (90d - 0.5d * (double)(ilat + 1))));
                tmp = submatrix[lnd - 1].Get(i, j);
            }
            return tmp;
        }
    }
}
