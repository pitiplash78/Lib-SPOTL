using System;
using System.Collections;
using System.IO;
using System.IO.Compression;
using System.Xml.Serialization;

namespace SPOTL
{
    /// <summary>
    /// Internal class of directory and file paths used in library
    /// </summary>
    internal class FilePath
    {
        /// <summary>
        /// Directory path for the Earth models.
        /// </summary>
        internal string EarthModels = null;
        internal const string _EarthModels = "EarthModels"; 

        /// <summary>
        /// Directory path for the tidal models.
        /// </summary>
        internal string OceanModels = null;
        internal const string _OceanModels = "OceanModels";

        /// <summary>
        /// Directory path for the local tidal models.
        /// </summary>
        internal string LocalModels = null;
        internal const string _LocalModels = "LocalModels";

        /// <summary>
        /// File path for the Waterdensity.
        /// </summary>
        internal string WaterDensity = null;
        internal const string _WaterDensityDirectory = "WaterDensity";
        internal const string _WaterDensityFileName = "seadens.gz";

        /// <summary>
        /// File path for the land sea matrix.
        /// </summary>
        internal string LndSeaMatrix = null;
        internal const string _LndSeaMatrixDirectory = "LandSeaMatrix";
        internal const string _LndSeaMatrixFileName = "landsea.matrix.gz";
    }

    /// <summary>
    /// Class for the needed parameter for ocean loading calculation.
    /// </summary>
    [XmlRoot("OceanLoadingProperties")]
    public class OceanLoadingProperties
    {
        /// <summary>
        /// Location of the site where the ocean loading is to be calculated.
        /// </summary>
        public class Location
        {
            /// <summary>
            /// Name of the location.
            /// </summary>
            internal string Name = null;

            /// <summary>
            /// Longitude in degree of the location.
            /// </summary>
            internal double Longitude = double.NaN;
            
            /// <summary>
            /// Latitude in degree of the location.
            /// </summary>
            internal double Latitude = double.NaN;

            /// <summary>
            /// Heigth in meter of the location.
            /// </summary>
            internal double Height = double.NaN;

            /// <summary>
            /// Constructor: Location of the site where the ocean loading is to be calculated.
            /// </summary>
            /// <param name="Name">Name of the location.</param>
            /// <param name="Longitude">Longitude in degree of the location.</param>
            /// <param name="Latitude">Latitude in degree of the location.</param>
            /// <param name="Height">Heigth in meter of the location.</param>
            public Location(string Name, double Longitude, double Latitude, double Height)
            {
                this.Name = Name;
                this.Longitude = Longitude;
                this.Latitude = Latitude;
                this.Height = Height;
            }
        }

        /// <summary>
        /// Constructor for the needed parameter for ocean loading calculation.
        /// </summary>
        public OceanLoadingProperties()
        {
            OceanModelData = new ArrayList();
            GreensFunctionData = new ArrayList();
            LocalOceanModelData = new ArrayList();
        }

        #region Greens Function
        [XmlElement("GreenFunctionData")]
        private ArrayList GreensFunctionData;

        /// <summary>
        /// Array of usable Green's function
        /// </summary>
        [XmlElement("GreenFunctionItems")]
        public GreenFunctionEntry[] GreenFunctionItems
        {
            get
            {
                GreenFunctionEntry[] items = new GreenFunctionEntry[GreensFunctionData.Count];
                GreensFunctionData.CopyTo(items);
                return items;
            }
            set
            {
                if (value == null) return;
                GreenFunctionEntry[] items = (GreenFunctionEntry[])value;
                GreensFunctionData.Clear();
                foreach (GreenFunctionEntry item in items)
                    GreensFunctionData.Add(item);
            }
        }

        public int AddGreenFunctionEntry(GreenFunctionEntry item)
        {
            return GreensFunctionData.Add(item);
        }

        /// <summary>
        /// Class containing all needed meta information of a Green's function
        /// </summary>
        public class GreenFunctionEntry
        {
            [XmlAttribute("Description")]
            public string Description = null;

            [XmlAttribute("FileName")]
            public string FileName = null;

            [XmlAttribute("EarthModel")]
            public string EarthModel = null;

            [XmlAttribute("Source")]
            public string Source = null;

            [XmlAttribute("Pattern")]
            public string Pattern = null;

            [XmlAttribute("Frame")]
            public string Frame = null;

            /// <summary>
            /// Constructor for the meta information of a Green's function
            /// </summary>
            public GreenFunctionEntry()
            { }

            /// <summary>
            /// Constructor for the metha information of a Green's function
            /// </summary>
            /// <param name="Entry">Green's function meta information.</param>
            public GreenFunctionEntry(GreenFunctionEntry Entry)
            {
                this.Description = Entry.Description;
                this.FileName = Entry.FileName;
                this.EarthModel = Entry.EarthModel;
                this.Source = Entry.Source;
                this.Pattern = Entry.Pattern;
                this.Frame = Entry.Frame;
            }

            /// <summary>
            /// Constructor for the metha information of a Green's function.
            /// </summary>
            /// <param name="Description">Description of the Green's function.</param>
            /// <param name="FileName">File name of the Green's function.</param>
            /// <param name="EarthModel">Earth model used for the calculation of the Green's function.</param>
            /// <param name="Source">Source (creator/developer) of the Green's function</param>
            /// <param name="Pattern">Pattern (fine/coarse) of the Green's function.</param>
            /// <param name="Frame">Frame of the Green's function</param>
            public GreenFunctionEntry(string Description,
                                     string FileName,
                                     string EarthModel,
                                     string Source,
                                     string Pattern,
                                     string Frame)
            {
                this.Description = Description;
                this.FileName = FileName;
                this.EarthModel = EarthModel;
                this.Source = Source;
                this.Pattern = Pattern;
                this.Frame = Frame;
            }

            /// <summary>
            /// Compares two Green's function.
            /// </summary>
            /// <param name="Entry">Green's function to be compared.</param>
            /// <returns>Equal: TRUE; not equal: FALSE.</returns>
            public bool Equal(GreenFunctionEntry Entry)
            {
                if (this != null && Entry != null)
                {
                    if (this.Description == Entry.Description &&
                        this.FileName == Entry.FileName &&
                        this.EarthModel == Entry.EarthModel &&
                        this.Source == Entry.Source &&
                        this.Pattern == Entry.Pattern &&
                        this.Frame == Entry.Frame)
                        return true;
                    else
                        return false;
                }
                else
                    return false;
            }
        }
        #endregion

        /// <summary>
        /// Names of waves of available global ocean models.
        /// </summary>
        [XmlArray("WaveSequenceOceanModel")]
        public string[] WaveSequenceOceanModel = new string[] { "M2", "N2", "S2", "K2", "O1", "K1", "P1", "Q1", "S1", "M4", "Mf", "Mm", "Ssa" };

        #region Ocean Model
        private ArrayList OceanModelData;

        [XmlElement("OceanModelItems")]
        public OceanModelEntry[] OceanModelItems
        {
            get
            {
                OceanModelEntry[] items = new OceanModelEntry[OceanModelData.Count];
                OceanModelData.CopyTo(items);
                return items;
            }
            set
            {
                if (value == null) return;
                OceanModelEntry[] items = (OceanModelEntry[])value;
                OceanModelData.Clear();
                foreach (OceanModelEntry item in items)
                    OceanModelData.Add(item);
            }
        }

        public int AddOceanModelEntryEntry(OceanModelEntry item)
        {
            int rValue = OceanModelData.Add(item);
            OceanModelData.Sort();
            return rValue;
        }

        public class OceanModelEntry : IComparable
        {
            [XmlAttribute("Description")]
            public string Description = null;

            [XmlAttribute("ModelName")]
            public string ModelName = null;

            [XmlAttribute("FilePath")]
            public string FilePath = null;

            [XmlAttribute("Year")]
            public int Year;

            [XmlArrayItem("Waves")]
            public string[] DarwinSymbol = null;

            public OceanModelEntry()
            { }

            public OceanModelEntry(string ModelName, string FilePath, int Year, string[] DarwinSymbol)
            {
                this.ModelName = ModelName;
                this.FilePath = FilePath;
                this.Year = Year;
                this.DarwinSymbol = DarwinSymbol;
            }

            public int CompareTo(Object obj)
            {
                if (obj is OceanModelEntry)
                    return (obj as OceanModelEntry).Year.CompareTo(this.Year);
                return 0;
            }

        }
        #endregion


        /// <summary>
        /// Names of waves of available local ocean models.
        /// </summary>
        [XmlArray("WaveSequenceLocalModel")]
        public string[] WaveSequenceLocalModel = new string[] { "M2", "N2", "S2", "K2", "O1", "K1", "P1", "Q1", "M4", "Mf", "Mm" };

        #region Local Ocean Model
        private ArrayList LocalOceanModelData;

        [XmlElement("LocalOceanModelItem")]
        public LocalOceanModelEntry[] LocalOceanModelItem
        {
            get
            {
                LocalOceanModelEntry[] items = new LocalOceanModelEntry[LocalOceanModelData.Count];
                LocalOceanModelData.CopyTo(items);
                return items;
            }
            set
            {
                if (value == null) return;
                LocalOceanModelEntry[] items = (LocalOceanModelEntry[])value;
                LocalOceanModelData.Clear();
                foreach (LocalOceanModelEntry item in items)
                    LocalOceanModelData.Add(item);
            }
        }

        public int AddLocalTideModelEntry(LocalOceanModelEntry item)
        {
             int rValue = LocalOceanModelData.Add(item);
             LocalOceanModelData.Sort();
             return rValue;
        }

        public class LocalOceanModelEntry : IComparable
        {
            [XmlAttribute("Number")]
            public int Number;

            [XmlAttribute("Location")]
            public string Location = null;

            [XmlAttribute("Description")]
            public string Description = null;

            [XmlAttribute("ModelName")]
            public string ModelName = null;

            [XmlAttribute("FilePath")]
            public string FilePath = null;

            [XmlAttribute("Year")]
            public int Year;

            [XmlArrayItem("Waves")]
            public string[] DarwinSymbol = null;

            public LocalOceanModelEntry()
            { }

            public LocalOceanModelEntry(int Number, string Location, string ModelName, string FilePath, int Year, string[] DarwinSymbol)
            {
                this.ModelName = ModelName;
                this.FilePath = FilePath;
                this.Year = Year;
                this.DarwinSymbol = DarwinSymbol;
            }

            public int CompareTo(Object obj)
            {
                if (obj is LocalOceanModelEntry)
                    return this.Number.CompareTo((obj as LocalOceanModelEntry).Number);
                return 0;
            }
        }
        #endregion  

        /// <summary>
        /// Get the all possible usable Green's function, global ocean models, and local ocean models.
        /// </summary>
        /// <param name="path">Directory base path for determination of usable Green's function, global ocean models, and local ocean models.</param>
        public void GetParameter(string path)
        {
            // GreensFuntion
            getGreensFunction(path);

            // OceanModels
            getOceanModels(path);

            // LocalModels
            getLocalModels(path);
        }

        #region internal functions
        internal void getGreensFunction(string path)
        {
            string[] dirs = Directory.GetFiles(path + FilePath._EarthModels + Path.DirectorySeparatorChar);

            for (int i = 0; i < dirs.Length; i++)
            {
                GreenFunctionEntry gfn = new GreenFunctionEntry();

                gfn.FileName = Path.GetFileName(dirs[i]);

                string[] parts = Path.GetFileName(dirs[i]).Split(new char[] { '.' }, StringSplitOptions.RemoveEmptyEntries);

                switch (parts[1])
                {
                    case "gbaver": gfn.EarthModel = "gbaver"; break;
                    case "gbcont": gfn.EarthModel = "gbcont"; break;
                    case "gbocen": gfn.EarthModel = "gbcont"; break;
                }

                switch (parts[2])
                {
                    case "wef": gfn.Source = "W.E. Farrel"; break;
                }

                switch (parts[3])
                {
                    case "p01": gfn.Pattern = "coarse"; break;
                    case "p02": gfn.Pattern = "fine"; break;
                }

                switch (parts[4])
                {
                    case "ce": gfn.Frame = "Common"; break;
                    case "cm": gfn.Frame = "GPS"; break;
                }

                FileStream fs = new FileStream(dirs[i], FileMode.Open);
                GZipStream decompressingStream = new GZipStream(fs, CompressionMode.Decompress);

                using (var sr = new StreamReader(decompressingStream))
                {
                    string tmp = sr.ReadLine().TrimEnd();
                    tmp = tmp.Remove(0, 3).TrimStart();
                    gfn.Description = tmp;
                }

                AddGreenFunctionEntry(gfn);
            }

        }

        internal void getOceanModels(string path)
        {
            string[] dirs = Directory.GetDirectories(path + FilePath._OceanModels + Path.DirectorySeparatorChar);

            for (int i = 0; i < dirs.Length; i++)
            {
                OceanModelEntry tm = new OceanModelEntry();
                tm.FilePath = Path.GetFileName(dirs[i]);

                string[] waves = Directory.GetFiles(dirs[i]);

                tm.DarwinSymbol = new string[waves.Length];
                string[] parts = null;
                for (int j = 0; j < waves.Length; j++)
                {
                    parts = Path.GetFileName(waves[j]).Split(new char[] { '.' }, StringSplitOptions.None);
                    tm.DarwinSymbol[j] = parts[0];
                }
                tm.Year = int.Parse(parts[parts.Length - 3]);
                if (parts.Length == 5)
                    tm.ModelName = parts[1];
                if (parts.Length == 6)
                    tm.ModelName = parts[1] + "." + parts[2];

                FileStream fs = new FileStream(waves[0], FileMode.Open);
                GZipStream decompressingStream = new GZipStream(fs, CompressionMode.Decompress);

                using (var sr = new StreamReader(decompressingStream))
                {
                    string tmp = null;
                    for (int j = 0; j < 8; j++)
                        tmp = sr.ReadLine().Trim();
                    tm.Description = tmp;
                }
                AddOceanModelEntryEntry(tm);
            }
        }

        internal string[,] localModelInformation = new string[,] 
        {
            
            {"1","osu.bering.2010","Bering Sea"},
            {"2","osu.hawaii.2010","Pacific Ocean around Hawaii"},
            {"3","sfbay.1984","San Francisco Bay"},
            {"4","osu.usawest.2010","West coast of United States and British Columbia"},
            {"5","cortez.1976","Sea of Cortez/Gulf of California "},
            {"6","osu.gulfmex.2010","Gulf of Mexico"},
            {"7","osu.hudson.2010","Hudson Bay and surrounding waters"},
            {"8","osu.namereast.2010","East coast of North America, Maryland to Labrador"},
            {"9","osu.patagonia.2010","Patagonian shelf "},
            {"10","osu.amazon.2010","Off the mouth of the Amazon"},
            {"11","osu.europeshelf.2008","NW European shelf"},
            {"12","osu.mediterranean.2011","Mediterranean and Black Seas "},
            {"13","osu.redsea.2010","Red Sea"},
            {"14","osu.persian.2010","Arabian Sea and Persian Gulf "},
            {"15","osu.bengal.2010","Bay of Bengal "},
            {"16","osu.chinasea.2010","East China Sea, South China Sea"},
            {"17","osu.northaustral.2009","North of Australia, Indian Ocean to Tasman Sea"},
            {"18","osu.tasmania.2010","Bass Strait and parts of the Tasman Sea and Great Australian Bight"},
            {"19","osu.okhotsk.2010","Seas of Okhotsk and Japan, NE Pacific "},
            {"20","naoregional.1999","Sea of Japan area (used in GOTIC2 package)"},
            {"21","esr.aotim5.2004","Arctic Ocean, and part of the North Atlantic"},
            {"22","esr.cats.2008","Southern Ocean and Antarctic waters"}
        };


        internal void getLocalModels(string path)
        {
            string[] dirs = Directory.GetDirectories(path + "LocalModels" + Path.DirectorySeparatorChar);

            for (int i = 0; i < dirs.Length; i++)
            {
                LocalOceanModelEntry tm = new LocalOceanModelEntry();
                tm.FilePath = Path.GetFileName(dirs[i]);
                
                string[] waves = Directory.GetFiles(dirs[i], "*gz");
                
                if (waves.Length > 0)
                {
                    tm.DarwinSymbol = new string[waves.Length];
                    string[] parts = null;
                    for (int j = 0; j < waves.Length; j++)
                    {
                        parts = Path.GetFileName(waves[j]).Split(new char[] { '.' }, StringSplitOptions.None);
                        tm.DarwinSymbol[j] = parts[0];
                    }

                    for (int j = 0; j < localModelInformation.GetLength(0); j++)
                    {
                        string li = localModelInformation[j, 1];
                        string lm = Path.GetFileName(dirs[i]);

                        if (li == lm)
                        {
                            tm.Number = int.Parse(localModelInformation[j, 0]);
                            tm.Location = localModelInformation[j, 2];
                        }
                    }

                    tm.Year = int.Parse(parts[parts.Length - 3]);
                    if (parts.Length == 5)
                        tm.ModelName = parts[1];
                    if (parts.Length == 6)
                        tm.ModelName = parts[1] + "." + parts[2];

                    FileStream fs = new FileStream(waves[0], FileMode.Open);
                    GZipStream decompressingStream = new GZipStream(fs, CompressionMode.Decompress);

                    using (var sr = new StreamReader(decompressingStream))
                    {
                        string tmp = null;
                        for (int j = 0; j < 8; j++)
                            tmp = sr.ReadLine().Trim();
                        tm.Description = tmp;
                    }
                    AddLocalTideModelEntry(tm);
                }
            }
        }
        #endregion

        public static void serialisieren(OceanLoadingProperties oceanloadingproperties, string path)
        {
            XmlSerializer s = new XmlSerializer(typeof(OceanLoadingProperties));
            TextWriter w = new StreamWriter(path);
            s.Serialize(w, oceanloadingproperties);
            w.Close();
        } // Serialization
        public static OceanLoadingProperties deserialisieren(string path)
        {
            XmlSerializer s = new XmlSerializer(typeof(OceanLoadingProperties));
            OceanLoadingProperties oceanloadingproperties;
            TextReader r = new StreamReader(path);
            oceanloadingproperties = (OceanLoadingProperties)s.Deserialize(r);
            r.Close();
            return oceanloadingproperties;
        } // Deserialization
    }
}