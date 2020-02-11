using System;
using System.Xml.Serialization;

namespace SPOTL
{
    public class Location
    {
        [XmlAttribute("Name")]
        public string Name = null;

        [XmlElement("Longitude")]
        public double Longitude = double.NaN;

        [XmlElement("Latitude")]
        public double Latitude = double.NaN;

        [XmlElement("Height")]
        public double Height = double.NaN;

        public Location()
        {

        }

        public Location(string Name, double Longitude, double Latitude, double Height)
        {
            this.Name = Name;
            this.Longitude = Longitude;
            this.Latitude = Latitude;
            this.Height = Height;
        }
    }

}
