using System.Drawing;

namespace TerrainGen
{
    public class MapParams
    {
        public Size extent;
        public int npts { get; set; }
        public int ncities { get; set; }
        public int nterrs { get; set; }
        public FontSizes fontsizes { get; set; }

        public delegate void Generator(MapParams paramObjects);
    }
}