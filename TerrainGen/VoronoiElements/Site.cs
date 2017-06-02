using System.Collections.Generic;

namespace TerrainGen
{
    public class Site
    {
        public Point coord;
        public int sitenbr;
        public double z;

        public Site()
        {
            coord = new Point();
            points = new List<Point>();
        }

        public List<Point> points { get; set; }
    }
}