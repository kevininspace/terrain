using System.Collections.Generic;

namespace TerrainGen
{
    public class SiteSorterYX : IComparer<Site>
    {
        public int Compare(Site p1, Site p2)
        {
            Point s1 = p1.coord;
            Point s2 = p2.coord;
            if (s1.y < s2.y) return -1;
            if (s1.y > s2.y) return 1;
            if (s1.x < s2.x) return -1;
            if (s1.x > s2.x) return 1;
            return 0;
        }
    }
}