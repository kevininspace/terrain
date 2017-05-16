using System.Collections.Generic;
using System.Drawing;

namespace TerrainGen
{
    public class RenderObject
    {
        public List<PointF> rivers;
        public HeightMap h { get; set; }

        public RenderObject()
        {
            h = new HeightMap();
            rivers = new List<PointF>();
        }
    }
}