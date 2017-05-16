using System;
using System.Collections.Generic;
using System.Drawing;

namespace TerrainGen
{
    public class Mesh
    {
        public List<Point> pts { get; set; }
        public List<GraphEdge> vor { get; set; }
        public List<Point> vxs { get; set; }

        public List<AdjacentSites> adj { get; set; }
        public List<Triangle> tris { get; set; }
        public List<GraphEdge> edges { get; set; }
        public Size extent { get; set; }

        public double Map(Func<double> f)
        {
            double fxnResult = f.Invoke();
            return fxnResult;
            //var mapped = vxs.map(f);
            //mapped.mesh = mesh;
            //return mapped;

        }
    }
}