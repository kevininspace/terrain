using System.Collections.Generic;
using System.Linq;

namespace TerrainGen
{
    public class HeightMap
    {
        public List<Site> Heights { get; set; }
        public Mesh mesh { get; set; }
        public object[] downhill = new object[5];

        public double[] Zero(Mesh mesh)
        {
            List<double> z = new List<double>();
            //z.ForEach(zero => zero = 0);

            IEnumerable<double> zz = mesh.vxs.Zip(z, (number, zeroes) => zeroes = 0);

            //for (var i = 0; i < mesh.vxs.Count; i++)
            //{
            //    z.Add = 0;
            //}
            //z.mesh = mesh;
            return zz.ToArray();
        }

        public void AddSlope(double[] zeroes, double[] direction)
        {
            //return mesh.Map(SlopeFunction(new[] { 2.0, 3.2 }, direction));
        }

        public double SlopeFunction(double[] x, double[] direction)
        {
            return x[0] * direction[0] + x[1] * direction[1];
        }

        public void AddCone(Mesh mesh, double slope)
        {
            //return mesh.map(function(x) {
            //    return Math.pow(x[0] * x[0] + x[1] * x[1], 0.5) * slope;
            //});
        }

        public void AddMountains(Mesh mesh, int i)
        {
            //private void mountains(Mesh mesh, double n, double r)

            //    r = r || 0.05;
            //    var mounts = [];
            //    for (var i = 0; i < n; i++)
            //    {
            //        mounts.push([mesh.extent.Width * (Math.random() - 0.5), mesh.extent.Height * (Math.random() - 0.5)]);
            //}
            //var newvals = zero(double mesh);
            //    for (var i = 0; i<mesh.vxs.length; i++) {
            //        var p = mesh.vxs[i];
            //        for (var j = 0; j<n; j++) {
            //            var m = mounts[j];
            //newvals[i] += Math.pow(Math.exp(-((p[0] - m[0]) * (p[0] - m[0]) + (p[1] - m[1]) * (p[1] - m[1])) / (2 * r * r)), 2);
            //        }
            //    }
            //    return newvals;
        }

        //For the length of the slope, first zero it, then adjust each by adding a slope, then cone, then mountain
        //private void add()
        //{
        //    var n = arguments[0].length; //length of slope
        //    var newvals = zero(arguments[0].mesh);
        //    for (var i = 0; i < n; i++)
        //    {
        //        for (var j = 0; j < arguments.length; j++)
        //        {
        //            newvals[i] += arguments[j][i];
        //        }
        //    }
        //    return newvals;
        //}
    }
}