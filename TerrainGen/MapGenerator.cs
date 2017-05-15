using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Svg;
using Svg.Pathing;

namespace TerrainGen
{

    internal class MapGenerator
    {
        public Size defaultExtent = new Size()
        {
            Height = 1000,
            Width = 1000
        };

        public List<Point> pts = new List<Point>();
        public List<GraphEdge> voronoiGraphEdges = new List<GraphEdge>();
        public Mesh mainMesh = new Mesh();
        public Random intRandom = new Random();
        private SvgDocument svg = new SvgDocument();

        public MapParams Params { get; set; }

        public void DoMap()
        {

            SvgUnit width = Params.extent.Width; //svg.Width;
            svg.Height = width * Params.extent.Height / Params.extent.Width;
            svg.ViewBox = new SvgViewBox(
                -1 * Params.extent.Width / 2,
                -1 * Params.extent.Height / 2,
                1 * Params.extent.Width,
                1 * Params.extent.Height);
            //svg.selectAll().remove();
            RenderObject render = new RenderObject();
            generateCoast(Params);
            //render.h =  paramObjects.generator(paramObjects);
            placeCities();
            drawMap(svg, render);
        }

        internal void drawMap(SvgDocument svg, RenderObject render)
        {
            //render.rivers = getRivers(render.h, 0.01);
            //render.coasts = contour(render.h, 0);
            //render.terr = getTerritories(render);
            //render.borders = getBorders(render);
            drawPaths(svg, "river", render.rivers);
            //drawPaths(svg, 'coast', render.coasts);
            //drawPaths(svg, 'border', render.borders);
            //visualizeSlopes(svg, render);
            //visualizeCities(svg, render);
            //drawLabels(svg, render);
        }

        private int runif(int lo, int hi)
        {
            return lo + intRandom.Next() * (hi - lo);
        }

        private double? rnorm()
        {
            double? z2 = null;
            if (z2 != null)
            {
                double? tmp = z2;
                z2 = null;
                return tmp;
            }
            var x1 = 0;
            var x2 = 0;
            var w = 2.0;
            while (w >= 1)
            {
                x1 = runif(-1, 1);
                x2 = runif(-1, 1);
                w = x1 * x1 + x2 * x2;
            }
            w = Math.Sqrt(-2 * Math.Log(w) / w);
            z2 = x2 * w;
            return x1 * w;
        }


        private double?[] randomVector(double? scale)
        {
            return new[] { scale * rnorm(), scale * rnorm() };
        }

        private List<Point> generatePoints(double n, Size extent)
        {
            extent = extent.IsEmpty ? defaultExtent : extent;
            for (int i = 0; i < n; i++)
            {
                Point p = new Point();
                p.setPoint((intRandom.NextDouble() - 0.5) * extent.Width, (intRandom.NextDouble() - 0.5) * extent.Height);
                pts.Add(p);
            }

            return pts;
        }

        private Point centroid(List<Point> pts)
        {
            double x = 0;
            double y = 0;
            for (int i = 0; i < pts.Count; i++)
            {
                x += pts[i].x;
                y += pts[i].y;
            }
            Point p = new Point();
            p.setPoint(x / pts.Count, y / pts.Count);
            return p;
        }

        private void improvePoints(double? n, Size extent)
        {
            n = n ?? 1;
            extent = extent.IsEmpty ? defaultExtent : extent;
            for (int i = 0; i < n; i++)
            {
                voronoi(pts, extent);

                //Polygon polygon = new Polygon(graphEdges);

                //pts.Select(p => centroid(p));
                //pts = voronoi(pts, extent)
                //    .polygons(pts)
                //    .map(centroid);
            }

        }

        private void generateGoodPoints(double n, Size extent)
        {
            extent = extent.IsEmpty ? defaultExtent : extent;
            pts = generatePoints(n, extent);
            pts.Sort(PointComparer);
            //pts = pts.sort(function(a, b) {
            //    return a[0] - b[0];
            //});
            improvePoints(1, extent);
        }

        private int PointComparer(Point a, Point b)
        {
            return (int)a.x - (int)b.x;
        }

        private void voronoi(List<Point> pts, Size extent)
        {
            extent = extent.IsEmpty ? defaultExtent : extent;
            double w = extent.Width / 2.0;
            double h = extent.Height / 2.0;

            Voronoi v = new Voronoi(5);
            double[] xs = pts.Select(x => x.x).ToArray();
            double[] ys = pts.Select(y => y.y).ToArray();
            voronoiGraphEdges = v.generateVoronoi(xs, ys, -w, w, -h, h);


            //return d3.voronoi().extent([[-w, -h], [w, h]])(pts);
        }

        private void makeMesh(Size extent)
        {
            extent = extent.IsEmpty ? defaultExtent : extent;

            List<Point> verticies = new List<Point>();
            List<VertexId> vertexIds = new List<VertexId>();
            List<AdjacentSites> adj = new List<AdjacentSites>();
            List<GraphEdge> edges = new List<GraphEdge>();
            List<Triangle> triangles = new List<Triangle>();
            foreach (GraphEdge graphEdge in voronoiGraphEdges)
            {
                if (graphEdge == null)
                {
                    continue;
                }
                //vertexIds[graphEdge.site1.sitenbr]
                //Point e0 = vertexIds.[graphEdge.site1];
                //        Point e1 = vertexIds[graphEdge.site2];
                //if (e0 == null)
                //{
                //    e0 = verticies.Count;
                //vertexIds[graphEdge.site1.sitenbr] = graphEdge.site1.sitenbr;
                verticies.Add(graphEdge.site1.coord);
                //}
                //if (e1 == null)
                //{
                //    e1 = verticies.Count;
                //vertexIds[graphEdge.site2.sitenbr] = graphEdge.site2.sitenbr;
                verticies.Add(graphEdge.site2.coord);
                //}
                //adj[e0] = adj[(int)e0]; // ||  double[];
                if (!adj.Contains(graphEdge.site1))
                {
                    adj.Add(new AdjacentSites()
                    {
                        coord = graphEdge.site1.coord,
                        nearbySites = new List<Site>() { graphEdge.site2 },
                        sitenbr = graphEdge.site1.sitenbr
                    });
                }
                if (!adj.Contains(graphEdge.site2))
                {
                    adj.Add(new AdjacentSites()
                    {
                        coord = graphEdge.site2.coord,
                        nearbySites = new List<Site>() { graphEdge.site1 },
                        sitenbr = graphEdge.site2.sitenbr
                    });
                }

                edges.Add(graphEdge);

                if (!triangles.Contains(graphEdge.site1))
                {
                    triangles.Add(new Triangle()
                    {
                        coord = graphEdge.site1.coord,
                        triangleEdges = new List<GraphEdge>() { graphEdge },
                        sitenbr = graphEdge.site1.sitenbr
                    });
                }

                //if (!triangles.Any(tri => tri.sitenbr == graphEdge.site2.sitenbr).Contains(graphEdge.site1))
                //{
                //    triangles.Add(new Triangle()
                //    {
                //        coord = graphEdge.site1.coord,
                //        triangleEdges = new List<GraphEdge>() { graphEdge },
                //        sitenbr = graphEdge.site1.sitenbr
                //    });
                //}

                //if (graphEdge.site2 && !triangles[e0].includes(graphEdge.site2))
                //    triangles[e0].push(graphEdge.site2);
                //triangles[e1] = triangles[e1]; // || []
                //;
                //if (!triangles[e1].includes(e.left))
                //    triangles[e1].push(e.left);
                //if (e.right && !triangles[e1].includes(e.right))
                //    triangles[e1].push(e.right);

            }

            mainMesh.pts = pts;
            //mainMesh.adj = adj;
            // mainMesh.edges = edges;
            //mainMesh.extent = extent;
            //mainMesh.tris = triangles;
            mainMesh.vor = voronoiGraphEdges;
            mainMesh.vxs = verticies;

        }






        private void generateGoodMesh(double n, Size extent)
        {
            extent = extent.IsEmpty ? defaultExtent : extent;
            generateGoodPoints(n, extent);
            makeMesh(extent);
        }

        private bool isedge(Mesh mesh, int i)
        {
            //return (mesh.adj[i].Length < 3);
            return false;
        }

        private bool isnearedge(Mesh mesh, int i)
        {
            var x = mesh.vxs[i].x;
            var y = mesh.vxs[i].y;
            var w = mesh.extent.Width;
            var h = mesh.extent.Height;
            return x < -0.45 * w || x > 0.45 * w || y < -0.45 * h || y > 0.45 * h;
        }

        private double[] neighbours(Mesh mesh, int i)
        {
            //double onbs = mesh.adj[i];
            //        double[] nbs = new double[] {};
            //for (var j = 0; j < onbs.Length; j++)
            //{
            //    nbs.push(onbs[j]);
            //}
            //return nbs;
            return new double[] { };
        }

        private double distance(Mesh mesh, int i, int j)
        {
            Point p = mesh.vxs[i];
            Point q = mesh.vxs[j];
            return Math.Sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y));
        }

        private void quantile(double h, double q)
        {
            //var sortedh = [];
            //for (var i = 0; i < h.length; i++)
            //{
            //    sortedh[i] = h[i];
            //}
            //sortedh.sort(d3.ascending);
            //return d3.quantile(sortedh, q);
        }

        private void zero(double mesh)
        {
            //var z = [];
            //for (var i = 0; i < mesh.vxs.length; i++)
            //{
            //    z[i] = 0;
            //}
            //z.mesh = mesh;
            //return z;
        }

        private void slope(double mesh, double direction)
        {
            //return mesh.map(function(x) {
            //    return x[0] * direction[0] + x[1] * direction[1];
            //});
        }

        private void cone(double mesh, double slope)
        {
            //return mesh.map(function(x) {
            //    return Math.pow(x[0] * x[0] + x[1] * x[1], 0.5) * slope;
            //});
        }

        private void map(double h, double f)
        {
            //var newh = h.map(f);
            //newh.mesh = h.mesh;
            //return newh;
        }

        private void normalize(double h)
        {
            //var lo = d3.min(h);
            //var hi = d3.max(h);
            //return map(h, function(x) { return (x - lo) / (hi - lo)});
        }

        private void peaky(double h)
        {
            //return map(normalize(h), Math.sqrt);
        }

        private void add()
        {
            //var n = arguments[0].length;
            //var newvals = zero(arguments[0].mesh);
            //for (var i = 0; i < n; i++)
            //{
            //    for (var j = 0; j < arguments.length; j++)
            //    {
            //        newvals[i] += arguments[j][i];
            //    }
            //}
            //return newvals;
        }

        private void mountains(double mesh, double n, double r)
        {
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

        private void relax(double h)
        {
            //    var newh = zero(h.mesh);
            //    for (var i = 0; i < h.length; i++)
            //    {
            //        var nbs = neighbours(h.mesh, i);
            //        if (nbs.length < 3)
            //        {
            //            newh[i] = 0;
            //            continue;
            //        }
            //        newh[i] = d3.mean(nbs.map(function(j) { return h[j]}));
            //}
            //    return newh;
        }

        private void downhill(double h)
        {
            //if (h.downhill) return h.downhill;
            //private void downfrom(i) {
            //    if (isedge(h.mesh, i)) return -2;
            //    var best = -1;
            //    var besth = h[i];
            //    var nbs = neighbours(h.mesh, i);
            //    for (var j = 0; j < nbs.length; j++)
            //    {
            //        if (h[nbs[j]] < besth)
            //        {
            //            besth = h[nbs[j]];
            //            best = nbs[j];
            //        }
            //    }
            //    return best;
            //}
            //var downs = [];
            //for (var i = 0; i < h.length; i++)
            //{
            //    downs[i] = downfrom(i);
            //}
            //h.downhill = downs;
            //return downs;
        }

        private void findSinks(double h)
        {
            //var dh = downhill(h);
            //var sinks = [];
            //for (var i = 0; i < dh.length; i++)
            //{
            //    var node = i;
            //    while (true)
            //    {
            //        if (isedge(h.mesh, node))
            //        {
            //            sinks[i] = -2;
            //            break;
            //        }
            //        if (dh[node] == -1)
            //        {
            //            sinks[i] = node;
            //            break;
            //        }
            //        node = dh[node];
            //    }
            //}
        }

        private void fillSinks(double h, double epsilon)
        {
            //epsilon = epsilon || 1e-5;
            //var infinity = 999999;
            //var newh = zero(h.mesh);
            //for (var i = 0; i < h.length; i++)
            //{
            //    if (isnearedge(h.mesh, i))
            //    {
            //        newh[i] = h[i];
            //    }
            //    else
            //    {
            //        newh[i] = infinity;
            //    }
            //}
            //while (true)
            //{
            //    var changed = false;
            //    for (var i = 0; i < h.length; i++)
            //    {
            //        if (newh[i] == h[i]) continue;
            //        var nbs = neighbours(h.mesh, i);
            //        for (var j = 0; j < nbs.length; j++)
            //        {
            //            if (h[i] >= newh[nbs[j]] + epsilon)
            //            {
            //                newh[i] = h[i];
            //                changed = true;
            //                break;
            //            }
            //            var oh = newh[nbs[j]] + epsilon;
            //            if ((newh[i] > oh) && (oh > h[i]))
            //            {
            //                newh[i] = oh;
            //                changed = true;
            //            }
            //        }
            //    }
            //    if (!changed) return newh;
            //}
        }

        private void getFlux(double h)
        {
            //var dh = downhill(h);
            //var idxs = [];
            //var flux = zero(h.mesh);
            //for (var i = 0; i < h.length; i++)
            //{
            //    idxs[i] = i;
            //    flux[i] = 1 / h.length;
            //}
            //idxs.sort(function(a, b) {
            //    return h[b] - h[a];
            //});
            //for (var i = 0; i < h.length; i++)
            //{
            //    var j = idxs[i];
            //    if (dh[j] >= 0)
            //    {
            //        flux[dh[j]] += flux[j];
            //    }
            //}
            //return flux;
        }

        private void getSlope(double h)
        {
            //var dh = downhill(h);
            //var slope = zero(h.mesh);
            //for (var i = 0; i < h.length; i++)
            //{
            //    var s = trislope(h, i);
            //    slope[i] = Math.sqrt(s[0] * s[0] + s[1] * s[1]);
            //    continue;
            //    if (dh[i] < 0)
            //    {
            //        slope[i] = 0;
            //    }
            //    else
            //    {
            //        slope[i] = (h[i] - h[dh[i]]) / distance(h.mesh, i, dh[i]);
            //    }
            //}
            //return slope;
        }

        private void erosionRate(double h)
        {
            //var flux = getFlux(h);
            //var slope = getSlope(h);
            //var newh = zero(h.mesh);
            //for (var i = 0; i < h.length; i++)
            //{
            //    var river = Math.sqrt(flux[i]) * slope[i];
            //    var creep = slope[i] * slope[i];
            //    var total = 1000 * river + creep;
            //    total = total > 200 ? 200 : total;
            //    newh[i] = total;
            //}
            //return newh;
        }

        private void erode(double h, double amount)
        {
            //var er = erosionRate(h);
            //var newh = zero(h.mesh);
            //var maxr = d3.max(er);
            //for (var i = 0; i < h.length; i++)
            //{
            //    newh[i] = h[i] - amount * (er[i] / maxr);
            //}
            //return newh;
        }

        private void doErosion(double h, double amount, double n)
        {
            //n = n || 1;
            //h = fillSinks(h);
            //for (var i = 0; i < n; i++)
            //{
            //    h = erode(h, amount);
            //    h = fillSinks(h);
            //}
            //return h;
        }

        private void setSeaLevel(double h, double q)
        {
            //var newh = zero(h.mesh);
            //var delta = quantile(h, q);
            //for (var i = 0; i < h.length; i++)
            //{
            //    newh[i] = h[i] - delta;
            //}
            //return newh;
        }

        private void cleanCoast(double h, double iters)
        {
            //for (var iter = 0; iter < iters; iter++)
            //{
            //    var changed = 0;
            //    var newh = zero(h.mesh);
            //    for (var i = 0; i < h.length; i++)
            //    {
            //        newh[i] = h[i];
            //        var nbs = neighbours(h.mesh, i);
            //        if (h[i] <= 0 || nbs.length != 3) continue;
            //        var count = 0;
            //        var best = -999999;
            //        for (var j = 0; j < nbs.length; j++)
            //        {
            //            if (h[nbs[j]] > 0)
            //            {
            //                count++;
            //            }
            //            else if (h[nbs[j]] > best)
            //            {
            //                best = h[nbs[j]];
            //            }
            //        }
            //        if (count > 1) continue;
            //        newh[i] = best / 2;
            //        changed++;
            //    }
            //    h = newh;
            //    newh = zero(h.mesh);
            //    for (var i = 0; i < h.length; i++)
            //    {
            //        newh[i] = h[i];
            //        var nbs = neighbours(h.mesh, i);
            //        if (h[i] > 0 || nbs.length != 3) continue;
            //        var count = 0;
            //        var best = 999999;
            //        for (var j = 0; j < nbs.length; j++)
            //        {
            //            if (h[nbs[j]] <= 0)
            //            {
            //                count++;
            //            }
            //            else if (h[nbs[j]] < best)
            //            {
            //                best = h[nbs[j]];
            //            }
            //        }
            //        if (count > 1) continue;
            //        newh[i] = best / 2;
            //        changed++;
            //    }
            //    h = newh;
            //}
            //return h;
        }

        private void trislope(double h, double i)
        {
            //var nbs = neighbours(h.mesh, i);
            //if (nbs.length != 3) return [0, 0];
            //var p0 = h.mesh.vxs[nbs[0]];
            //var p1 = h.mesh.vxs[nbs[1]];
            //var p2 = h.mesh.vxs[nbs[2]];

            //var x1 = p1[0] - p0[0];
            //var x2 = p2[0] - p0[0];
            //var y1 = p1[1] - p0[1];
            //var y2 = p2[1] - p0[1];

            //var det = x1 * y2 - x2 * y1;
            //var h1 = h[nbs[1]] - h[nbs[0]];
            //var h2 = h[nbs[2]] - h[nbs[0]];

            //return [(y2 * h1 - y1 * h2) / det,
            //        (-x2 * h1 + x1 * h2) / det];
        }

        private void cityScore(double h, double cities)
        {
            //var score = map(getFlux(h), Math.sqrt);
            //for (var i = 0; i < h.length; i++)
            //{
            //    if (h[i] <= 0 || isnearedge(h.mesh, i))
            //    {
            //        score[i] = -999999;
            //        continue;
            //    }
            //    score[i] += 0.01 / (1e-9 + Math.abs(h.mesh.vxs[i][0]) - h.mesh.extent.Width / 2)
            //    score[i] += 0.01 / (1e-9 + Math.abs(h.mesh.vxs[i][1]) - h.mesh.extent.Height / 2)
            //    for (var j = 0; j < cities.length; j++)
            //    {
            //        score[i] -= 0.02 / (distance(h.mesh, cities[j], i) + 1e-9);
            //    }
            //}
            //return score;
        }

        private void placeCity(double render)
        {
            //render.cities = render.cities || [];
            //var score = cityScore(render.h, render.cities);
            //var newcity = d3.scan(score, d3.descending);
            //render.cities.push(newcity);
        }

        private void placeCities()
        {
            //var params = render.params;
            //var h = render.h;
            //var n = params.ncities;
            //for (var i = 0; i < n; i++)
            //{
            //    placeCity(render);
            //}
        }

        private void contour(double h, double level)
        {
            //    level = level || 0;
            //    var edges = [];
            //    for (var i = 0; i < h.mesh.edges.length; i++)
            //    {
            //        var e = h.mesh.edges[i];
            //        if (e[3] == undefined) continue;
            //        if (isnearedge(h.mesh, e[0]) || isnearedge(h.mesh, e[1])) continue;
            //        if ((h[e[0]] > level && h[e[1]] <= level) ||
            //            (h[e[1]] > level && h[e[0]] <= level))
            //        {
            //            edges.push([e[2], e[3]]);
            //}
            //    }
            //    return mergeSegments(edges);
        }

        private void getRivers(double h, double limit)
        {
            //    var dh = downhill(h);
            //    var flux = getFlux(h);
            //    var links = [];
            //    var above = 0;
            //    for (var i = 0; i < h.length; i++)
            //    {
            //        if (h[i] > 0) above++;
            //    }
            //    limit *= above / h.length;
            //    for (var i = 0; i < dh.length; i++)
            //    {
            //        if (isnearedge(h.mesh, i)) continue;
            //        if (flux[i] > limit && h[i] > 0 && dh[i] >= 0)
            //        {
            //            var up = h.mesh.vxs[i];
            //            var down = h.mesh.vxs[dh[i]];
            //            if (h[dh[i]] > 0)
            //            {
            //                links.push([up, down]);
            //} else {
            //                links.push([up, [(up[0] + down[0])/2, (up[1] + down[1])/2]]);
            //            }
            //        }
            //    }
            //    return mergeSegments(links).map(relaxPath);
        }

        private void getTerritories(double render)
        {
            //    var h = render.h;
            //    var cities = render.cities;
            //    var n = render.params.nterrs;
            //    if (n > render.cities.length) n = render.cities.length;
            //    var flux = getFlux(h);
            //    var terr = [];
            //    var queue = new PriorityQueue({ comparator: function(a, b) { return a.score - b.score}
            //});
            //    private void weight(double u, double v)
            //{
            //    var horiz = distance(h.mesh, u, v);
            //    var vert = h[v] - h[u];
            //    if (vert > 0) vert /= 10;
            //    var diff = 1 + 0.25 * Math.pow(vert / horiz, 2);
            //    diff += 100 * Math.sqrt(flux[u]);
            //    if (h[u] <= 0) diff = 100;
            //    if ((h[u] > 0) != (h[v] > 0)) return 1000;
            //    return horiz * diff;
            //}
            //    for (var i = 0; i<n; i++) {
            //        terr[cities[i]] = cities[i];
            //        var nbs = neighbours(h.mesh, cities[i]);
            //        for (var j = 0; j<nbs.length; j++) {
            //            queue.queue({
            //                score: weight(cities[i], nbs[j]),
            //                city: cities[i],
            //                vx: nbs[j]
            //            });
            //        }
            //    }
            //    while (queue.length) {
            //        var u = queue.dequeue();
            //        if (terr[u.vx] != undefined) continue;
            //        terr[u.vx] = u.city;
            //        var nbs = neighbours(h.mesh, u.vx);
            //        for (var i = 0; i<nbs.length; i++) {
            //            var v = nbs[i];
            //            if (terr[v] != undefined) continue;
            //            var newdist = weight(u.vx, v);
            //queue.queue({
            //                score: u.score + newdist,
            //                city: u.city,
            //                vx: v
            //            });
            //        }
            //    }
            //    terr.mesh = h.mesh;
            //    return terr;
        }

        private void getBorders(double render)
        {
            //    var terr = render.terr;
            //    var h = render.h;
            //    var edges = [];
            //    for (var i = 0; i < terr.mesh.edges.length; i++)
            //    {
            //        var e = terr.mesh.edges[i];
            //        if (e[3] == undefined) continue;
            //        if (isnearedge(terr.mesh, e[0]) || isnearedge(terr.mesh, e[1])) continue;
            //        if (h[e[0]] < 0 || h[e[1]] < 0) continue;
            //        if (terr[e[0]] != terr[e[1]])
            //        {
            //            edges.push([e[2], e[3]]);
            //}
            //    }
            //    return mergeSegments(edges).map(relaxPath);
        }

        private void mergeSegments(double segs)
        {
            //var adj = { };
            //for (var i = 0; i < segs.length; i++)
            //{
            //    var seg = segs[i];
            //    var a0 = adj[seg[0]] || [];
            //    var a1 = adj[seg[1]] || [];
            //    a0.push(seg[1]);
            //    a1.push(seg[0]);
            //    adj[seg[0]] = a0;
            //    adj[seg[1]] = a1;
            //}
            //var done = [];
            //var paths = [];
            //var path = null;
            //while (true)
            //{
            //    if (path == null)
            //    {
            //        for (var i = 0; i < segs.length; i++)
            //        {
            //            if (done[i]) continue;
            //            done[i] = true;
            //            path = [segs[i][0], segs[i][1]];
            //            break;
            //        }
            //        if (path == null) break;
            //    }
            //    var changed = false;
            //    for (var i = 0; i < segs.length; i++)
            //    {
            //        if (done[i]) continue;
            //        if (adj[path[0]].length == 2 && segs[i][0] == path[0])
            //        {
            //            path.unshift(segs[i][1]);
            //        }
            //        else if (adj[path[0]].length == 2 && segs[i][1] == path[0])
            //        {
            //            path.unshift(segs[i][0]);
            //        }
            //        else if (adj[path[path.length - 1]].length == 2 && segs[i][0] == path[path.length - 1])
            //        {
            //            path.push(segs[i][1]);
            //        }
            //        else if (adj[path[path.length - 1]].length == 2 && segs[i][1] == path[path.length - 1])
            //        {
            //            path.push(segs[i][0]);
            //        }
            //        else
            //        {
            //            continue;
            //        }
            //        done[i] = true;
            //        changed = true;
            //        break;
            //    }
            //    if (!changed)
            //    {
            //        paths.push(path);
            //        path = null;
            //    }
            //}
            //return paths;
        }

        private void relaxPath(double path)
        {
            //var newpath = [path[0]];
            //for (var i = 1; i < path.length - 1; i++)
            //{
            //    var newpt = [0.25 * path[i - 1][0] + 0.5 * path[i][0] + 0.25 * path[i + 1][0],
            //                 0.25 * path[i - 1][1] + 0.5 * path[i][1] + 0.25 * path[i + 1][1]];
            //    newpath.push(newpt);
            //}
            //newpath.push(path[path.length - 1]);
            //return newpath;
        }

        private void visualizePoints(double svg, double pts)
        {
            //var circle = svg.selectAll('circle').data(pts);
            //circle.enter()
            //    .append('circle');
            //circle.exit().remove();
            //d3.selectAll('circle')
            //    .attr('cx', function(d) { return 1000 * d[0]})
            //    .attr('cy', function(d) { return 1000 * d[1]})
            //    .attr('r', 100 / Math.sqrt(pts.length));
        }

        private string makeD3Path(double path)
        {
            SvgGroup group = new SvgGroup();
            svg.Children.Add(group);

            group.Children.Add(new SvgCircle
            {
                Radius = 100,
                Fill = new SvgColourServer(Color.Red),
                Stroke = new SvgColourServer(Color.Black),
                StrokeWidth = 2
            });

            //SvgPath p = new SvgPath();
            SvgPathSegmentList pathlist = new SvgPathSegmentList();
            //SvgMoveToSegment
            SvgPathSegment spSegment = new SvgLineSegment(new PointF(25, 63), new PointF(789, 854));
            pathlist.Add(spSegment);
            //p.PathData = pathlist;

            group.Children.Add(new SvgPath
            {
                PathData = pathlist,
                Stroke = new SvgColourServer(Color.Red),
                StrokeWidth = 3
            });

            //var p = d3.path();
            //p.moveTo(1000 * path[0][0], 1000 * path[0][1]);
            //for (var i = 1; i < path.length; i++)
            //{
            //    p.lineTo(1000 * path[i][0], 1000 * path[i][1]);
            //}
            //return p.toString();
            return "";
        }

        private void visualizeVoronoi(double svg, double field, double lo, double hi)
        {
            //if (hi == undefined) hi = d3.max(field) + 1e-9;
            //if (lo == undefined) lo = d3.min(field) - 1e-9;
            //var mappedvals = field.map(function(x) { return x > hi ? 1 : x < lo ? 0 : (x - lo) / (hi - lo)});
            //var tris = svg.selectAll('path.field').data(field.mesh.tris)
            //tris.enter()
            //    .append('path')
            //    .classed('field', true);

            //tris.exit()
            //    .remove();

            //svg.selectAll('path.field')
            //    .attr('d', makeD3Path)
            //    .style('fill', function(d, i) {
            //    return d3.interpolateViridis(mappedvals[i]);
            //});
        }

        private void visualizeDownhill(double h)
        {
            //var links = getRivers(h, 0.01);
            //drawPaths('river', links);
        }

        private void drawPaths(SvgDocument svg, string cls, double paths)
        {
            //var paths = svg.selectAll('path.' + cls).data(paths)
            //paths.enter()
            //        .append('path')
            //        .classed(cls, true)
            //paths.exit()
            //        .remove();




            string content = svg.Content;
            svg.Write("kevinsvg.svg");

            Bitmap bitmap = new Bitmap(1000, 1000);
            svg.Draw(bitmap);
            bitmap.Save("kevintest.bmp");


            //svg.selectAll('path.' + cls)
            //    .attr('d', makeD3Path);
        }

        private void visualizeSlopes(double svg, double render)
        {
            //    var h = render.h;
            //    var strokes = [];
            //    var r = 0.25 / Math.sqrt(h.length);
            //    for (var i = 0; i < h.length; i++)
            //    {
            //        if (h[i] <= 0 || isnearedge(h.mesh, i)) continue;
            //        var nbs = neighbours(h.mesh, i);
            //        nbs.push(i);
            //        var s = 0;
            //        var s2 = 0;
            //        for (var j = 0; j < nbs.length; j++)
            //        {
            //            var slopes = trislope(h, nbs[j]);
            //            s += slopes[0] / 10;
            //            s2 += slopes[1];
            //        }
            //        s /= nbs.length;
            //        s2 /= nbs.length;
            //        if (Math.abs(s) < runif(0.1, 0.4)) continue;
            //        var l = r * runif(1, 2) * (1 - 0.2 * Math.pow(Math.atan(s), 2)) * Math.exp(s2 / 100);
            //        var x = h.mesh.vxs[i][0];
            //        var y = h.mesh.vxs[i][1];
            //        if (Math.abs(l * s) > 2 * r)
            //        {
            //            var n = Math.floor(Math.abs(l * s / r));
            //            l /= n;
            //            if (n > 4) n = 4;
            //            for (var j = 0; j < n; j++)
            //            {
            //                var u = rnorm() * r;
            //                var v = rnorm() * r;
            //                strokes.push([[x + u - l, y + v + l * s], [x+u+l, y+v-l*s]]);
            //            }
            //        } else {
            //            strokes.push([[x-l, y+l*s], [x+l, y-l*s]]);
            //        }
            //    }
            //    var lines = svg.selectAll('line.slope').data(strokes)
            //    lines.enter()
            //            .append('line')
            //            .classed('slope', true);
            //lines.exit()
            //            .remove();
            //svg.selectAll('line.slope')
            //        .attr('x1', function (d) { return 1000 * d[0][0]})
            //        .attr('y1', function (d) { return 1000 * d[0][1]})
            //        .attr('x2', function (d) { return 1000 * d[1][0]})
            //        .attr('y2', function (d) { return 1000 * d[1][1]})
        }

        private void visualizeContour(double h, double level)
        {
            //level = level || 0;
            //var links = contour(h, level);
            //drawPaths('coast', links);
        }

        private void visualizeBorders(double h, double cities, double n)
        {
            //var links = getBorders(h, getTerritories(h, cities, n));
            //drawPaths('border', links);
        }

        private void visualizeCities(double svg, double render)
        {
            //var cities = render.cities;
            //var h = render.h;
            //var n = render.params.nterrs;

            //var circs = svg.selectAll('circle.city').data(cities);
            //circs.enter()
            //        .append('circle')
            //        .classed('city', true);
            //circs.exit()
            //        .remove();
            //svg.selectAll('circle.city')
            //    .attr('cx', function(d) { return 1000 * h.mesh.vxs[d][0]})
            //    .attr('cy', function(d) { return 1000 * h.mesh.vxs[d][1]})
            //    .attr('r', function(d, i) { return i >= n ? 4 : 10})
            //    .style('fill', 'white')
            //    .style('stroke-Width', 5)
            //    .style('stroke-linecap', 'round')
            //    .style('stroke', 'black')
            //    .raise();
        }

        private void dropEdge(double h, double p)
        {
            //p = p || 4
            //var newh = zero(h.mesh);
            //for (var i = 0; i < h.length; i++)
            //{
            //    var v = h.mesh.vxs[i];
            //    var x = 2.4 * v[0] / h.mesh.extent.Width;
            //    var y = 2.4 * v[1] / h.mesh.extent.Height;
            //    newh[i] = h[i] - Math.exp(10 * (Math.pow(Math.pow(x, p) + Math.pow(y, p), 1 / p) - 1));
            //}
            //return newh;
        }

        private void generateCoast(MapParams mapParams)
        {
            generateGoodMesh(mapParams.npts, mapParams.extent);
            //var h = add(
            //        slope(mesh, randomVector(4)),
            //        cone(mesh, runif(-1, -1)),
            //        mountains(mesh, 50)
            //        );
            //for (var i = 0; i < 10; i++)
            //{
            //    h = relax(h);
            //}
            //h = peaky(h);
            //h = doErosion(h, runif(0, 0.1), 5);
            //h = setSeaLevel(h, runif(0.2, 0.6));
            //h = fillSinks(h);
            //h = cleanCoast(h, 3);
            //return h;
            return;
        }

        private void terrCenter(double h, double terr, double city, double landOnly)
        {
            //var x = 0;
            //var y = 0;
            //var n = 0;
            //for (var i = 0; i < terr.length; i++)
            //{
            //    if (terr[i] != city) continue;
            //    if (landOnly && h[i] <= 0) continue;
            //    x += terr.mesh.vxs[i][0];
            //    y += terr.mesh.vxs[i][1];
            //    n++;
            //}
            //return [x / n, y / n];
        }

        private void drawLabels(double svg, double render)
        {
            //    var params = render.params;
            //    var h = render.h;
            //    var terr = render.terr;
            //    var cities = render.cities;
            //    var nterrs = render.params.nterrs;
            //    var avoids = [render.rivers, render.coasts, render.borders];
            //    var lang = makeRandomLanguage();
            //    var citylabels = [];
            //    function penalty(label) {
            //        var pen = 0;
            //        if (label.x0 < -0.45 * h.mesh.extent.Width) pen += 100;
            //        if (label.x1 > 0.45 * h.mesh.extent.Width) pen += 100;
            //        if (label.y0 < -0.45 * h.mesh.extent.Height) pen += 100;
            //        if (label.y1 > 0.45 * h.mesh.extent.Height) pen += 100;
            //        for (var i = 0; i < citylabels.length; i++)
            //        {
            //            var olabel = citylabels[i];
            //            if (label.x0 < olabel.x1 && label.x1 > olabel.x0 &&
            //                label.y0 < olabel.y1 && label.y1 > olabel.y0)
            //            {
            //                pen += 100;
            //            }
            //        }

            //        for (var i = 0; i < cities.length; i++)
            //        {
            //            var c = h.mesh.vxs[cities[i]];
            //            if (label.x0 < c[0] && label.x1 > c[0] && label.y0 < c[1] && label.y1 > c[1])
            //            {
            //                pen += 100;
            //            }
            //        }
            //        for (var i = 0; i < avoids.length; i++)
            //        {
            //            var avoid = avoids[i];
            //            for (var j = 0; j < avoid.length; j++)
            //            {
            //                var avpath = avoid[j];
            //                for (var k = 0; k < avpath.length; k++)
            //                {
            //                    var pt = avpath[k];
            //                    if (pt[0] > label.x0 && pt[0] < label.x1 && pt[1] > label.y0 && pt[1] < label.y1)
            //                    {
            //                        pen++;
            //                    }
            //                }
            //            }
            //        }
            //        return pen;
            //    }
            //    for (var i = 0; i < cities.length; i++)
            //    {
            //        var x = h.mesh.vxs[cities[i]][0];
            //        var y = h.mesh.vxs[cities[i]][1];
            //        var text = makeName(lang, 'city');
            //        var size = i < nterrs ? params.fontsizes.city : params.fontsizes.town;
            //        var sx = 0.65 * size / 1000 * text.length;
            //        var sy = size / 1000;
            //        var posslabels = [
            //        {
            //            x: x + 0.8 * sy,
            //            y: y + 0.3 * sy,
            //            align: 'start',
            //            x0: x + 0.7 * sy,
            //            y0: y - 0.6 * sy,
            //            x1: x + 0.7 * sy + sx,
            //            y1: y + 0.6 * sy
            //        },
            //        {
            //            x: x - 0.8 * sy,
            //            y: y + 0.3 * sy,
            //            align: 'end',
            //            x0: x - 0.9 * sy - sx,
            //            y0: y - 0.7 * sy,
            //            x1: x - 0.9 * sy,
            //            y1: y + 0.7 * sy
            //        },
            //        {
            //            x: x,
            //            y: y - 0.8 * sy,
            //            align: 'middle',
            //            x0: x - sx / 2,
            //            y0: y - 1.9 * sy,
            //            x1: x + sx / 2,
            //            y1: y - 0.7 * sy
            //        },
            //        {
            //            x: x,
            //            y: y + 1.2 * sy,
            //            align: 'middle',
            //            x0: x - sx / 2,
            //            y0: y + 0.1 * sy,
            //            x1: x + sx / 2,
            //            y1: y + 1.3 * sy
            //        }
            //        ];
            //        var label = posslabels[d3.scan(posslabels, function(a, b) { return penalty(a) - penalty(b)})];
            //    label.text = text;
            //    label.size = size;
            //    citylabels.push(label);
            //}
            //var texts = svg.selectAll('text.city').data(citylabels);
            //texts.enter()
            //        .append('text')
            //        .classed('city', true);
            //texts.exit()
            //        .remove();
            //svg.selectAll('text.city')
            //        .attr('x', function (d) { return 1000 * d.x})
            //        .attr('y', function (d) { return 1000 * d.y})
            //        .style('font-size', function (d) { return d.size})
            //        .style('text-anchor', function (d) { return d.align})
            //        .text(function (d) { return d.text})
            //        .raise();

            //var reglabels = [];
            //    for (var i = 0; i<nterrs; i++) {
            //        var city = cities[i];
            //var text = makeName(lang, 'region');
            //var sy = params.fontsizes.region / 1000;
            //        var sx = 0.6 * text.length * sy;
            //var lc = terrCenter(h, terr, city, true);
            //var oc = terrCenter(h, terr, city, false);
            //var best = 0;
            //var bestscore = -999999;
            //        for (var j = 0; j<h.length; j++) {
            //            var score = 0;
            //var v = h.mesh.vxs[j];
            //score -= 3000 * Math.sqrt((v[0] - lc[0]) * (v[0] - lc[0]) + (v[1] - lc[1]) * (v[1] - lc[1]));
            //            score -= 1000 * Math.sqrt((v[0] - oc[0]) * (v[0] - oc[0]) + (v[1] - oc[1]) * (v[1] - oc[1]));
            //            if (terr[j] != city) score -= 3000;
            //            for (var k = 0; k<cities.length; k++) {
            //                var u = h.mesh.vxs[cities[k]];
            //                if (Math.abs(v[0] - u[0]) < sx && 
            //                    Math.abs(v[1] - sy/2 - u[1]) < sy) {
            //                    score -= k<nterrs? 4000 : 500;
            //                }
            //                if (v[0] - sx/2 < citylabels[k].x1 &&
            //                    v[0] + sx/2 > citylabels[k].x0 &&
            //                    v[1] - sy<citylabels[k].y1 &&
            //                    v[1]> citylabels[k].y0) {
            //                    score -= 5000;
            //                }
            //            }
            //            for (var k = 0; k<reglabels.length; k++) {
            //                var label = reglabels[k];
            //                if (v[0] - sx/2 < label.x + label.Width/2 &&
            //                    v[0] + sx/2 > label.x - label.Width/2 &&
            //                    v[1] - sy<label.y &&
            //                    v[1]> label.y - label.size)
            //{
            //    score -= 20000;
            //}
            //            }
            //            if (h[j] <= 0) score -= 500;
            //            if (v[0] + sx/2 > 0.5 * h.mesh.extent.Width) score -= 50000;
            //            if (v[0] - sx/2 < -0.5 * h.mesh.extent.Width) score -= 50000;
            //            if (v[1] > 0.5 * h.mesh.extent.Height) score -= 50000;
            //            if (v[1] - sy< -0.5 * h.mesh.extent.Height) score -= 50000;
            //            if (score > bestscore) {
            //                bestscore = score;
            //                best = j;
            //            }
            //        }
            //        reglabels.push({
            //            text: text, 
            //            x: h.mesh.vxs[best][0], 
            //            y: h.mesh.vxs[best][1], 
            //            size:sy, 
            //            Width:sx
            //        });
            //    }
            //    texts = svg.selectAll('text.region').data(reglabels);
            //texts.enter()
            //        .append('text')
            //        .classed('region', true);
            //texts.exit()
            //        .remove();
            //svg.selectAll('text.region')
            //        .attr('x', function (d) { return 1000 * d.x})
            //        .attr('y', function (d) { return 1000 * d.y})
            //        .style('font-size', function (d) { return 1000 * d.size})
            //        .style('text-anchor', 'middle')
            //        .text(function (d) { return d.text})
            //        .raise();

        }
    }

    public class MapParams
    {
        public Size extent;
        public int npts { get; set; }
        public int ncities { get; set; }
        public int nterrs { get; set; }
        public FontSizes fontsizes { get; set; }

        public delegate void Generator(MapParams paramObjects);
    }

    internal class RenderObject
    {
        public double rivers;
        public object h { get; set; }
    }

    public class FontSizes
    {
        public int region { get; set; }
        public int city { get; set; }
        public int town { get; set; }
    }
}