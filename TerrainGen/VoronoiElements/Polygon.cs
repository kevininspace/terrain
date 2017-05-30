using System.Collections.Generic;

namespace TerrainGen
{
    public class Polygon
    {
        public List<GraphEdge> ge;

        public Polygon(List<GraphEdge> ge)
        {
            this.ge = ge;
        }
        //    Edge edges = this.edges;

        //return this.cells.map(function(cell)
        //    {
        //        var polygon = cell.halfedges.map(function(i) { return cellHalfedgeStart(cell, edges[i]); });
        //        polygon.data = cell.site.data;
        //        return polygon;
        //    });
    }
}