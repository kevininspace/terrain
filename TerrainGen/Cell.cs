using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TerrainGen
{
    public class Cell
    {
        public object[] site { get; set; }
        public Halfedge halfedges { get; set; }
        private Cell cell;
        private int cellnb; // = new int[] {};
        //        import {createBorderEdge
        //    }
        //    from "./Edge";
        //import {cells, edges, epsilon
        //}
        //from "./Diagram";
        //public Cell()
        //{

        //}



        public Cell createCell(object[] site)
        {
            //Cell c = new Cell
            //{
            //    site = site,
            //    halfedges = new Halfedge(),
            //    cellnb = site.edgenbr

            //};
            //return c;
            ////return cells[site.edgenbr] = new Cell()
            ////    site: site,
            ////halfedges: []
            ////};
            return null;
        }

        public double cellHalfedgeAngle(Cell cell, Edge edge)
        {
            //object[] site = cell.site;
            //object[] va = edge.ep;
            //object[] vb = edge.reg;
            ////var va = edge.Left;
            ////        var vb = edge.Right;
            //if (site.Equals(vb))
            //{
            //    vb = va;
            //    va = site;
            //}
            //if (vb != null)
            //{
            //    return Math.Atan2(vb.coord.y - va.coord.y, vb.coord.x - va.coord.x);
            //}
            //if (site.Equals(va))
            //{
            //    //va = edge[1];
            //    // vb = edge[0];
            //}
            //else
            //{
            //    //va = edge[0];
            //    //vb = edge[1];
            //}
            //return Math.Atan2(va.coord.x - vb.coord.x, vb.coord.y - va.coord.y);
            return 0.0;
        }

        public void cellHalfedgeStart(Cell cell, Edge edge)
        {
            //return edge[+(edge.left !== cell.site)];
        }

        public void cellHalfedgeEnd(Cell cell, Edge edge)
        {
            //return edge[+(edge.left === cell.site)];
        }

        public void sortCellHalfedges()
        {
            //    for (var i = 0, n = cells.length, cell, halfedges, j, m; i < n; ++i)
            //    {
            //        if ((cell = cells[i]) && (m = (halfedges = cell.halfedges).length))
            //        {
            //            var index = new Array(m),
            //                array = new Array(m);
            //            for (j = 0; j < m; ++j) index[j] = j, array[j] = cellHalfedgeAngle(cell, edges[halfedges[j]]);
            //            index.sort(function(i, j) { return array[j] - array[i]; });
            //    for (j = 0; j < m; ++j) array[j] = halfedges[index[j]];
            //    for (j = 0; j < m; ++j) halfedges[j] = array[j];
            //}
            //  }
        }

        public void clipCells(int x0, int y0, int x1, int y1)
        {
            //var nCells = cells.length;
            //    iCell,
            //    cell,
            //    site,
            //    iHalfedge,
            //    halfedges,
            //    nHalfedges,
            //    start,
            //    startX,
            //    startY,
            //    end,
            //    endX,
            //    endY,
            //    cover = true;

            //for (iCell = 0; iCell < nCells; ++iCell)
            //{
            //    if (cell = cells[iCell])
            //    {
            //        site = cell.site;
            //        halfedges = cell.halfedges;
            //        iHalfedge = halfedges.length;

            //        // Remove any dangling clipped edges.
            //        while (iHalfedge--)
            //        {
            //            if (!edges[halfedges[iHalfedge]])
            //            {
            //                halfedges.splice(iHalfedge, 1);
            //            }
            //        }

            //        // Insert any border edges as necessary.
            //        iHalfedge = 0, nHalfedges = halfedges.length;
            //        while (iHalfedge < nHalfedges)
            //        {
            //            end = cellHalfedgeEnd(cell, edges[halfedges[iHalfedge]]), endX = end[0], endY = end[1];
            //            start = cellHalfedgeStart(cell, edges[halfedges[++iHalfedge % nHalfedges]]), startX = start[0], startY = start[1];
            //            if (Math.abs(endX - startX) > epsilon || Math.abs(endY - startY) > epsilon)
            //            {
            //                halfedges.splice(iHalfedge, 0, edges.push(createBorderEdge(site, end,
            //                    Math.abs(endX - x0) < epsilon && y1 - endY > epsilon?[x0, Math.abs(startX - x0) < epsilon ? startY : y1]
            //          : Math.abs(endY - y1) < epsilon && x1 - endX > epsilon?[Math.abs(startY - y1) < epsilon ? startX : x1, y1]
            //          : Math.abs(endX - x1) < epsilon && endY - y0 > epsilon?[x1, Math.abs(startX - x1) < epsilon ? startY : y0]
            //          : Math.abs(endY - y0) < epsilon && endX - x0 > epsilon?[Math.abs(startY - y0) < epsilon ? startX : x0, y0]
            //          : null)) - 1);
            //                ++nHalfedges;
            //            }
            //        }

            //        if (nHalfedges) cover = false;
            //    }
            //}

            //// If there weren’t any edges, have the closest site cover the extent.
            //// It doesn’t matter which corner of the extent we measure!
            //if (cover)
            //{
            //    var dx, dy, d2, dc = Infinity;

            //    for (iCell = 0, cover = null; iCell < nCells; ++iCell)
            //    {
            //        if (cell = cells[iCell])
            //        {
            //            site = cell.site;
            //            dx = site[0] - x0;
            //            dy = site[1] - y0;
            //            d2 = dx * dx + dy * dy;
            //            if (d2 < dc) dc = d2, cover = cell;
            //        }
            //    }

            //    if (cover)
            //    {
            //        var v00 = [x0, y0], v01 = [x0, y1], v11 = [x1, y1], v10 = [x1, y0];
            //        cover.halfedges.push(
            //          edges.push(createBorderEdge(site = cover.site, v00, v01)) - 1,
            //          edges.push(createBorderEdge(site, v01, v11)) - 1,
            //          edges.push(createBorderEdge(site, v11, v10)) - 1,
            //          edges.push(createBorderEdge(site, v10, v00)) - 1
            //        );
            //    }
            //}

            //// Lastly delete any cells with no edges; these were entirely clipped.
            //for (iCell = 0; iCell < nCells; ++iCell)
            //{
            //    if (cell = cells[iCell])
            //    {
            //        if (!cell.halfedges.length)
            //        {
            //            delete cells[iCell];
            //        }
            //    }
            //}
        }
    }
}
