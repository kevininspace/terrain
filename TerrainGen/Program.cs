using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Threading;
using Svg;
using Svg.Pathing;

namespace TerrainGen
{
    class Program
    {
        static Size defaultExtent = new Size()
        {
            Height = 400,
            Width = 400
        };



        static void Main(string[] args)
        {
            MapParams param = new MapParams()
            {
                extent = defaultExtent,
                npts = 20, //16384,
                ncities = 15,
                nterrs = 5,
                fontsizes = new FontSizes()
                {
                    region = 40,
                    city = 25,
                    town = 20
                }

            };
            
            TerrainMap tm = new TerrainMap(param);
            tm.doMap();
            //param.Generator = new Func<Delegate>(generateCoast);
            //doMap(new SvgDocument(), param);
        }

    }
}
