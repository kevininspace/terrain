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
        static readonly Size DefaultExtent = new Size()
        {
            Height = 1000,
            Width = 1000
        };

        static void Main(string[] args)
        {
            MapParams param = new MapParams()
            {
                extent = DefaultExtent,
                npts = 5, //16384,
                ncities = 15,
                nterrs = 5,
                fontsizes = new FontSizes()
                {
                    region = 40,
                    city = 25,
                    town = 20
                }

            };

            MapGenerator mapGenerator = new MapGenerator
            {
                Params = param
            };
            mapGenerator.DoMap();
        }





    }





}
