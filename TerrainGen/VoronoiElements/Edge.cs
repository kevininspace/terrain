namespace TerrainGen
{
    public class Edge
    {
        public double a = 0, b = 0, c = 0;
        public Site[] ep;
        public Site[] reg;
        public int edgenbr;

        public Edge()
        {
            ep = new Site[2];
            reg = new Site[2];
        }
    }
}