namespace TerrainGen
{
    public class Halfedge
    {
        public Halfedge ELleft, ELright;
        public Edge ELedge;
        public bool deleted;
        public int ELpm;
        public Site vertex;
        public double ystar;
        public Halfedge PQnext;

        public Halfedge()
        {
            PQnext = null;
        }
    }
}