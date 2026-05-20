Start Project 2 : 

## Lab 1 : 13/05/2026


int main() {

    //show a polygon
    Polygon p;
    p.vertices.push_back(Vector(0.1, 0.2));
    p.vertices.push_back(Vector(0.6, 0.4));
    p.vertices.push_back(Vector(0.5, 0.7));
    p.vertices.push_back(Vector(0.2, 0.5));

    std::vector<Polygon> s;
    s.push_back(p);

    /*
    save_frame(s, "toto");
    save_svg(s, "toto.svg");
    */

    //VORONOI DIARGAM
    VoronoiDiagram v;
    v.points.push_back(Vector(0.1, 0.2));
    v.points.push_back(Vector(0.6, 0.4));
    v.points.push_back(Vector(0.5, 0.7));
    v.points.push_back(Vector(0.2, 0.5));

    v.compute();

    save_svg(v.cells, "voronoi.svg", &v.points, "none");

    return 0;
}


time : 
real    0m0.043s
user    0m0.036s
sys     0m0.004s

## Lab 2 : 20/05/2026

