#define _CRT_SECURE_NO_WARNINGS 1

#include <iostream>
#include <sstream>

#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "lbfgs.h"

/*
Voronoi diagram:
many points in a square --> 
for every point, which region of the square is closer to this point than to any other point?
point -> its region (voronoi cell)

Polygon clipping:
for one point, we don't know its final cell
we start with the whole square
- we compare Pi with another point Pj, there is a line in the middle btwn Pi and Pj (bisector)
--> separate the plan in 2 sides (closer to Pi / closer to Pj)

Sutherland-Hodgman
algorithm to cut a polygon by a line (half-space)

*/




double sqr(double x) { return x * x; };

class Vector {
public:
    explicit Vector(double x = 0, double y = 0) {
        data[0] = x;
        data[1] = y;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[2];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a * b[0], a * b[1]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0] * b, a[1] * b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1];
}
double cross(const Vector& a, const Vector& b){
    return a[0]*b[1] - a[1]*b[0];
}


class Polygon {
public:

    double area() {
        if (vertices.size() < 3) return 0;
        // TODO Lab 2
        // Compute the area of the polygon

        double A = 0.0;

        for (int i = 0; i < vertices.size(); i++){
            Vector vi = vertices[i];
            //i forgot the last segment so i need to add one step with vertices[-1]
            Vector vip1 = vertices[(i+1) % vertices.size()]; //don't forget last segment vertices[-1] -- vertices[0]
            A+= vi[0]*vip1[1] - vip1[0]*vi[1];
        }

        A = 0.5 * std::abs(A);
        return A;
    }

    Vector centroid() {
        if (vertices.size() < 3) return Vector(0, 0);
        // TODO Lab 3
        // Compute the centroid of the polygon

        return Vector(-111,-111);
    }

    double integral_square_distance(const Vector& Pi) {
        if (vertices.size() < 3) return 0;

        // TODO Lab 2
        // Compute the integral of ||x-Pi||^2 over the polygon

        double integral = 0.0;

        //decompose the polygon into triangles
        for (int i = 0; i < vertices.size() -1; i++){
            Vector triangle[3] = {vertices[0], vertices[i], vertices[i+1]};
            
            double T = 0.5 * std::abs(cross(triangle[1] - triangle[0], triangle[2] - triangle[0]));
            //compute sum for this triangle
            
            double sum = 0.0;
            for (int k = 0; k < 3; k++){
                for (int l = 0; l < 3; l++){
                    sum += dot(triangle[k] - Pi, triangle[l] - Pi);
                }

            }
            integral += (T/6.0) * sum;
        }


        // sum all integrals of the triangles 
        return integral;    
    }

    std::vector<Vector> vertices;
};


void save_frame(const std::vector<Polygon>& cells, std::string filename, int frameid = 0) {
    constexpr int W = 800, H = 800;
    constexpr double edge_width = 2.0;
    constexpr double edge_width2 = edge_width * edge_width;

    std::vector<unsigned char> inside(W * H, 0), edge(W * H, 0);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)cells.size(); ++i) {
        const auto& V = cells[i].vertices;
        const int n = (int)V.size();
        if (n < 3) continue;

        std::vector<double> xs(n), ys(n);
        double xmin = 1e30, ymin = 1e30, xmax = -1e30, ymax = -1e30;
        for (int j = 0; j < n; ++j) {
            xs[j] = V[j][0] * W;
            ys[j] = V[j][1] * H;
            xmin = std::min(xmin, xs[j]);
            ymin = std::min(ymin, ys[j]);
            xmax = std::max(xmax, xs[j]);
            ymax = std::max(ymax, ys[j]);
        }

        int x0 = std::max(0, (int)std::floor(xmin - edge_width));
        int y0 = std::max(0, (int)std::floor(ymin - edge_width));
        int x1 = std::min(W - 1, (int)std::ceil(xmax + edge_width));
        int y1 = std::min(H - 1, (int)std::ceil(ymax + edge_width));
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                const double px = x + 0.5, py = y + 0.5;

                int prev_sign = 0;
                bool isInside = true;
                bool isEdge = false;

                for (int j = 0; j < n; ++j) {
                    int k = (j + 1) % n;

                    double ax = xs[j], ay = ys[j];
                    double bx = xs[k], by = ys[k];
                    double dx = bx - ax, dy = by - ay;
                    double qx = px - ax, qy = py - ay;

                    double det = qx * dy - qy * dx;
                    int s = (det > 1e-12) - (det < -1e-12);

                    if (s != 0) {
                        if (prev_sign != 0 && s != prev_sign) {
                            isInside = false;
                            break;
                        }
                        prev_sign = s;
                    }

                    double len2 = dx * dx + dy * dy;
                    double dot = qx * dx + qy * dy;
                    if (dot >= 0.0 && dot <= len2 && det * det <= edge_width2 * len2)
                        isEdge = true;
                }

                if (isInside) {
                    int id = (H - 1 - y) * W + x;
                    inside[id] = 1;
                    if (isEdge) edge[id] = 1;
                }
            }
        }
    }

    std::vector<unsigned char> image(W * H * 3, 255);

#pragma omp parallel for
    for (int i = 0; i < W * H; ++i) {
        if (edge[i]) {
            image[3 * i + 0] = 0;
            image[3 * i + 1] = 0;
            image[3 * i + 2] = 0;
        }
        else if (inside[i]) {
            image[3 * i + 0] = 0;
            image[3 * i + 1] = 0;
            image[3 * i + 2] = 255;
        }
    }

    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, image.data(), W * 3);
}


class VoronoiDiagram {

public:

    VoronoiDiagram() {
    };


    void compute() {

        // TODO Lab 1 (Voronoi)
        // For all sites Pi (in parallel) :
        //      Start with a unit square
        //      For all other sites Pj (optionally, only k nearest neighbors) :
        //          Clip it with bisector of [Pi,Pj]

        //size weights
        weights.resize(points.size());

        //size of cells 
        cells.resize(points.size());

    #pragma omp parallel for schedule(dynamic)
        for (int i= 0; i < points.size(); i++){
            //current point Pi
            Vector Pi = points[i];

            //unit square [0,1] * [0,1]
            Polygon cell;
            cell.vertices.push_back(Vector(0,0));
            cell.vertices.push_back(Vector(1,0));
            cell.vertices.push_back(Vector(1,1));
            cell.vertices.push_back(Vector(0,1));

            //lab1 : weights = 0
            double wi = weights[i];

            //clip cell Pi with all other bisectors [Pi, Pj]
            for (int j = 0; j < points.size(); j++){
                //i different from j
                if (i == j){
                    continue;
                }

                Vector Pj = points[j];
                double wj = weights[j];

                //keep only side closer to Pi than to Pj
                cell = clip_by_bisector(cell, Pi, Pj, wi, wj);

                //if cell empty no need to continue
                if (cell.vertices.empty()){
                    break;
                }
            }

            //store voronoi cell of Pi
            cells[i] = cell; 

        }




        //      (Lab 3, fluids) : also clip it by a disk of radius sqrt(w_i - w_air) centered at Pi
    }


    static Polygon clip_by_edge(const Polygon& V, const Vector& u, const Vector& v) {

        // TODO Lab 3 (fluids)
        // Clip a polygon by an edge defined by vertices u and v
        // Will be used to clip a polygon (a cell) by all the edges of a (discretized) disk

        Polygon result;

        return result;
    }

    static Polygon clip_by_bisector(const Polygon& V, const Vector& P0, const Vector& Pi, double w0, double wi) {

        // TODO Lab 1 (Voronoi) : in Lab 1, we assume w0 = w1 = 0
        // Clip a polygon by the bisector of the segment defined by P0 (the current site of the Voronoi cell being computed) and Pi (another site)
        
        /*
        we want to build the cell of P0
        For a point Pi, we keep only points x such that : distance(x, P0) <= distance(x, Pi)
        */
        
        Polygon outPolygon;

        //wi = w0 = 0 in lab 1
        //weights associated to points

        if (V.vertices.empty()){
            return outPolygon;
        }

        for (int i = 0; i < V.vertices.size(); i++ ){
            //for each vertex of the subject polygon

            //test the subject polygon edge with vertices (i-1, i)
            Vector curVertex = V.vertices[i];
            Vector prevVertex = V.vertices[(i>0)?(i-1):(V.vertices.size()-1)];

            //check if curVertex is inside 
            // ||curVertex - P0|^2 - w0 <= ||curVErtex - Pi||^2 - wi 
            double curDistance0 = (curVertex - P0).norm2() - w0;
            double curDistancei = (curVertex - Pi).norm2() - wi; 
            bool curInside = curDistance0 <= curDistancei;

            //same for prevVertex
            double prevDistance0 = (prevVertex - P0).norm2() - w0;
            double prevDistancei = (prevVertex - Pi).norm2() - wi; 
            bool prevInside = prevDistance0 <= prevDistancei;

            //find center of P0, Pi
            Vector M = (P0 + Pi) / 2.0;
            Vector M_prime = M + (Pi - P0) * ((w0 - wi) / (2 * (P0 - Pi).norm2()));

            //find direction of mediator
            Vector normal = Pi - P0;

            //find intersection btwn [prevVertex, curVertex] and bisector of [P0, Pi]
            // P = A + t * (B-A)
            Vector intersection;

            double denominator = dot(curVertex - prevVertex, normal);

            if (std::abs(denominator) > 1e-10){
                double t = dot(M_prime - prevVertex, normal) / denominator;
                intersection = prevVertex + t *(curVertex - prevVertex);
            }

            if (curInside){
                if(!prevInside){
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            }
            else if (prevInside){
                outPolygon.vertices.push_back(intersection);
            }

        }

        return outPolygon;


        



       



        // TODO Lab 2 (Semi-Discrete Optimal Transport) : extend to Laguerre cells, i.e., w0 != w1

    }

    std::vector<Vector> points;    // Lab 1 (Voronoi) : the sites to consider

    std::vector<double> weights;   // Lab 2 (OT) : the weight associated to each site (the Laguerre weight, i.e. the dual optimal transport variables to be optimized)
    
    std::vector<Polygon> cells;   // Lab 1 : the polygons representing each individual cell

};


// Lab 2 
class OptimalTransport {

public:
    OptimalTransport() {};

    void optimize();

    VoronoiDiagram vor;
};


// Labs 2 and 3
static lbfgsfloatval_t evaluate(
    void* instance,
    const lbfgsfloatval_t* x,
    lbfgsfloatval_t* g,
    const int n,
    const lbfgsfloatval_t step
)
{
    OptimalTransport* ot = (OptimalTransport*)(instance);

    // first compute the Voronoi diagram at the current optimization step
    memcpy(&ot->vor.weights[0], x, n * sizeof(x[0]));
    ot->vor.compute();
  
   
    // Lab 2 (Optimal transport) : compute the function to be minimized (fx) and its gradient (g[i], i=0..n-1)
    // Lab 3 (fluid) : adapt these functions to support partial optimal transport (now "n" has been increased by 1 to account for the air variable)
    
    lbfgsfloatval_t fx = 0.0;
    // g[i] = ...
    // fx = ...

    
    for (int i = 0; i < ot->vor.cells.size(); i++){
        double lambda = 1.0 / n;
        double area = ot->vor.cells[i].area();
        double integral = ot->vor.cells[i].integral_square_distance(ot->vor.points[i]);
        double wi = ot->vor.weights[i];

        double gi = integral - wi * area + lambda * wi;
        fx -= gi; //lbfgs minimizes --> we want to maximize
        g[i] = area - lambda; //gradient
    }
    

    return fx;
}

// Labs 2 and 3 : you may use this function to print debugging info.
static int progress(
    void* instance, const lbfgsfloatval_t* x, const lbfgsfloatval_t* g, const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
    int n, int k, int ls) {
    printf("Iteration %d:\n", k);
    printf("  fx = %f\n", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}


// Lab 2
void OptimalTransport::optimize() {

    lbfgsfloatval_t fx;
    std::vector<double> weights(vor.weights);

    lbfgs_parameter_t param;
    // Initialize the parameters for the L-BFGS optimization.
    lbfgs_parameter_init(&param);

    // run the LBFGS optimizer
    int ret = lbfgs(weights.size(), &weights[0], &fx, evaluate, progress, (void*)this, &param);

    // copy the result back to the voronoi structure
    vor.weights = weights;

    // finally recompute the Voronoi diagram with the final optimized weights
    vor.compute();
}


// Lab 3 (fluids)
class Fluid {
public:
    Fluid(int N_particles = 1000) : N_particles(N_particles) {
    }

    // Lab 3 : advance the simulation dt in time
    void time_step(double dt) {

        double epsilon2 = 0.004 * 0.004;
        Vector g(0, -9.81);
        double m_i = 200;

        // TODO Lab 3 : 
        // Compute semi-discrete partial optimal transport
        // for all particles, add gravity and spring force towards cell centroid, integrate acceleration->velocity and velocity->position
    }

    // just run the full simulation
    void run_simulation() {
        double dt = 0.002;
        for (int i = 0; i < 1000; i++) {
            time_step(dt);
            save_frame(ot.vor.cells, "test", i);
        }
    }

    int N_particles;

    OptimalTransport ot;
    std::vector<Vector> particles;  // the position of all particles
    std::vector<Vector> velocities; // the velocities of all particles
    double fluid_volume; // you decide the fraction of the unit square occupied by the fluid
};

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon>& polygons, std::string filename, const std::vector<Vector>* points = NULL, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }

    if (points) {
        fprintf(f, "<g>\n");
        for (int i = 0; i < points->size(); i++) {
            fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"3\" />\n", (*points)[i][0] * 1000., 1000. - (*points)[i][1] * 1000);
        }
        fprintf(f, "</g>\n");

    }

    fprintf(f, "</svg>\n");
    fclose(f);
}



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
    /*
    save_svg(v.cells, "voronoi.svg", &v.points, "none");
    */

    OptimalTransport ot;

    ot.vor.points.push_back(Vector(0.1, 0.2));
    ot.vor.points.push_back(Vector(0.6, 0.4));
    ot.vor.points.push_back(Vector(0.5, 0.7));
    ot.vor.points.push_back(Vector(0.2, 0.5));

    ot.vor.weights.resize(ot.vor.points.size(), 0.0);

    ot.optimize();
    save_frame(ot.vor.cells, "lab2");
    save_svg(ot.vor.cells, "lab2.svg", &ot.vor.points, "white");
    return 0;
}


