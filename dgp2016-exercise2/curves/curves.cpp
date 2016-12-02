#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>

using namespace OpenGP;

struct MainWindow : public TrackballWindow {
    PointsRenderer render_points = PointsRenderer ();
    SegmentsRenderer render_segments = SegmentsRenderer ();

    MatMxN points;
    MatMxN points_3d_render;
    int num_points;
    SegmentsRenderer::Segments segments;

    // ============================================================================
    // Exercise 2 : fill the 2 functions below (see PDF for instructions)
    // To test your implementation, use the S key for laplacian smoothing and the
    // C key for the osculating circle.
    // Hint : try to play with epsilon
    // ============================================================================
    // time step for smoothing
    double epsilon = 0.05;

    // **** Utilities **** //
    //Computes the pow2
    float pow2(float f){
        return f*f;
    }
    //Computes the distance between 2 vec2
    float distance(Vec2 a, Vec2 b){
        return std::sqrt( pow2(a(0) - b(0)) + pow2(a(1) - b(1)) );
    }
    //Compute the length of the curve
    float computeCurveLength(){
        float length = 0;
        for (int i = 0; i<=num_points-2; ++i){
            length += distance(points.col(i), points.col(i+1));
        }
        length += distance(points.col(num_points-1), points.col(0));
        return length;
    }
    //Computes the center of the Circumcircle given the 3 points of the triangle
    Vec2 computeCenter(Vec2 a, Vec2 b, Vec2 c){
        float dA, dB, dC, aux1, aux2, div;
        dA = pow2(a(0)) + pow2(a(1));
        dB = pow2(b(0)) + pow2(b(1));
        dC = pow2(c(0)) + pow2(c(1));
        aux1 = dA*(c(1) - b(1)) + dB*(a(1) - c(1)) + dC*(b(1) - a(1));
        aux2 = dA*(c(0) - b(0)) + dB*(a(0) - c(0)) + dC*(b(0) - a(0));
        aux2 = -aux2;
        div = 2*(a(0)*(c(1)-b(1)) + b(0)*(a(1)-c(1)) + c(0)*(b(1)-a(1)));
        Vec2 center = Vec2(aux1/div, aux2/div);
        return center;
    }
    //Compute a "center" of the curve by computing the circum center of 3 distanced points
    Vec2 computeCurveCenter(){
        int third = ceil(num_points/3.0);
        Vec2 center = computeCenter(points.col(0), points.col(third), points.col(2*third));
        return center;
    }
    //Shift all the points of the points matrix by the input value
    void shiftPoints(Vec2 shift){
        for(int i = 0; i<points.cols(); ++i){
            points.col(i) += shift;
        }
    }
    // Creates a copy of the points matrix with 1 cyclic repetition at each end.
    // I.e. [a, b, c, d, e] ===> [e, a, b, c, d, e, a]
    MatMxN createPointsCopy(){
        int nb = points.cols();
        MatMxN copy = MatMxN::Zero(2, nb + 2);

        for (int i = 1; i <= nb; ++i){
            copy.col(i) = points.col(i-1);
        }
        copy.col(0) = copy.col(nb);
        copy.col(nb+1) = copy.col(1);
        return copy;
    }



    // **** Curve Smoothing Functions **** //
    void laplacianSmoothing() {
        float old_length = computeCurveLength();
        Vec2 oldCenter = computeCurveCenter();

        MatMxN pointsCopy = createPointsCopy();
        Vec2 previous, current, next;

        for(int i = 0; i<num_points; ++i){
            previous = pointsCopy.col(i);
            current = pointsCopy.col(i+1);
            next = pointsCopy.col(i+2);
            points.col(i) = (1.0 - epsilon)*current + epsilon * ((next + previous)/2.0);
        }
        points = (old_length /  computeCurveLength()) * points;
        shiftPoints(oldCenter - computeCurveCenter());
    }

    void osculatingCircle() {
        float old_length = computeCurveLength();
        Vec2 oldCenter = computeCurveCenter();

        MatMxN pointsCopy = createPointsCopy();
        Vec2 previous, current, next;

        for(int i = 0; i<num_points; ++i){
            previous = pointsCopy.col(i);
            current = pointsCopy.col(i+1);
            next = pointsCopy.col(i+2);
            Vec2 center = computeCenter(previous, current, next);
            points.col(i) = current + epsilon*((center - current)/pow2(distance(center, current)));
        }

        points = (old_length / computeCurveLength()) * points;
        shiftPoints(oldCenter - computeCurveCenter());
    }


    // ============================================================================
    // END OF Exercise 2 (do not thouch the rest of the code)
    // ============================================================================

    void generateRandomizedClosedPolyline() {
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0., 5*3e-2);

        Vec2 center(3e-2, 2e-3);
        const double radius = 0.3;

        points = MatMxN::Zero(2, num_points);
        for (int i = 0; i < num_points; ++i)
        {
            double frac = static_cast<double>(i) / static_cast<double>(num_points);
            points(0, i) = center(0) + radius * cos (2. * M_PI * frac) + distribution(generator);
            points(1, i) = center(1) + radius * sin (2. * M_PI * frac) + distribution(generator);
        }
    }

    void render () {

        // Prepare the render points
        points_3d_render = MatMxN::Zero(3, points.cols());
        points_3d_render.block(0, 0, 2, points.cols()) = points;

        // Rebuild the segments
        segments.clear();
        for (int i = 0; i < points_3d_render.cols(); ++i) {
            segments.push_back({ points_3d_render.col(i), points_3d_render.col((i+1) % points_3d_render.cols()) });
        }
        render_points.init_data(points_3d_render);
        render_segments.init_data(segments);
    }

    MainWindow(int argc, char** argv) : TrackballWindow("2D Viewer", 640, 480) {
        num_points = 50;
        generateRandomizedClosedPolyline();

        this->scene.add(render_points);
        this->scene.add(render_segments);

        render();
    }

    bool key_callback(int key, int scancode, int action, int mods) override {
        TrackballWindow::key_callback(key, scancode, action, mods);
        if (key == GLFW_KEY_S && action == GLFW_RELEASE)
        {
            laplacianSmoothing();
        }
        else if (key == GLFW_KEY_C && action == GLFW_RELEASE)
        {
            osculatingCircle();
        }
        render();
        return true;
    }
};


int main(int argc, char** argv)
{
    MainWindow window(argc, argv);
    return window.run();


}
