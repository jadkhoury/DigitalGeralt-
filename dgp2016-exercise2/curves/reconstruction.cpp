#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>

using namespace OpenGP;

struct MainWindow : public TrackballWindow {
    PointsRenderer render_points = PointsRenderer();
    SegmentsRenderer render_segments = SegmentsRenderer();

    double epsilon;
    MatMxN points;
    MatMxN points_3d_render;
    Vec2 tangent;
    SegmentsRenderer::Segments segments;


    // ============================================================================
    // Exercise 3 : fill the function below (see PDF for instructions)
    // To test your implementation, use the SPACE key
    // ============================================================================
    float pow2(float f) {return f*f;}
    float pow3(float f) {return f*f*f;}

    void reconstructCurveStep() {
        //every time you press space, a new point of the curve is added, so the number of points increases
        points.conservativeResize(2, points.cols() + 1);
        int n = points.cols() - 1;
        Vec2 t; float phi;

        //phi = n*epsilon; //first curve
        phi = pow2(n*epsilon)/2.0; //second curve
        //phi = pow3(n*epsilon)/3.0 - 2.19*n*epsilon; //third curve
        if(n == 1) t = tangent;
        else t = Vec2(cos(phi), sin(phi));
        Vec2 u = points.col(n-1) + epsilon * t;
        points.col(n) = u;
    }
    // ============================================================================
    // END OF Exercise 2 (do not thouch the rest of the code except to uncomment
    // the part where you can change the function)
    // ============================================================================

    void render() {
        // Prepare the render points
        points_3d_render = MatMxN::Zero(3, points.cols());
        points_3d_render.block(0, 0, 2, points.cols()) = points * 1e-1;

        // Rebuild the segments
        segments.clear();

        for (int i = 1; i < points_3d_render.cols (); ++i) {
            segments.push_back({ points_3d_render.col(i), points_3d_render.col(i-1)});
        }

        render_points.init_data(points_3d_render);
        render_segments.init_data(segments);
    }

    MainWindow(int argc, char** argv) : TrackballWindow("2D Viewer", 640, 480) {
        points = MatMxN(2, 1);

        // First, second, third curve
        //points.col(0) = Vec2 (0, -1.0); // first
        points.col (0) = Vec2 (0, 0.0); // second
        // points.col (0) = Vec2 (0, 0.0); // third

        tangent = Vec2 (1.0, 0.0);

        epsilon = 0.1;

        this->scene.add(render_points);
        this->scene.add(render_segments);

        render ();
    }

    bool key_callback(int key, int scancode, int action, int mods) override {
        TrackballWindow::key_callback(key, scancode, action, mods);
        if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
        {
            reconstructCurveStep();
        }

        render ();
        return true;
    }
};


int main(int argc, char** argv)
{
    MainWindow window(argc, argv);
    return window.run();
}
