//=================================================================================
//   Code framework for the lecture "Digital 3D Geometry Processing"
//   Copyright Gaspard Zoss (C) 2016 Computer Graphics and Geometry Laboratory, EPFL
//----------------------------------------------------------------------------------
#include <nanogui/opengl.h>
#include <nanogui/glutil.h>
#include <nanogui/screen.h>
#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/label.h>
#include <nanogui/button.h>
#include <nanogui/tabwidget.h>
#include <surface_mesh/Surface_Mesh.h>

#include "trackball.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
#if defined(_WIN32)
#  pragma warning(push)
#  pragma warning(disable: 4457 4456 4005 4312)
#endif

#if defined(_WIN32)
#  pragma warning(pop)
#endif
#if defined(_WIN32)
#  if defined(APIENTRY)
#    undef APIENTRY
#  endif
#  include <windows.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::to_string;
using namespace surface_mesh;

class Viewer : public nanogui::Screen {
public:
    Viewer(surface_mesh::Surface_mesh *mesh) :
        nanogui::Screen(Eigen::Vector2i(1024, 768), "DGP Viewer") {
        using namespace nanogui;

        // Setting up the Window and the GUI
        window = new Window(this, "Controls");
        window->setPosition(Vector2i(15, 15));
        window->setLayout(new GroupLayout());

        Button *b = new Button(window, "wireframe");
        b->setFlags(Button::ToggleButton);
        b->setChangeCallback([this](bool wireframe) {
            this->wireframe =! this->wireframe;
        });

        b = new Button(window, "valence display");
        b->setFlags(Button::ToggleButton);
        b->setChangeCallback([this](bool valence_coloring) {
            this->valence_coloring =! this->valence_coloring;
        });

        b = new Button(window, "display normals");
        b->setFlags(Button::ToggleButton);
        b->setChangeCallback([this](bool normals) {
            this->normals =! this->normals;
        });

        TabWidget* tabWidget = window->add<TabWidget>();
        Widget* layer = tabWidget->createTab("normals");

        layer->setLayout(new GroupLayout());
        b = layer->add<Button>();
        b->setCaption("constant weights");
        b->setFlags(Button::RadioButton);
        b->setPushed(true);
        b->setChangeCallback([this](int normals_computation) {
            this->normals_computation = 0;
        });
        b = layer->add<Button>();
        b->setCaption("area weights");
        b->setFlags(Button::RadioButton);
        b->setChangeCallback([this](int normals_computation) {
            this->normals_computation = 1;
        });
        b = layer->add<Button>();
        b->setCaption("angle weights");
        b->setFlags(Button::RadioButton);
        b->setChangeCallback([this](int normals_computation) {
            this->normals_computation = 2;
        });

        performLayout();
        // end of window and GUI setup

        // trackball initial parameters
        trackball_matrix = Matrix4f::Identity();
        old_trackball_matrix = Matrix4f::Identity();
        view_matrix = lookAt(Vector3f(0.0, 0.0, 3.0),  // camera position
                             Vector3f(0.0, 0.0, 0.0),  // looking where ?
                             Vector3f(0.0, 1.0, 0.0)); // up vector
        project_matrix = perspective(45.0, (float)width()/height(), 0.01, 100.0);

        // scale and move the mesh
        model_matrix << 10.0f, 0.0f, 0.0f, 0.2f,
                        0.0f, 10.0f, 0.0f, -1.0f,
                        0.0f, 0.0f, 10.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 1.0f;

        // Shaders
        mShader.init(
            "a_simple_shader",

            /* Vertex shader */
            "#version 330\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"
            "uniform int normal_selector;\n"

            "in vec3 position;\n"
            "in vec3 vcolor;\n"
            "in vec3 n_cste_weights;\n"
            "in vec3 n_area_weights;\n"
            "in vec3 n_angle_weights;\n"

            "out vec3 fcolor;\n"
            "out vec3 fnormal;\n"
            "out vec3 view_dir;\n"
            "out vec3 light_dir;\n"

            "void main() {\n"
            "    vec4 vpoint_mv = MV * vec4(position, 1.0);\n"
            "    gl_Position = P * vpoint_mv;\n"
            "    fcolor = vcolor;\n"
            "    if (normal_selector == 0) {\n"
            "       fnormal = mat3(transpose(inverse(MV))) * n_cste_weights;\n"
            "    } else if(normal_selector == 1) {\n"
            "       fnormal = mat3(transpose(inverse(MV))) * n_area_weights;\n"
            "    } else {\n"
            "       fnormal = mat3(transpose(inverse(MV))) * n_angle_weights;\n"
            "    }\n"
            "    light_dir = vec3(0.0, 3.0, 3.0) - vpoint_mv.xyz;\n"
            "    view_dir = -vpoint_mv.xyz;\n"
            "}",

            /* Fragment shader */
            "#version 330\n"
            "uniform vec3 intensity;\n"
            "uniform int valenceColor;\n"

            "in vec3 fcolor;\n"
            "in vec3 fnormal;\n"
            "in vec3 view_dir;\n"
            "in vec3 light_dir;\n"

            "out vec4 color;\n"

            "void main() {\n"
            "    vec3 c = vec3(0.0);\n"
            "    if (valenceColor == 1) {\n"
            "        c = fcolor;\n"
            "    }\n"
            "    else {\n"
            "        c += vec3(1.0)*vec3(0.18, 0.1, 0.1);\n"
            "        vec3 n = normalize(fnormal);\n"
            "        vec3 v = normalize(view_dir);\n"
            "        vec3 l = normalize(light_dir);\n"
            "        float lambert = dot(n,l);\n"
            "        if(lambert > 0.0) {\n"
            "            c += vec3(1.0)*vec3(0.9, 0.5, 0.5)*lambert;\n"
            "            vec3 v = normalize(view_dir);\n"
            "            vec3 r = reflect(-l,n);\n"
            "            c += vec3(1.0)*vec3(0.8, 0.8, 0.8)*pow(max(dot(r,v), 0.0), 90.0);\n"
            "        }\n"
            "        c *= intensity;\n"
            "    }\n"
            "    if (intensity == vec3(0.0)) {\n"
            "        c = intensity;\n"
            "    }\n"
            "    color = vec4(c, 1.0);\n"
            "}"
        );

        mShaderNormals.init(
            "normal_shader",
            /* Vertex shader */
            "#version 330\n\n"
            "in vec3 position;\n"
            "in vec3 n_cste_weights;\n"
            "in vec3 n_area_weights;\n"
            "in vec3 n_angle_weights;\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"
            "uniform int normal_selector;\n"
            "out VS_OUT {\n"
            "    mat3 normal_mat;\n"
            "    vec3 normal;\n"
            "} vs_out;\n"
            "void main() {\n"
            "  gl_Position = vec4(position, 1.0);\n"
            "    if (normal_selector == 0) {\n"
            "       vs_out.normal = n_cste_weights;\n"
            "    } else if(normal_selector == 1) {\n"
            "       vs_out.normal = n_area_weights;\n"
            "    } else {\n"
            "       vs_out.normal = n_angle_weights;\n"
            "    }\n"
            "    vs_out.normal_mat = mat3(transpose(inverse(MV)));\n"
            "}",
            /* Fragment shader */
            "#version 330\n\n"
            "out vec4 frag_color;\n"
            "void main() {\n"
            "   frag_color = vec4(0.0, 1.0, 0.0, 1.0);\n"
            "}",
            /* Geometry shader */
            "#version 330\n\n"
            "layout (triangles) in;\n"
            "layout (line_strip, max_vertices = 6) out;\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"
            "in VS_OUT {\n"
            "    mat3 normal_mat;\n"
            "    vec3 normal;\n"
            "} gs_in[];\n"
            "void createline(int index) {\n"
            "   gl_Position = P * MV * gl_in[index].gl_Position;\n"
            "   EmitVertex();\n"
            "   gl_Position = P * (MV * gl_in[index].gl_Position "
            "+ vec4(normalize(gs_in[index].normal_mat * "
            "gs_in[index].normal), 0.0f)*0.01f);\n"
            "   EmitVertex();\n"
            "   EndPrimitive();\n"
            "}\n"
            "void main() {\n"
            "   createline(0);\n"
            "   createline(1);\n"
            "   createline(2);\n"
            "}"
        );

        cout << "Mesh loaded :" << endl;
        n_vertices = mesh->n_vertices();
        cout << "# of vertices : " << n_vertices << endl;
        n_faces = mesh->n_faces();
        cout << "# of faces : " << n_faces << endl;
        n_edges = mesh->n_edges();
        cout << "# of edges : " << n_edges << endl;


        Surface_mesh::Vertex_property<Point> vertex_color =
                colorsFromValence(mesh);
        Surface_mesh::Vertex_property<Point> v_cste_weights_n =
                mesh->vertex_property<Point>("v:cste_weights_n");
        Surface_mesh::Vertex_property<Point> v_area_weights_n =
                mesh->vertex_property<Point>("v:area_weight_n");
        Surface_mesh::Vertex_property<Point> v_angle_weights_n =
                mesh->vertex_property<Point>("v:angle_weight_n");


        int j = 0;
        MatrixXf mesh_points(3, n_vertices);
        MatrixXu indices(3, n_faces);

        for(auto f: mesh->faces()) {
            vector<float> vv(3);
            int k = 0;
            for (auto v: mesh->vertices(f)) {
                vv[k] = v.idx();
                ++k;
            }
            indices.col(j) << vv[0], vv[1], vv[2];
            ++j;
        }

        // Create big matrices to send the data to the GPU with the required
        // format
        MatrixXf color_attrib(3, n_vertices);
        MatrixXf normal_cste_weights_attrib(3, n_vertices);
        MatrixXf normal_area_weights_attrib(3, n_vertices);
        MatrixXf normal_angle_weights_attrib(3, n_vertices);

        j = 0;
        for (auto v: mesh->vertices()) {
            mesh_points.col(j) << mesh->position(v).x,
                                  mesh->position(v).y,
                                  mesh->position(v).z;

            color_attrib.col(j) << vertex_color[v].x,
                                   vertex_color[v].y,
                                   vertex_color[v].z;

            normal_cste_weights_attrib.col(j) << v_cste_weights_n[v].x,
                                                 v_cste_weights_n[v].y,
                                                 v_cste_weights_n[v].z;

            normal_area_weights_attrib.col(j) << v_area_weights_n[v].x,
                                                 v_area_weights_n[v].y,
                                                 v_area_weights_n[v].z;

            normal_angle_weights_attrib.col(j) << v_angle_weights_n[v].x,
                                                  v_angle_weights_n[v].y,
                                                  v_angle_weights_n[v].z;
            ++j;
        }

        mShader.bind();
        mShader.uploadIndices(indices);
        mShader.uploadAttrib("position", mesh_points);
        mShader.uploadAttrib("vcolor", color_attrib);
        mShader.uploadAttrib("n_cste_weights", normal_cste_weights_attrib);
        mShader.uploadAttrib("n_area_weights", normal_area_weights_attrib);
        mShader.uploadAttrib("n_angle_weights", normal_angle_weights_attrib);
        mShader.setUniform("valenceColor", (valence_coloring)? 1 : 0);
        mShader.setUniform("intensity", Vector3f(0.98, 0.59, 0.04));

        mShaderNormals.bind();
        mShaderNormals.uploadIndices(indices);
        mShaderNormals.uploadAttrib("position", mesh_points);
        mShaderNormals.uploadAttrib("n_cste_weights", normal_cste_weights_attrib);
        mShaderNormals.uploadAttrib("n_area_weights", normal_area_weights_attrib);
        mShaderNormals.uploadAttrib("n_angle_weights", normal_angle_weights_attrib);
    }

    Surface_mesh::Vertex_property<Point> colorsFromValence(Surface_mesh* mesh) {

        Surface_mesh::Vertex_property<unsigned int> vertex_valence =
                mesh->vertex_property<unsigned int>("v:valence");

        Point default_color(1.0, 1.0, 1.0);
        Surface_mesh::Vertex_property<Point> vertex_color =
                mesh->vertex_property<Point>("v:color", default_color);

        // put all valence values into one array
        vector<unsigned int> valence_values = vertex_valence.vector();

        unsigned int n = valence_values.size()-1;
        unsigned int i = valence_values.size() / 100.0;
        std::sort(valence_values.begin(), valence_values.end());
        double min_valence = valence_values[i];
        double max_valence = valence_values[n-1-i];

        // define uniform color intervalls [v0,v1,v2,v3,v4]
        double v0, v1, v2, v3, v4;
        v0 = min_valence + 0.0/4.0 * (max_valence - min_valence);
        v1 = min_valence + 1.0/4.0 * (max_valence - min_valence);
        v2 = min_valence + 2.0/4.0 * (max_valence - min_valence);
        v3 = min_valence + 3.0/4.0 * (max_valence - min_valence);
        v4 = min_valence + 4.0/4.0 * (max_valence - min_valence);

        // map valences to colors
        for (auto v: mesh->vertices())
        {
            double valen = vertex_valence[v];
            Vector<float,3> col(1.0,1.0,1.0);

            if (valen < v0) {
                col = Vector<float,3>(0, 0, 1.0);
            } else if (valen > v4) {
                col = Vector<float,3>(1.0, 0, 0);
            } else if (valen <= v2) {
                if (valen <= v1) { // [v0, v1]
                    double u = (valen - v0) / (v1 - v0);
                    col = Vector<float,3>(0, u, 1.0);
                } else { // ]v1, v2]
                    double u = (valen - v1) / (v2 - v1);
                    col = Vector<float,3>(0, 1.0, 1.0 -u);
                }
            } else {
                if (valen <= v3) { // ]v2, v3]
                    double u = (valen - v2) / (v3 - v2);
                    col = Vector<float,3>(u, 1.0, 0);
                } else { // ]v3, v4]
                    double u = (valen - v3) / (v4 - v3);
                    col = Vector<float,3>(1.0, 1.0-u, 0);
                }
            }
            vertex_color[v] = col;
        }
        return vertex_color;
    }

    ~Viewer() {
        mShader.free();
        mShaderNormals.free();
    }

    virtual bool keyboardEvent(int key, int scancode, int action, int modifiers) {
        if (Screen::keyboardEvent(key, scancode, action, modifiers)) {
            return true;
        }
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
            setVisible(false);
            return true;
        }
        if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS && !left_control_pan) {
            left_control_pan = true;
            return true;
        }
        if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE && left_control_pan) {
            left_control_pan = false;
            return true;
        }
        return false;
    }

    virtual void draw(NVGcontext *ctx) {
        /* Draw the user interface */
        Screen::draw(ctx);
    }

    Eigen::Matrix4f perspective(double fovy, double aspect,
                                double zNear, double zFar) {
       assert(aspect > 0);
       assert(zFar > zNear);
       double radf = fovy*M_PI/180.0;
       double tanHalfFovy = tan(radf / 2.0);
       Matrix4f res = Matrix4f::Zero();
       res(0,0) = 1.0 / (aspect * tanHalfFovy);
       res(1,1) = 1.0 / (tanHalfFovy);
       res(2,2) = - (zFar + zNear) / (zFar - zNear);
       res(3,2) = - 1.0;
       res(2,3) = - (2.0 * zFar * zNear) / (zFar - zNear);
       return res;
   }

    Vector2f getScreenCoord() {
        Vector2i pos = mousePos();
        return Vector2f(2.0f * (float)pos.x() / width() - 1.0f,
                        1.0f - 2.0f * (float)pos.y() / height());
    }

    void trackballUpdate() {
        if (!pressing_button && mMouseState == 1 && !left_control_pan &&
            !window->contains(mousePos())) {
            pressing_button = true;
            Vector2f p = getScreenCoord();
            trackball.BeingDrag(p.x(), p.y());
            old_trackball_matrix = trackball_matrix;
        }

        if (!pressing_button && mMouseState == 2 &&
            !window->contains(mousePos())) {
            pressing_button = true;
            Vector2f p = getScreenCoord();
            trackball.beginZoom(p.x(), p.y());
        }

        if (!pressing_button && ((mMouseState == 1 && left_control_pan)
            || mMouseState == 4) && !window->contains(mousePos())) {
            pressing_button = true;
            Vector2f p = getScreenCoord();
            trackball.beginPan(p.x(), p.y());
        }

        if (pressing_button && ((mMouseState == 1 && left_control_pan)
            || mMouseState == 4)) {
            Vector2f p = getScreenCoord();
            view_matrix = trackball.pan(p.x(), p.y())*view_matrix;
        }

        if (pressing_button && mMouseState == 2) {
            Vector2f p = getScreenCoord();
            view_matrix = trackball.zoom(p.x(), p.y())*view_matrix;
        }

        if (pressing_button && mMouseState == 1 && !left_control_pan) {
            Vector2f p = getScreenCoord();
            trackball_matrix = trackball.Drag(p.x(), p.y()) * old_trackball_matrix;
        }

        if (pressing_button && mMouseState == 0) {
            pressing_button = false;
        }

    }

    virtual void drawContents() {
        using namespace nanogui;

        /* Update the trackball states */
        trackballUpdate();

        /* Draw the window contents using OpenGL */
        mShader.bind();

        /* MVP uniforms */
        Matrix4f mv = view_matrix*trackball_matrix*model_matrix;
        Matrix4f p = project_matrix;
        mShader.setUniform("MV", mv);
        mShader.setUniform("P", p);

        /* Setup OpenGL (making sure the GUI doesn't disable these */
        glEnable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);

        /* Render everything */
        if (wireframe) {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0, 1.0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        Vector3f colors(0.98, 0.59, 0.04);
        mShader.setUniform("intensity", colors);
        mShader.setUniform("valenceColor", (valence_coloring)? 1 : 0);
        mShader.setUniform("normal_selector", normals_computation);
        mShader.drawIndexed(GL_TRIANGLES, 0, n_faces);

        if (wireframe) {
            glDisable(GL_POLYGON_OFFSET_FILL);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            colors << 0.0, 0.0, 0.0;
            mShader.setUniform("intensity", colors);
            mShader.drawIndexed(GL_TRIANGLES, 0, n_faces);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        if (normals) {
            mShaderNormals.bind();
            mShaderNormals.setUniform("MV", mv);
            mShaderNormals.setUniform("P", p);
            mShaderNormals.setUniform("normal_selector", normals_computation);
            mShaderNormals.drawIndexed(GL_TRIANGLES, 0, n_faces);
        }
    }
private:
    // Variables for the viewer
    nanogui::GLShader mShader;
    nanogui::GLShader mShaderNormals;
    nanogui::Window *window;

    // Variables for the trackball
    Eigen::Matrix4f trackball_matrix;
    Eigen::Matrix4f old_trackball_matrix;
    bool pressing_button = false;
    bool left_control_pan = false;
    Trackball trackball;

    // the usual MVP
    Eigen::Matrix4f view_matrix;
    Eigen::Matrix4f model_matrix;
    Eigen::Matrix4f project_matrix;

    // Boolean for the viewer
    bool wireframe = false;
    bool normals = false;
    int normals_computation = 0;
    bool valence_coloring = false;

    // Mesh informations
    int n_vertices = 0;
    int n_faces = 0;
    int n_edges = 0;
};
