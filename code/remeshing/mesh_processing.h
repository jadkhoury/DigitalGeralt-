//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#ifndef MESH_PROCESSING_H
#define MESH_PROCESSING_H

#include <surface_mesh/Surface_mesh.h>
#include <Eigen/Sparse>

typedef surface_mesh::Surface_mesh Mesh;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;

namespace mesh_processing {

using std::string;
using std::vector;
enum REMESHING_TYPE : int { AVERAGE = 0, CURV = 1, DIRECTIONAL = 2, DIRECTIONAL_QX = 3};

class MeshProcessing {

public:
    MeshProcessing(const string& filename);
    ~MeshProcessing();

    const surface_mesh::Point get_mesh_center() { return mesh_center_; }
    const float get_dist_max() { return dist_max_; }
    const Eigen::MatrixXf* get_points() { return &points_; }
    const MatrixXu* get_indices() { return &indices_; }
    const Eigen::MatrixXf* get_normals() { return &normals_; }
    const Eigen::MatrixXf* get_colors_valence() { return &color_valence_; }
    const Eigen::MatrixXf* get_colors_unicurvature() { return &color_unicurvature_; }
    const Eigen::MatrixXf* get_colors_gaussian_curv() { return &color_gaussian_curv_; }
    const Eigen::MatrixXf* get_color_curvature() { return &color_curvature_; }
    const unsigned int get_number_of_face() { return mesh_.n_faces(); }


    void remesh (const REMESHING_TYPE &remeshing_type, const int &num_iterations);
    void calc_target_length (const REMESHING_TYPE &remeshing_type);
    void split_long_edges ();
    void collapse_short_edges ();
    void equalize_valences ();
    void tangential_relaxation ();
    //Line added here
    void smooth_target_length();
    vector<Mesh::Vertex> findCommonNeighbours(Mesh::Vertex v1, Mesh::Vertex v2);

    void write_mesh(const string &filename);
    void convertToWireframe();
    std::pair<float, float> computeAllFaceHeight();

    void load_mesh(const string& filename);
    void compute_mesh_properties();

    void calc_mean_curvature();
    void calc_uniform_mean_curvature();
    void calc_gauss_curvature();

    void stars();
    void flowers();
    void add_anti_spearhead( Mesh &new_mesh,Mesh::Vertex middle, Mesh::Vertex  v1, Mesh::Vertex  v2,  Mesh::Vertex  v1r, Mesh::Vertex  v2r, Mesh::Vertex  v1cv2 );
    void add_hexagon(Mesh &new_mesh, Mesh::Vertex vertices[], int size,  Mesh::Vertex origin);
    void add_kite(Mesh &new_mesh, Mesh::Vertex bv0v1, Mesh::Vertex bv1v0, Mesh::Vertex bv1v2,
            Mesh::Vertex bv2v1, Mesh::Vertex bv2v0, Mesh::Vertex bv0v2, Mesh::Vertex middle);

    void give_thickness(float thickness);


        private:
    void calc_weights();
    void calc_edges_weights();
    void calc_vertices_weights();

private:
    Mesh mesh_;
    Mesh mesh_init_;
    surface_mesh::Point mesh_center_ = surface_mesh::Point(0.0f, 0.0f, 0.0f);
    float dist_max_ = 0.0f;

    Eigen::MatrixXf points_;
    MatrixXu indices_;
    Eigen::MatrixXf normals_;
    Eigen::MatrixXf color_valence_;
    Eigen::MatrixXf color_unicurvature_;
    Eigen::MatrixXf color_gaussian_curv_;
    Eigen::MatrixXf color_curvature_;

    void color_coding(Mesh::Vertex_property<surface_mesh::Scalar> prop,
                      Mesh *mesh,
                      Mesh::Vertex_property<surface_mesh::Color> color_prop,
                      surface_mesh::Scalar min_value = 0.0,
                      surface_mesh::Scalar max_value = 0.0, int bound = 20);
    void set_color(Mesh::Vertex v, const surface_mesh::Color& col,
                   Mesh::Vertex_property<surface_mesh::Color> color_prop);
    surface_mesh::Color value_to_color(surface_mesh::Scalar value,
                                       surface_mesh::Scalar min_value,
                                       surface_mesh::Scalar max_value);
    };

}

#endif // MESH_PROCESSING_H
