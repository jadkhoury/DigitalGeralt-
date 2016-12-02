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
#include "mesh_processing.h"
#include <set>

namespace mesh_processing {

using surface_mesh::Point;
using surface_mesh::Scalar;
using surface_mesh::Color;
using std::min;
using std::max;
using std::cout;
using std::endl;

MeshProcessing::MeshProcessing(const string &filename) {
    load_mesh(filename);
}

MeshProcessing::~MeshProcessing() {
}

void MeshProcessing::remesh(const REMESHING_TYPE &remeshing_type,
                            const int &num_iterations) {
    calc_weights();
    calc_mean_curvature();
    calc_uniform_mean_curvature();
    calc_gauss_curvature();
    calc_target_length(remeshing_type);

    // main remeshing loop
    for (int i = 0; i < 5; ++i) {
        split_long_edges();
        collapse_short_edges();
        equalize_valences();
        // tangential_relaxation ();
    }
}

void MeshProcessing::calc_target_length(const REMESHING_TYPE &remeshing_type) {

    Mesh::Vertex_iterator v_it, v_begin, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator vv_c, vv_end;
    Scalar length;
    Scalar mean_length;
    Scalar H;
    Scalar K;

    Mesh::Vertex_property <Scalar> curvature = mesh_.vertex_property<Scalar>("v:mean_curvature", 0);
    Mesh::Vertex_property <Scalar> gauss_curvature = mesh_.vertex_property<Scalar>("v:gauss_curvature", 0);
    Mesh::Vertex_property <Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
    Mesh::Vertex_property <Scalar> target_new_length = mesh_.vertex_property<Scalar>("v:newlength", 0);


    // we caluclate length summing up over all edges
    length = 0.;

    Mesh::Edge_iterator e_it, e_end, e_begin;
    e_end = mesh_.edges_end();
    e_begin = mesh_.edges_begin();
    for (e_it = e_begin; e_it != e_end; ++e_it) {
        length = length + mesh_.edge_length(*e_it);
    }

    //We comppute the mean and define the max and min length in function of the mean
    mean_length = length / mesh_.n_edges();
    float max_length = 1.5 * mean_length;
    float min_length = 0.5 * mean_length;


    if (remeshing_type == AVERAGE) {
        //we simply set the target length to be the mean
        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
            target_length[*v_it] = mean_length;
        }

        float factor = 1.0;
        for (auto v:mesh_.vertices()) {
            target_length[v] *= factor;
        }
    } else if (remeshing_type == CURV) {
        // get max min mean curvature
        v_begin = mesh_.vertices_begin();
        float min_k = fabs(gauss_curvature[*v_begin]);
        float max_k = fabs(gauss_curvature[*v_begin]);
        float mean_k = 0;
        float temp_k = 0;
        // compute mean_k
        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
            temp_k = fabs(gauss_curvature[*v_it]);
            max_k = (temp_k > max_k) ? temp_k : max_k;
            min_k = (temp_k < min_k) ? temp_k : min_k;
            mean_k += temp_k;
        }
        mean_k = mean_k / mesh_.vertices_size();

        // compute target , just put min if k is bigger than the average and max if k is smaller than the average
        float temp_target_length = 0;

        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
            temp_k = fabs(gauss_curvature[*v_it]);
            if( temp_k < mean_k){
                temp_target_length = max_length;
            }
            else {
                temp_target_length = min_length;
            }
            target_length[*v_it] = temp_target_length;

        }

        // smooth desired length
        float total_mean;
        int neigh_ctr;
        Mesh::Vertex current_neighbour;

        for (int i = 0; i < 1; i++) {
            for (auto v: mesh_.vertices()) {

                total_mean = 0.0;
                neigh_ctr = 0;
                vv_c = mesh_.vertices(v);
                vv_end = vv_c;
                do {
                    current_neighbour = *vv_c;
                    total_mean += target_length[current_neighbour];
                    neigh_ctr++;

                } while (++vv_c != vv_end);
                target_new_length[v] = total_mean / float(neigh_ctr);

            }

            for (auto v: mesh_.vertices()) {
                if (!mesh_.is_boundary(v)) {
                    target_length[v] = target_new_length[v];
                    target_new_length[v] = 0.0;
                }
            }

        }



    } else if (remeshing_type == HEIGHT) {
        float max_h = 0.0;
        float min_h = 999999999.9;
        float current_h, normalized_h;
        for(auto v : mesh_.vertices()){
            auto position = mesh_.position(v);
            current_h = position[1];
            max_h = (current_h > max_h) ? current_h : max_h;
            min_h = (current_h < min_h) ? current_h : min_h;
        }
        for (auto v: mesh_.vertices()){
            auto position = mesh_.position(v);
            current_h = position[1];
            normalized_h = (current_h - min_h)/(max_h - min_h);
            target_length[v] = min_length + (max_length - min_length) * normalized_h;
        }
    }
}

void MeshProcessing::split_long_edges() {
    Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
    Mesh::Vertex v0, v1, v;
    Mesh::Halfedge h0, h1;
    bool finished;
    int i;
    double target_L;
    double L;
    Mesh::Vertex_property <Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property <Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    for (finished = false, i = 0; !finished && i < 100; ++i) {
        finished = true;
        // iteration sur les edges
        // No matter with boundaries here
        for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
            h0 = mesh_.halfedge(*e_it, 0);
            h1 = mesh_.halfedge(*e_it, 1);
            v0 = mesh_.to_vertex(h0);
            v1 = mesh_.to_vertex(h1);
            target_L = (target_length[v0] + target_length[v1]) / 2.0;
            L = norm(mesh_.position(v1) - mesh_.position(v0));

            if (L > (4. / 3.) * target_L) {
                finished = false;
                v = mesh_.add_vertex((mesh_.position(v0) + mesh_.position(v1)) / 2.0);
                // comment calculer la normale ????
                // Pour l'instant moyenne des 2 normales
                normals[v] = normalize(normals[v0] + normals[v1]);
                target_length[v] = target_L;
                mesh_.split(*e_it, v);
            }
        }
    }
    mesh_.garbage_collection();
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
}

void MeshProcessing::collapse_short_edges() {
    Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
    Mesh::Vertex v0, v1;
    Mesh::Halfedge h0, h1;
    bool finished, b0, b1;
    int i, val0, val1;
    bool hcol01, hcol10;
    double target_L;
    double L;

    Mesh::Vertex_property <Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

    for (finished = false, i = 0; !finished && i < 100; ++i) {
        finished = true;

        for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
            if (!mesh_.is_deleted(*e_it)){
                h1 = mesh_.halfedge(*e_it, 0);
                h0 = mesh_.halfedge(*e_it, 1);
                v0 = mesh_.to_vertex(h0);
                v1 = mesh_.to_vertex(h1);

                // Are we on boundaries ?
                bool v0_bound = mesh_.is_boundary(v0);
                bool v1_bound = mesh_.is_boundary(v1);


                target_L = (target_length[v0] + target_length[v1]) / 2;
                L = norm(mesh_.position(v1) - mesh_.position(v0));

                if (L < (4. / 5.) * target_L) {
                    // If we can collapse both
                    if (mesh_.is_collapse_ok(h0) && mesh_.is_collapse_ok(h1)) {
                        if (!v0_bound && !v1_bound) {
                            val0 = mesh_.valence(v0);
                            val1 = mesh_.valence(v1);
                            if (val0 >= val1) {
                                mesh_.collapse(h0);
                                finished = false;
                            } else {
                                mesh_.collapse(h1);
                                finished = false;
                            }
                        } else if (!v0_bound && v1_bound) {
                            mesh_.collapse(h1);
                            finished = false;
                        } else if (v0_bound && !v1_bound) {
                            mesh_.collapse(h0);
                            finished = false;
                        }
                    }

                    // If we can collapse only h0
                    // Check if v1 is not on the boundary
                    else if (mesh_.is_collapse_ok(h0) && !mesh_.is_collapse_ok(h1) && !v1_bound) {
                        finished = false;
                        mesh_.collapse(h0);
                    }

                    // If we can collapse only h1
                    // Check if v0 is not on the boundary
                    else if (!mesh_.is_collapse_ok(h0) && mesh_.is_collapse_ok(h1) && !v0_bound) {
                        finished = false;
                        mesh_.collapse(h1);
                    }
                }
            }

        }
    }

    mesh_.garbage_collection();
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();

    if (i == 100) std::cerr << "collapse break\n";
}

void MeshProcessing::equalize_valences() {
    Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
    Mesh::Vertex v0, v1, v2, v3;
    Mesh::Halfedge h, h0, h1, h2, h3;
    int val0, val1, val2, val3;
    int val_opt0, val_opt1, val_opt2, val_opt3;
    int ve0, ve1, ve2, ve3, ve_before, ve_after;
    bool finished;
    int i;


    // flip all edges
    for (finished = false, i = 0; !finished && i < 100; ++i) {
        finished = true;
        for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
            if (!mesh_.is_boundary(*e_it)) {
                // we access to the 4 vertices by halfedges
                h0 = mesh_.halfedge(*e_it, 0);
                h1 = mesh_.halfedge(*e_it, 1);
                v0 = mesh_.to_vertex(h0);
                v1 = mesh_.to_vertex(h1);
                h2 = mesh_.next_halfedge(h0);
                h3 = mesh_.next_halfedge(h1);
                v2 = mesh_.to_vertex(h2);
                v3 = mesh_.to_vertex(h3);

                // we compute the valences for 4 vertices
                val0 = mesh_.valence(v0);
                val1 = mesh_.valence(v1);
                val2 = mesh_.valence(v2);
                val3 = mesh_.valence(v3);

                // we compute optimal valences for 4 verices
                // for v0
                if (mesh_.is_boundary(v0)) {
                    val_opt0 = 4;
                } else {
                    val_opt0 = 6;
                }

                // for v1
                if (mesh_.is_boundary(v1)) {
                    val_opt1 = 4;
                } else {
                    val_opt1 = 6;
                }

                // for v2
                if (mesh_.is_boundary(v2)) {
                    val_opt2 = 4;
                } else {
                    val_opt2 = 6;
                }

                // for v3
                if (mesh_.is_boundary(v3)) {
                    val_opt3 = 4;
                } else {
                    val_opt3 = 6;
                }

                // we compute valence deviation for 4 vertices
                ve0 = abs(val0 - val_opt0);
                ve1 = abs(val0 - val_opt1);
                ve2 = abs(val0 - val_opt2);
                ve3 = abs(val0 - val_opt3);

                // sum of square of deviations for 2 edges
                ve_before = ve0 * ve0 + ve1 * ve1;
                ve_after = ve2 * ve2 + ve3 * ve3;


                if (ve_after < ve_before && mesh_.is_flip_ok(*e_it)) {
                    mesh_.flip(*e_it);
                    finished = false;
                }
            }
        }
    }

    mesh_.garbage_collection();
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
    if (i == 100) std::cerr << "flip break\n";
}

void MeshProcessing::tangential_relaxation() {
    Mesh::Vertex_iterator v_it, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator vv_c, vv_end;
    int valence;
    Point u, n;
    Point laplace;

    Mesh::Vertex_property <Point> normals = mesh_.vertex_property<Point>("v:normal");
    Mesh::Vertex_property <Point> update = mesh_.vertex_property<Point>("v:update");


    // smooth
    for (int iters = 0; iters < 10; ++iters) {
        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
            if (!mesh_.is_boundary(*v_it)) {
            }
        }

        for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
            if (!mesh_.is_boundary(*v_it))
                mesh_.position(*v_it) += update[*v_it];
    }
}
// ========================================================================
// EXERCISE 1.1
// ========================================================================



void MeshProcessing::calc_uniform_mean_curvature() {
    Mesh::Vertex_property <Scalar> v_unicurvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0);
    // ------------- IMPLEMENT HERE ---------
    // Define vertex iterator
    Mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();
    // Define vertex circulator
    Mesh::Vertex_around_vertex_circulator vc, vc_end;
    // Current vertex value
    Mesh::Vertex current_v;
    int neighbors_counter;
    surface_mesh::Vec3 approx, final_approx;
    // Iterating over all vertices
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        neighbors_counter = 0;
        approx = surface_mesh::Vec3(0.0, 0.0, 0.0);
        current_v = *v_it;
        vc = mesh_.vertices(current_v);
        vc_end = vc;
        do {
            neighbors_counter++;
            approx = approx + mesh_.position(*vc) - mesh_.position(current_v);
        } while (++vc != vc_end);
        final_approx = approx / neighbors_counter;
        // Storing the computing v_unicurvature
        v_unicurvature[current_v] = norm(final_approx) / 2;
    }
    // ------------- IMPLEMENT HERE ---------
}

// ========================================================================
// EXERCISE 1.2
// ========================================================================


float computeAngle(surface_mesh::Vec3 origin, surface_mesh::Vec3 p1, surface_mesh::Vec3 p2) {
    surface_mesh::Vec3 v1 = normalize(p1 - origin);
    surface_mesh::Vec3 v2 = normalize(p2 - origin);
    return acos(dot(v1, v2));
}

void MeshProcessing::calc_mean_curvature() {
    Mesh::Vertex_property <Scalar> v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property <Scalar> e_weight =
            mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property <Scalar> v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    // Define Half-edge circulator ...
    Mesh::Halfedge_around_vertex_circulator h, h_end;
    Mesh::Halfedge current_h;
    Mesh::Edge current_e;

    // Define vertex circulator
    // Surface_mesh::Vertex_around_vertex_circulator vc, vc_end;

    surface_mesh::Vec3 approx;
    surface_mesh::Vec3 final_approx;
    float wi, w;

    // Iterating over all vertices
    for (auto current_v: mesh_.vertices()) {

        approx = surface_mesh::Vec3(0.0, 0.0, 0.0);

        wi = 0.;

        // initialisant le Halfedge circulator
        h = h_end = mesh_.halfedges(current_v);

        do {
            current_h = *h;
            current_e = mesh_.edge(current_h);
            wi = e_weight[current_e];
            approx = approx + wi * (mesh_.position(mesh_.to_vertex(current_h)) - mesh_.position(current_v));

        } while (++h != h_end);

        w = v_weight[current_v];
        final_approx = w * approx;

        // Storing the computing v_curvature
        v_curvature[current_v] = norm(final_approx) / 2;
    }

}
// ========================================================================
// EXERCISE 1.3
// ========================================================================

void MeshProcessing::calc_gauss_curvature() {
    Mesh::Vertex_property <Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property <Scalar> v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    // Vertex iterator to iterate the current vertex over all vertices
    Mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();

    // Vertex circulators, one added to access to consecutive neighbours
    Mesh::Vertex_around_vertex_circulator vc, vc_next, vc_end;

    //Current vertex variables
    Mesh::Vertex current_v;
    float angle;
    float sum_angles;
    float gauss_approx;
    surface_mesh::Vec3 current_pos, n1, n2;

    // Iterating current_v over all vertices
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        // Update current vertex related variables
        current_v = *v_it;
        current_pos = mesh_.position(current_v);
        vc = mesh_.vertices(current_v);
        vc_next = vc;
        vc_end = vc;
        sum_angles = 0.;

        //Iterating over all current_v neighbors
        do {

            //vc_next will always be one step ahead vc to be its sucessive neighbour
            ++vc_next;
            //Computing angle
            n1 = mesh_.position(*vc);
            n2 = mesh_.position(*vc_next);
            angle = computeAngle(current_pos, n1, n2);

            // computing sum of angles
            sum_angles = sum_angles + angle;
        } while (++vc != vc_end);

        gauss_approx = (2 * M_PI - sum_angles) * 2 * v_weight[current_v];

        // Store the gaussian curvature as a property
        v_gauss_curvature[current_v] = gauss_approx;
    }

}


void MeshProcessing::calc_weights() {
    calc_edges_weights();
    calc_vertices_weights();
}

void MeshProcessing::calc_edges_weights() {
    auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    auto points = mesh_.vertex_property<Point>("v:point");

    Mesh::Halfedge h0, h1, h2;
    Point p0, p1, p2, d0, d1;

    for (auto e: mesh_.edges()) {
        e_weight[e] = 0.0;

        h0 = mesh_.halfedge(e, 0);
        p0 = points[mesh_.to_vertex(h0)];

        h1 = mesh_.halfedge(e, 1);
        p1 = points[mesh_.to_vertex(h1)];

        if (!mesh_.is_boundary(h0)) {
            h2 = mesh_.next_halfedge(h0);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
        }

        if (!mesh_.is_boundary(h1)) {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
        }
    }
}

void MeshProcessing::calc_vertices_weights() {
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    for (auto v: mesh_.vertices()) {
        area = 0.0;
        vf_c = mesh_.faces(v);

        if (!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh_.vertices(*vf_c);

            const Point &P = mesh_.position(*fv_c);
            ++fv_c;
            const Point &Q = mesh_.position(*fv_c);
            ++fv_c;
            const Point &R = mesh_.position(*fv_c);

            area += norm(cross(Q - P, R - P)) * 0.5f * 0.3333f;

        } while (++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
}

void MeshProcessing::load_mesh(const string &filename) {
    if (!mesh_.read(filename)) {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

    cout << "Mesh " << filename << " loaded." << endl;
    cout << "# of vertices : " << mesh_.n_vertices() << endl;
    cout << "# of faces : " << mesh_.n_faces() << endl;
    cout << "# of edges : " << mesh_.n_edges() << endl;

    // Compute the center of the mesh
    mesh_center_ = Point(0.0f, 0.0f, 0.0f);
    for (auto v: mesh_.vertices()) {
        mesh_center_ += mesh_.position(v);
    }
    mesh_center_ /= mesh_.n_vertices();

    // Compute the maximum distance from all points in the mesh and the center
    dist_max_ = 0.0f;
    for (auto v: mesh_.vertices()) {
        if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
            dist_max_ = distance(mesh_center_, mesh_.position(v));
        }
    }

    compute_mesh_properties();

    // Store the original mesh, this might be useful for some computations
    mesh_init_ = mesh_;
}

void MeshProcessing::compute_mesh_properties() {
    Mesh::Vertex_property <Point> vertex_normal =
            mesh_.vertex_property<Point>("v:normal");
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
    Mesh::Vertex_property <Color> v_color_valence =
            mesh_.vertex_property<Color>("v:color_valence",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property <Color> v_color_unicurvature =
            mesh_.vertex_property<Color>("v:color_unicurvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property <Color> v_color_curvature =
            mesh_.vertex_property<Color>("v:color_curvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property <Color> v_color_gaussian_curv =
            mesh_.vertex_property<Color>("v:color_gaussian_curv",
                                         Color(1.0f, 1.0f, 1.0f));

    Mesh::Vertex_property <Scalar> vertex_valence =
            mesh_.vertex_property<Scalar>("v:valence", 0.0f);
    for (auto v: mesh_.vertices()) {
        vertex_valence[v] = mesh_.valence(v);
    }

    Mesh::Vertex_property <Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_property <Scalar> v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Vertex_property <Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);

    calc_weights();
    calc_uniform_mean_curvature();
    calc_mean_curvature();
    calc_gauss_curvature();
    color_coding(vertex_valence, &mesh_, v_color_valence, 3 /* min */,
                 8 /* max */);
    color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
    color_coding(v_curvature, &mesh_, v_color_curvature);
    color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);

    // get the mesh attributes and upload them to the GPU
    int j = 0;
    unsigned int n_vertices(mesh_.n_vertices());

    // Create big matrices to send the data to the GPU with the required
    // format
    color_valence_ = Eigen::MatrixXf(3, n_vertices);
    color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
    color_curvature_ = Eigen::MatrixXf(3, n_vertices);
    color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
    normals_ = Eigen::MatrixXf(3, n_vertices);
    points_ = Eigen::MatrixXf(3, n_vertices);
    indices_ = MatrixXu(3, mesh_.n_faces());

    for (auto f: mesh_.faces()) {
        std::vector<float> vv(3);
        int k = 0;
        for (auto v: mesh_.vertices(f)) {
            vv[k] = v.idx();
            ++k;
        }
        indices_.col(j) << vv[0], vv[1], vv[2];
        ++j;
    }

    j = 0;
    for (auto v: mesh_.vertices()) {
        points_.col(j) << mesh_.position(v).x,
                mesh_.position(v).y,
                mesh_.position(v).z;

        normals_.col(j) << vertex_normal[v].x,
                vertex_normal[v].y,
                vertex_normal[v].z;

        color_valence_.col(j) << v_color_valence[v].x,
                v_color_valence[v].y,
                v_color_valence[v].z;

        color_unicurvature_.col(j) << v_color_unicurvature[v].x,
                v_color_unicurvature[v].y,
                v_color_unicurvature[v].z;

        color_curvature_.col(j) << v_color_curvature[v].x,
                v_color_curvature[v].y,
                v_color_curvature[v].z;

        color_gaussian_curv_.col(j) << v_color_gaussian_curv[v].x,
                v_color_gaussian_curv[v].y,
                v_color_gaussian_curv[v].z;
        ++j;
    }
}

void MeshProcessing::color_coding(Mesh::Vertex_property <Scalar> prop, Mesh *mesh,
                                  Mesh::Vertex_property <Color> color_prop, Scalar min_value,
                                  Scalar max_value, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    if (min_value == 0.0 && max_value == 0.0) {
        // discard upper and lower bound
        unsigned int n = values.size() - 1;
        unsigned int i = n / bound;
        std::sort(values.begin(), values.end());
        min_value = values[i];
        max_value = values[n - 1 - i];
    }

    // map values to colors
    for (auto v: mesh->vertices()) {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color &col,
                               Mesh::Vertex_property <Color> color_prop) {
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0 / 4.0 * (max_value - min_value);
    v1 = min_value + 1.0 / 4.0 * (max_value - min_value);
    v2 = min_value + 2.0 / 4.0 * (max_value - min_value);
    v3 = min_value + 3.0 / 4.0 * (max_value - min_value);
    v4 = min_value + 4.0 / 4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0) {
        col = Color(0, 0, 1);
    } else if (value > v4) {
        col = Color(1, 0, 0);
    } else if (value <= v2) {
        if (value <= v1) { // [v0, v1]
            Scalar u = (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        } else { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1 - u);
        }
    } else {
        if (value <= v3) { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        } else { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1 - u, 0);
        }
    }
    return col;
}


}


