//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss, Alexandru Ichim
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



// ============================================================================
// EXERCISE 5.2
// ============================================================================
    void MeshProcessing::give_thickness() {

        const float thickness = 1;

        Mesh::Vertex_property <Point> vertex_normal =
                mesh_.vertex_property<Point>("v:normal");


        // Add propreties to generate tickness

        // why in not reconize the type ???
        Mesh::Vertex_property <surface_mesh::Surface_mesh::Vertex> associated_vertex =
                mesh_.add_vertex_property<surface_mesh::Surface_mesh::Vertex>("v:associated_vertex");

        // determine if this vertex was generate to give a thickness to the surface
        Mesh::Vertex_property<bool> is_thickness =
                mesh_.add_vertex_property<bool>("v:is_thickness");







        /*
         * Iterate over all vertex and create a new vertex in the oposite direction of the normal
         * We call this new vertex the assiociated vertex.
         */
        Mesh::Vertex_iterator v_it, v_begin, v_end;

        v_begin = mesh_.vertices_begin();
        v_end = mesh_.vertices_end();
        Point p, normal;
        Mesh::Vertex assiociated_v;


        for (v_it = v_begin; v_it != v_end; ++v_it) {
            p = mesh_.position(*v_it);
            normal = vertex_normal[*v_it];

            assiociated_v = mesh_.add_vertex(p - ( thickness * normal) );
            associated_vertex[*v_it] = assiociated_v;


        }





        // the new vertex should not have faces

        Mesh::Face_iterator f_it, f_begin, f_end;
        f_begin = mesh_.faces_begin();
        f_end = mesh_.faces_end();

        Mesh::Vertex_around_face_circulator vc, vc_end;


        for (f_it = f_begin; f_it != f_end; ++f_it) {

            Mesh::Face f = *f_it;

            vc = mesh_.vertices(f);
            vc_end = vc;


            Mesh::Vertex associated_vertices[3];

            int associated_v;

            int i = 0;
            do {

                associated_vertices[i] = associated_vertex[*vc];
                ++i;


            } while (++vc != vc_end);

            mesh_.add_triangle(associated_vertices[0], associated_vertices[1], associated_vertices[2] );

        }


        cout << "# of vertices : " << mesh_.n_vertices() << endl;
        cout << "# of faces : " << mesh_.n_faces() << endl;
        cout << "# of edges : " << mesh_.n_edges() << endl;


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

    }


    MeshProcessing::~MeshProcessing() {}
}
