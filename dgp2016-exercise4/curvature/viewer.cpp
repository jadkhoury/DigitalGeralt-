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
#include "viewer.h"
#include <surface_mesh/Surface_mesh.h>

using std::min;
using std::max;
using namespace surface_mesh;

typedef Surface_mesh Mesh;

// ========================================================================
// NOTE : We've only included the functions you need to implement (or the
//        ones that you will need to use) in the cpp file. This is not the
//        best practice as you normaly would have all the implementation of
//        the functions here and only the declaration in the header file
//        but it allows you to have all the things you need here.
// ========================================================================

void print(std::string s){
    std::cout << s << std::endl;
}

float computeAngle(Vec3 origin, Vec3 p1, Vec3 p2){
    Vec3 v1 = normalize(p1 - origin);
    Vec3 v2 = normalize(p2 - origin);
    return acos(dot(v1, v2));
}


// ========================================================================
// EXERCISE 1.1
// ========================================================================
void Viewer::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature = mesh.vertex_property<Scalar>("v:unicurvature", 0);
    // ------------- IMPLEMENT HERE ---------
    // Define vertex iterator
    Surface_mesh :: Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    // Define vertex circulator
    Surface_mesh :: Vertex_around_vertex_circulator vc, vc_end;

    // Current vertex value
    Surface_mesh :: Vertex current_v;
    int neighbors_counter;
    Vec3 approx, final_approx ;
    // Iterating over all vertices
    for(v_it = v_begin; v_it != v_end; ++ v_it){
        neighbors_counter = 0;
        approx = Vec3(0.0, 0.0, 0.0);
        current_v = *v_it;
        vc = mesh.vertices(current_v);
        vc_end = vc;
        do {
            neighbors_counter++;
            approx = approx + mesh.position(*vc) - mesh.position(current_v);
        } while(++vc != vc_end);
        final_approx = approx/neighbors_counter;
        // Storing the computing v_unicurvature
        v_unicurvature[current_v] = norm(final_approx)/2;
    }
    // ------------- IMPLEMENT HERE ---------
}
// ========================================================================
// EXERCISE 1.2
// ========================================================================
void Viewer::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature = mesh.vertex_property<Scalar>("v:curvature", 0);
    Mesh::Edge_property<Scalar> e_weight = mesh.edge_property<Scalar>("e:weight", 0);
    Mesh::Vertex_property<Scalar>  v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    // ------------- IMPLEMENT HERE ---------
    Surface_mesh::Halfedge_around_vertex_circulator h, h_end;
    Surface_mesh::Halfedge current_h;
    Surface_mesh::Edge current_e;
    float w, wi;
    Vec3 tmp;
    for (auto current_v: mesh.vertices()) {
        tmp = Vec3(0.0, 0.0, 0.0);
        w = v_weight[current_v];
        h = h_end = mesh.halfedges(current_v);
        do{
            current_h = *h;
            current_e = mesh.edge(current_h);
            wi = e_weight[current_e];
            tmp += wi * (mesh.position(mesh.to_vertex(current_h)) - mesh.position(current_v));
        }while(++h != h_end);
        tmp *= w;
        v_curvature[current_v] = 0.5*norm(tmp);
    }
    // ------------- IMPLEMENT HERE ---------
}
// ========================================================================
// EXERCISE 1.3
// ========================================================================
void Viewer::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh.vertex_property<Scalar>("v:gauss_curvature", 0);
    Mesh::Vertex_property<Scalar> v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    // ------------- IMPLEMENT HERE ---------
    Surface_mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();
    // Vertex circulators, one added to access to consecutive neighbours
    Surface_mesh::Vertex_around_vertex_circulator vc, vc_next, vc_end;
    //Current vertex and others useful variables
    Surface_mesh::Vertex current_v;
    float angle, sum_angles;
    float gauss_approx;
    Vec3 current_pos, n1, n2;
    // Iterating current_v over all vertices
    for(v_it = v_begin; v_it != v_end; ++v_it){
        // Update current vertex related variables
        current_v = *v_it;
        current_pos = mesh.position(current_v);
        vc = mesh.vertices(current_v);
        vc_next = vc;
        vc_end = vc;
        sum_angles = 0.0;
        //Iterating over all current_v neighbors
        do{
            //vc_next will always be one step ahead vc to be its sucessive neighbour
            ++vc_next;
            //Computing angle
            n1 =  mesh.position(*vc);
            n2 =  mesh.position(*vc_next);
            angle = computeAngle(current_pos, n1, n2);
            // computing sum of angles
            sum_angles=sum_angles+angle;
        } while(++vc != vc_end);
        gauss_approx = (2*M_PI-sum_angles)*2*v_weight[current_v];
        // Store the gaussian curvature as a property
        v_gauss_curvature[current_v] = gauss_approx;
    }
    // ------------- IMPLEMENT HERE ---------
}

// ========================================================================
// EXERCISE 2.1
// ========================================================================

void Viewer::uniform_smooth(unsigned int n_iters) {

    for (unsigned int iter=0; iter<n_iters; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, update its position according to the uniform Laplacian operator
        calc_edges_weights();
        calc_vertices_weights();
        Vec3 Lu, current_pos;
        int n = 0;
        Surface_mesh::Vertex_around_vertex_circulator vc, vc_end;
        //iterating over the vertices
        for(auto v : mesh.vertices()){
            if(mesh.is_boundary(v) == false){
                n = 0;
                Lu = Vec3(0.0, 0.0, 0.0);
                vc = mesh.vertices(v);
                vc_end = vc;
                current_pos = mesh.position(v);
                //Iterating over the neighbours
                do{
                    n++;
                    Lu  = Lu + mesh.position(*vc);
                }while(++vc != vc_end);
                Lu = (Lu/n) - current_pos;
                mesh.position(v) = mesh.position(v) + 0.5 * Lu;
            }
        }
    }
    // ------------- IMPLEMENT HERE ---------

    // update face and vertex normals
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
// ========================================================================
// EXERCISE 2.2
// ========================================================================
void Viewer::smooth(unsigned int n_iters) {
    for (unsigned int iter=0; iter<n_iters; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        calc_edges_weights();
        Mesh::Edge_property<Scalar> e_weight = mesh.edge_property<Scalar>("e:weight", 0);
        // Define Half-edge circulator ...
        Surface_mesh::Halfedge_around_vertex_circulator h, h_end;
        Surface_mesh::Halfedge current_h;
        Surface_mesh::Edge current_e;
        // Defining usefull variables
        Vec3 approx;
        float wi;
        float sum_wi;
        Vec3 nLb;
        // Iterating over all vertices
        for (auto current_v: mesh.vertices()){
            if(mesh.is_boundary(current_v) == false){
                approx = Vec3(0.0, 0.0, 0.0);
                wi=0.;
                sum_wi=0.;
                // Initializing the Halfedge circulator
                h = mesh.halfedges(current_v);
                h_end = h;
                do {
                    current_h = *h;
                    current_e = mesh.edge(current_h);
                    wi = e_weight[current_e];
                    sum_wi = sum_wi + wi;
                    approx = approx + wi*(mesh.position(mesh.to_vertex(current_h)) - mesh.position(current_v)) ;
                } while(++h != h_end);
                nLb = approx / sum_wi;
                // computing the new position of each vertex;
                mesh.position(current_v) = mesh.position(current_v) +  nLb*0.5;
            }
        }
        // ------------- IMPLEMENT HERE ---------
        // update face and vertex normals
        mesh.update_face_normals();
        mesh.update_vertex_normals();
    }
}
// ========================================================================
// EXERCISE 3
// ========================================================================
void Viewer::uniform_laplacian_enhance_feature(int enhancement_smoothing_iterations,
                                               float enhancement_coef) {
    // ------------- IMPLEMENT HERE ---------
    // Creating new property to store old position if it doesn't already exists
    Surface_mesh::Vertex_property<Point> old_pos;
    if(mesh.vertex_property<Point>("v:old")){
        print("v:old already exists");
        old_pos = mesh.vertex_property<Point>("v:old");
    } else {
        print("v:old is not yet defined");
        old_pos = mesh.add_vertex_property<Point>("v:old");
    }
    // Storing old position in the old pos property
    for(auto v: mesh.vertices()){
        old_pos[v] = mesh.position(v);
    }
    //Applying the smoothing
    Viewer::uniform_smooth(enhancement_smoothing_iterations);
    //Computing new position
    Vec3 v_in, v_out;
    for(auto v: mesh.vertices()){
        v_in = old_pos[v];
        v_out = mesh.position(v);
        mesh.position(v) = v_out + enhancement_coef * (v_in-v_out);
    }
    // ------------- IMPLEMENT HERE ---------
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
// ========================================================================
// EXERCISE 3
// ========================================================================
void Viewer::laplace_beltrami_enhance_feature(int enhancement_smoothing_iterations,
                                              float enhancement_coef) {
    // ------------- IMPLEMENT HERE ---------
    //We do exactly the same thing with smooth instead of uniform_smooth
    Surface_mesh::Vertex_property<Point> old_pos;
    if(mesh.vertex_property<Point>("v:old")){
        print("v:old already exists");
        old_pos = mesh.vertex_property<Point>("v:old");
    } else {
        print("v:old is not yet defined");
        old_pos = mesh.add_vertex_property<Point>("v:old");
    }
    for(auto v: mesh.vertices()){
        old_pos[v] = mesh.position(v);
    }
    Viewer::smooth(enhancement_smoothing_iterations);
    Vec3 v_in, v_out;
    for(auto v: mesh.vertices()){
        v_in = old_pos[v];
        v_out = mesh.position(v);
        mesh.position(v) = v_out + enhancement_coef * (v_in-v_out);
    }
    // ------------- IMPLEMENT HERE ---------
    mesh.update_face_normals();
    mesh.update_vertex_normals();
}
