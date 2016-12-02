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

MeshProcessing::MeshProcessing(const string& filename) {
    load_mesh(filename);
}

// ============================================================================
// EXERCISE 5.1
// ============================================================================
void print(std::string s){
    std::cout << s << std::endl;
}

float computeAngle(surface_mesh::Vec3 origin, surface_mesh::Vec3 p1, surface_mesh::Vec3 p2){
    surface_mesh::Vec3 v1 = normalize(p1 - origin);
    surface_mesh::Vec3 v2 = normalize(p2 - origin);
    return acos(dot(v1, v2));
}


void MeshProcessing::implicit_smoothing(const double timestep) {

    const int n = mesh_.n_vertices();

    // get vertex position
    auto points = mesh_.vertex_property<Point>("v:point");

    // compute cotan edge weights and vertex areas
    calc_weights ();
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

    // A*X = B
    Eigen::SparseMatrix<double> A(n,n);
    Eigen::MatrixXd B(n,3);

    // nonzero elements of A as triplets: (row, column, value)
    std::vector< Eigen::Triplet<double> > triplets;

    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 5.1 HERE


    // We add an edge property to keep track of the edges we already used to compute M
    Mesh::Edge_property<bool> edge_done;
    if(mesh_.edge_property<bool>("e:done")){
        print("e:done already exists");
        edge_done = mesh_.edge_property<bool>("e:done", false);
    } else {
        print("e:done is not yet defined");
        edge_done = mesh_.add_edge_property<bool>("e:done", false);
    }
    //Declare halfedge circulator
    Mesh::Halfedge_around_vertex_circulator h, h_end;
    //Declare current iteration variables
    Mesh::Halfedge current_h;
    Mesh::Edge current_e;
    float Ei;
    float v_w;
    int  idx, neighbour_idx;
    float lambda = timestep / 50.0;
    Eigen::Triplet<double> value;
    surface_mesh::Vec3 Bi;

    //We iterate over all vertices
    for(auto current_v: mesh_.vertices()){
        /// We start by filling A
        //We initialize the current vertex variables
        Ei = 0;
        v_w = area_inv[current_v];
        h = h_end = mesh_.halfedges(current_v);
        idx = current_v.idx();
        //We iterate over all the halfedges from the current vertex
        do{
            current_h = *h;
            current_e = mesh_.edge(current_h);
            //Summing over the outgoing edges weight for Ei
            Ei += cotan[current_e];
            //If we didn't work with the current edge yet
            if(!edge_done[current_e]){
                //Get the index of the neighbour pointed by the halfedge
                neighbour_idx = mesh_.to_vertex(current_h).idx();
                //Create the new triplets and push them
                value = Eigen::Triplet<double>(idx, neighbour_idx, -lambda * cotan[current_e]);
                triplets.push_back(value);
                value = Eigen::Triplet<double>(neighbour_idx, idx, -lambda * cotan[current_e]);
                triplets.push_back(value);
                //Mark the current edge as done so we don't add twice its weights to M
                edge_done[current_e] = true;
            }
        }while(++h != h_end);
        //D is diagonal so D-1 is D with all value inversed
        value =  Eigen::Triplet<double>(idx, idx, (1.0/v_w) + lambda*Ei);
        triplets.push_back(value);

        ///Then we fill B
        Bi = mesh_.position(current_v) / v_w;
        B(idx, 0) = Bi[0];
        B(idx, 1) = Bi[1];
        B(idx, 2) = Bi[2];
    }

    //When we are done with computing A and B, we set the edge_done to false for all edges for the next call of the function
    for(auto e: mesh_.edges()){
        edge_done[e] = false;
    }

    // ========================================================================

    // build sparse matrix from triplets
    A.setFromTriplets(triplets.begin(), triplets.end());

    // solve A*X = B
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver(A);
    Eigen::MatrixXd X = solver.solve(B);

    // copy solution
    for (int i = 0; i < n; ++i)
    {
        Mesh::Vertex v(i);
        for (int dim = 0; dim < 3; ++dim)
            points[v][dim] = X(i, dim);
    }

    // clean-up
    mesh_.remove_vertex_property(area_inv);
    mesh_.remove_edge_property(cotan);
}

// ============================================================================
// EXERCISE 5.2
// ============================================================================
void MeshProcessing::minimal_surface() {

    const int n = mesh_.n_vertices();

    // get vertex position
    auto points = mesh_.vertex_property<Point>("v:point");
    auto points_init = mesh_init_.vertex_property<Point>("v:point");

    // compute cotan edge weights and vertex areas
    calc_weights ();
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

    // A*X = B
    Eigen::SparseMatrix<double> L (n, n);
    Eigen::MatrixXd rhs (Eigen::MatrixXd::Zero (n, 3));

    // nonzero elements of A as triplets: (row, column, value)
    std::vector< Eigen::Triplet<double> > triplets_L;

    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 5.2 HERE

    //Declare halfedge circulator
    Mesh::Halfedge_around_vertex_circulator h, h_end, h1, h1_end;

    //Declare current iteration variables
    Mesh::Halfedge current_h;
    Mesh::Edge current_e;

    float Ei;
    float v_w;
    float v_w2;
    int  idx, neighbour_idx;
    Eigen::Triplet<double> value;

    //We iterate over all vertices
    for(auto current_v: mesh_.vertices()){
        Ei = 0;
        v_w = area_inv[current_v];
        h = h_end = mesh_.halfedges(current_v);
        idx = current_v.idx();



        //We iterate over all the halfedges from the current vertex
        do{
            current_h = *h;
            current_e = mesh_.edge(current_h);

            //Summing over the outgoing edges weight for Ei
            Ei += cotan[current_e];

            if(mesh_.is_boundary(current_v)==false && mesh_.is_boundary(mesh_.to_vertex(current_h))==false ){

                //Get the index of the neighbour pointed by the halfedge
                neighbour_idx = mesh_.to_vertex(current_h).idx();

                //Create the new triplets and push them
                value = Eigen::Triplet<double>(idx, neighbour_idx, cotan[current_e]*v_w);
                triplets_L.push_back(value);
            }


        }while(++h != h_end);

        if(mesh_.is_boundary(current_v) == false){
            value = Eigen::Triplet<double>(idx, idx, -Ei*v_w);
            triplets_L.push_back(value);
        }

        if(mesh_.is_boundary(current_v) == true){
            value = Eigen::Triplet<double>(idx, idx, 1);
            triplets_L.push_back(value);
            rhs(idx,0)=rhs(idx,0)+mesh_.position(current_v)[0];
            rhs(idx,1)=rhs(idx,1)+mesh_.position(current_v)[1];
            rhs(idx,2)=rhs(idx,2)+mesh_.position(current_v)[2];
            h1 = h1_end = mesh_.halfedges(current_v);

            do{
                current_h = *h1;
                current_e = mesh_.edge(current_h);
                //Get the index of the neighbour pointed by the halfedge
                neighbour_idx = mesh_.to_vertex(current_h).idx();

                v_w2 = area_inv[mesh_.to_vertex(current_h)];
                if(mesh_.is_boundary(mesh_.to_vertex(current_h)) == false){
                    rhs(neighbour_idx,0)=rhs(neighbour_idx,0)-cotan[current_e]*v_w2*mesh_.position(current_v)[0];
                    rhs(neighbour_idx,1)=rhs(neighbour_idx,1)-cotan[current_e]*v_w2*mesh_.position(current_v)[1];
                    rhs(neighbour_idx,2)=rhs(neighbour_idx,2)-cotan[current_e]*v_w2*mesh_.position(current_v)[2];
                }
            }while(++h1 != h1_end);

        }
    }




    // ========================================================================

    L.setFromTriplets (triplets_L.begin (), triplets_L.end ());

    // solve A*X = B
    Eigen::SparseLU< Eigen::SparseMatrix<double> > solver(L);
    if (solver.info () != Eigen::Success) {
        printf("linear solver init failed.\n");
    }

    Eigen::MatrixXd X = solver.solve(rhs);
    if (solver.info () != Eigen::Success) {
        printf("linear solver failed.\n");
    }

    // copy solution
    for (int i = 0; i < n; ++i) {
        Mesh::Vertex v(i);
        for (int dim = 0; dim < 3; ++dim) {
            points[v][dim] += 1. * (X(i, dim) - points[v][dim]);
        }
    }

    // clean-up
    mesh_.remove_vertex_property(area_inv);
    mesh_.remove_edge_property(cotan);
}

void MeshProcessing::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0);
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
    surface_mesh::Vec3 approx, final_approx ;
    // Iterating over all vertices
    for(v_it = v_begin; v_it != v_end; ++ v_it){
        neighbors_counter = 0;
        approx = surface_mesh::Vec3(0.0, 0.0, 0.0);
        current_v = *v_it;
        vc = mesh_.vertices(current_v);
        vc_end = vc;
        do {
            neighbors_counter++;
            approx = approx + mesh_.position(*vc) - mesh_.position(current_v);
        } while(++vc != vc_end);
        final_approx = approx/neighbors_counter;
        // Storing the computing v_unicurvature
        v_unicurvature[current_v] = norm(final_approx)/2;
    }
    // ------------- IMPLEMENT HERE ---------
}

void MeshProcessing::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight =
            mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar>  v_weight =
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
    float wi,w;

    // Iterating over all vertices
    for (auto current_v: mesh_.vertices()){

        approx = surface_mesh::Vec3(0.0, 0.0, 0.0);

        wi=0.;

        // initialisant le Halfedge circulator
        h = h_end = mesh_.halfedges(current_v);

        do {
            current_h = *h;
            current_e = mesh_.edge(current_h);
            wi = e_weight[current_e];
            approx = approx + wi*(mesh_.position(mesh_.to_vertex(current_h)) - mesh_.position(current_v)) ;

        } while(++h != h_end);

        w = v_weight[current_v];
        final_approx = w*approx;

        // Storing the computing v_curvature
        v_curvature[current_v] = norm(final_approx)/2;
    }

}

void MeshProcessing::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
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
    for(v_it = v_begin; v_it != v_end; ++v_it){
        // Update current vertex related variables
        current_v = *v_it;
        current_pos = mesh_.position(current_v);
        vc = mesh_.vertices(current_v);
        vc_next = vc;
        vc_end = vc;
        sum_angles = 0.;

        //Iterating over all current_v neighbors
        do{

            //vc_next will always be one step ahead vc to be its sucessive neighbour
            ++vc_next;
            //Computing angle
            n1 =  mesh_.position(*vc);
            n2 =  mesh_.position(*vc_next);
            angle = computeAngle(current_pos, n1, n2);

            // computing sum of angles
            sum_angles=sum_angles+angle;
        } while(++vc != vc_end);

        gauss_approx = (2*M_PI-sum_angles)*2*v_weight[current_v];

        // Store the gaussian curvature as a property
        v_gauss_curvature[current_v] = gauss_approx;
    }

}

void MeshProcessing::uniform_smooth(const unsigned int iterations) {

    Mesh::Vertex_property<Scalar> v_unicurvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0);
    for (unsigned int iter=0; iter<iterations; ++iter) {
        // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------

        // Define vertex iterator
        Mesh::Vertex_iterator v_it, v_begin, v_end;
        v_begin = mesh_.vertices_begin();
        v_end = mesh_.vertices_end();

        // Define vertex circulator
        Mesh::Vertex_around_vertex_circulator vc, vc_end;

        // Define current vertex value
        Mesh::Vertex current_v;
        int neighbors_counter;
        surface_mesh::Vec3 Lu;

        // Iterating over all vertices
        for(v_it = v_begin; v_it != v_end; ++ v_it){
            current_v = *v_it;

            // Check if it is on the boundary or not
            if(mesh_.is_boundary(current_v) == false){

                neighbors_counter = 0;
                surface_mesh::Vec3 approx;
                approx = surface_mesh::Vec3(0.,0.,0.);

                vc = mesh_.vertices(current_v);
                vc_end = vc;

                do {
                    neighbors_counter++;
                    approx = approx + mesh_.position(*vc) - mesh_.position(current_v);
                } while(++vc != vc_end);

                Lu = approx/neighbors_counter;

                mesh_.position(current_v) =  mesh_.position(current_v) +  Lu/2;
            }
        }

    }
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
}

void MeshProcessing::smooth(const unsigned int iterations) {

    Mesh::Edge_property<Scalar> e_weight = mesh_.edge_property<Scalar>("e:weight", 0);
    for (unsigned int iter=0; iter<iterations; ++iter) {
        // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------

        // Perform Laplace-Beltrami smoothing:
        // 1) precompute edge weights using calc_edge_weights()
        calc_edges_weights();

        // 2) for each non-boundary vertex, update its position using the normalized Laplace-Beltrami operator
        //    (Hint: use the precomputed edge weights in the edge property "e:weight")

        // Define vertex iterator
        Mesh::Vertex_iterator v_it, v_begin, v_end;
        v_begin = mesh_.vertices_begin();
        v_end = mesh_.vertices_end();

        // Define Half-edge circulator ...
        Mesh::Halfedge_around_vertex_circulator h, h_end;
        Mesh::Halfedge current_h;
        Mesh::Edge current_e;

        // Define current vertex value
        Mesh::Vertex current_v;

        // Defining usefull variables
        surface_mesh::Vec3 approx;
        float wi;
        float sum_wi;
        surface_mesh::Vec3 nLb;

        // Iterating over all vertices
        for (auto current_v: mesh_.vertices()){

            if(mesh_.is_boundary(current_v) == false){

                approx = surface_mesh::Vec3(0.0, 0.0, 0.0);
                wi=0.;
                sum_wi=0.;

                // initialisant le Halfedge circulator
                h = h_end = mesh_.halfedges(current_v);

                do {
                    current_h = *h;
                    current_e = mesh_.edge(current_h);
                    wi = e_weight[current_e];
                    sum_wi = sum_wi + wi;
                    approx = approx + wi*(mesh_.position(mesh_.to_vertex(current_h)) - mesh_.position(current_v)) ;
                } while(++h != h_end);

                nLb = approx / sum_wi;
                // computin the new position of each vertex;
                mesh_.position(current_v) = mesh_.position(current_v) +  nLb/2;

            }
        }
        mesh_.update_face_normals();
        mesh_.update_vertex_normals();
    }
}

void MeshProcessing::uniform_laplacian_enhance_feature(const unsigned int iterations,
                                                       const unsigned int coefficient) {
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    // Creating new property to store old position if it doesn't already exists
    Mesh::Vertex_property<Point> old_pos;
    if(mesh_.vertex_property<Point>("v:old")){
        std::cout <<"v:old already exists" << std::endl;
        old_pos = mesh_.vertex_property<Point>("v:old");
    } else {
        std::cout <<"v:old is not yet defined"<<std::endl;
        old_pos = mesh_.add_vertex_property<Point>("v:old");
    }
    // Storing old position in the old pos property
    for(auto v: mesh_.vertices()){
        old_pos[v] = mesh_.position(v);
    }

    // 1) perform uniform Laplacian smoothing for enhancement_smoothing_iterations iterations
    //Applying the smoothing
    MeshProcessing::uniform_smooth(iterations);

    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    surface_mesh::Vec3 v_in, v_out;
    for(auto v: mesh_.vertices()){
        v_in = old_pos[v];
        v_out = mesh_.position(v);
        mesh_.position(v) = v_out + coefficient * (v_in-v_out);
    }
}

void MeshProcessing::laplace_beltrami_enhance_feature(const unsigned int iterations,
                                                      const unsigned int coefficient) {
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    //We do exactly the same thing with smooth instead of uniform_smooth
    Mesh::Vertex_property<Point> old_pos;
    if(mesh_.vertex_property<Point>("v:old")){
        std::cout << "v:old already exists"<< std::endl;
        old_pos = mesh_.vertex_property<Point>("v:old");
    } else {
        std::cout << "v:old is not yet defined" << std::endl;
        old_pos = mesh_.add_vertex_property<Point>("v:old");
    }
    for(auto v: mesh_.vertices()){
        old_pos[v] = mesh_.position(v);
    }
    // 1) perform Laplace-Beltrami smoothing for enhancement_smoothing_iterations iterations
    MeshProcessing::smooth(iterations);

    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    surface_mesh::Vec3 v_in, v_out;
    for(auto v: mesh_.vertices()){
        v_in = old_pos[v];
        v_out = mesh_.position(v);
        mesh_.position(v) = v_out + coefficient * (v_in-v_out);
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

    for (auto e: mesh_.edges())
    {
        e_weight[e] = 0.0;

        h0 = mesh_.halfedge(e, 0);
        p0 = points[mesh_.to_vertex(h0)];

        h1 = mesh_.halfedge(e, 1);
        p1 = points[mesh_.to_vertex(h1)];

        if (!mesh_.is_boundary(h0))
        {
            h2 = mesh_.next_halfedge(h0);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }

        if (!mesh_.is_boundary(h1))
        {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
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

        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh_.vertices(*vf_c);

            const Point& P = mesh_.position(*fv_c);  ++fv_c;
            const Point& Q = mesh_.position(*fv_c);  ++fv_c;
            const Point& R = mesh_.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        } while(++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
}

void MeshProcessing::load_mesh(const string &filename) {
    if (!mesh_.read(filename)) {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

    cout << "Mesh "<< filename << " loaded." << endl;
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
    Mesh::Vertex_property<Point> vertex_normal =
            mesh_.vertex_property<Point>("v:normal");
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
    Mesh::Vertex_property<Color> v_color_valence =
            mesh_.vertex_property<Color>("v:color_valence",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_unicurvature =
            mesh_.vertex_property<Color>("v:color_unicurvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_curvature =
            mesh_.vertex_property<Color>("v:color_curvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_gaussian_curv =
            mesh_.vertex_property<Color>("v:color_gaussian_curv",
                                         Color(1.0f, 1.0f, 1.0f));

    Mesh::Vertex_property<Scalar> vertex_valence =
            mesh_.vertex_property<Scalar>("v:valence", 0.0f);
    for (auto v: mesh_.vertices()) {
        vertex_valence[v] = mesh_.valence(v);
    }

    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);

    calc_weights();
    calc_uniform_mean_curvature();
    calc_mean_curvature();
    calc_gauss_curvature();
    color_coding(vertex_valence, &mesh_, v_color_valence, 100 /* bound */);
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

    for(auto f: mesh_.faces()) {
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

void MeshProcessing::color_coding(Mesh::Vertex_property<Scalar> prop, Mesh *mesh,
                                  Mesh::Vertex_property<Color> color_prop, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    // discard upper and lower bound
    unsigned int n = values.size()-1;
    unsigned int i = n / bound;
    std::sort(values.begin(), values.end());
    Scalar min_value = values[i], max_value = values[n-1-i];

    // map values to colors
    for (auto v: mesh->vertices())
    {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
                               Mesh::Vertex_property<Color> color_prop)
{
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0/4.0 * (max_value - min_value);
    v1 = min_value + 1.0/4.0 * (max_value - min_value);
    v2 = min_value + 2.0/4.0 * (max_value - min_value);
    v3 = min_value + 3.0/4.0 * (max_value - min_value);
    v4 = min_value + 4.0/4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0) {
        col = Color(0, 0, 1);
    } else if (value > v4) {
        col = Color(1, 0, 0);
    } else if (value <= v2) {
        if (value <= v1) { // [v0, v1]
            Scalar u =  (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        } else { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1-u);
        }
    } else {
        if (value <= v3) { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        } else { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1-u, 0);
        }
    }
    return col;
}

MeshProcessing::~MeshProcessing() {}
}
