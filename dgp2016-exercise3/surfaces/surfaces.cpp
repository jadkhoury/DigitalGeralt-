//=================================================================================
//   Code framework for the lecture "Digital 3D Geometry Processing"
//   Copyright Gaspard Zoss (C) 2016 Computer Graphics and Geometry Laboratory, EPFL
//----------------------------------------------------------------------------------
#include "viewer.h"

using namespace surface_mesh;

float computeAngle(Vec3 origin, Vec3 p1, Vec3 p2){
    Vec3 v1 = normalize(p1 - origin);
    Vec3 v2 = normalize(p2 - origin);
    return acos(dot(v1, v2));
}



void computeValence(Surface_mesh *mesh) {

    Surface_mesh::Vertex_property<unsigned int> vertex_valence =
            mesh->vertex_property<unsigned int>("v:valence", 0);

    // Vertex Iterator to iterate over all vertices
    Surface_mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh->vertices_begin();
    v_end = mesh->vertices_end();
    // Vertex circulator to iterate over the neighbours of the current vertex
    Surface_mesh::Vertex_around_vertex_circulator vc, vc_end;
    //Current vertex related variables
    Surface_mesh::Vertex current_v;
    int neighbors_counter;

    //Iterating over all the vertices
    for(v_it = v_begin; v_it != v_end; ++v_it){
        neighbors_counter = 0;
        current_v = *v_it;
        vc = mesh->vertices(current_v);
        vc_end = vc;
        //Iterating over all current_v neighbors
        do{
            neighbors_counter++;
        } while(++vc != vc_end);
        //Storing the computed valence
        vertex_valence[current_v] = neighbors_counter;
    }


}

//!\\ ONE LINE ADDED IN MAIN FUNCTION
void computeAllFacesNormals(Surface_mesh *mesh){

    // Creating the new property
    Surface_mesh::Face_property<Point> f_normal =
            mesh->add_face_property<Point>("f:normal");
    // Face iterator to iterate over all faces
    Surface_mesh::Face_iterator f_it, f_begin, f_end;
    f_begin = mesh->faces_begin();
    f_end = mesh->faces_end();
    // Vertex Circulator to iterate over the connected vertices to the current face
    Surface_mesh::Vertex_around_face_circulator vc, vc_end;
    // Currrent face variables
    Surface_mesh::Face current_f;
    Vec3 v[3];
    Vec3 v1, v2, normal;

    // Iterating over all faces
    for(f_it = f_begin; f_it != f_end; ++f_it){
        current_f = *f_it;
        vc = mesh->vertices(current_f);
        vc_end = vc;
        int ctr = 0;
        // Creating a small array of all the vertices arount the current face
        do{
            v[ctr] = mesh->position(*vc);
            ctr++;
        }while(++vc != vc_end);
        // Computing face normal, -1 needed for some reason
        v1 = v[1] - v[0];
        v2 = v[2] - v[0];
        normal = -cross(v2, v1);
        // Storing the face normal in the created property
        // !!!!! Normals are not normalized !!!! (on purpose)
        f_normal[current_f] = normal;
    }

}


void computeNormalsWithConstantWeights(Surface_mesh *mesh) {

    Point default_normal(0.0, 1.0, 0.0);
    // Preparing access to the vertex property and the faces normals
    Surface_mesh::Vertex_property<Point> v_cste_weights_n =
            mesh->vertex_property<Point>("v:cste_weights_n", default_normal);
    Surface_mesh::Face_property<Point> f_normal =
            mesh->face_property<Point>("f:normal", default_normal);
    // Vertex Iterator to iterate the current vertex over all vertices
    Surface_mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh->vertices_begin();
    v_end = mesh->vertices_end();
    // Face circulator to iterate over all connected faces to the current vertex
    Surface_mesh::Face_around_vertex_circulator fc, fc_end;
    Surface_mesh::Vertex current_v;
    // Current vertex variables
    Vec3 normal_vec;
    // Iterating over all vertices
    for(v_it = v_begin; v_it!=v_end; ++v_it){
        current_v = *v_it;
        // For each new vertex we set everything to 0
        normal_vec = Vec3(0.0, 0.0, 0.0);
        fc = mesh->faces(current_v);
        fc_end = fc;
        //Iterating over all connected faces to the current vertex
        do{
            // Summing the connected faces normals
            normal_vec += normalize(f_normal[*fc]);
        } while(++fc != fc_end);
        // Averaging the connected faces normals
        v_cste_weights_n[current_v] = normalize(normal_vec);
    }
    //Displaying time info

}

void computeNormalsByAreaWeights(Surface_mesh *mesh) {
    //We do exactly the same thing except we normalize the vectors after the average instead
    //of before such that the faces normal norms are taken into account

    Point default_normal(0.0, 1.0, 0.0);
    Surface_mesh::Vertex_property<Point> v_area_weights_n =
            mesh->vertex_property<Point>("v:area_weight_n", default_normal);
    Surface_mesh::Face_property<Point> f_normal =
            mesh->face_property<Point>("f:normal", default_normal);
    Surface_mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh->vertices_begin();
    v_end = mesh->vertices_end();
    Surface_mesh::Face_around_vertex_circulator fc, fc_end;
    Surface_mesh::Vertex current_v;
    Vec3 normal_vec;
    for(v_it = v_begin; v_it!=v_end; ++v_it){
        current_v = *v_it;
        normal_vec = Vec3(0.0, 0.0, 0.0);
        fc = mesh->faces(current_v);
        fc_end = fc;
        do{
            // change here: we don't normalize to have area information in the average
            normal_vec += f_normal[*fc];
        } while(++fc != fc_end);
        // normalizing after averaging
        v_area_weights_n[current_v] = normalize(normal_vec);
    }
    //Displaying time info

}

void computeNormalsWithAngleWeights(Surface_mesh *mesh) {

    Point default_normal(0.0, 1.0, 0.0);
    // Preparing vertex properties for writing
    Surface_mesh::Vertex_property<Point> v_angle_weights_n =
            mesh->vertex_property<Point>("v:angle_weight_n", default_normal);
    // Vertex iterator to iterate the current vertex over all vertices
    Surface_mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh->vertices_begin();
    v_end = mesh->vertices_end();
    // Vertex circulators, one added to access to consecutive neighbours
    Surface_mesh::Vertex_around_vertex_circulator vc, vc_next, vc_end;
    //Current vertex variables
    Surface_mesh::Vertex current_v;
    float angle;
    Vec3 current_pos, n1, n2, normal;
    // Iterating current_v over all vertices
    for(v_it = v_begin; v_it != v_end; ++v_it){
        // Update current vertex related variables
        normal = Vec3(0.0, 0.0, 0.0);
        current_v = *v_it;
        current_pos = mesh->position(current_v);
        vc = mesh->vertices(current_v);
        vc_next = vc;
        vc_end = vc;
        //Iterating over all current_v neighbors
        do{
            //vc_next will always be one step ahead vc to be its sucessive neighbour
            ++vc_next;
            //Computing angle
            n1 =  mesh->position(*vc);
            n2 =  mesh->position(*vc_next);
            angle = computeAngle(current_pos, n1, n2);
            //Computing the normal with the angles as weight of normalize norm, no??
            normal += angle * normalize(cross(n1-current_pos, n2-current_pos)) ;
        } while(++vc != vc_end);
        v_angle_weights_n[current_v] = normalize(normal);
    }

}

// #############################################################################
int main(int /* argc */, char ** /* argv */) {
    try {
        nanogui::init();
        {
            // Load the Mesh
            Surface_mesh mesh;
            if (!mesh.read("../data/bunny.off")) {
                std::cerr << "Mesh not found, exiting." << std::endl;
                return -1;
            }
            //!\\ LINE ADDED HERE //!\\

            computeValence(&mesh);
            computeAllFacesNormals(&mesh);
            computeNormalsWithConstantWeights(&mesh);
            computeNormalsByAreaWeights(&mesh);
            computeNormalsWithAngleWeights(&mesh);

            nanogui::ref<Viewer> app = new Viewer(&mesh);
            app->drawAll();
            app->setVisible(true);
            nanogui::mainloop();
        }
        nanogui::shutdown();
    } catch (const std::runtime_error &e) {
        std::string error_msg = std::string("Caught a fatal error: ") + std::string(e.what());
#if defined(_WIN32)
        MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
        std::cerr << error_msg << endl;
#endif
        return -1;
    }

    return 0;
}
