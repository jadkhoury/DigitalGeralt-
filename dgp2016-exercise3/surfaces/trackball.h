//=================================================================================
//   Code framework for the lecture "Digital 3D Geometry Processing"
//   Copyright Gaspard Zoss (C) 2016 Computer Graphics and Geometry Laboratory, EPFL
//----------------------------------------------------------------------------------
#pragma once

#include <Eigen/Dense>

using namespace Eigen;

class Trackball {
public:
    Trackball() : radius_(1.0f) {}

    void BeingDrag(float x, float y) {
        anchor_pos_ = Vector3f(x, y, 0.0f);
        ProjectOntoSurface(anchor_pos_);
    }

    void beginZoom(float x, float y) {
        old_y = y;
    }

    void beginPan(float x, float y) {
        old_x = x;
        old_y = y;
    }

    Matrix4f pan(float x, float y) {
        const float dx = x - old_x;
        const float dy = y - old_y;
        const float scale = 0.5f;
        Matrix4f translation = Eigen::Affine3f(Eigen::Translation3f(scale * dx, scale * dy, 0.0f)).matrix();
        old_x = x;
        old_y = y;
        return translation;
    }

    Matrix4f zoom(float x, float y) {
        const float dy = y - old_y;
        const float scale = 0.5f;
        Matrix4f translation = Eigen::Affine3f(Eigen::Translation3f(0.0f, 0.0f, -scale * dy)).matrix();
        old_y = y;
        return translation;
    }

    Matrix4f Drag(float x, float y) {
        Vector3f current_pos(x, y, 0.0f);
        ProjectOntoSurface(current_pos);

        Matrix4f rotation = Matrix4f::Identity();
        if (current_pos == anchor_pos_) {
            return rotation;
        }

        const float angle_boost = 2.0f;
        Vector3f axis = anchor_pos_.cross(current_pos);
        float angle = angle_boost * atan2(axis.norm(), anchor_pos_.dot(current_pos));
        rotation = Affine3f(AngleAxisf(angle, axis.normalized())).matrix();

        return rotation;
    }

private:
    void ProjectOntoSurface(Vector3f& p) const {
        const float rad2 = radius_ * radius_;
        const float length2 = p.x() * p.x() + p.y() * p.y();

        if (length2 <= rad2 * 0.5f) {
            p.z() = std::sqrt(rad2 - length2);
        } else {
            p.z() = rad2 / (2.0f * std::sqrt(length2));
        }
        float length = std::sqrt(length2 + p.z() * p.z());
        p /= length;
    }

    float radius_;
    Vector3f anchor_pos_;
    Matrix4f rotation_;
    float old_x;
    float old_y;
};
