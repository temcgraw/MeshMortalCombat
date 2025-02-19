#ifndef GEOMETRY_PREPROCESS_H
#define GEOMETRY_PREPROCESS_H

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <glm/glm.hpp>

// a struct to store the position and normal of a vertex
struct pos_norm{
    glm::vec3 pos;
    glm::vec3 norm;
};


// bilinear interpolation, interpolate both position and normal
// normal is normalized first before interpolation
pos_norm Lerp(const pos_norm& v0, const pos_norm& v1, float t) {
    pos_norm result;
    result.pos = v0.pos + (v1.pos - v0.pos) * t;
    result.norm = glm::normalize(v0.norm + (v1.norm - v0.norm) * t);
    return result;
}

// define two function types for clipping
typedef bool (*InsideFunc)(const pos_norm& v);// tell if a vertex is inside the clip plane
typedef pos_norm (*IntersectFunc)(const pos_norm& v1, const pos_norm& v2);// given two vertices that are on different sides of the clip plane, return the intersection point

// Sutherland–Hodgman 算法：用一个裁剪平面对多边形进行裁剪
void ClipPolygon(std::vector<pos_norm>& poly, InsideFunc inside, IntersectFunc intersect) {
    std::vector<pos_norm> input = poly;
    poly.clear();
    if (input.empty())
        return;
    pos_norm S = input.back();
    bool S_inside = inside(S);
    for (const pos_norm& E : input) {
        bool E_inside = inside(E);
        if (E_inside) {
            if (!S_inside)
                poly.push_back(intersect(S, E));
            poly.push_back(E);
        } else {
            if (S_inside)
                poly.push_back(intersect(S, E));
        }
        S = E;
        S_inside = E_inside;
    }
}

// define inside and intersect functions for the 6 clip planes

// x = 0 plane: inside condition is v.pos.x >= 0
bool inside_xmin(const pos_norm& v) { return v.pos.x >= 0.0f; }
pos_norm intersect_xmin(const pos_norm& v1, const pos_norm& v2) {
    float t = (0.0f - v1.pos.x) / (v2.pos.x - v1.pos.x);
    return Lerp(v1, v2, t);
}

// x = 1 plane: inside condition is v.pos.x <= 1
bool inside_xmax(const pos_norm& v) { return v.pos.x <= 1.0f; }
pos_norm intersect_xmax(const pos_norm& v1, const pos_norm& v2) {
    float t = (1.0f - v1.pos.x) / (v2.pos.x - v1.pos.x);
    return Lerp(v1, v2, t);
}

// y = 0 plane: inside condition is v.pos.y >= 0
bool inside_ymin(const pos_norm& v) { return v.pos.y >= 0.0f; }
pos_norm intersect_ymin(const pos_norm& v1, const pos_norm& v2) {
    float t = (0.0f - v1.pos.y) / (v2.pos.y - v1.pos.y);
    return Lerp(v1, v2, t);
}

// y = 1 plane: inside condition is v.pos.y <= 1
bool inside_ymax(const pos_norm& v) { return v.pos.y <= 1.0f; }
pos_norm intersect_ymax(const pos_norm& v1, const pos_norm& v2) {
    float t = (1.0f - v1.pos.y) / (v2.pos.y - v1.pos.y);
    return Lerp(v1, v2, t);
}

// z = 0 plane: inside condition is v.pos.z >= 0
bool inside_zmin(const pos_norm& v) { return v.pos.z >= 0.0f; }
pos_norm intersect_zmin(const pos_norm& v1, const pos_norm& v2) {
    float t = (0.0f - v1.pos.z) / (v2.pos.z - v1.pos.z);
    return Lerp(v1, v2, t);
}

// z = 1 plane: inside condition is v.pos.z <= 1
bool inside_zmax(const pos_norm& v) { return v.pos.z <= 1.0f; }
pos_norm intersect_zmax(const pos_norm& v1, const pos_norm& v2) {
    float t = (1.0f - v1.pos.z) / (v2.pos.z - v1.pos.z);
    return Lerp(v1, v2, t);
}

// fan triangulation of the polygon
// this function assumes the input polygon is convex
// since the input triangle is convex, the output polygon is also convex, thus no worry about concave polygon
std::vector<std::array<pos_norm, 3>> triangulatePolygon(const std::vector<pos_norm>& poly) {
    std::vector<std::array<pos_norm, 3>> triangles;
    if (poly.size() < 3)
        return triangles;
    // split the polygon into triangles by connecting the first vertex with each pair of adjacent vertices
    for (size_t i = 1; i < poly.size() - 1; ++i) {
        triangles.push_back({ poly[0], poly[i], poly[i+1] });
    }
    return triangles;
}

// clip a triangle to the unit cube [0, 1]^3
std::vector<std::array<pos_norm, 3>> clipTriangleToUnitCube(const pos_norm& a, const pos_norm& b, const pos_norm& c) {
    std::vector<pos_norm> poly = { a, b, c };

    // sequentially clip the triangle to the 6 clip planes
    ClipPolygon(poly, inside_xmin, intersect_xmin);
    ClipPolygon(poly, inside_xmax, intersect_xmax);
    ClipPolygon(poly, inside_ymin, intersect_ymin);
    ClipPolygon(poly, inside_ymax, intersect_ymax);
    ClipPolygon(poly, inside_zmin, intersect_zmin);
    ClipPolygon(poly, inside_zmax, intersect_zmax);

    // if the polygon is empty, return empty triangles
    // could happen if the triangle is completely outside the unit cube and only some vertices are on the clip planes
    if (poly.size() < 3)
        return std::vector<std::array<pos_norm, 3>>();

    // fan triangulation of the polygon
    return triangulatePolygon(poly);
}



#endif // GEOMETRY_PREPROCESS_H
