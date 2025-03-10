#ifndef GEOMETRY_PREPROCESS_H
#define GEOMETRY_PREPROCESS_H

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <queue>
#include <glm/glm.hpp>


// ------------------------------------utility functions for voxelization------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
#include "vector3d.hpp" // probably rewrite this later...

// flood fill algorithm
// it will set external empty voxels to 0.0f
// and internal voxels to 1.0f
// and boundary voxels to 2.0f
vector3d<float> imfill(const vector3d<float>& im)
{
	vector3d<float> filled(im.size());
	for (int i = 0; i < filled.mVec.size(); i++)
	{
		filled.mVec[i] = 1.0f;
	}

    // adding only one seed is not enough, at least test all the boundary voxels
	// const glm::ivec3 seed(0, 0, 0);
	// std::queue<glm::ivec3> q;
	// q.push(seed);
    std::queue<glm::ivec3> q;
    for(int i = 0; i < im.size().x; i++)
    for(int j = 0; j < im.size().y; j++)
    {
        q.push(glm::ivec3(i, j, 0));
        q.push(glm::ivec3(i, j, im.size().z-1));
    }
    for(int i = 0; i < im.size().x; i++)
    for(int k = 0; k < im.size().z; k++)
    {
        q.push(glm::ivec3(i, 0, k));
        q.push(glm::ivec3(i, im.size().y-1, k));
    }
    for(int j = 0; j < im.size().y; j++)
    for(int k = 0; k < im.size().z; k++)
    {
        q.push(glm::ivec3(0, j, k));
        q.push(glm::ivec3(im.size().x-1, j, k));
    }

	while (!q.empty())
	{
		glm::ivec3 p = q.front();
		q.pop();

		//neighbors
		const glm::ivec3 n[6] = { p + glm::ivec3(-1, 0, 0), p + glm::ivec3(0, -1, 0), p + glm::ivec3(0, 0, -1),
											p + glm::ivec3(+1, 0, 0), p + glm::ivec3(0, +1, 0), p + glm::ivec3(0, 0, +1) };

		if (im.get(p) == 0.0f && filled.get(p) == 1.0f)
		{
			filled.set(p, 0.0f);
			for (int i = 0; i < 6; i++)
			{
				if (filled.valid_index(n[i]))
				{
					q.push(n[i]);
				}
			}
		}
        else if(im.get(p) == 1.0f && filled.get(p) == 1.0f)
        {
            filled.set(p, 2.0f);
        }
	}
	return filled;
}

//triangle/voxel intersection (https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt)

bool planeBoxOverlap(glm::vec3 normal, glm::vec3 vert, glm::vec3 maxbox)
{
   float v;
   glm::vec3 vmin, vmax;
   for (int q = 0; q < 3; q++)
   {
      v = vert[q];					

      if (normal[q] > 0.0f)
      {
         vmin[q] = -maxbox[q] - v;	
         vmax[q] = maxbox[q] - v;	
      }
      else
      {
         vmin[q] = maxbox[q] - v;	
         vmax[q] = -maxbox[q] - v;
      }
   }
   if(glm::dot(normal, vmin) > 0.0f) return false;	
   if(glm::dot(normal, vmax) >= 0.0f) return true;

   return false;
}

/*======================== X-tests ========================*/

bool AXISTEST_X01(float a, float b, float fa, float fb, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 boxhalfsize)
{
   float p0 = a * v0.y - b * v0.z;			       	   
   float p2 = a * v2.y - b * v2.z;		
	       	   
   float min, max;
   if (p0 < p2) { min = p0; max = p2; }
   else { min = p2; max = p0; } 

   float rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   
   if (min > rad || max < -rad) return false;

   return true;
}


bool AXISTEST_X2(float a, float b, float fa, float fb, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 boxhalfsize)
{
   float p0 = a * v0.y - b * v0.z;			           
   float p1 = a * v1.y - b * v1.z;		
	       	   
   float min, max;
   if (p0 < p1) { min = p0; max = p1; }
   else { min = p1; max = p0; } 

   float rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   
   if (min > rad || max < -rad) return false;

   return true;
}


/*======================== Y-tests ========================*/

bool AXISTEST_Y02(float a, float b, float fa, float fb, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 boxhalfsize)
{
   float p0 = -a * v0.x + b * v0.z;		      	   
   float p2 = -a * v2.x + b * v2.z;	   
    	       	   
   float min, max;
   if (p0 < p2) { min = p0; max = p2; }
   else { min = p2; max = p0; } \

   float rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   
   if (min > rad || max < -rad) return false;

   return true;
}


bool AXISTEST_Y1(float a, float b, float fa, float fb, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 boxhalfsize)
{
   float p0 = -a * v0.x + b * v0.z;		      	   
   float p1 = -a * v1.x + b * v1.z;	     	       	   

   float min, max;
   if (p0 < p1) { min = p0; max = p1; }
   else { min = p1; max = p0; } \

   float rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   
   if (min > rad || max < -rad) return false;
   
   return true;
}


/*======================== Z-tests ========================*/


bool AXISTEST_Z12(float a, float b, float fa, float fb, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 boxhalfsize)
{		  
   float p1 = a * v1.x - b * v1.y;			          
   float p2 = a * v2.x - b * v2.y;		
	      
   float min, max;
   if (p2 < p1) { min = p2; max = p1; }
   else { min = p1; max = p2; } 

   float rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   
   if (min > rad || max < -rad) return false;

   return true;
}


bool AXISTEST_Z0(float a, float b, float fa, float fb, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 boxhalfsize)
{		   
   float p0 = a * v0.x - b * v0.y;				   
   float p1 = a * v1.x - b * v1.y;			           

   float min, max;
   if (p0 < p1) { min = p0; max = p1; }
   else { min = p1; max = p0; } 

   float rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   
   if (min > rad || max < -rad) return false;

   return true;
}

void FINDMINMAX(float x0, float x1, float x2, float& min, float& max) 
{
   min = max = x0;   
   if (x1 < min) min = x1; 
   if (x1 > max) max = x1; 
   if (x2 < min) min = x2; 
   if (x2 > max) max = x2;
}


bool triBoxOverlap(glm::vec3 boxcenter, glm::vec3 boxhalfsize, glm::vec3 triverts[3])
{

   /*    use separating axis theorem to test overlap between triangle and box */
   /*    need to test for overlap in these directions: */
   /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
   /*       we do not even need to test these) */
   /*    2) normal of the triangle */
   /*    3) crossproduct(edge from tri, {x,y,z}-direction) */
   /*       this gives 3x3=9 more tests */

   glm::vec3 v0, v1, v2;
   float min, max, fex, fey, fez;		// -NJMP- "d" local variable removed
   glm::vec3 normal, e0, e1, e2;

   //move everything so that the boxcenter is at (0,0,0) 
   v0 = triverts[0]-boxcenter;
   v1 = triverts[1]-boxcenter;
   v2 = triverts[2]-boxcenter;

   /* compute triangle edges */
   e0=v1-v0;      /* tri edge 0 */
   e1=v2-v1;      /* tri edge 1 */
   e2=v0-v2;      /* tri edge 2 */

   /* Bullet 3:  */

   //  test the 9 tests first (this was faster) 

   fex = fabsf(e0.x);
   fey = fabsf(e0.y);
   fez = fabsf(e0.z);

   if(AXISTEST_X01(e0.z, e0.y, fez, fey, v0, v1, v2, boxhalfsize) == false) return false;
   if(AXISTEST_Y02(e0.z, e0.x, fez, fex, v0, v1, v2, boxhalfsize) == false) return false;
   if(AXISTEST_Z12(e0.y, e0.x, fey, fex, v0, v1, v2, boxhalfsize) == false) return false;


   fex = fabsf(e1.x);
   fey = fabsf(e1.y);
   fez = fabsf(e1.z);

   if(AXISTEST_X01(e1.z, e1.y, fez, fey, v0, v1, v2, boxhalfsize) == false) return false;
   if(AXISTEST_Y02(e1.z, e1.x, fez, fex, v0, v1, v2, boxhalfsize) == false) return false;
   if(AXISTEST_Z0(e1.y, e1.x, fey, fex, v0, v1, v2, boxhalfsize) == false) return false;


   fex = fabsf(e2.x);
   fey = fabsf(e2.y);
   fez = fabsf(e2.z);

   if(AXISTEST_X2(e2.z, e2.y, fez, fey, v0, v1, v2, boxhalfsize) == false) return false;
   if(AXISTEST_Y1(e2.z, e2.x, fez, fex, v0, v1, v2, boxhalfsize) == false) return false;
   if(AXISTEST_Z12(e2.y, e2.x, fey, fex, v0, v1, v2, boxhalfsize) == false) return false;


   /* Bullet 1: */
   /*  first test overlap in the {x,y,z}-directions */
   /*  find min, max of the triangle each direction, and test for overlap in */
   /*  that direction -- this is equivalent to testing a minimal AABB around */
   /*  the triangle against the AABB */

   /* test in X-direction */

   FINDMINMAX(v0.x, v1.x, v2.x, min, max);
   if (min > boxhalfsize.x || max < -boxhalfsize.x) return false;

   /* test in Y-direction */
   FINDMINMAX(v0.y, v1.y, v2.y, min, max);
   if (min > boxhalfsize.y || max < -boxhalfsize.y) return false;

   /* test in Z-direction */
   FINDMINMAX(v0.z, v1.z, v2.z, min, max);
   if (min > boxhalfsize.z || max < -boxhalfsize.z) return false;



   /* Bullet 2: */
   /*  test if the box intersects the plane of the triangle */
   /*  compute plane equation of triangle: normal*x+d=0 */

   normal = glm::cross(e0, e1);

   if (!planeBoxOverlap(normal, v0, boxhalfsize)) return false;

   return true;   /* box and triangle overlaps */

}

// -----------------------------utility functions for generating embedded mesh-------------------------------
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

// ---------------------------------necessary data structures------------------------------------------------
// a struct to store the position and normal of a vertex
struct pos_norm{
    glm::vec3 pos;
    glm::vec3 norm;
};
typedef std::array<pos_norm, 3> pos_norm_Triangle;
// a struct for edge definition, used in topological data structure
struct Edge {
    glm::vec3 start;
    glm::vec3 end;
}; 
// simple ray definition
// the direction is assumed to be normalized
struct Ray {
    glm::vec3 origin;
    glm::vec3 direction;
};
// ----------------------------------------------------------------------------------------------------------

// ---------------------------------for clipping triangles---------------------------------------------------
// a simple function to compare two float numbers with torlerance
bool almostEqual(const glm::vec3& a, const glm::vec3& b, float epsilon = 1e-5f) {
    return glm::length(a - b) < epsilon;
}

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

// Sutherland–Hodgman algorithm：use a clip plane to clip a polygon
// the input function is the inside function and intersect function for the clip plane
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
std::vector<pos_norm_Triangle> triangulatePolygon(const std::vector<pos_norm>& poly) {
    std::vector<pos_norm_Triangle> triangles;
    if (poly.size() < 3)
        return triangles;
    // split the polygon into triangles by connecting the first vertex with each pair of adjacent vertices
    for (size_t i = 1; i < poly.size() - 1; ++i) {
        triangles.push_back({ poly[0], poly[i], poly[i+1] });
    }
    return triangles;
}

// clip a triangle to the unit cube [0, 1]^3
std::vector<pos_norm_Triangle> clipTriangleToUnitCube(const pos_norm& a, const pos_norm& b, const pos_norm& c) {
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
        return std::vector<pos_norm_Triangle>();

    // fan triangulation of the polygon
    return triangulatePolygon(poly);
}



// ---------------------------------for clipping voxels------------------------------------------------------

// classic Möller–Trumbore algorithm for ray-triangle intersection
// input：
//   ray      - the ray to be tested
//   triangle - triangle structure std::array<pos_norm, 3> TODO: probably need uv attribute in the future
//   t        - the distance from the ray origin to the intersection point（output parameter）
//   u, v     - the interpolation weights of the intersection point on the triangle uv. Since we don't have uv attribute, we don't use it here
//   normal   - the normal of the triangle at the intersection point. Still, not needed for this project
// return true if the ray intersects the triangle, false otherwise
bool rayIntersectsTriangle(const Ray& ray,
                           const pos_norm_Triangle& triangle,
                           float &t, float &u, float &v, glm::vec3& normal)
{
    const float EPSILON = 1e-6f;
    glm::vec3 v0 = triangle[0].pos;
    glm::vec3 v1 = triangle[1].pos;
    glm::vec3 v2 = triangle[2].pos;
    glm::vec3 n0 = triangle[0].norm;
    glm::vec3 n1 = triangle[1].norm;
    glm::vec3 n2 = triangle[2].norm;

    glm::vec3 E1, E2, S1, S, S2;
    glm::vec3 origin = ray.origin;
    glm::vec3 D = ray.direction;
    float S1E1, inv_S1E1, b1, b2;

    E1 = v1 - v0;
    E2 = v2 - v0;
    S1 = glm::cross(D, E2);
    S1E1 = glm::dot(E1, S1);
    //std::cout << "hit triangle" << std::endl;
    if (S1E1 > -EPSILON && S1E1 < EPSILON)
        return false;    // ray is parallel to the triangle

    inv_S1E1 = 1.0f / S1E1;
    S = origin - v0;
    b1 = inv_S1E1 * glm::dot(S, S1);

    if (b1 < 0.0 || b1 > 1.0)
        return false;

    S2 = glm::cross(S, E1);
    b2 = inv_S1E1 * glm::dot(D, S2);

    if (b2 < 0.0 || b1 + b2 > 1.0)
        return false;

    // compute t to find intersection point on the triangle
    t = inv_S1E1 * glm::dot(E2, S2);
    
    // determine whether use point normal or face normal
    float b0 = 1 - b1 - b2;
    normal = glm::normalize(b0 * n0 + b1 * n1 + b2 * n2);


    
    // triangle's UV coordinate (currently not used because no UV attribute yet)
    glm::vec2 uv0, uv1, uv2;
    u = b0 * uv0.x + b1 * uv1.x + b2 * uv2.x;
    v = b0 * uv0.y + b1 * uv1.y + b2 * uv2.y;

    return true;
}


// input a set of triangles and a ray, return the intersection points of the ray with all triangles
// note that the intersection points are not sorted, so you may need to sort them by distance to the ray origin if needed
std::vector<glm::vec3> raycastTriangles(const std::vector<pos_norm_Triangle>& triangles,
                                        const Ray& ray)
{
    std::vector<glm::vec3> intersections;
    float t, u, v;
    glm::vec3 normal;
    for (const auto& tri : triangles) {
        if (rayIntersectsTriangle(ray, tri, t, u, v, normal)) {
            // calculate the intersection point
            glm::vec3 intersectionPoint = ray.origin + ray.direction * t;
            intersections.push_back(intersectionPoint);
        }
    }
    return intersections;
}

// 1D line clipping algorithm
// segments: a list of values representing line segments, each two values represent a line segment
// the length of segments should be even, if not, the last value will be duplicated as a degenerated line segment
// clipMin, clipMax: the range to clip the line segments into
// return a list of values representing the clipped line segments
std::vector<float> clipLineSegments(const std::vector<float>& segments, float clipMin, float clipMax) {
    std::vector<float> clipped;
    size_t n = segments.size();
    
    // deal with each line segment (every two points)
    for (size_t i = 0; i + 1 < n; i += 2) {
        float segStart = segments[i];
        float segEnd   = segments[i + 1];
        // calculate the intersection of the line segment with the clip range
        float newStart = std::max(segStart, clipMin);
        float newEnd   = std::min(segEnd, clipMax);
        if (newStart < newEnd) {
            clipped.push_back(newStart);
            clipped.push_back(newEnd);
        }
    }
    
    // if the number of points is odd, then the last point will be duplicated as a degenerated line segment
    if (n % 2 == 1) {
        float lastPoint = segments[n - 1];
        float clamped = std::min(std::max(lastPoint, clipMin), clipMax);
        if (clamped >= clipMin && clamped <= clipMax) {
            clipped.push_back(clamped);
            clipped.push_back(clamped);
        }
    }
    
    return clipped;
}



// generate a closed polygon from a set of edges
// return a list of polygons, and each polygon is closed, meaning the last point is the same as the first point
// each polygon looks like: [p0, p1, p2, ..., pN, p0]
std::vector<std::vector<glm::vec3>> buildAllClosedPolygons(const std::vector<Edge>& inputEdges) {
    // make a mutable copy
    std::vector<Edge> remainingEdges = inputEdges;
    std::vector<std::vector<glm::vec3>> polygons;
    
    // it will by the end process all edges and try best to form closed polygons as many as possible
    while (!remainingEdges.empty()) {
        // select the first edge in the remaining edges as the starting edge of the polygon
        Edge currentEdge = remainingEdges.front();
        remainingEdges.erase(remainingEdges.begin());
        
        std::vector<glm::vec3> poly;
        poly.push_back(currentEdge.start);
        poly.push_back(currentEdge.end);
        glm::vec3 currentPoint = currentEdge.end;
        
        bool closed = false;
        bool progress = true;
        
        // expand the polygon by finding the edge that connects to the current point
        while (progress && !closed) {
            progress = false;
            // iterate remainingEdges to find edge that connect to currentPoint
            for (auto it = remainingEdges.begin(); it != remainingEdges.end(); ++it) {
                if (almostEqual(it->start, currentPoint)) {
                    currentPoint = it->end;
                    poly.push_back(currentPoint);
                    remainingEdges.erase(it);
                    progress = true;
                    break;
                }
                // check the other side of the edge
                else if (almostEqual(it->end, currentPoint)) {
                    currentPoint = it->start;
                    poly.push_back(currentPoint);
                    remainingEdges.erase(it);
                    progress = true;
                    break;
                }
            }
            // if currentPoint is the same as the first point, then the polygon is closed
            if (almostEqual(poly.front(), currentPoint)) {
                closed = true;
                // add the first point to the end to make it closed (actually not necessary)
                poly.back() = poly.front();
            }
        }
        
        if (closed) {
            polygons.push_back(poly);
        }
        else {
            // deal with open chain
            // either ignore it
            // or force close it
            // poly.push_back(poly.front());
            // polygons.push_back(poly);
            // we don't want such thing happen, but unfortunately it happens a lot, and outputing a warning is not helpful
            // std::cerr << "Warning: found an open chain that does not form a closed polygon." << std::endl;
        }
    }
    
    return polygons;
}


// determine if a point is inside a triangle in 2D
// utility function of ear clipping algorithm
bool pointInTriangle2D(const glm::vec2& pt, const glm::vec2& A, const glm::vec2& B, const glm::vec2& C) {
    glm::vec2 v0 = C - A;
    glm::vec2 v1 = B - A;
    glm::vec2 v2 = pt - A;
    float dot00 = glm::dot(v0, v0);
    float dot01 = glm::dot(v0, v1);
    float dot02 = glm::dot(v0, v2);
    float dot11 = glm::dot(v1, v1);
    float dot12 = glm::dot(v1, v2);
    float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    return (u >= 0.0f) && (v >= 0.0f) && (u + v <= 1.0f);
}

// according to the faceDirection, project a 3D polygon to 2D
// faceDirection = 0: projection discards x, return (y,z)
// faceDirection = 1: projection discards y, return (x,z)
// faceDirection = 2: projection discards z, return (x,y)
std::vector<glm::vec2> projectPolygonTo2D_fixed(const std::vector<glm::vec3>& poly, int faceDirection) {
    std::vector<glm::vec2> poly2d;
    poly2d.resize(poly.size());
    for (size_t i = 0; i < poly.size(); i++) {
        if (faceDirection == 0)
            poly2d[i] = glm::vec2(poly[i].y, poly[i].z);
        else if (faceDirection == 1)
            poly2d[i] = glm::vec2(poly[i].x, poly[i].z);
        else // faceDirection == 2
            poly2d[i] = glm::vec2(poly[i].x, poly[i].y);
    }
    return poly2d;
}

// calculate the oriented area of a 2D polygon (positive value means the vertices are in counter-clockwise order)
float polygonArea2D(const std::vector<glm::vec2>& poly2d) {
    float area = 0.0f;
    int n = poly2d.size();
    for (int i = 0; i < n; i++) {
        int next = (i + 1) % n;
        area += poly2d[i].x * poly2d[next].y - poly2d[next].x * poly2d[i].y;
    }
    return area * 0.5f;
}



// correct the orientation of a triangle to make sure its computed normal is in the same direction as the desired normal
// if not, swap the second and third vertices
// this is needed because our preprocess may generate triangles in different orientations
// while rendering usually requires the triangles to be in counter-clockwise order
std::array<glm::vec3, 3> correctTriangleOrientation(
    const std::array<glm::vec3, 3>& tri, 
    const glm::vec3& desiredNormal)
{
    // calculate the computed normal of the triangle
    glm::vec3 computedNormal = glm::normalize(glm::cross(tri[1] - tri[0], tri[2] - tri[0]));
    glm::vec3 normDesired = glm::normalize(desiredNormal);
    
    // if the computed normal is in the opposite direction of the desired normal, swap the second and third vertices
    if(glm::dot(computedNormal, normDesired) < 0.0f) {
        return { tri[0], tri[2], tri[1] };
    }
    return tri;
}





// using ear clipping algorithm to triangulate a polygon that lies on a fixed axis plane
// the input polygon should be in 3D, and all vertices should lie on the same plane defined by faceDirection
// the faceDirection is 0, 1, or 2, indicating the fixed axis of the plane (e.g., x plane if faceDirection==0)
// return a list of triangles, each triangle is a std::array<glm::vec3,3>
std::vector<std::array<glm::vec3, 3>> triangulatePolygonFixed(std::vector<glm::vec3>& polygon, int faceDirection) {
    std::vector<std::array<glm::vec3, 3>> triangles;
    int n = polygon.size();
    if (n < 3)
        return triangles;
    
    // project the 3D polygon to 2D (assuming the polygon is counter-clockwise oriented)
    std::vector<glm::vec2> poly2d = projectPolygonTo2D_fixed(polygon, faceDirection);
    // the ear clipping algorithm requires the polygon to be counter-clockwise oriented
    // if you don't do so, the result may be terrible!!!
    if (polygonArea2D(poly2d) < 0.0f) {
        std::reverse(poly2d.begin(), poly2d.end());
        std::reverse(polygon.begin(), polygon.end());
    }
    
    // use indices to represent the current polygon
    std::vector<int> indices(n);
    for (int i = 0; i < n; i++) {
        indices[i] = i;
    }
    
    auto isConvex = [&](int i0, int i1, int i2) -> bool {
        const glm::vec2 &A = poly2d[i0], &B = poly2d[i1], &C = poly2d[i2];
        float cross = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);
        return cross > 0;
    };
    
    while (indices.size() > 3) {
        bool earFound = false;
        int m = indices.size();
        for (int i = 0; i < m; i++) {
            int prev = indices[(i + m - 1) % m];
            int curr = indices[i];
            int next = indices[(i + 1) % m];
            
            if (!isConvex(prev, curr, next))
                continue;
            
            bool hasInside = false;
            for (int j = 0; j < m; j++) {
                if (j == (i + m - 1) % m || j == i || j == (i + 1) % m)
                    continue;
                int idx = indices[j];
                if (pointInTriangle2D(poly2d[idx], poly2d[prev], poly2d[curr], poly2d[next])) {
                    hasInside = true;
                    break;
                }
            }
            if (hasInside)
                continue;
            
            // find an ear, add the triangle (using the original 3D vertices)
            triangles.push_back({ polygon[prev], polygon[curr], polygon[next] });
            indices.erase(indices.begin() + i);
            earFound = true;
            break;
        }
        if (!earFound)
            break;
    }
    
    // if there are only 3 vertices left, add the last triangle
    // the remaining 3 vertices form the last triangle
    if (indices.size() == 3) {
        triangles.push_back({ polygon[indices[0]], polygon[indices[1]], polygon[indices[2]] });
    }
    
    return triangles;
}






#endif // GEOMETRY_PREPROCESS_H
