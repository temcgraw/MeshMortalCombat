#ifndef VOXELGENERATOR_H
#define VOXELGENERATOR_H


#include "GraphicObject.h"
#include "CommonSceneObject.h"



#include <unordered_set>

#include "GeometryPreprocess.h"

#include "vector3d.hpp"



// one important thing: the voxelization resolution should be the same in all three dimensions
// because we want the voxels to be cubic instead of rectangular parallelepiped
// and actually, making the resolution different in different dimensions is meaningless
class VoxelGeneratorComponent : public IComponent {
public:
   // !!! the grid size refers to the voxelization resolution
   // it has nothing to do with the scene's data such as uniform grid size
   VoxelGeneratorComponent(GModel * _model,glm::ivec3 _grid_size = glm::ivec3(32,32,32)) {
      if (_model == nullptr) {
         std::cerr << "[VoxelGeneratorComponent]: model is nullptr" << std::endl;
         return;
      }
      model = _model;
      grid_size = _grid_size;
   }
    // it won't tick, because it does voxelization only during initialization
    // or when user orders to do so (should be implemented via callbacks)
   void Tick(float deltaTime) override {
      return;
   }
   void GenerateVoxelData() {
      if(model->meshes.size() == 0) {
         std::cerr << "[VoxelGeneratorComponent]: no mesh detected" << std::endl;
         return;
      }
      // check the resolution, make sure three dimensions are the same
      if (grid_size.x != grid_size.y || grid_size.x != grid_size.z) {
         std::cerr << "[VoxelGeneratorComponent]: grid size is not cubic" << std::endl;
         return;
      }
      
      // generate voxel data
      // need to first get the bounding box of the mesh
      // then generate the voxel data from a normalized version of the mesh
      glm::vec3 boundingBox_min = model->AA;
      glm::vec3 boundingBox_max = model->BB;
      glm::vec3 boundingBox_size = boundingBox_max - boundingBox_min;
      glm::vec3 boundingBox_center = (boundingBox_max + boundingBox_min) * 0.5f;
      // !!! we use the max axis to scale the voxel data, because in the end we want
      // voxels be cubic instead of rectangular parallelepiped
      float maxAxisLen = std::max(std::max(boundingBox_size.x, boundingBox_size.y), boundingBox_size.z);
      glm::vec3 scale = glm::vec3(maxAxisLen);

      //std::cout << "bounding box min: " << boundingBox_min.x << " " << boundingBox_min.y << " " << boundingBox_min.z << std::endl;
      //std::cout << "bounding box max: " << boundingBox_max.x << " " << boundingBox_max.y << " " << boundingBox_max.z << std::endl;
      //std::cout << "bounding box size: " << boundingBox_size.x << " " << boundingBox_size.y << " " << boundingBox_size.z << std::endl;
      //std::cout << "bounding box center: " << boundingBox_center.x << " " << boundingBox_center.y << " " << boundingBox_center.z << std::endl;

      mVoxels.resize(grid_size);
      mVoxels.fill(0.0f);

      // for each mesh:
      for(int m = 0; m < model->meshes.size(); m++){
         GMesh *mesh = model->meshes[m];
         
         
         int num_triangles = mesh->indices.size() / 3;


         int currentVoxelNum = 0;
         //for each triangle in mesh
         for (int t = 0; t < num_triangles; t++)
         {
            glm::uvec3 ix = glm::uvec3(mesh->indices[t * 3], mesh->indices[t * 3 + 1], mesh->indices[t * 3 + 2]);
            glm::vec4 po[3]; // the normalized vertex position between 0-1
            glm::vec3 pvox[3]; // vertex in voxel space, between 0 to grid_size

            glm::vec3 bbmin = glm::vec3(1e10f);
            glm::vec3 bbmax = glm::vec3(-1e10f);

            //for each vertex in triangle
            for (int v = 0; v < 3; v++)
            {
               po[v] = glm::vec4(mesh->positions[ix[v]], 1.0f);
               // normalize vertex position into 0-1 range
               po[v] = po[v] - glm::vec4(boundingBox_min, 0.0f);
               po[v] = po[v] / glm::vec4(scale, 1.0f);

               // convert normalized vertex position into voxel space where each dimension is an integer between 0 to grid_size-1
               // each voxel is a cubic 1x1x1 unit
               pvox[v] = glm::vec3(po[v]) * glm::vec3(grid_size);

               //update bbox
               bbmin = glm::min(bbmin, pvox[v]);
               bbmax = glm::max(bbmax, pvox[v]);

            
            }
            //for each voxel overlapped by bounding box
            for (int i = bbmin.x; i <= bbmax.x; i++)
            for (int j = bbmin.y; j <= bbmax.y; j++)
            for (int k = bbmin.z; k <= bbmax.z; k++)
            {
               if(mVoxels.valid_index(glm::vec3(i,j,k)))
               {
                  const float eps = 0.001f;
                  glm::vec3 cen = glm::vec3(i,j,k)+glm::vec3(0.5f);
                  const glm::vec3 half(0.5f + eps);
                  if(triBoxOverlap(cen, half, pvox))
                  {
                     mVoxels.set(i, j, k, 1.0f);
                  }
               }
            }
         }
      }
      //filled voxelization
      mVoxels = imfill(mVoxels);

   }

   // the most painful part:
   /*  
    It is not very easy to render the boundary voxels since they intersect with the model. 
    We not only need to clip the model's surface triangles within each voxel block, 
    but also clip each voxel block’s six faces against these surface triangles. 
    This involves traditional geometric processing, which can be quite verbose to implement.
   */
   void GenerateEmbeddedMeshSurfaceData(){
      // assume the voxel data has been generated
      if(model->meshes.size() == 0) {
         std::cerr << "[VoxelGeneratorComponent]: no mesh detected" << std::endl;
         return;
      }
      // check the resolution, make sure three dimensions are the same
      if (grid_size.x != grid_size.y || grid_size.x != grid_size.z) {
         std::cerr << "[VoxelGeneratorComponent]: grid size is not cubic" << std::endl;
         return;
      }

      
      glm::vec3 boundingBox_min = model->AA;
      glm::vec3 boundingBox_max = model->BB;
      glm::vec3 boundingBox_size = boundingBox_max - boundingBox_min;
      glm::vec3 boundingBox_center = (boundingBox_max + boundingBox_min) * 0.5f;
      float maxAxisLen = std::max(std::max(boundingBox_size.x, boundingBox_size.y), boundingBox_size.z);
      glm::vec3 scale = glm::vec3(maxAxisLen);

      // check whether the voxel data has been generated
      if(mVoxels.size().x == 0) {
         std::cerr << "[VoxelGeneratorComponent]: voxel data has not been generated" << std::endl;
         return;
      }
      embeddedSurfaceMeshData.resize(grid_size);
      embeddedVoxelMeshData.resize(grid_size);

      // for each mesh:
      for(int m = 0; m < model->meshes.size(); m++){
         GMesh *mesh = model->meshes[m];
         int num_triangles = mesh->indices.size() / 3;
         // ----------surface mesh part: clip the triangles to the voxel's bounding box----------
         // for each triangle in mesh
         for (int t = 0; t < num_triangles; t++)
         {
            glm::uvec3 ix = glm::uvec3(mesh->indices[t * 3], mesh->indices[t * 3 + 1], mesh->indices[t * 3 + 2]);
            glm::vec4 po[3]; // the normalized vertex position between 0-1
            glm::vec3 pn[3]; // the vertex normal, needed for shading so we need to store it
            glm::vec3 pvox[3]; // vertex in voxel space, between 0 to grid_size

            glm::vec3 bbmin = glm::vec3(1e10f);
            glm::vec3 bbmax = glm::vec3(-1e10f);

            //for each vertex in triangle
            for (int v = 0; v < 3; v++)
            {
               po[v] = glm::vec4(mesh->positions[ix[v]], 1.0f);
               // normalize vertex position into 0-1 range
               po[v] = po[v] - glm::vec4(boundingBox_min, 0.0f);
               po[v] = po[v] / glm::vec4(scale, 1.0f);

               // convert normalized vertex position into voxel space where each dimension is an integer between 0 to grid_size-1
               // each voxel is a cubic 1x1x1 unit
               pvox[v] = glm::vec3(po[v]) * glm::vec3(grid_size);

               pn[v] = mesh->normals[ix[v]];

               //update bbox
               bbmin = glm::min(bbmin, pvox[v]);
               bbmax = glm::max(bbmax, pvox[v]);

            
            }

            //for each voxel overlapped by bounding box
            for (int i = bbmin.x; i <= bbmax.x; i++)
            for (int j = bbmin.y; j <= bbmax.y; j++)
            for (int k = bbmin.z; k <= bbmax.z; k++)
            {
               if(mVoxels.valid_index(glm::vec3(i,j,k)))
               {
                  const float eps = 0.001f;
                  glm::vec3 cen = glm::vec3(i,j,k)+glm::vec3(0.5f);
                  const glm::vec3 half(0.5f + eps);
                  if(triBoxOverlap(cen, half, pvox))
                  {
                     // 1. convert the triangle to the target voxel's local space
                     // which ranging from 0,0,0 to 1,1,1
                     glm::vec3 localTriPos[3];
                     for (int v = 0; v < 3; v++) {
                        localTriPos[v] = pvox[v] - glm::vec3(i, j, k);
                     }

                     pos_norm localTri[3];
                     for (int v = 0; v < 3; v++) {
                        localTri[v] = {localTriPos[v], pn[v]};
                     }

                     std::vector<pos_norm> embeddedTriangleVertices;
                     // 2. clipping the triangle to the voxel's bounding box
                     std::vector<pos_norm_Triangle> clippedTriangles = clipTriangleToUnitCube(localTri[0], localTri[1], localTri[2]);
                     for(const auto& clippedTri : clippedTriangles) {
                        embeddedTriangleVertices.push_back(clippedTri[0]);
                        embeddedTriangleVertices.push_back(clippedTri[1]);
                        embeddedTriangleVertices.push_back(clippedTri[2]);
                     }
                     // or, only store the original triangle's vertices
                     // embeddedTriangleVertices.push_back(localTri[0]);
                     // embeddedTriangleVertices.push_back(localTri[1]);
                     // embeddedTriangleVertices.push_back(localTri[2]);

                     // 3. store the embedded triangle data into the voxel data structure
                     if(embeddedSurfaceMeshData.valid_index(glm::vec3(i,j,k))) {
                        // add the embedded triangle data to the corresponding voxel
                        auto& voxelData = embeddedSurfaceMeshData.getRef(i, j, k);
                        for (const auto& vertex : embeddedTriangleVertices) {
                           voxelData.push_back(vertex);
                        }
                     }
                     
                     


                  }
               }
            }
         }
         // ----------boundary voxel part: clip the voxels to the surface of the mesh----------
         // the hardest part...
         // 1. build a triangle uniform grid for fast triangle intersection test
         std::vector<pos_norm_Triangle> meshTrianglesData; // the flattened triangle raw data in voxel space
         for(int t = 0; t < num_triangles; t++) {
            glm::uvec3 ix = glm::uvec3(mesh->indices[t * 3], mesh->indices[t * 3 + 1], mesh->indices[t * 3 + 2]);
            glm::vec4 po[3]; // the normalized vertex position between 0-1
            glm::vec3 pn[3]; // the vertex normal, needed for shading so we need to store it
            glm::vec3 pvox[3]; // vertex in voxel space, between 0 to grid_size

            //for each vertex in triangle
            for (int v = 0; v < 3; v++)
            {
               po[v] = glm::vec4(mesh->positions[ix[v]], 1.0f);
               // normalize vertex position into 0-1 range
               po[v] = po[v] - glm::vec4(boundingBox_min, 0.0f);
               po[v] = po[v] / glm::vec4(scale, 1.0f);

               // convert normalized vertex position into voxel space where each dimension is an integer between 0 to grid_size-1
               // each voxel is a cubic 1x1x1 unit
               pvox[v] = glm::vec3(po[v]) * glm::vec3(grid_size);

               pn[v] = mesh->normals[ix[v]];

               
            }
            pos_norm_Triangle tri = { { {pvox[0], pn[0]}, {pvox[1], pn[1]}, {pvox[2], pn[2]} } };
            meshTrianglesData.push_back(tri);
         }

         // the triangle ugrid, each voxel stores the triangle indices in the meshTrianglesData
         vector3d<std::vector<int>> triangleGrid(grid_size); 
         for (int t = 0; t < num_triangles; t++)
         {
            glm::vec3 pvox[3] = {meshTrianglesData[t][0].pos, meshTrianglesData[t][1].pos, meshTrianglesData[t][2].pos};
            glm::vec3 pn[3] = {meshTrianglesData[t][0].norm, meshTrianglesData[t][1].norm, meshTrianglesData[t][2].norm};

            glm::vec3 bbmin = glm::vec3(1e10f);
            glm::vec3 bbmax = glm::vec3(-1e10f);

            //for each vertex in triangle
            for (int v = 0; v < 3; v++)
            {
            
               //update bbox
               bbmin = glm::min(bbmin, pvox[v]);
               bbmax = glm::max(bbmax, pvox[v]);

            }

            for (int i = bbmin.x; i <= bbmax.x; i++)
            for (int j = bbmin.y; j <= bbmax.y; j++)
            for (int k = bbmin.z; k <= bbmax.z; k++)
            {
               if(mVoxels.valid_index(glm::vec3(i,j,k)))
               {
                  const float eps = 0.001f;
                  glm::vec3 cen = glm::vec3(i,j,k)+glm::vec3(0.5f);
                  const glm::vec3 half(0.5f + eps);
                  if(triBoxOverlap(cen, half, pvox))
                  {
                     // add the triangle to the grid
                     triangleGrid.getRef(i, j, k).push_back(t);
                     
                  }
               }
            }
         }
         // 2. for each possible voxel edge direction, do raycasting to find all the intersection points
         // and then connect the intersection points to form the boundary mesh's edges

         // initialize the rays for all three directions, it covered all the possible voxel edges
         std::vector<Ray> raysXdir((grid_size.y+1)*(grid_size.z+1));
         std::vector<Ray> raysYdir((grid_size.x+1)*(grid_size.z+1));
         std::vector<Ray> raysZdir((grid_size.x+1)*(grid_size.y+1));
         // helper function to find all the voxels that a ray passes through
         auto findRayRelatedVoxels = [&](const Ray& ray, std::vector<glm::ivec3>& relatedVoxelIndices) {
            relatedVoxelIndices.clear();
            glm::vec3 shift[4];
            int axis = -1;
            // find the axis of the ray (X/Y/Z) and the possible shift of the ray
            if (ray.direction.x != 0.0f) {// x direction
               axis = 0;
               shift[0] = glm::vec3(0,0,0);shift[1] = glm::vec3(0,0,-1);shift[2] = glm::vec3(0,-1,0);shift[3] = glm::vec3(0,-1,-1);                                   
            }
            else if (ray.direction.y != 0.0f) { // y direction
               axis = 1;
               shift[0] = glm::vec3(0,0,0);shift[1] = glm::vec3(0,0,-1);shift[2] = glm::vec3(-1,0,0);shift[3] = glm::vec3(-1,0,-1);         
            }
            else if (ray.direction.z != 0.0f) {// z direction
               axis = 2;
               shift[0] = glm::vec3(0,0,0);shift[1] = glm::vec3(0,-1,0);shift[2] = glm::vec3(-1,0,0);shift[3] = glm::vec3(-1,-1,0);
            }
            else{
               std::cerr << "[VoxelGeneratorComponent]: ray direction is not along x, y or z axis" << std::endl;
               return;
            }
            // ensure the ray is at four corners of the every four voxels it passes through
            if(std::fmod(ray.origin[(axis+1)%3],1.0) != 0.0f || std::fmod(ray.origin[(axis+2)%3],1.0) != 0.0f) {
               std::cerr << "[VoxelGeneratorComponent]: ray origin is not at the corner of the voxel" << std::endl;
               return;
            }
            // ensure direction is unit vector and is positive 1
            if(ray.direction[axis] != 1.0f || abs(glm::length(ray.direction)-1) > 1e-6) {
               std::cerr << "[VoxelGeneratorComponent]: ray direction is not unit positive vector" << std::endl;
               return;
            }
            for(int i = 0; i < grid_size[axis]; i++) {
               for(int j = 0; j < 4; j++) {
                  glm::vec3 p = ray.origin + ray.direction * (float)i + shift[j];
                  glm::ivec3 idx = glm::ivec3(p);
                  if (idx.x >= 0 && idx.x < grid_size.x && idx.y >= 0 && idx.y < grid_size.y && idx.z >= 0 && idx.z < grid_size.z) {
                     relatedVoxelIndices.push_back(idx);
                  }
               }
            }
         };

         // initialize the rays, for three directions and for each face, every ray is at one of the corner of the voxels
         for (int i = 0; i < grid_size.x+1; i++)
         for (int j = 0; j < grid_size.y+1; j++)
         {
            int idx = i*(grid_size.y+1)+j;
            raysZdir[idx] = {glm::vec3(i, j, 0), glm::vec3(0, 0, 1)};
         }
         for (int i = 0; i < grid_size.x+1; i++)
         for (int j = 0; j < grid_size.z+1; j++)
         {
            int idx = i*(grid_size.z+1)+j;
            raysYdir[idx] = {glm::vec3(i, 0, j), glm::vec3(0, 1, 0)};
         }
         for (int i = 0; i < grid_size.y+1; i++)
         for (int j = 0; j < grid_size.z+1; j++)
         {
            int idx = i*(grid_size.z+1)+j;
            raysXdir[idx] = {glm::vec3(0, i, j), glm::vec3(1, 0, 0)};
         }

         // for each ray, do raycasting on related triangles to find all the intersection points
         std::vector<std::vector<glm::vec3>> intersectionPointsXdir((grid_size.y+1)*(grid_size.z+1));
         std::vector<std::vector<glm::vec3>> intersectionPointsYdir((grid_size.x+1)*(grid_size.z+1));
         std::vector<std::vector<glm::vec3>> intersectionPointsZdir((grid_size.x+1)*(grid_size.y+1));
         for (int i = 0; i < grid_size.y+1; i++)
         for (int j = 0; j < grid_size.z+1; j++)
         {
            int idx = i*(grid_size.z+1)+j;
            // for each ray, it might intersect with triangles in multiple voxels
            std::vector<glm::ivec3> relatedVoxelIndices;
            findRayRelatedVoxels(raysXdir[idx], relatedVoxelIndices);
            std::vector<pos_norm_Triangle> triangles(0);
            std::unordered_set<int> candidateTriangleIndices; // need to use set to avoid duplicate triangles
            for(const auto& voxelIdx : relatedVoxelIndices) {
               for(const auto& triIdx : triangleGrid.get(voxelIdx)) {
                  candidateTriangleIndices.insert(triIdx);
               }
            }
            for(const auto& triIdx : candidateTriangleIndices) {
               triangles.push_back(meshTrianglesData[triIdx]);
            }
            std::vector<glm::vec3> intersectionPoints =  raycastTriangles(triangles, raysXdir[idx]);
            intersectionPointsXdir[idx].insert(intersectionPointsXdir[idx].end(), intersectionPoints.begin(), intersectionPoints.end());
         }
         for (int i = 0; i < grid_size.x+1; i++)
         for (int j = 0; j < grid_size.z+1; j++)
         {
            int idx = i*(grid_size.z+1)+j;
            // for each ray, it might intersect with triangles in multiple voxels
            std::vector<glm::ivec3> relatedVoxelIndices;
            findRayRelatedVoxels(raysYdir[idx], relatedVoxelIndices);
            std::vector<pos_norm_Triangle> triangles(0);
            std::unordered_set<int> candidateTriangleIndices; // need to use set to avoid duplicate triangles
            for(const auto& voxelIdx : relatedVoxelIndices) {
               for(const auto& triIdx : triangleGrid.get(voxelIdx)) {
                  candidateTriangleIndices.insert(triIdx);
               }
            }
            for(const auto& triIdx : candidateTriangleIndices) {
               triangles.push_back(meshTrianglesData[triIdx]);
            }
            std::vector<glm::vec3> intersectionPoints =  raycastTriangles(triangles, raysYdir[idx]);
            intersectionPointsYdir[idx].insert(intersectionPointsYdir[idx].end(), intersectionPoints.begin(), intersectionPoints.end());
         }
         for (int i = 0; i < grid_size.x+1; i++)
         for (int j = 0; j < grid_size.y+1; j++)
         {
            int idx = i*(grid_size.y+1)+j;
            // for each ray, it might intersect with triangles in multiple voxels
            std::vector<glm::ivec3> relatedVoxelIndices;
            findRayRelatedVoxels(raysZdir[idx], relatedVoxelIndices);
            std::vector<pos_norm_Triangle> triangles(0);
            std::unordered_set<int> candidateTriangleIndices; // need to use set to avoid duplicate triangles
            for(const auto& voxelIdx : relatedVoxelIndices) {
               for(const auto& triIdx : triangleGrid.get(voxelIdx)) {
                  candidateTriangleIndices.insert(triIdx);
               }
            }
            for(const auto& triIdx : candidateTriangleIndices) {
               triangles.push_back(meshTrianglesData[triIdx]);
            }
            std::vector<glm::vec3> intersectionPoints =  raycastTriangles(triangles, raysZdir[idx]);
            intersectionPointsZdir[idx].insert(intersectionPointsZdir[idx].end(), intersectionPoints.begin(), intersectionPoints.end());
         }


         std::vector<std::vector<std::vector<glm::vec3>>> intersectionPoints(3);
         intersectionPoints[0] = intersectionPointsXdir;
         intersectionPoints[1] = intersectionPointsYdir;
         intersectionPoints[2] = intersectionPointsZdir;


         // 3. for each voxel, find all of its intersection edges and in the end form the boundary voxel mesh triangles

         // for simplicity define parameter set before the loop
         // for each voxel, we need to process its six faces
         // and for each face, we need to process four edges
         // so now define these parameters:
         std::vector<std::vector<std::vector<int>>> faceEdges(6);// outer vector is the face(6), inner vector is the edge(4), each edge is a vector <x_shift, y_shift, z_shift, direction>
         faceEdges[0] = {{0, 0, 0, 1}, {0, 0, 1, 1}, {0, 0, 0, 2}, {0, 1, 0, 2}};// x_min face: y0, y1, z0, z1
         faceEdges[1] = {{1, 0, 0, 1}, {1, 0, 1, 1}, {1, 0, 0, 2}, {1, 1, 0, 2}};// x_max face: y0, y1, z0, z1
         faceEdges[2] = {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 2}, {1, 0, 0, 2}};// y_min face: x0, x1, z0, z1
         faceEdges[3] = {{0, 1, 0, 0}, {0, 1, 1, 0}, {0, 1, 0, 2}, {1, 1, 0, 2}};// y_max face: x0, x1, z0, z1
         faceEdges[4] = {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}};// z_min face: x0, x1, y0, y1
         faceEdges[5] = {{0, 0, 1, 0}, {0, 1, 1, 0}, {0, 0, 1, 1}, {1, 0, 1, 1}};// z_max face: x0, x1, y0, y1
         std::vector<glm::vec3> faceNormals = {{-1.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {0.0f, -1.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 0.0f, 1.0f}};
         // raw index: {i,j,k}
         // idx of A direction =( index[indexCalculation[A].first]+shift[indexCalculation[A].first]) * (grid_size[indexCalculation[A].second]+1) + (index[indexCalculation[A].second]+shift[indexCalculation[A].second])
         // e.g.: x direction (direction = 0): int idx = j*(grid_size.z+1)+k;
         std::vector<std::pair<int, int>> indexCalculation = {{1, 2}, {0, 2}, {0, 1}};
         

         for (int i = 0; i < grid_size.x; i++)
         for (int j = 0; j < grid_size.y; j++)
         for (int k = 0; k < grid_size.z; k++)
         {
            // ensure this is a surface voxel
            if(embeddedSurfaceMeshData.get(i,j,k).size()==0){
               continue;
            }
            auto& voxelData = embeddedVoxelMeshData.getRef(i, j, k);
            std::vector<pos_norm> embeddedTriangleVertices[6];
            std::vector<Edge> embeddedEdges[6];
            std::vector<int> rawIndex = {i, j, k};
            for(int face = 0; face < 6; face++) {
               // first find the intersection edges from four boundary rays of the face
               for(int edge = 0; edge < 4; edge++) {
                  std::vector<int> xyz_shift = faceEdges[face][edge];// {x_shift, y_shift, z_shift, direction} we only need the first three
                  std::vector<int> shiftedIndex = {rawIndex[0]+xyz_shift[0], rawIndex[1]+xyz_shift[1], rawIndex[2]+xyz_shift[2]};
                  int direction = faceEdges[face][edge][3];
                  auto [rowId, colId] = indexCalculation[direction];
                  int rayIdx = (shiftedIndex[rowId])*(grid_size[colId]+1)+(shiftedIndex[colId]);
                  std::vector<float> inputSegments;
                  for(int i = 0; i < intersectionPoints[direction][rayIdx].size(); i++) {
                     inputSegments.push_back(intersectionPoints[direction][rayIdx][i][direction]);
                  }
                  std::sort(inputSegments.begin(), inputSegments.end());//!!!!!!!!!!very important
                  std::vector<float> outputSegments = clipLineSegments(inputSegments, rawIndex[direction], rawIndex[direction]+1);
                  for(int id = 0; id < outputSegments.size(); id+=2) {
                     glm::vec3 p0 = glm::vec3(shiftedIndex[0], shiftedIndex[1], shiftedIndex[2]);
                     p0[direction] = outputSegments[id];
                     p0 = p0 - glm::vec3(i, j, k);
                     glm::vec3 p1 = glm::vec3(shiftedIndex[0], shiftedIndex[1], shiftedIndex[2]);
                     p1[direction] = outputSegments[id+1];
                     p1 = p1 - glm::vec3(i, j, k);
                     embeddedEdges[face].push_back({p0, p1});
                  }             
               }
               // then find the intersection edges from the intersection triangles on the face
               auto& voxelSurfaceData = embeddedSurfaceMeshData.getRef(i, j, k);
               int faceDirection = face/2; // face 0,1 is x direction, face 2,3 is y direction, face 4,5 is z direction
               float faceDirectionAxisOffsetValue = (face%2==0)?0.0f:1.0f; // face 0,2,4 is min face, face 1,3,5 is max face
               // find all the intersection edges
               for(int i=0;i<voxelSurfaceData.size();i+=3){
                  glm::vec3 p[3] = {voxelSurfaceData[i].pos, voxelSurfaceData[i+1].pos, voxelSurfaceData[i+2].pos};
                  // check if any edge is on the face
                  for(int e = 0; e < 3; e++) {
                     Edge edge = {p[e], p[(e+1)%3]};
                     // if both vertices are on the face
                     if(abs(edge.start[faceDirection]-faceDirectionAxisOffsetValue)<1e-6 && abs(edge.end[faceDirection]-faceDirectionAxisOffsetValue)<1e-6) {
                        embeddedEdges[face].push_back(edge);
                     }
                  }
               }
               // connect the edges to form the triangles
               std::vector<std::vector<glm::vec3>> closedPolygons = buildAllClosedPolygons(embeddedEdges[face]);
               // triangulate the polygons
               for(auto& polygon : closedPolygons) {
                  // remove the last point if it is the same as the first point
                  // not doing so is not a problem, but it is better to remove it
                  if(polygon.back()==polygon.front()) {
                     polygon.pop_back();                
                  }
                  std::vector<std::array<glm::vec3, 3>> triangles = triangulatePolygonFixed(polygon, faceDirection);
                  for(const auto& tri : triangles) {
                     // ensure the triangle is in the correct direction
                     auto correctedTri = correctTriangleOrientation(tri, faceNormals[face]);
                     embeddedTriangleVertices[face].push_back({correctedTri[0], faceNormals[face]});
                     embeddedTriangleVertices[face].push_back({correctedTri[1], faceNormals[face]});
                     embeddedTriangleVertices[face].push_back({correctedTri[2], faceNormals[face]});
                  }
               }
               for(int id = 0; id < embeddedTriangleVertices[face].size(); id++) {
                  voxelData.push_back(embeddedTriangleVertices[face][id]);
               }
            }
         }
      }
   }
   

   // get a copy of the voxel data
   vector3d<float> getVoxels() {
      return mVoxels;
   }

   // get the embedded mesh data
   vector3d<std::vector<pos_norm>> getEmbeddedSurfaceMeshData() {
      return embeddedSurfaceMeshData;
   }

   vector3d<std::vector<pos_norm>> getEmbeddedVoxelMeshData() {
      return embeddedVoxelMeshData;
   }

private:
GModel * model;// the model to be voxelized, currently we only voxelized the model's first mesh
vector3d<float> mVoxels;
glm::ivec3 grid_size;
vector3d<std::vector<pos_norm>> embeddedSurfaceMeshData;
vector3d<std::vector<pos_norm>> embeddedVoxelMeshData; // temporary data structure, in future we should merge it with embeddedSurfaceMeshData


};




#endif // VOXELGENERATOR_H
