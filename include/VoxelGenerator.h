#ifndef VOXELGENERATOR_H
#define VOXELGENERATOR_H


#include "GraphicObject.h"
#include "CommonSceneObject.h"

#include "vector3d.hpp" // probably rewrite this later...

#include <queue>
vector3d<float> imfill(const vector3d<float>& im)
{
	vector3d<float> filled(im.size());
	for (int i = 0; i < filled.mVec.size(); i++)
	{
		filled.mVec[i] = 1.0f;
	}

	const glm::ivec3 seed(0, 0, 0);
	std::queue<glm::ivec3> q;
	q.push(seed);

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
      mesh = model->meshes[0];
      grid_size = _grid_size;
      GenerateVoxelData();
   }
    // it won't tick, because it does voxelization only during initialization
    // or when user orders to do so (should be implemented via callbacks)
   void Tick(float deltaTime) override {
      return;
   }
   void GenerateVoxelData() {
      if(mesh == nullptr) {
         std::cerr << "[VoxelGeneratorComponent]: mesh is nullptr" << std::endl;
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
      glm::vec3 boundingBox_min = mesh->AA;
      glm::vec3 boundingBox_max = mesh->BB;
      glm::vec3 boundingBox_size = boundingBox_max - boundingBox_min;
      glm::vec3 boundingBox_center = (boundingBox_max + boundingBox_min) * 0.5f;
      // !!! we use the max axis to scale the voxel data, because in the end we want
      // voxels be cubic instead of rectangular parallelepiped
      float maxAxis = std::max(std::max(boundingBox_size.x, boundingBox_size.y), boundingBox_size.z);
      glm::vec3 scale = glm::vec3(maxAxis);

      //std::cout << "bounding box min: " << boundingBox_min.x << " " << boundingBox_min.y << " " << boundingBox_min.z << std::endl;
      //std::cout << "bounding box max: " << boundingBox_max.x << " " << boundingBox_max.y << " " << boundingBox_max.z << std::endl;
      //std::cout << "bounding box size: " << boundingBox_size.x << " " << boundingBox_size.y << " " << boundingBox_size.z << std::endl;
      //std::cout << "bounding box center: " << boundingBox_center.x << " " << boundingBox_center.y << " " << boundingBox_center.z << std::endl;

      mVoxels.resize(grid_size);
      mVoxels.fill(0.0f);

      int num_triangles = mesh->indices.size() / 3;



      //for each triangle in mesh
      for (int t = 0; t < num_triangles; t++)
      {
         glm::uvec3 ix = glm::uvec3(mesh->indices[t * 3], mesh->indices[t * 3 + 1], mesh->indices[t * 3 + 2]);
         glm::vec4 po[3];
         glm::vec4 pw[3];
         glm::vec4 pn[3];
         glm::vec3 pvox[3];

         glm::vec3 bbmin = glm::vec3(1e10f);
         glm::vec3 bbmax = glm::vec3(-1e10f);

         //for each vertex in triangle
         for (int v = 0; v < 3; v++)
         {
               po[v] = glm::vec4(mesh->positions[ix[v]], 1.0f);
               // normalize vertex position into 0-1 range
               po[v] = po[v] - glm::vec4(boundingBox_min, 0.0f);
               po[v] = po[v] / glm::vec4(scale, 1.0f);

              
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
      //filled voxelization
      mVoxels = imfill(mVoxels);

   }
   // get a copy of the voxel data
   vector3d<float> getVoxels() {
      return mVoxels;
   }

private:
GModel * model;// the model to be voxelized, currently we only voxelized the model's first mesh
GMesh * mesh;// the mesh to be voxelized
vector3d<float> mVoxels;
glm::ivec3 grid_size;


};



















#endif // VOXELGENERATOR_H
