#ifndef DestructiveCSSceneObject_H
#define DestructiveCSSceneObject_H



#include "CommonSceneObject.h"
#include "VoxelGenerator.h"

#define particleCollisionMode_UPDATE_V_AND_X 0
#define particleCollisionMode_COLLISION 1

#define particleUniformGridMode_COMPUTE_COUNT 0
#define particleUniformGridMode_PREFIX_SUM 1 // to compute start buffer, not used for now since we use the prefix sum shader
#define particleUniformGridMode_INSERT_POINTS 2

#define particlePrefixSumMode_UPSWEEP 0
#define particlePrefixSumMode_DOWNSWEEP 1

#define local_workgroup_size_x_common 1024
#define local_workgroup_size_x_prefix_sum 256
#define local_workgroup_size_x_vgs_voxel 256
#define local_workgroup_size_x_vgs_face 64



// this is the particle struct for the destructive particle system
// it is not only used for computing, but also for rendering
struct alignas(16) Particle
{
    glm::vec4 curPos;   //pos.w = radius
    glm::vec4 prevPos;  //xprev.w = 1.0/mass
    glm::vec4 vel;    //xvel.w = mu (coefficient of friction)
    glm::uvec4 flags = glm::uvec4(0); //flags.y == lod, flags.z = type, flags.w = rgba8 color

   Particle(glm::vec3 pos, float radius, float w)
    {
        curPos = glm::vec4(pos, radius);
        prevPos = glm::vec4(pos, w);
        vel = glm::vec4(0.0f);
        flags = glm::uvec4(0, 0, 0, 0);
    }
};

// not valid for now
struct alignas(16) FaceConstraint
{
    glm::ivec2 indices; // the two voxel indices which the face constraint is applied to
    glm::ivec2 voxelFaceIndices; // each voxel has 6 faces, the face index of the voxel, x = voxel 1, y = voxel 2, for each voxel, 0 = -x, 1 = +x, 2 = -y, 3 = +y, 4 = -z, 5 = +z
    glm::vec2 strainLimit; // the strain limit of the face constraint: x = compression, y = stretch
    FaceConstraint(std::vector<int> indices, std::vector<int> voxelFaceIndices, std::vector<float> strainLimit) {
        if (indices.size() != 2 || voxelFaceIndices.size() != 2 || strainLimit.size() != 2) {
            std::cerr << "Error: [FaceConstraint] indices size should be 2, voxelFaceIndices size should be 2, strainLimit size should be 2" << std::endl;
            return;
        }
        for (int i = 0; i < 2; i++) {
            this->indices[i] = indices[i];
            this->voxelFaceIndices[i] = voxelFaceIndices[i];
            this->strainLimit[i] = strainLimit[i];
        }
    }
};

struct alignas(16) VoxelConstraint
{
    glm::ivec4 indices_0123; // the lower 4 particle indices of the voxel constraint
    glm::ivec4 indices_4567; // the upper 4 particle indices of the voxel constraint
    glm::ivec4 flags;        // flags.x = whether surface voxel (for rendering, 0->no, 1->yes), yzw not used
    VoxelConstraint(std::vector<int> indices, std::vector<int> _flags = {0, 0, 0, 0}) {
        if (indices.size() != 8) {
            std::cerr << "Error: [VoxelConstraint] indices size should be 8" << std::endl;
            return;
        }
        for (int i = 0; i < 8; i++) {
            if (i < 4) {
                indices_0123[i] = indices[i];
            } else {
                indices_4567[i-4] = indices[i];
            }
        }
        if (_flags.size() != 4) {
            std::cerr << "Error: [VoxelConstraint] flags size should be 4" << std::endl;
            return;
        }
        for (int i = 0; i < 4; i++) {
            this->flags[i] = _flags[i];
        }
        
    }
};


// for simplicity, we use a big UBO to store the uniform grid layout and 
// some other settings for the uniform grid and particle system
// all of the compute shaders will use this UBO to get the settings
// use alignas to ensure the UBO is aligned to 16 bytes
// so no worries about the alignment issue in the shader
struct alignas(16) DestructiveSystemUBO{
    alignas(4) int numParticles;
    alignas(4) int numCells;
    alignas(16) glm::vec3 gridMin;          // the min AABB corner of the uniform grid
    alignas(16) glm::vec3 gridMax;          // the max AABB corner of the uniform grid
    alignas(16) glm::ivec3 gridRes;         // the resolution of the uniform grid
    alignas(16) glm::vec3 cellSize;         // the size of each cell in the uniform grid, = (gridMax - gridMin) / gridRes
    alignas(4) int maxCellsPerElement;     // =1 for now, we might use it for the future extension
    // wall constraints, consists of 6 planes, each plane is a vec4
    alignas(16) glm::vec4 wallConstraints[6]; 
    alignas(16) glm::vec4 g;                // gravity vector
    alignas(4) float dt;                   // time step
    alignas(4) float c;                    // damping
    alignas(4) float omega_collision;      // collision omega
};


// a special structure representing embedded mesh data
// which not only contains the vertex data of the embedded mesh in voxel local space
// but also has the voxel index of the embedded mesh
struct EmbeddedVertex {
    glm::vec3 position;    // position of the vertex in voxel local space
    int voxelIndex;        // the voxel index of the embedded triangle
    glm::vec3 normal;      // normal of the vertex
    int padding;           // reserve for future use
};



class DestructiveCSComponent : public ComputeComponent {
public:
    DestructiveCSComponent() {
        //// ------------------------------------------------------------------------------------------------
        //// 1. initialize particle buffers (ping-pong buffer)
        //// simple case: 8 particles in a cube
        //particles.push_back(Particle(glm::vec3(0.0f, 0.0f, 0.0f), 0.01781f, 1.0f));
        //particles.push_back(Particle(glm::vec3(0.2f, 0.0f, 0.0f), 0.01781f, 1.0f));
        //particles.push_back(Particle(glm::vec3(0.0f, 0.2f, 0.0f), 0.01781f, 1.0f));
        //particles.push_back(Particle(glm::vec3(0.2f, 0.2f, 0.0f), 0.01781f, 1.0f));
        //particles.push_back(Particle(glm::vec3(0.0f, 0.0f, 0.2f), 0.01781f, 1.0f));
        //particles.push_back(Particle(glm::vec3(0.2f, 0.0f, 0.2f), 0.01781f, 1.0f));
        //particles.push_back(Particle(glm::vec3(0.0f, 0.2f, 0.2f), 0.01781f, 1.0f));
        //particles.push_back(Particle(glm::vec3(0.2f, 0.2f, 0.2f), 0.01781f, 1.0f));
//
//
        //// add 100 random particles
        //for (int i = 0; i < 35000; i++) {
        //    float x = (rand() % 100) / 100.0f;
        //    float y = (rand() % 100) / 100.0f;
        //    float z = (rand() % 100) / 100.0f;
        //    particles.push_back(Particle(glm::vec3(x, y, z), 0.01781f, 1.0f));
        //}


    }
    void initializeVoxels(const vector3d<float> &_voxels, const vector3d<int> &_voxelID, const vector3d<std::vector<pos_norm>> & _embeddedMesh, const vector3d<std::vector<pos_norm>> & _embeddedVoxelMesh, glm::mat4 _modelMatrix = glm::mat4(1.0f)) {
        particles.clear();
        particles.shrink_to_fit();
        // for face constraints, we need to know the index of each voxel in voxel buffer
        // so we define a voxelID buffer to temporarily store the index of each voxel
        vector3d<int> voxelID(_voxels.size(),-1); 
        // ------------------------------------------------------------------------------------------------
        // 1. initialize particle buffers (ping-pong buffer) and corresponding voxel constraint buffer
        // generate particles via mesh voxelization
        // we use the voxel data to generate particles
        // the voxel data is a 3D array, each element is a float value, 0.0f means no voxel, 1.0f means voxel
        // we will generate particles at eight corners of each voxel

        // now I don't want to handle the situation when the model get non-uniformly scaled
        // so let check if the model matrix is uniformly scaled
        if(_modelMatrix[0][0] != _modelMatrix[1][1] || _modelMatrix[0][0] != _modelMatrix[2][2]){
            std::cerr << "Error: [DestructiveCSComponent] the model matrix is not uniformly scaled" << std::endl;
            return;
        }

        // get scale ratio from model matrix
        float scaleRatio = _modelMatrix[0][0];
        float voxelSize = (1.0f / _voxels.size().x) * scaleRatio;
        float halfVoxelSize = voxelSize * 0.5f;
        float particleRadius = voxelSize * 0.25f;
        for (int i = 0; i < _voxels.size().x; i++) {
        for (int j = 0; j < _voxels.size().y; j++) {
        for (int k = 0; k < _voxels.size().z; k++) {
            if (_voxels.get(i, j, k) > 0.0f) {
                voxelID.set(i, j, k, particles.size() / 8);
                // pos normalize to vec3[0, 1]
                glm::vec3 pos_center = (glm::vec3(i, j, k) + glm::vec3(0.5,0.5,0.5));
                for(int v = 0; v < 8; v++){
                    glm::vec3 pos = pos_center;
                    static const glm::vec3 offsets[8] = {
                        glm::vec3(-0.25f, -0.25f, -0.25f),
                        glm::vec3(0.25f, -0.25f, -0.25f),
                        glm::vec3(-0.25f, 0.25f, -0.25f),
                        glm::vec3(0.25f, 0.25f, -0.25f),
                        glm::vec3(-0.25f, -0.25f, 0.25f),
                        glm::vec3(0.25f, -0.25f, 0.25f),
                        glm::vec3(-0.25f, 0.25f, 0.25f),
                        glm::vec3(0.25f, 0.25f, 0.25f)
                    };
                    pos += offsets[v];
                    pos /= glm::vec3(_voxels.size());
                    // then transform to match the bounding box of the mesh
                    // we let DestructiveCSSceneObject to handle the transformation
                    // and this component only handles the particle generation by the given voxel data and model matrix
                    pos = glm::vec3(_modelMatrix * glm::vec4(pos, 1.0f));
                    //std::cout << "pos: " << pos.x << " " << pos.y << " " << pos.z << std::endl;
                    particles.push_back(Particle(pos, particleRadius, 1.0f));
                }
                
                // voxel constraints: each voxel has 8 particles
                int baseIndex = particles.size() - 8;
                std::vector<int> indices;
                for(int j=0;j<8;j++){
                    indices.push_back(baseIndex+j);
                }
                std::vector<int> flags({0, 0, 0, 0});
                // check if the voxel is a surface voxel
                if(_voxels.valid_index(glm::vec3(i,j,k))){
                    // surface voxel refers to the voxel which has surface mesh triangles to render
                    // so in shader we don't want to render them unlike the inner voxels
                    // the ideal way to identify surface voxels is to check if the voxel has surface mesh triangles
                    // but sometimes some boundary voxels may not have surface mesh triangles after triangle culling
                    // so we also check whether the voxel is a boundary voxel
                    // surface voxels must be boundary voxels, but not all boundary voxels are surface voxels!
                    if(_voxels.get(i,j,k) == 2.0f||_embeddedMesh.get(i,j,k).size() > 0){
                        flags[0] = 1;// indicate that this voxel is a surface voxel
                    }
                }

                voxelConstraints.push_back(VoxelConstraint(indices, flags));

            }
        }
        }
        }
        // ------------------------------------------------------------------------------------------------
        // 2. initialize the face constraints for particles
        // face constraints: each voxel has 6 face constraints, and every two adjacent voxels share a face constraint
        // we need to generate face constraints for each pair of adjacent voxels
        const glm::ivec3 grid_size = _voxels.size();
        //std::vector<float> strainLimit({0.5f, -0.5f}); // strain limit for face constraints
        std::vector<float> strainLimit({0.0f, -0.0f}); // strain limit for face constraints
        //Create face constraints
        glm::ivec3 vox_ix;
        for (vox_ix.z = 0; vox_ix.z < grid_size.z; vox_ix.z++)
        for (vox_ix.y = 0; vox_ix.y < grid_size.y; vox_ix.y++)
        for (vox_ix.x = 0; vox_ix.x < grid_size.x; vox_ix.x++)
        {
                int id0 = voxelID.get(vox_ix);
                
                //x face con
                if (vox_ix.x < grid_size.x - 1)
                {
                    int id1 = voxelID.get(vox_ix + glm::ivec3(1, 0, 0));
                    if (id0 >= 0 && id1 >= 0)
                    {
                        FaceConstraint face_con({8 * id0, 8 * id1}, {1,-1}, strainLimit);
                        faceConstraints[0].push_back(face_con);
                    

                    }
                }
                //y face con
                if (vox_ix.y < grid_size.y - 1)
                {
                    int id1 = voxelID.get(vox_ix + glm::ivec3(0, 1, 0));
                    if (id0 >= 0 && id1 >= 0)
                    {
                        FaceConstraint face_con({8 * id0, 8 * id1}, {3,-1}, strainLimit);
                        faceConstraints[1].push_back(face_con);
                    

                    }
                }
                //z face con
                if (vox_ix.z < grid_size.z - 1)
                {
                    int id1 = voxelID.get(vox_ix + glm::ivec3(0, 0, 1));
                    if (id0 >= 0 && id1 >= 0)
                    {
                        FaceConstraint face_con({8 * id0, 8 * id1}, {5,-1}, strainLimit);
                        faceConstraints[2].push_back(face_con);
                    }
                }
        }
        // ------------------------------------------------------------------------------------------------
        // 3. initialize the embedded mesh data (for rendering only)
        embeddedMesh.clear();
        embeddedMesh.shrink_to_fit();
        //EmbeddedTriangle
        for (int i = 0; i < _voxels.size().x; i++) {
        for (int j = 0; j < _voxels.size().y; j++) {
        for (int k = 0; k < _voxels.size().z; k++) {
            //if (_voxels.get(i, j, k) > 0.0f) {
                for(const auto& vertex : _embeddedMesh.get(i,j,k)){
                    EmbeddedVertex et;
                    et.position = vertex.pos;
                    et.normal = vertex.norm;
                    et.voxelIndex = voxelID.get(i,j,k);
                    if(et.voxelIndex==-1){
                        std::cerr << "Error: [DestructiveCSComponent] voxel index is -1, I wrote bugs" << std::endl;
                    }
                    embeddedMesh.push_back(et);
                }
            //}
        }
        }
        }
        //std::cout<<"[DestructiveCSComponent] embedded mesh data size: "<<embeddedMesh.size()<<std::endl;
        // !!! temp code for debugging
        embeddedVoxelMesh.clear();
        embeddedVoxelMesh.shrink_to_fit();
        for (int i = 0; i < _voxels.size().x; i++) {
        for (int j = 0; j < _voxels.size().y; j++) {
        for (int k = 0; k < _voxels.size().z; k++) {
            //if (_voxels.get(i, j, k) > 0.0f) {
                for(const auto& vertex : _embeddedVoxelMesh.get(i,j,k)){
                    EmbeddedVertex et;
                    et.position = vertex.pos;
                    et.normal = vertex.norm;
                    et.voxelIndex = voxelID.get(i,j,k);
                    if(et.voxelIndex==-1){
                        std::cerr << "Error: [DestructiveCSComponent] voxel index is -1, I wrote bugs" << std::endl;
                    }
                    embeddedVoxelMesh.push_back(et);
                }
            //}
        }
        }
        }

    }

    void initializeBuffersAndShaders(){

        // initialize the ping-pong buffer for particles
        particleBuffers.resize(2);
        // create the buffer 1
        // use glNamedBufferStorage 
        glCreateBuffers(1, &particleBuffers[0]); 
        glNamedBufferStorage(particleBuffers[0], particles.size() * sizeof(Particle), particles.data(), GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT | GL_MAP_WRITE_BIT);
        // create the buffer 2
        glCreateBuffers(1, &particleBuffers[1]);
        glNamedBufferStorage(particleBuffers[1], particles.size() * sizeof(Particle), particles.data(), GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT | GL_MAP_WRITE_BIT);
        // ------------------------------------------------------------------------------------------------
        // 1. define the uniform grid layout and upload it to the GPU via a UBO
        SystemUBO.numParticles = particles.size();
        SystemUBO.numCells = 32 * 32 * 32;
        SystemUBO.gridMin = glm::vec3(-1.0f, -1.0f, -1.0f);
        SystemUBO.gridMax = glm::vec3(1.0f, 1.0f, 1.0f);
        SystemUBO.gridRes = glm::ivec3(32, 32, 32);
        SystemUBO.cellSize = (SystemUBO.gridMax - SystemUBO.gridMin) / glm::vec3(SystemUBO.gridRes);
        SystemUBO.maxCellsPerElement = 1;
        glm::vec4 walls[6] = {
            glm::vec4(1.0f, 0.0f, 0.0f, -1.0f),
            glm::vec4(-1.0f, 0.0f, 0.0f, -1.0f),
            glm::vec4(0.0f, 1.0f, 0.0f, -1.0f),
            glm::vec4(0.0f, -1.0f, 0.0f, -1.0f),
            glm::vec4(0.0f, 0.0f, 1.0f, -1.0f),
            glm::vec4(0.0f, 0.0f, -1.0f, -1.0f)
        };
        for (int i = 0; i < 6; i++) {
            SystemUBO.wallConstraints[i] = walls[i];
        }
        SystemUBO.g = glm::vec4(0.0f, -9.8f, 0.0f, 0.0f);
        SystemUBO.dt = 0.0015f;
        SystemUBO.c = 0.99f;
        SystemUBO.omega_collision = 0.90f;
        // create the UBO
        glGenBuffers(1, &systemUBOBuffer);
        glBindBuffer(GL_UNIFORM_BUFFER, systemUBOBuffer);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(DestructiveSystemUBO), &SystemUBO, GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_UNIFORM_BUFFER, systemUBOBindingPoint, systemUBOBuffer); // we should bind it just before the compute shader is executed, but do it here is also fine
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        

        // ------------------------------------------------------------------------------------------------
        // 2. initialize ugrid's count buffer and start buffer, and compact content buffer
        // we need to create 3 buffers for the uniform grid, they are all read-write buffers
        // 2.1. count buffer
        glGenBuffers(1, &ugridCountBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ugridCountBuffer);
        glBufferData(GL_SHADER_STORAGE_BUFFER, SystemUBO.numCells * sizeof(int), nullptr, GL_DYNAMIC_COPY);// dynamic copy means we will write to it in the compute shader, CPU will not read or write to it
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, countBindingPoint, ugridCountBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        // 2.2. start buffer
        glGenBuffers(1, &ugridStartBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ugridStartBuffer);
        glBufferData(GL_SHADER_STORAGE_BUFFER, SystemUBO.numCells * sizeof(int), nullptr, GL_DYNAMIC_COPY);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, startBindingPoint, ugridStartBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        // 2.3. compact content buffer
        glGenBuffers(1, &ugridCompactContentBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ugridCompactContentBuffer);
        glBufferData(GL_SHADER_STORAGE_BUFFER, SystemUBO.numParticles * sizeof(Particle), nullptr, GL_DYNAMIC_COPY);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, compactContentBindingPoint, ugridCompactContentBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        // ------------------------------------------------------------------------------------------------
        // 3. initialize the voxel constraints buffer and face constraints buffer
        // std::cout << "Size of FaceConstraint: " << sizeof(FaceConstraint) << std::endl;
        // std::cout << "Alignment of FaceConstraint: " << alignof(FaceConstraint) << std::endl;
        // voxel constraints buffer
        glGenBuffers(1, &voxelConstraintsBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, voxelConstraintsBuffer);
        glBufferData(GL_SHADER_STORAGE_BUFFER, voxelConstraints.size() * sizeof(VoxelConstraint), voxelConstraints.data(), GL_DYNAMIC_COPY);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, voxelConstraintsBindingPoint, voxelConstraintsBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        // face constraints buffer
        glGenBuffers(3, faceConstraintsBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, faceConstraintsBuffer[0]);
        glBufferData(GL_SHADER_STORAGE_BUFFER, faceConstraints[0].size() * sizeof(FaceConstraint), faceConstraints[0].data(), GL_DYNAMIC_COPY);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, faceConstraintsBuffer[1]);
        glBufferData(GL_SHADER_STORAGE_BUFFER, faceConstraints[1].size() * sizeof(FaceConstraint), faceConstraints[1].data(), GL_DYNAMIC_COPY);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, faceConstraintsBuffer[2]);
        glBufferData(GL_SHADER_STORAGE_BUFFER, faceConstraints[2].size() * sizeof(FaceConstraint), faceConstraints[2].data(), GL_DYNAMIC_COPY);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, faceConstraintsBindingPoint, faceConstraintsBuffer[0]);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);



        // ------------------------------------------------------------------------------------------------
        // 4. initialize the particle collision shader
        particleCollisionShader = new ComputeShader("shaders/particle_collision.comp");
        

        // ------------------------------------------------------------------------------------------------
        // 5. initialize the count/start(UniformGrid) shader and the prefix sum shader
        particleUniformGridShader = new ComputeShader("shaders/particle_uniform_grid.comp");
        particlePrefixSumShader = new ComputeShader("shaders/particle_prefix_sum.comp");


        // ------------------------------------------------------------------------------------------------
        // 6. initialize the VGS shader (both face and voxel)
        particleVGSFaceShader = new ComputeShader("shaders/particle_vgs_face.comp");
        particleVGSVoxelShader = new ComputeShader("shaders/particle_vgs_voxel.comp");


        // ------------------------------------------------------------------------------------------------
        // 7. initialize embedded mesh buffer
        glGenBuffers(1, &embeddedMeshBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, embeddedMeshBuffer);
        glBufferData(GL_SHADER_STORAGE_BUFFER, embeddedMesh.size() * sizeof(EmbeddedVertex), embeddedMesh.data(), GL_DYNAMIC_COPY);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, embeddedMeshBufferBindingPoint, embeddedMeshBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        // temporal voxel shit buffer
        glGenBuffers(1, &embeddedVoxelMeshBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, embeddedVoxelMeshBuffer);
        glBufferData(GL_SHADER_STORAGE_BUFFER, embeddedVoxelMesh.size() * sizeof(EmbeddedVertex), embeddedVoxelMesh.data(), GL_DYNAMIC_COPY);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, embeddedMeshBufferBindingPoint, embeddedVoxelMeshBuffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    }

    // ugrid related buffers do not need to be reuploaded for reinitialization
    // but particle buffers, voxel constraints buffer, and face constraints buffer need to be reuploaded
    void reuploadBuffer(){
        // reupload particle buffers
        glNamedBufferSubData(particleBuffers[0], 0, particles.size() * sizeof(Particle), particles.data());
        glNamedBufferSubData(particleBuffers[1], 0, particles.size() * sizeof(Particle), particles.data());
        // reupload voxel constraints buffer
        glNamedBufferSubData(voxelConstraintsBuffer, 0, voxelConstraints.size() * sizeof(VoxelConstraint), voxelConstraints.data());
        // reupload face constraints buffer
        glNamedBufferSubData(faceConstraintsBuffer[0], 0, faceConstraints[0].size() * sizeof(FaceConstraint), faceConstraints[0].data());
        glNamedBufferSubData(faceConstraintsBuffer[1], 0, faceConstraints[1].size() * sizeof(FaceConstraint), faceConstraints[1].data());
        glNamedBufferSubData(faceConstraintsBuffer[2], 0, faceConstraints[2].size() * sizeof(FaceConstraint), faceConstraints[2].data());
    }




    void Compute() override {
        //return; // don't do anything here

        for(int i = 0; i < 3; i++){
            // the system needs multiple iterations to generate feasible results
            // ------------------------------------------------------------------------------------------------
            // ------------------------------------------------------------------------------------------------
            // 1. update particle position and velocity
            particleCollisionShader->use();
            // bind the particle buffer as SSBO for instance's positions
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, inputParticleBindingPoint, particleBuffers[0]);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, outputParticleBindingPoint, particleBuffers[1]);
            // set the shader as MODE_UPDATE_V_AND_X mode
            particleCollisionShader->setInt("uMode", particleCollisionMode_UPDATE_V_AND_X);
            // dispatch compute shader
            glDispatchCompute((SystemUBO.numParticles / local_workgroup_size_x_common)+1, 1, 1);
            // wait for the compute shader to finish
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            // ping-pong the particle buffer
            pingPongParticleBuffer();


            // ------------------------------------------------------------------------------------------------
            // 2. uniform grid update
            particleUniformGridShader->use();
            // bind the particle buffer as SSBO for instance's positions
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, inputParticleBindingPoint, particleBuffers[0]);
            // bind the uniform grid related buffers
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, countBindingPoint, ugridCountBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, startBindingPoint, ugridStartBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, compactContentBindingPoint, ugridCompactContentBuffer);    

            // 2.1. compute count of particles in each cell
            // first, clear the count buffer to zero
            clearBufferToZero(ugridCountBuffer);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            // set the shader as MODE_COMPUTE_COUNT mode
            particleUniformGridShader->setInt("uMode", particleUniformGridMode_COMPUTE_COUNT);
            // dispatch compute shader
            glDispatchCompute((SystemUBO.numParticles / local_workgroup_size_x_common)+1, 1, 1);
            // wait for the compute shader to finish
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            // 2.2. compute prefix sum of the count (so that we reserve the space for each cell's particles)
            // we need to do prefix sum on the count buffer
            // but since we want to output the prefix sum to the start buffer
            // we copy the count buffer to the start buffer first
            glCopyNamedBufferSubData(ugridCountBuffer, ugridStartBuffer, 0, 0, SystemUBO.numCells * sizeof(int));
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            // bind the start buffer and let shader calculate the prefix sum
            particlePrefixSumShader->use();
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, inputParticleBindingPoint, ugridStartBuffer);
            particlePrefixSumShader->setInt("uNumElements", SystemUBO.numCells);
            //upsweep is a reduction
            particlePrefixSumShader->setInt("uMode", particlePrefixSumMode_UPSWEEP);
            int n = SystemUBO.numCells;
            assert(((n & (n - 1)) == 0)); //n must be a power of 2
            n = n / 2;
            int pass = 0;
            bool mBarrierEnabled = true;
            for (;;)
            {
                int stride = 2 << pass;
                //particlePrefixSumShader->SetUniform(mStrideLoc, stride); //STRIDE = 2, 4, 8, ...
                //particlePrefixSumShader->SetGridSize(glm::ivec3(n, 1, 1));
                //particlePrefixSumShader->Dispatch();
                particlePrefixSumShader->setInt("stride", stride); //STRIDE = 2, 4, 8, ...
                glDispatchCompute((n / local_workgroup_size_x_prefix_sum)+1, 1, 1);
                if (mBarrierEnabled) glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

                if (n == 1) break;
                n = n / 2;
                pass = pass + 1;
            }

            //downsweep
            particlePrefixSumShader->setInt("uMode", particlePrefixSumMode_DOWNSWEEP);
            for (;;)
            {
                int stride = 2 << pass;
                particlePrefixSumShader->setInt("stride", stride); //STRIDE = 2, 4, 8, ...
                glDispatchCompute((n / local_workgroup_size_x_prefix_sum)+1, 1, 1);
               
                if (mBarrierEnabled) glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

                if (n == SystemUBO.numCells / 2) break;
                n = n * 2;
                pass = pass - 1;
            }

            // 2.3. reorder the particles based on the prefix sum (inside particle IDs into each cell)
            // again, clear the count buffer to zero
            clearBufferToZero(ugridCountBuffer);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            particleUniformGridShader->use();
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, inputParticleBindingPoint, particleBuffers[0]);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, countBindingPoint, ugridCountBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, startBindingPoint, ugridStartBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, compactContentBindingPoint, ugridCompactContentBuffer);    
            // set the shader as MODE_INSERT_POINTS mode
            particleUniformGridShader->setInt("uMode", particleUniformGridMode_INSERT_POINTS);
            // dispatch compute shader
            glDispatchCompute((SystemUBO.numParticles / local_workgroup_size_x_common)+1, 1, 1);
            // wait for the compute shader to finish
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);


            // ------------------------------------------------------------------------------------------------
            // 3. particle collision
            particleCollisionShader->use();
            // bind the particle buffer as SSBO for instance's positions
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, inputParticleBindingPoint, particleBuffers[0]);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, outputParticleBindingPoint, particleBuffers[1]);
            // // bind the uniform grid related buffers
            // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, countBindingPoint, ugridCountBuffer);
            // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, startBindingPoint, ugridStartBuffer);
            // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, compactContentBindingPoint, ugridCompactContentBuffer);
            // set the shader as MODE_COLLISION mode
            particleCollisionShader->setInt("uMode", particleCollisionMode_COLLISION);
            // dispatch compute shader
            glDispatchCompute((SystemUBO.numParticles / local_workgroup_size_x_common)+1, 1, 1);
            // wait for the compute shader to finish
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            // ping-pong the particle buffer
            pingPongParticleBuffer();
            
            


            // ------------------------------------------------------------------------------------------------
            // 4. Voxel-Gram-Schmidt orthogonalization for shape correction (core part of the algorithm)
            

            // 4.1 VGS on voxels
            particleVGSVoxelShader->use();
            // bind the particle buffer as SSBO for instance's positions
            // no ping-pong buffer here since no read-write conflict, only mapping computation
            // directly modify the particle buffer 0 (input buffer)
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, inputParticleBindingPoint, particleBuffers[0]);
            // bind the uniform grid related buffers, not necessary probably
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, countBindingPoint, ugridCountBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, startBindingPoint, ugridStartBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, compactContentBindingPoint, ugridCompactContentBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, voxelConstraintsBindingPoint, voxelConstraintsBuffer);
            // set the shader as MODE_COLLISION mode
            particleVGSVoxelShader->setInt("uMode", 5);
            // dispatch compute shader, for each voxel, we need 8 particles, but in shader, each thread handles 8 particles, so we need to divide the size by 8 and multiply by 8
            glDispatchCompute(((voxelConstraints.size()/8*8) / local_workgroup_size_x_vgs_voxel)+1, 1, 1);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            
            //std::cout<<"voxelConstraints.size(): "<<voxelConstraints.size()<<std::endl;
            //std::cout<<"faceConstraints[0].size(): "<<faceConstraints[0].size()<<std::endl;
            //std::cout<<"faceConstraints[1].size(): "<<faceConstraints[1].size()<<std::endl;
            //std::cout<<"faceConstraints[2].size(): "<<faceConstraints[2].size()<<std::endl;
            //std::cout<<"particles.size(): "<<particles.size()<<std::endl;



            // 4.2 VGS on faces
            particleVGSFaceShader->use();
            // bind the particle buffer 0 as inout buffer
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, inputParticleBindingPoint, particleBuffers[0]);
            // bind the uniform grid related buffers, not necessary probably
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, countBindingPoint, ugridCountBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, startBindingPoint, ugridStartBuffer);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, compactContentBindingPoint, ugridCompactContentBuffer);
            
            // set the shader as MODE_COLLISION mode
            particleVGSFaceShader->setInt("uMode", 5);
            for(int i = 0; i < 3; i++){
                glBindBufferBase(GL_SHADER_STORAGE_BUFFER, faceConstraintsBindingPoint, faceConstraintsBuffer[i]);
                // dispatch compute shader
                glDispatchCompute((faceConstraints[i].size() / local_workgroup_size_x_vgs_face)+1, 1, 1);
                glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            }


            // ------------------------------------------------------------------------------------------------
        }
        

        return;
    }
    const std::vector<GLuint> & getParticleBuffer() const{
        return particleBuffers;
    }
    const GLuint & getVoxelConstraintsBuffer() const{
        return voxelConstraintsBuffer;
    }
    int getNumParticles() const {
        return particles.size();
    }
    int getNumVoxels() const {
        return voxelConstraints.size();
    }
    int getNumEmbeddedVertices() const {
        return embeddedMesh.size();
    }
    int getNumEmbeddedVoxelVertices() const {
        return embeddedVoxelMesh.size();
    }
    GLuint getSystemUBOBuffer() const {
        return systemUBOBuffer;
    }
    GLuint getEmbeddedMeshBuffer() const {
        return embeddedMeshBuffer;
    }
    GLuint getEmbeddedMeshVoxelBuffer() const {
        return embeddedVoxelMeshBuffer;
    }

private:
    // some data
    // CPU side initial particle data
    std::vector<Particle> particles;
    // GPU side particle data, should be binded as SSBO
    // need ping-pong buffer for updating
    std::vector<GLuint> particleBuffers;// size should be 2
    // CPU side voxel constraints and face constraints data
    std::vector<VoxelConstraint> voxelConstraints;
    std::vector<FaceConstraint> faceConstraints[3]; // faces have X, Y, Z three directions, so we need 3 buffers for partition
    // CPU side embedded mesh data
    std::vector<EmbeddedVertex> embeddedMesh;// x1,y1,z1,x2,y2,z2,... for now
    std::vector<EmbeddedVertex> embeddedVoxelMesh;// a temorary vertex buffer, in the future, we may want to merge it with embeddedMesh
    // GPU side embedded mesh data, should be binded as SSBO
    GLuint embeddedMeshBuffer;
    GLuint embeddedVoxelMeshBuffer;
    // GPU side voxel constraints and face constraints data, should be binded as SSBO
    GLuint voxelConstraintsBuffer;
    GLuint faceConstraintsBuffer[3]; // faces have X, Y, Z three directions, so we need 3 buffers for partition
    // CPU side uniform grid layout and settings
    DestructiveSystemUBO SystemUBO;
    // GPU side uniform grid layout and settings, should be binded as UBO
    GLuint systemUBOBuffer;
    // Uniform Grid related buffers
    GLuint ugridCountBuffer;
    GLuint ugridStartBuffer;
    GLuint ugridCompactContentBuffer;
    // we need several compute shaders to update the particle data
    ComputeShader * particleCollisionShader;
    ComputeShader * particleVGSShader;
    ComputeShader * particlePrefixSumShader;
    ComputeShader * particleUniformGridShader;
    ComputeShader * particleVGSFaceShader;
    ComputeShader * particleVGSVoxelShader;

    void pingPongParticleBuffer() {
        std::swap(particleBuffers[0], particleBuffers[1]);
    }
        // just a helper function to clear the COUNT buffers to zero
    void clearBufferToZero(GLuint buffer) {
        static int zero = 0;
        glClearNamedBufferData(buffer, GL_R32I, GL_RED_INTEGER, GL_INT, &zero);
    }

public:
    // binding points for the compute shaders
    const GLuint systemUBOBindingPoint = 1; // 0 has been used by the camera UBO, so we use 1 here
    const GLuint inputParticleBindingPoint = 2; // 0 has been used by the camera UBO, 1 has been used by the DestructiveSystemUBO, so we use 2 here
    const GLuint outputParticleBindingPoint = 3;
    const GLuint countBindingPoint = 4;
    const GLuint startBindingPoint = 5;
    const GLuint compactContentBindingPoint = 6;
    const GLuint voxelConstraintsBindingPoint = 7;
    const GLuint faceConstraintsBindingPoint = 8;
    const GLuint embeddedMeshBufferBindingPoint = 9; // I hope I can get rid of this buffer in the future

};


class DestructiveRenderComponent : public RenderComponent {
public:
    DestructiveRenderComponent(RenderContext * _context = nullptr, enum renderQueue _renderPriority = OPAQUE) {
        this->renderPriority = _renderPriority;
        // initialize the particle shader
        particleShader = new Shader("shaders/particle.vert", "shaders/particle.frag");
        // initialize the VAO, VBO, EBO for particle rendering--a sphere
        generateBufferResource();
        // initialize the voxel shader
        // common cubic voxel shader
        voxelShader = new Shader("shaders/voxel.vert", "shaders/voxel.frag");
        // the voxel plus embedded mesh shader, used for rendering surface voxels of the mesh
        skinShader = new Shader("shaders/voxel_skin.vert", "shaders/voxel_skin.frag");


        debugShader = new Shader("shaders/debug_AABB.vert", "shaders/debug_AABB.frag");
        glGenVertexArrays(1, &debug_VAO); // openGL sometimes requires a VAO to draw something, even if we don't use it
        targetMesh = nullptr;
    }
    void Render() override {
        // if the compute component is not set, we cannot render the particles
        if (computeComponent.expired()) {
            std::cout<<"ERROR: [DestructiveRenderComponent] compute component is expired"<<std::endl;
            return;
        }
        // --------------------particle rendering--------------------------------
        if(bRenderParticles){
            particleShader->use();
		    glBindVertexArray(particle_VAO);
            // bind the particle buffer as SSBO for instance's positions
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, computeComponent.lock()->inputParticleBindingPoint, computeComponent.lock()->getParticleBuffer()[0]);
            // instance rendering, use glDrawArraysInstanced
            glDrawElementsInstanced(GL_TRIANGLES, 1200, GL_UNSIGNED_INT, 0, computeComponent.lock()->getNumParticles());
		    glBindVertexArray(0);
        }
        // use debug shader to render the boundary box of the system
        debugShader->use();
        glBindVertexArray(debug_VAO);
        glm::vec3 AABB_min = glm::vec3(-1.0f, -1.0f, -1.0f);
        glm::vec3 AABB_max = glm::vec3(1.0f, 1.0f, 1.0f);
        glm::mat4 model = glm::mat4(1.0f);
        debugShader->setMat4("model", model);
        debugShader->setVec3("AABB_min", AABB_min);
        debugShader->setVec3("AABB_max", AABB_max);
        debugShader->setVec4("lineColor", glm::vec4(0.3, 1, 1, 1));
        glDrawArrays(GL_LINES, 0, 24);
        glBindVertexArray(0);
        // -------------------------------------------------------------------
        // --------------------voxel rendering--------------------------------
        if(bRenderVoxels){
            voxelShader->use();
            glBindVertexArray(debug_VAO);// nothing to do with the VAO, just to avoid the OpenGL error
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, computeComponent.lock()->inputParticleBindingPoint, computeComponent.lock()->getParticleBuffer()[0]);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, computeComponent.lock()->voxelConstraintsBindingPoint, computeComponent.lock()->getVoxelConstraintsBuffer());
            // no instance rendering, just draw the voxel mesh by drawing a lot of vertices which are generated by reading the voxel data and particle data
            glDrawArrays(GL_TRIANGLES, 0, computeComponent.lock()->getNumVoxels() * 36);
            glBindVertexArray(0);
        }
        


        // -------------------------------------------------------------------
        // --------------------skin mesh rendering----------------------------
        if(bRenderSkinMesh){
            skinShader->use();
            
            glBindVertexArray(debug_VAO);// nothing to do with the VAO, just to avoid the OpenGL error
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, computeComponent.lock()->inputParticleBindingPoint, computeComponent.lock()->getParticleBuffer()[0]);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, computeComponent.lock()->voxelConstraintsBindingPoint, computeComponent.lock()->getVoxelConstraintsBuffer());
            //glDisable(GL_DEPTH_TEST);
            //glDisable(GL_CULL_FACE); // disable culling for now, have some bug with the ccw and cw embedding voxels
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, computeComponent.lock()->embeddedMeshBufferBindingPoint, computeComponent.lock()->getEmbeddedMeshVoxelBuffer());
            skinShader->setVec4("color", glm::vec4(1.0f, 0.5f, 0.8f, 1.0f));
            glDrawArrays(GL_TRIANGLES, 0, computeComponent.lock()->getNumEmbeddedVoxelVertices());
            //glEnable(GL_DEPTH_TEST);
            //glEnable(GL_CULL_FACE); // enable culling
            // no instance rendering, just draw the voxel mesh by drawing a lot of vertices which are generated by reading the voxel data and particle data
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, computeComponent.lock()->embeddedMeshBufferBindingPoint, computeComponent.lock()->getEmbeddedMeshBuffer());
            skinShader->setVec4("color", glm::vec4(0.6f, 0.5f, 0.7f, 1.0f));
            glDrawArrays(GL_TRIANGLES, 0, computeComponent.lock()->getNumEmbeddedVertices());
            
            glBindVertexArray(0);
        }




        // --------------------original mesh rendering------------------------
        if(targetMesh != nullptr && brenderOriginalMesh){
            // render the target mesh
            // TODO: get rid of uniform context but use ubo context instead
            // if (context != nullptr) {
            //     targetMesh->prepareDraw(context);
            // }
            glDisable(GL_DEPTH_TEST);
            //targetMesh->draw();
            glEnable(GL_DEPTH_TEST);
        }
        // -------------------------------------------------------------------
    }
    void bindComputeComponent(std::shared_ptr<DestructiveCSComponent> _computeComponent) {
        this->computeComponent = _computeComponent;
    }
    void setTargetMesh(GMVPObject * _targetMesh) {
        this->targetMesh = _targetMesh;
    }

    void setRenderParticles(bool b) {
        bRenderParticles = b;
    }
    void setRenderVoxels(bool b) {
        bRenderVoxels = b;
    }
    void setRenderSkinMesh(bool b) {
        bRenderSkinMesh = b;
    }
    
    
private:
    // it needs some data which belongs to the compute component for rendering
    // such as the particle buffer and the number of particles
    std::weak_ptr<DestructiveCSComponent> computeComponent;
    // we use instancing to render the particles
    // but we still need to use a VAO, VBO, EBO to store the particle mesh data
    GLuint particle_VAO, particle_VBO, particle_EBO;
    // we also need a shader to render the particles
    Shader * particleShader;
    // and a shader to optionally render the voxels via particle data
    Shader * voxelShader;
    // and a shader to render the skin mesh triangles
    Shader * skinShader;
    // also a debug shader to render the AABB stuff ... maybe not necessary
    Shader * debugShader;
    GMVPObject * targetMesh;
    bool brenderOriginalMesh = true;
    // though we don't need any VBO or EBO for debug rendering, I found it is necessary to create a VAO
    // otherwise, the draw call will be automatically ignored by the OpenGL...
    GLuint debug_VAO; 

    bool bRenderVoxels = false;
    bool bRenderParticles = true;
    bool bRenderSkinMesh = false;

    void generateBufferResource(){
		std::vector<GLfloat> sphereVertices;
		int sectors = 20;
		int stacks = 10;
		int num = (stacks * 2) * sectors * 3;
        float radius = 1.0f;
        glm::vec3 center = glm::vec3(0.0f, 0.0f, 0.0f);
		for (int i = 0; i <= stacks; ++i) {
			float stackAngle = glm::pi<float>() / 2.0f - i * glm::pi<float>() / stacks;
			float y = radius * sin(stackAngle);
			for (int j = 0; j <= sectors; ++j) {
				float sectorAngle = 2.0f * glm::pi<float>() * j / sectors;
				float x = radius * cos(stackAngle) * cos(sectorAngle);
				float z = radius * cos(stackAngle) * sin(sectorAngle);

				// position
				sphereVertices.push_back(x+center.x);
				sphereVertices.push_back(y+center.y);
				sphereVertices.push_back(z+center.z);

				// normal
				glm::vec3 normal = glm::normalize(glm::vec3(x, y, z));
				sphereVertices.push_back(normal.x);
				sphereVertices.push_back(normal.y);
				sphereVertices.push_back(normal.z);
				

				// uv
				sphereVertices.push_back((float)j / sectors);
				sphereVertices.push_back((float)i / stacks);
			}
		}
		std::vector<GLuint> sphereIndices;
		for (int i = 0; i < stacks; ++i) {
			for (int j = 0; j < sectors; ++j) {
				int top = i * (sectors + 1) + j;
				int bottom = top + sectors + 1;

				sphereIndices.push_back(top);
				sphereIndices.push_back(top + 1);
				sphereIndices.push_back(bottom);

				sphereIndices.push_back(bottom);
				sphereIndices.push_back(top + 1);
				sphereIndices.push_back(bottom + 1);
			}
		}
		glGenBuffers(1, &particle_EBO);
		glGenVertexArrays(1, &particle_VAO);
		glGenBuffers(1, &particle_VBO);
		glBindVertexArray(particle_VAO);
		glBindBuffer(GL_ARRAY_BUFFER, particle_VBO); 
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * sphereVertices.size(), sphereVertices.data(), GL_STATIC_DRAW); 
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0); // set vertex attribute pointers: position
		glEnableVertexAttribArray(0); // activate vertex attribute
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float))); // set vertex attribute pointers: normal
		glEnableVertexAttribArray(1); // activate vertex attribute
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float))); // set vertex attribute pointers: uv
		glEnableVertexAttribArray(2); // activate vertex attribute
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, particle_EBO); // bind EBO
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * sphereIndices.size(), sphereIndices.data(), GL_STATIC_DRAW); // copy the index data to EBO
		glBindVertexArray(0);
	}
};


// DestructiveCSSceneObject is a Compute Shader-based Scene Object
// it has a pipeline that uses compute shader to update the object's data
// and also in the end it uses the data to render the object via vs and fs
class DestructiveCSSceneObject : public RenderableSceneObject, public ComputeSceneObject {
public:
    DestructiveCSSceneObject(std::shared_ptr<GModel> _targetMesh, enum renderQueue _renderPriority = OPAQUE, RenderContext * _context = nullptr) : RenderableSceneObject() {
        this->targetMesh = _targetMesh;
        if(targetMesh){
            this->voxelGeneratorComponent = std::make_shared<VoxelGeneratorComponent>(targetMesh.get(),glm::ivec3(32,32,32));
            voxelGeneratorComponent->GenerateVoxelData();
            voxelGeneratorComponent->GenerateEmbeddedMeshSurfaceData();
            // get necessary data from the voxel generator component
            mVoxels = voxelGeneratorComponent->getVoxels();
            voxelID = voxelGeneratorComponent->getVoxelID();
            embeddedMeshData = voxelGeneratorComponent->getEmbeddedMeshData();
            embeddedVoxelMeshData = voxelGeneratorComponent->getEmbeddedVoxelMeshData();
            
        }
        else{
            std::cerr << "Error: [DestructiveCSSceneObject] failed to cast targetMesh" << std::endl;
        }
        
        // first we initialize the compute component
        this->computeComponent = std::make_shared<DestructiveCSComponent>();
        // then we initialize the render component
        this->renderComponent = std::make_shared<DestructiveRenderComponent>(_context, _renderPriority);
        // computeComponent initializes the particle data
        // and the renderComponent uses the particle data to render the object
        auto renderComponentDerived = std::dynamic_pointer_cast<DestructiveRenderComponent>(this->renderComponent);
        auto computeComponentDerived = std::dynamic_pointer_cast<DestructiveCSComponent>(this->computeComponent);
        if(computeComponentDerived){
            // calculate the transformation matrix of the voxel data(normalized voxel data is 0-1)
            glm::mat4 modelMatrix = glm::mat4(1.0f);
            if(targetMesh){
                // 1. scale the voxel data to the target mesh's bounding box, translate the voxel data to the target mesh's center
                glm::vec3 min = targetMesh->meshes[0]->AA;
                glm::vec3 max = targetMesh->meshes[0]->BB;
                glm::vec3 size = max - min;
                // we use the largest axis to scale the voxel data, must be uniform scale because I don't want to handle non-uniform scale...
                float maxAxis = std::max(std::max(size.x, size.y), size.z);
                glm::vec3 scale = glm::vec3(maxAxis);
                glm::mat4 scaleMatrix = glm::scale(glm::mat4(1.0f), scale);
                glm::mat4 translateMatrix = glm::translate(glm::mat4(1.0f), min);
                modelMatrix = translateMatrix * scaleMatrix;
                std::cout<<"modelMatrix: "<<modelMatrix[0][0]<<" "<<modelMatrix[0][1]<<" "<<modelMatrix[0][2]<<" "<<modelMatrix[0][3]<<std::endl;
                std::cout<<"modelMatrix: "<<modelMatrix[1][0]<<" "<<modelMatrix[1][1]<<" "<<modelMatrix[1][2]<<" "<<modelMatrix[1][3]<<std::endl;
                std::cout<<"modelMatrix: "<<modelMatrix[2][0]<<" "<<modelMatrix[2][1]<<" "<<modelMatrix[2][2]<<" "<<modelMatrix[2][3]<<std::endl;
                std::cout<<"modelMatrix: "<<modelMatrix[3][0]<<" "<<modelMatrix[3][1]<<" "<<modelMatrix[3][2]<<" "<<modelMatrix[3][3]<<std::endl;
                // 2. apply mesh's model matrix to the voxel data
                modelMatrix = targetMesh->meshes[0]->model * modelMatrix;
            }
           
            

            
            
            computeComponentDerived->initializeVoxels(mVoxels,voxelID, embeddedMeshData, embeddedVoxelMeshData, modelMatrix);
            computeComponentDerived->initializeBuffersAndShaders();
        }
        else{
            std::cerr << "Error: [DestructiveCSSceneObject] failed to cast computeComponent" << std::endl;
        }
        if (renderComponentDerived && computeComponentDerived) {
            renderComponentDerived->bindComputeComponent(computeComponentDerived);
        } 
        else {
            std::cerr << "Error: [DestructiveCSSceneObject] failed to cast renderComponent or computeComponent" << std::endl;
        }
        
        if(renderComponentDerived && targetMesh){
            renderComponentDerived->setTargetMesh(targetMesh.get());
        }
        else{
            std::cerr << "Error: [DestructiveCSSceneObject] failed to set target mesh" << std::endl;
        }
    
    }
    ~DestructiveCSSceneObject() {
    }

    void reinitializeSystem() {
        auto computeComponentDerived = std::dynamic_pointer_cast<DestructiveCSComponent>(this->computeComponent);
        if(computeComponentDerived){
            computeComponentDerived->reuploadBuffer();
        }
    }


    void setRenderParticles(bool b) {
        auto renderComponentDerived = std::dynamic_pointer_cast<DestructiveRenderComponent>(this->renderComponent);
        if(renderComponentDerived){
            renderComponentDerived->setRenderParticles(b);
        }
    }
    void setRenderVoxels(bool b) {
        auto renderComponentDerived = std::dynamic_pointer_cast<DestructiveRenderComponent>(this->renderComponent);
        if(renderComponentDerived){
            renderComponentDerived->setRenderVoxels(b);
        }
    }
    void setRenderSkinMesh(bool b) {
        auto renderComponentDerived = std::dynamic_pointer_cast<DestructiveRenderComponent>(this->renderComponent);
        if(renderComponentDerived){
            renderComponentDerived->setRenderSkinMesh(b);
        }
    }
    void setComputeComponentActive(bool b) {
        this->computeComponent->active = b;
        
    }

protected:
    // it also has a render component, which is inherited from CommonSceneObject
    // and also a voxel render component, to initialize the voxel data and particle data
    std::shared_ptr<VoxelGeneratorComponent> voxelGeneratorComponent;
    std::shared_ptr<GModel> targetMesh;// we assume the target mesh is a GModel with only one mesh currently
    // get them from VoxelGeneratorComponent
    vector3d<float> mVoxels;
    vector3d<std::vector<pos_norm>> embeddedMeshData;
    vector3d<std::vector<pos_norm>> embeddedVoxelMeshData;
    vector3d<int> voxelID;
};



#include "UIManager.h"
// the UI for such system, I put it here
class DestructiveCSUI_left : public UIWindow {
public:
    DestructiveCSUI_left(DestructiveCSSceneObject * _systemRef) {
        systemRef = _systemRef;
        systemRef->setComputeComponentActive(bSystemComputeGPU);
        systemRef->setRenderSkinMesh(bDrawMesh);
        systemRef->setRenderParticles(bDrawParticle);
        systemRef->setRenderVoxels(bDrawVoxel);
    }
    void drawWindow() override {
        if (!display) return;
        const char* sceneNames[] = { "SCENE1", "SCENE2", "SCENE3", "SCENE4" };
        ImGui::SetNextWindowPos(ImVec2(5, 10), ImGuiCond_Always);
		ImGui::SetNextWindowSize(ImVec2(340, ImGui::GetIO().DisplaySize.y-20), ImGuiCond_Always);
        if (ImGui::Begin("DestructiveSystem", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize)) {
			ImGui::Text("XXX");
            if (ImGui::Button("Reinitialize Scene")) {
                // reinitialize the system
                if(systemRef){
                    systemRef->reinitializeSystem();
                }
                else{
                    std::cerr << "Error: [DestructiveCSUI_left] systemRef is nullptr" << std::endl;
                }
            }
            if (ImGui::Checkbox("GPU simulation", &bSystemComputeGPU)) {
                systemRef->setComputeComponentActive(bSystemComputeGPU);
            }
            if (ImGui::Combo("Scene Selection", &currentScene, sceneNames, IM_ARRAYSIZE(sceneNames))) {

            }
            if (ImGui::Checkbox("Draw Mesh", &bDrawMesh)) {
                systemRef->setRenderSkinMesh(bDrawMesh);
            }
            if (ImGui::Checkbox("Draw Particle", &bDrawParticle)) {
                systemRef->setRenderParticles(bDrawParticle);
            }
            if (ImGui::Checkbox("Draw Voxel", &bDrawVoxel)) {
                systemRef->setRenderVoxels(bDrawVoxel);
            }
            ImGui::End();
		}
    }
private:
    bool bSystemComputeGPU = false;
    bool bDrawMesh = true;
    bool bDrawParticle = false;
    bool bDrawVoxel = true;
    int currentScene = 0;
    DestructiveCSSceneObject * systemRef;
};





#endif // DestructiveCSSceneObject_H
