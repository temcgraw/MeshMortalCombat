#version 460 core

const int kPointsInBinding = 2;
const int kVoxelConstraintsBinding = 7;
const int kEmbeddedVerticesBinding = 9;

layout(std140, binding = 0) uniform RenderInfo{
    mat4 view;
    mat4 projection;
    vec3 cameraPos;
};

layout(std140, binding = 10) uniform RenderLightInfo{
    bool bUsePointLight;
    int PointLightCount;
    vec3 PointLightPos1;
    vec3 PointLightColor1;
    vec3 PointLightPos2;
    vec3 PointLightColor2;
    vec3 PointLightPos3;
    vec3 PointLightColor3;
    bool bUseDirectionalLight;
    vec3 DirectionalLightDir;
    vec3 DirectionalLightColor;
    vec3 AmbientLightColor;
};

// instance properties from the particle system SSBO
struct Particle {
    vec4 curPos;   //we need curPos.xyz to get the position of the particle, and w to get the size of the particle (radius)
    vec4 prevPos;  
    vec4 vel;     
    uvec4 flags;
};

struct VoxelConstraint{
	int indices[8]; // 0-7
	int isBoundary;
	int padding[3];
};

struct EmbeddedVertex{
    vec3 pos;
    int voxelIndex;
    vec3 normal;
    int padding;
};

uniform int meshIndexStart = 0;


layout(std430, binding = kPointsInBinding) restrict readonly buffer ParticleBuffer {
    Particle particles[];
};
layout (std430, binding = kVoxelConstraintsBinding) restrict readonly buffer VOXELS
{
	VoxelConstraint voxelConstraints[];
};

layout(std430, binding = kEmbeddedVerticesBinding) restrict readonly buffer EmbeddedMeshBuffer {
    EmbeddedVertex vertices[];
};

out vec3 Normal;

out vec3 FragPos;

out flat int vertexIndex;


mat4 translate(mat4 m, vec3 v) {
    mat4 translationMatrix = mat4(1.0);
    translationMatrix[3] = vec4(v, 1.0);
    return m * translationMatrix;
}

mat4 scale(mat4 m, vec3 v) {
    mat4 scaleMatrix = mat4(1.0);
    scaleMatrix[0][0] = v.x;
    scaleMatrix[1][1] = v.y;
    scaleMatrix[2][2] = v.z;
    return m * scaleMatrix;
}

const vec3 vertex_cube_coord[8] = vec3[8]( ivec3(0.0, 0.0, 0.0), ivec3(1.0, 0.0, 0.0),
                                          ivec3(0.0, 1.0, 0.0), ivec3(1.0, 1.0, 0.0),
                                          ivec3(0.0, 0.0, 1.0), ivec3(1.0, 0.0, 1.0),
                                          ivec3(0.0, 1.0, 1.0), ivec3(1.0, 1.0, 1.0));

//local indices of cube faces
const ivec3 face_indices[12] = ivec3[12]( ivec3(0,4,2), ivec3(4,6,2),   //-x face
                                          ivec3(1,3,5), ivec3(3,7,5),   //+x face
                                          ivec3(4,0,1), ivec3(4,1,5),   //-y face
                                          ivec3(3,2,6), ivec3(7,3,6),   //+y face
                                          ivec3(0,2,1), ivec3(2,3,1),   //-z face
                                          ivec3(4,5,6), ivec3(7,6,5));  //+z face

//cube corner offsets
const vec3 voffset[8] = vec3[8]( vec3(-1.0, -1.0, -1.0), vec3(+1.0, -1.0, -1.0),
                                 vec3(-1.0, +1.0, -1.0), vec3(+1.0, +1.0, -1.0), 
                                 vec3(-1.0, -1.0, +1.0), vec3(+1.0, -1.0, +1.0), 
                                 vec3(-1.0, +1.0, +1.0), vec3(+1.0, +1.0, +1.0));


void main()
{
    int VertexID = gl_VertexID + meshIndexStart;
    vec3 pos = vertices[VertexID].pos; // the position of the vertex in the current voxel's local space
    // now we need to find the voxel that contains this vertex
    // and convert the vertex position to world space via the voxel's transformation matrix
    int vox_id = vertices[VertexID].voxelIndex;
    //if(vox_id==-1){
    //    vox_id = 0;
    //}
    int voxel_vertices[8];
    for(int i=0; i<8; i++)
    {
        voxel_vertices[i] = voxelConstraints[vox_id].indices[i];
    }
    vec3 voxel_basis[3], voxel_basis_scaled[3];
    voxel_basis[0] = normalize(particles[voxel_vertices[1]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz);
    voxel_basis[1] = normalize(particles[voxel_vertices[2]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz);
    voxel_basis[2] = normalize(particles[voxel_vertices[4]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz);
    // maybe we should not normalize the basis vectors 
    // to keep the scale of the voxel
    voxel_basis_scaled[0] = (particles[voxel_vertices[1]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz)/(particles[voxel_vertices[0]].curPos.w*2);
    voxel_basis_scaled[1] = (particles[voxel_vertices[2]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz)/(particles[voxel_vertices[0]].curPos.w*2);
    voxel_basis_scaled[2] = (particles[voxel_vertices[4]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz)/(particles[voxel_vertices[0]].curPos.w*2);

    vec3 particle_base_pos = particles[voxel_vertices[0]].curPos.xyz; // the position of the first particle in the voxel, but still has some distance to the corner of the voxel
    vec3 voxel_base_pos = particle_base_pos - (voxel_basis[0] + voxel_basis[1] + voxel_basis[2]) * particles[voxel_vertices[0]].curPos.w;// notice this part we don't use voxel_basis_scaled
    vec3 vertex_offset = (voxel_basis_scaled[0] * pos.x + voxel_basis_scaled[1] * pos.y + voxel_basis_scaled[2] * pos.z) * particles[voxel_vertices[0]].curPos.w * 4;

    vec3 vertex_pos = voxel_base_pos + vertex_offset;
    gl_Position = projection * view * vec4(vertex_pos, 1.0);
    FragPos = vertex_pos;

    Normal = normalize(voxel_basis[0] * vertices[VertexID].normal.x + voxel_basis[1] * vertices[VertexID].normal.y + voxel_basis[2] * vertices[VertexID].normal.z);
    vertexIndex = VertexID;
}