#version 460 core

const int kPointsInBinding = 2;
const int kVoxelConstraintsBinding = 7;

layout(std140, binding = 0) uniform RenderInfo{
    mat4 view;
    mat4 projection;
    vec3 cameraPos;
};

// instance properties from the particle system SSBO
struct Particle {
    vec4 curPos;   //we need curPos.xyz to get the position of the particle, and w to get the size of the particle (radius)
    vec4 prevPos;  
    vec4 vel;     
    uvec4 flags;
};

layout(std430, binding = kPointsInBinding) restrict readonly buffer ParticleBuffer {
    Particle particles[];
};
layout (std430, binding = kVoxelConstraintsBinding) buffer VOXELS
{
	int voxelConstraints[];
};

out vec3 Normal;

out vec3 FragPos;

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
    // each voxel is a cube with 8 vertices
    // and 12 triangles
    // so each voxel has 36 vertices in drawcall
    int vox_id = int(gl_VertexID) / 36;
    if(vox_id >= voxelConstraints.length()/8) return;
    int voxel_vertex_base = vox_id * 8;
    int voxel_vertices[8];
    for(int i=0; i<8; i++)
    {
        voxel_vertices[i] = voxelConstraints[voxel_vertex_base + i];
    }
    // for each voxel's 36 vertices, determine the belonging triangle id among 12 triangles
    int triangle_id = (int(gl_VertexID) % 36) / 3;// 0~11
    ivec3 face_in_voxel = face_indices[triangle_id];// the three vertices of the triangle in local indices of the voxel cube
    ivec3 triangle_vertex_indices = ivec3(voxel_vertices[face_in_voxel.x], voxel_vertices[face_in_voxel.y], voxel_vertices[face_in_voxel.z]);// the three vertice, denoted by their global indices in the particle array
    // finally, for each triangle, determine the vertex this thread is responsible for
    int vertex_id = int(gl_VertexID) % 3;
    int vertex_index = triangle_vertex_indices[vertex_id];
    vec3 vertex = particles[vertex_index].curPos.xyz;
    
    // shift the vertex to the corner of the voxel
    // first, get the three basis vectors of the voxel
    vec3 voxel_basis[3];
    
    voxel_basis[0] = normalize(particles[voxel_vertices[1]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz);
    voxel_basis[1] = normalize(particles[voxel_vertices[2]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz);
    voxel_basis[2] = normalize(particles[voxel_vertices[4]].curPos.xyz - particles[voxel_vertices[0]].curPos.xyz);


    int vertex_id_in_voxel = face_indices[triangle_id][vertex_id];// 0~7, the local index of the vertex inside the voxel
    vec3 vertex_offset_voxel_space = voffset[vertex_id_in_voxel];// in each voxel's local space, the offset of the particle center from the corresponding corner of the voxel

    // then, get the offset of the vertex from the corner of the voxel
    vec3 vertex_offset = (vertex_offset_voxel_space.x * voxel_basis[0] + vertex_offset_voxel_space.y * voxel_basis[1] +vertex_offset_voxel_space.z * voxel_basis[2])*particles[voxel_vertices[0]].curPos.w*1.005;
    // then, shift the vertex to the corner of the voxel
    vertex += vertex_offset;
    


    gl_Position = projection * view * vec4(vertex, 1.0f);
    FragPos = vertex;
}