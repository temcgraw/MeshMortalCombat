#version 430 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;


layout(std140, binding = 0) uniform RenderInfo{
    mat4 view;
    mat4 projection;
    vec3 cameraPos;
};

// instance properties from the particle system SSBO
struct Particle {
    vec4 curPos;   //we need curPos.xyz to get the position of the particle, and w to get the size of the particle
    vec4 prevPos;  
    vec4 vel;     
    uvec4 flags;
};

layout(std430, binding = 2) buffer ParticleBuffer {
    Particle particles[];
};

out vec2 TexCoords;
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


void main()
{
    mat4 model = mat4(1.0);
    model = translate(model, particles[gl_InstanceID].curPos.xyz);
    model = scale(model, particles[gl_InstanceID].curPos.www);
    gl_Position = projection * view * model * vec4(aPos, 1.0f);
    TexCoords = vec2(aTexCoord.x, aTexCoord.y);
    Normal = mat3(transpose(inverse(model))) * aNormal;
    FragPos = vec3(model * vec4(aPos, 1.0));
}