#version 460 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;

uniform mat4 model;

layout(std140, binding = 0) uniform RenderCameraInfo{
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

out vec2 TexCoords;
out vec3 Normal;

out vec3 FragPos;


void main()
{
    gl_Position = projection * view * model * vec4(aPos, 1.0f);
    TexCoords = vec2(aTexCoord.x, aTexCoord.y);
    Normal = mat3(transpose(inverse(model))) * aNormal;
    FragPos = vec3(model * vec4(aPos, 1.0));
}