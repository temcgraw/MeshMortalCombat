#version 460 core
layout (location = 0) in vec3 aPos;

out vec3 TexCoords;

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

void main()
{
    //remove translation from view matrix
    mat3 rotationMatrix = mat3(view); // 移除位移，保留旋转和缩放

    TexCoords = aPos;
    vec4 pos = projection * mat4(rotationMatrix) * vec4(aPos, 1.0);
    gl_Position = pos.xyww;
}  