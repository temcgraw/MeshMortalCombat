#version 430 core


layout(std140, binding = 0) uniform RenderInfo{
    mat4 view;
    mat4 projection;
    vec3 cameraPos;
};


uniform vec3 AABB_min; // minimum vertex of the AABB
uniform vec3 AABB_max; // maximum vertex of the AABB
uniform mat4 model;

void main()
{
    const vec3 cubeVertices[24] = vec3[](
        // front face
        vec3(-1.0, -1.0, -1.0), vec3( 1.0, -1.0, -1.0),
        vec3( 1.0, -1.0, -1.0), vec3( 1.0,  1.0, -1.0),
        vec3( 1.0,  1.0, -1.0), vec3(-1.0,  1.0, -1.0),
        vec3(-1.0,  1.0, -1.0), vec3(-1.0, -1.0, -1.0),

        // back face
        vec3(-1.0, -1.0,  1.0), vec3( 1.0, -1.0,  1.0),
        vec3( 1.0, -1.0,  1.0), vec3( 1.0,  1.0,  1.0),
        vec3( 1.0,  1.0,  1.0), vec3(-1.0,  1.0,  1.0),
        vec3(-1.0,  1.0,  1.0), vec3(-1.0, -1.0,  1.0),

        // 4 connecting edges
        vec3(-1.0, -1.0, -1.0), vec3(-1.0, -1.0,  1.0),
        vec3( 1.0, -1.0, -1.0), vec3( 1.0, -1.0,  1.0),
        vec3( 1.0,  1.0, -1.0), vec3( 1.0,  1.0,  1.0),
        vec3(-1.0,  1.0, -1.0), vec3(-1.0,  1.0,  1.0)
    );

    // Map the vertex position from the [-1, 1] range to the [AABB_min, AABB_max] range
    vec3 position = 0.5 * (cubeVertices[gl_VertexID] + 1.0) * (AABB_max - AABB_min) + AABB_min;
    gl_Position = projection * view * model * vec4(position, 1.0);
}
