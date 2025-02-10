#version 430 core
out vec4 FragColor;

in vec2 TexCoords;
uniform bool flipYCoord = false;

uniform sampler2D texture1;

void main()
{
    // commonly we do not flip the y cords in opengl
    if (!flipYCoord){
        FragColor = texture(texture1, TexCoords);
    }
    // if such texture is a png or other image format that has the y axis flipped
    else{
        FragColor = texture(texture1, vec2(TexCoords.x, 1.0 - TexCoords.y));
    }

}