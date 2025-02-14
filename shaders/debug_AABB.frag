#version 460 core
out vec4 FragColor;

uniform vec4 lineColor = vec4(1.0, 1.0, 1.0, 1.0);

void main()
{
    FragColor = lineColor;
}
