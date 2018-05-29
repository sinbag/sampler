#version 330

uniform mat4 u_projectionMatrix;

in vec2 a_position;
in vec4 a_color;

out vertexData
{
	vec4 color;
} vertex;

void main(void)
{
	vertex.color = a_color;
	gl_Position = u_projectionMatrix*vec4(a_position, 0.0, 1.0);
}
