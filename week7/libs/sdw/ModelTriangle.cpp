#include "ModelTriangle.h"
#include <utility>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour) :
		vertices({{v0, v1, v2}}), vertex_brightness(), texturePoints(), colour(std::move(trigColour)), normal() {}

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &n0, const glm::vec3 &n1, const glm::vec3 &n2, Colour trigColour) :
		vertices({{v0, v1, v2}}), vertex_normals({{n0, n1, n2}}), vertex_brightness(), texturePoints(), colour(std::move(trigColour)), normal() {}		

std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle) {
	os << "(" << triangle.vertex_normals[0].x << ", " << triangle.vertex_normals[0].y << ", " << triangle.vertex_normals[0].z << ")\n";
	os << "(" << triangle.vertex_normals[1].x << ", " << triangle.vertex_normals[1].y << ", " << triangle.vertex_normals[1].z << ")\n";
	os << "(" << triangle.vertex_normals[2].x << ", " << triangle.vertex_normals[2].y << ", " << triangle.vertex_normals[2].z << ")\n";
	return os;
}
