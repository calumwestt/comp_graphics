#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <Utils.h>
#include <glm/glm.hpp> 
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>

#define WIDTH 400
#define HEIGHT 400

#define PI atan(1)*4

glm::vec3 camera(0.0, 0.0, 3.0);
glm::vec3 light(0.0, 1.0, 0.0);
glm::mat3 camera_orientation(
	glm::vec3(1.0, 0.0, 0.0),
	glm::vec3(0.0, 1.0, 0.0),
	glm::vec3(0.0, 0.0, 1.0)
);
float focal_length = 300;

float light_strength = 30;
bool incidence = false;

std::string drawing_method = "rasterise";
std::string shading_method = "gouraurd";
bool orbiting = false;

std::vector<std::vector<glm::vec3>> all_vertex_normals;

void draw_line(CanvasPoint from_point, CanvasPoint to_point, Colour line_colour, DrawingWindow &window) {
	float x_diff = to_point.x - from_point.x;
	float y_diff = to_point.y - from_point.y;
	float number_of_steps = std::max(abs(x_diff), abs(y_diff));
	float x_step_size = x_diff/number_of_steps;
	float y_step_size = y_diff/number_of_steps;
	uint32_t colour = (255 << 24) + (int(line_colour.red) << 16) + (int(line_colour.green) << 8) + int(line_colour.blue);
	for (float i=0; i<number_of_steps; i++) {
		float x = from_point.x + (x_step_size * i);
		float y = from_point.y + (y_step_size * i);
		window.setPixelColour(round(x), round(y), colour);
	}
}

void draw_unfilled_triangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window) {
	draw_line(triangle.v0(), triangle.v1(), colour, window);
	draw_line(triangle.v1(), triangle.v2(), colour, window);
	draw_line(triangle.v2(), triangle.v0(), colour, window);
}

void fill_triangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window) {
	CanvasPoint top = triangle.v0();
	CanvasPoint mid = triangle.v1();
	CanvasPoint bot = triangle.v2();

	float x1_diff = mid.x - top.x;
	float y1_diff = mid.y - top.y;
	float x2_diff = bot.x - top.x;
	float y2_diff = bot.y - top.y;

	float steps1 = std::max(abs(x1_diff),abs(y1_diff));
	float steps2 = std::max(abs(x2_diff),abs(y2_diff));

	float x1_step = x1_diff/steps1;
	float x2_step = x2_diff/steps2;
	float y1_step = y1_diff/steps1;
	float y2_step = y2_diff/steps2;

	for(float i = 0.0; i < steps1; i++) {
		for(float j = 0.0; j < steps2; j++) {
			float x1 = top.x + (x1_step * i);
			float y1 = top.y + (y1_step * i);
			float x2 = top.x + (x2_step * j);
			float y2 = top.y + (y2_step * j);
			draw_line(CanvasPoint(round(x1),round(y1)),CanvasPoint(round(x2),round(y2)), colour, window);
		}
	}
}

float find_mid_x_value(CanvasPoint mid_point, CanvasPoint from_point, CanvasPoint to_point) {
	float x = from_point.x + ((mid_point.y - from_point.y) / (to_point.y - from_point.y) * (to_point.x - from_point.x));
	return x;
}

float find_mid_depth_value(CanvasPoint mid_point, CanvasPoint from_point, CanvasPoint to_point) {
	float depth = from_point.depth + ((mid_point.y - from_point.y) / (to_point.y - from_point.y) * (to_point.depth - from_point.depth));
	return depth;
}

CanvasTriangle sort_vertices(CanvasTriangle triangle) {
	if(triangle.v1().y < triangle.v0().y) {
		std::swap(triangle.v0(), triangle.v1());
	}
	if(triangle.v2().y < triangle.v1().y) {
		std::swap(triangle.v1(), triangle.v2());
		if(triangle.v1().y < triangle.v0().y) {
			std::swap(triangle.v0(), triangle.v1());
		}
	}
	return triangle;
}

void draw_filled_triangle(CanvasTriangle triangle, Colour triangle_colour, DrawingWindow &window) {
	triangle = sort_vertices(triangle);

	fill_triangle(triangle, triangle_colour, window);
}

std::vector<CanvasPoint> interpolate_points(CanvasPoint start, CanvasPoint end, int steps) {
	std::vector<CanvasPoint> result;

	float x_step = (end.x - start.x)/(steps);
	float y_step = (end.y - start.y)/(steps);
	float depth_step = (-1/end.depth - -1/start.depth)/(steps);

	for (int i = 0; i < steps+1; i++) {
		CanvasPoint temp = CanvasPoint(start.x + (x_step * i), start.y + (y_step * i), -1/start.depth + (depth_step * i));
		result.push_back(temp);
	}

	return result;
}

std::vector<TexturePoint> interpolate_points(TexturePoint start, TexturePoint end, int steps) {
	std::vector<TexturePoint> result;
	
	for (size_t i = 0; i < steps + 1; i++) {
		TexturePoint p;
		if (steps == 1) steps = 2;
		p.x = start.x + ((end.x - start.x) * i / (steps));
		p.y = start.y + ((end.y - start.y) * i / (steps));
		result.push_back(p);
	}

	return result;
}

void fill_texture_top(CanvasTriangle triangle, TextureMap texture, std::vector<std::vector<float>> &depths, DrawingWindow &window) {
	CanvasPoint top = triangle.v0();
	CanvasPoint mid = triangle.v1();
	CanvasPoint bot = triangle.v2();

	std::vector<CanvasPoint> left = interpolate_points(top, mid, mid.y-top.y+1);
	std::vector<TexturePoint> left_texture = interpolate_points(top.texturePoint, mid.texturePoint, mid.y-top.y+1);

	std::vector<CanvasPoint> right = interpolate_points(top, bot, mid.y-top.y+1);
	std::vector<TexturePoint> right_texture = interpolate_points(top.texturePoint, bot.texturePoint, mid.y-top.y+1);

	for (float i = 0.0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);
				
		std::vector<CanvasPoint> points = interpolate_points(left[i], right[i], steps+1);
		std::vector<TexturePoint> texturePoints = interpolate_points(left_texture[i], right_texture[i], steps+1);

		for (int c = 0; c < points.size(); c++) {
			int rounded_x_texture = round(texturePoints.at(c).x);
			int rounded_y_texture = round(texturePoints.at(c).y);
			int rounded_x = round(points[c].x);
			int rounded_y = round(points[c].y);
			if(rounded_x >= 0 && rounded_x < window.width && rounded_y >= 0 && rounded_y < window.height) {
				if(-1/points[c].depth > depths[rounded_x][rounded_y]) {
					depths[rounded_x][rounded_y] = -1/points[c].depth;
					uint32_t colour = texture.pixels.at((rounded_y_texture*texture.width) + rounded_x_texture);
					window.setPixelColour(rounded_x, rounded_y, colour);
				}	
			}
		}

	}
}

void fill_texture_bottom(CanvasTriangle triangle, TextureMap texture, std::vector<std::vector<float>> &depths, DrawingWindow &window) {
	CanvasPoint mid = triangle.v0();
	CanvasPoint split = triangle.v1();
	CanvasPoint bot = triangle.v2();

	std::vector<CanvasPoint> left = interpolate_points(bot, mid, bot.y-mid.y+1);
	std::vector<TexturePoint> left_texture = interpolate_points(bot.texturePoint, mid.texturePoint, bot.y-mid.y+1);

	std::vector<CanvasPoint> right = interpolate_points(bot, split, bot.y-mid.y+1);
	std::vector<TexturePoint> right_texture = interpolate_points(bot.texturePoint, split.texturePoint, bot.y-mid.y+1);

	for (float i = 0.0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);
				
		std::vector<CanvasPoint> points = interpolate_points(left[i], right[i], steps+1);
		std::vector<TexturePoint> texturePoints = interpolate_points(left_texture[i], right_texture[i], steps+1);

		for (int c = 0; c < points.size(); c++) {
			int rounded_x_texture = round(texturePoints.at(c).x);
			int rounded_y_texture = round(texturePoints.at(c).y);
			int rounded_x = round(points[c].x);
			int rounded_y = round(points[c].y);
			if(rounded_x >= 0 && rounded_x < window.width && rounded_y >= 0 && rounded_y < window.height) {
				if(-1/points[c].depth > depths[rounded_x][rounded_y]) {
					depths[rounded_x][rounded_y] = -1/points[c].depth;
					uint32_t colour = texture.pixels.at((rounded_y_texture*texture.width) + rounded_x_texture);
					window.setPixelColour(rounded_x, rounded_y, colour);
				}	
			}
		}

	}
}

void map_texture(CanvasTriangle triangle, TextureMap texture, std::vector<std::vector<float>> &depths, DrawingWindow &window) {
	triangle = sort_vertices(triangle);

	CanvasPoint top = triangle.v0();
	CanvasPoint mid = triangle.v1();
	CanvasPoint bot = triangle.v2();

	CanvasPoint split;
	split.y = mid.y;

	split.x = round(find_mid_x_value(mid, top, bot));
	split.depth = find_mid_depth_value(mid, top, bot);

	float scale = (mid.y - top.y)/(bot.y-top.y);

	split.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	split.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

	CanvasTriangle top_triangle = CanvasTriangle(top, mid, split);
	fill_texture_top(top_triangle, texture, depths, window);

	CanvasTriangle bottom_triangle = CanvasTriangle(mid, split, bot);
	fill_texture_bottom(bottom_triangle, texture, depths, window);
}

void fill_cornell_top(CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depths, DrawingWindow &window) {
	CanvasPoint top = triangle.v0();
	CanvasPoint mid = triangle.v1();
	CanvasPoint split = triangle.v2();

	std::vector<CanvasPoint> left = interpolate_points(top, mid, mid.y-top.y+1);
	std::vector<CanvasPoint> right = interpolate_points(top, split, mid.y-top.y+1);

	for (float i = 0.0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);
		std::vector<CanvasPoint> points = interpolate_points(left[i], right[i], steps+1);
		for (int c = 0; c < points.size(); c++) {
			int rounded_x = round(points[c].x);
			int rounded_y = round(points[c].y);
			if(rounded_x >= 0 && rounded_x < window.width && rounded_y >= 0 && rounded_y < window.height) {
				if(-1/points[c].depth > depths[rounded_x][rounded_y]) {
					depths[rounded_x][rounded_y] = -1/points[c].depth;
					uint32_t tri_colour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
					window.setPixelColour(rounded_x, rounded_y, tri_colour);
				}	
			}
		}

	}
}

void fill_cornell_bottom(CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depths, DrawingWindow &window) {
	CanvasPoint mid = triangle.v0();
	CanvasPoint split = triangle.v1();
	CanvasPoint bot = triangle.v2();

	std::vector<CanvasPoint> left = interpolate_points(bot, mid, bot.y-mid.y+1);
	std::vector<CanvasPoint> right = interpolate_points(bot, split, bot.y-mid.y+1);

	for (float i = 0.0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);	
		std::vector<CanvasPoint> points = interpolate_points(left[i], right[i], steps+1);
		for (int c = 0; c < points.size(); c++) {
			int rounded_x = round(points[c].x);
			int rounded_y = round(points[c].y);
			if(rounded_x >= 0 && rounded_x < window.width && rounded_y >= 0 && rounded_y < window.height) {
				if(-1/points[c].depth > depths[round(points[c].x)][round(points[c].y)]) {
					depths[round(points[c].x)][round(points[c].y)] = -1/points[c].depth;
					uint32_t tri_colour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
					window.setPixelColour(round(points[c].x), round(points[c].y), tri_colour);
				}
			}
		}

	}
}

void fill_cornell(CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depths, DrawingWindow &window) {
	triangle = sort_vertices(triangle);

	CanvasPoint top = triangle.v0();
	CanvasPoint mid = triangle.v1();
	CanvasPoint bot = triangle.v2();

	CanvasPoint split;
	split.y = mid.y;

	split.x = round(find_mid_x_value(mid, top, bot));

	split.depth = find_mid_depth_value(mid, top, bot);

	CanvasTriangle top_triangle = CanvasTriangle(top, mid, split);
	fill_cornell_top(top_triangle, colour, depths, window);

	CanvasTriangle bottom_triangle = CanvasTriangle(mid, split, bot);
	fill_cornell_bottom(bottom_triangle, colour, depths, window);

}

void wireframe_draw_cornell(std::vector<ModelTriangle> triangles, DrawingWindow &window) {
	for(int i = 0; i < triangles.size(); i++) {
		CanvasTriangle triangle;
		for(int j = 0; j < 3; j++) {
			glm::vec3 camera_to_vertex = glm::vec3(triangles[i].vertices[j].x - camera.x, triangles[i].vertices[j].y - camera.y, triangles[i].vertices[j].z - camera.z);

			glm::vec3 adjusted_coords = camera_to_vertex * camera_orientation;

			int u = -(focal_length * (adjusted_coords.x)/(adjusted_coords.z)) + (window.width / 2);
			int v = (focal_length * (adjusted_coords.y)/(adjusted_coords.z)) + (window.height / 2);

			triangle.vertices[j] = CanvasPoint(u, v, adjusted_coords.z);
		}
		draw_unfilled_triangle(triangle, triangles[i].colour, window);
	}

	glm::vec3 cam_to_vertex = glm::vec3(light.x - camera.x, light.y - camera.y, light.z - camera.z);
    glm::vec3 adjusted_vector = cam_to_vertex * camera_orientation;

    int u = -(focal_length * (adjusted_vector.x)/(adjusted_vector.z)) + (window.width / 2);
    int v = (focal_length * (adjusted_vector.y)/(adjusted_vector.z)) + (window.height / 2);

    // prints red pixels to show light location
    window.setPixelColour(u,   v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
    window.setPixelColour(u+1, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
    window.setPixelColour(u, v+1, (255 << 24) + (255 << 16) + (0 << 8) + 0);
    window.setPixelColour(u-1, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
    window.setPixelColour(u, v-1, (255 << 24) + (255 << 16) + (0 << 8) + 0);
}

void rasterise_draw_cornell(std::vector<ModelTriangle> triangles, DrawingWindow &window) {
	std::vector<std::vector<float>> depths(window.width, std::vector<float> (window.height, -std::numeric_limits<float>::infinity()));

	for (int i = 0; i < triangles.size(); i++) {
		CanvasTriangle triangle;
		bool is_texture = false;
		TextureMap texture;

		if (triangles[i].colour.name != "") {
			texture = TextureMap(triangles[i].colour.name);
			is_texture = true;
		}
		for (int j = 0; j < 3; j++) {
			glm::vec3 camera_to_vertex = glm::vec3(triangles[i].vertices[j].x - camera.x, triangles[i].vertices[j].y - camera.y, triangles[i].vertices[j].z - camera.z);

			glm::vec3 adjusted_coords = camera_to_vertex * camera_orientation;

			int u = -(focal_length * (adjusted_coords.x)/(adjusted_coords.z)) + (window.width / 2);
			int v = (focal_length * (adjusted_coords.y)/(adjusted_coords.z)) + (window.height / 2);

			triangle.vertices[j] = CanvasPoint(u, v, adjusted_coords.z);

			if (is_texture) {
				triangle.vertices[j].texturePoint = triangles[i].texturePoints[j];
				triangle.vertices[j].texturePoint.x *= texture.width;
				triangle.vertices[j].texturePoint.y *= texture.height;

			}	
		}

		if (is_texture) map_texture(triangle, texture, depths, window);
		else fill_cornell(triangle, triangles[i].colour, depths, window);

		
		
 	}
}

bool is_shadow(RayTriangleIntersection intersection, std::vector<ModelTriangle> triangles) {	
	glm::vec3 shadow_ray = light - intersection.intersectionPoint;

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle triangle = triangles[i];

		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 sp_vector = intersection.intersectionPoint - triangle.vertices[0];
		glm::mat3 de_matrix(-glm::normalize(shadow_ray), e0, e1);
		glm::vec3 possible_solution = glm::inverse(de_matrix) * sp_vector;

		float t = possible_solution.x;
		float u = possible_solution.y;
		float v = possible_solution.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v <= 1.0)) {
			if(t < glm::length(shadow_ray) && t > 0.05 && i != intersection.triangleIndex) return true;
		}
	}
	return false;
}

double get_triangle_brightness(glm::vec3 intersection_point, glm::vec3 normal) {
	glm::vec3 ray = light - intersection_point;
	glm::vec3 viewing_angle = camera - intersection_point;
	glm::vec3 reflection = glm::normalize(ray) - ((2.0f * normal) * (glm::dot(glm::normalize(ray), normal)));

	double length = glm::length(ray);
	double brightness = light_strength / (4 * PI * (length * length));

	double angle_of_incidence = glm::dot(glm::normalize(ray), normal);

	double specular = glm::dot(glm::normalize(reflection), glm::normalize(viewing_angle));
	specular = std::pow(specular, 64);

	if (angle_of_incidence > 0) {
		brightness *= angle_of_incidence;
	}

	if (specular > 0) {
		brightness += specular;
	} 

	if (brightness > 1) {
		brightness = 1;
	}

	return brightness;
}

double gouraurd(RayTriangleIntersection intersection) {
	std::vector<float> brightnesses;
	std::vector<glm::vec3> reflections;

	ModelTriangle triangle = intersection.intersectedTriangle;

    glm::vec3 light_ray = light - intersection.intersectionPoint;
	glm::vec3 viewing_angle = camera - intersection.intersectionPoint;
    float length = glm::length(light_ray);

    for(int i = 0; i < 3; i++) {
        double temp_brightness = glm::dot(triangle.vertex_normals[i], glm::normalize(light_ray));
        brightnesses.push_back(temp_brightness);

		glm::vec3 temp_reflection = glm::normalize(light_ray) - ((2.0f * triangle.vertex_normals[i]) * glm::dot(glm::normalize(light_ray), triangle.vertex_normals[i]));
        reflections.push_back(temp_reflection);
    }

    double brightness = (1 - intersection.scale_u - intersection.scale_v) * brightnesses[0] + intersection.scale_u * brightnesses[1] + intersection.scale_v * brightnesses[2];

    glm::vec3 angle_of_reflection = (1 - intersection.scale_u - intersection.scale_v) * reflections[0] + intersection.scale_u * reflections[1] + intersection.scale_v * reflections[2];
    float specular = std::pow(glm::dot(glm::normalize(angle_of_reflection), glm::normalize(viewing_angle)), 256);

    if (specular) {
        brightness += specular;
    }

    if (brightness > 1) {
        brightness = 1;
    } 
    
    return brightness;
}

double phong(RayTriangleIntersection intersection) {
	ModelTriangle triangle = intersection.intersectedTriangle;

	glm::vec3 normal = (1 - intersection.scale_u - intersection.scale_v) * triangle.vertex_normals[0] + intersection.scale_u * triangle.vertex_normals[1] + intersection.scale_v * triangle.vertex_normals[2];
	glm::vec3 light_ray = light - intersection.intersectionPoint;
	glm::vec3 viewing_angle = glm::normalize(camera - intersection.intersectionPoint);
	glm::vec3 reflection_ray = glm::normalize(glm::normalize(light_ray) - (normal * 2.0f * glm::dot(glm::normalize(light_ray), normal)));
	double length = glm::length(light_ray);

	double angle_of_incidence = glm::dot(normal, glm::normalize(light_ray));
	double brightness = light_strength * angle_of_incidence / (4 * PI * (length * length));
	double specular = std::pow(glm::dot(reflection_ray, viewing_angle), 256);

	if (specular > 0) brightness += specular;
	
	if (brightness > 1) brightness = 1;

	return brightness;

}

RayTriangleIntersection get_closest_intersection(glm::vec3 direction, std::vector<ModelTriangle> triangles) {	
	RayTriangleIntersection intersection;
	intersection.distanceFromCamera = std::numeric_limits<float>::infinity();
	glm::vec3 ray_direction = camera_orientation * (camera - direction); 

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle triangle = triangles[i];

		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 sp_vector = camera - triangle.vertices[0];
		glm::mat3 de_matrix(-ray_direction, e0, e1);
		glm::vec3 possible_solution = glm::inverse(de_matrix) * sp_vector;

		float t = possible_solution.x;
		float u = possible_solution.y;
		float v = possible_solution.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v <= 1.0)) {
			if((intersection.distanceFromCamera > t) && (t > 0)) {
				intersection.distanceFromCamera = t;
				intersection.intersectedTriangle = triangle;
				intersection.triangleIndex = i;
				intersection.intersectionPoint = glm::vec3(triangle.vertices[0] + u * e0 + v * e1);
				intersection.scale_u = u;
				intersection.scale_v = v;
			}
		}
	}
	return intersection;
}

void raytrace_draw_cornell(std::vector<ModelTriangle> triangles, DrawingWindow &window) {
	for(int x = 0; x < window.width; x++) {
		for(int y = 0; y < window.height; y++) {
			glm::vec3 direction = glm::vec3(int(window.width / 2) - x, y - int(window.height / 2), focal_length);
			RayTriangleIntersection intersection = get_closest_intersection(direction, triangles);

			if(!isinf(intersection.distanceFromCamera)) {
				Colour colour = intersection.intersectedTriangle.colour;
				double brightness;
				if (shading_method == "gouraurd") brightness = gouraurd(intersection);
				else if (shading_method == "phong") brightness = phong(intersection);

				if(is_shadow(intersection, triangles)) {
					double shadow_brightness = brightness / 3;
					if (shadow_brightness < 0.15) shadow_brightness = 0.15;
					uint32_t shadow_colour = (255 << 24) + (int(colour.red * shadow_brightness) << 16) + (int(colour.green * shadow_brightness) << 8) + int(colour.blue * shadow_brightness);
					window.setPixelColour(x, y, shadow_colour);
				}
				else {
					if (brightness < 0.15) brightness = 0.15;
					uint32_t tri_colour = (255 << 24) + (int(colour.red * brightness) << 16) + (int(colour.green * brightness) << 8) + int(colour.blue * brightness);
					window.setPixelColour(x, y, tri_colour);
				}
			}
		}
	}
}

void rotate_x(float theta) {
	glm::mat3 rotation_matrix = glm::mat3(
				1, 0, 0,
				0, cos(theta), -sin(theta),
				0, sin(theta), cos(theta)
			);
	camera = camera * rotation_matrix;
}

void rotate_y(float theta) {
	glm::mat3 rotation_matrix = glm::mat3(
				cos(theta), 0, -sin(theta),
				0, 1, 0,
				sin(theta), 0, cos(theta)
			);
	camera = camera * rotation_matrix;
}

void tilt(float theta) {
	glm::mat3 rotation_matrix(
				glm::vec3(1, 0, 0),
				glm::vec3(0, cos(theta), -sin(theta)),
				glm::vec3(0, sin(theta), cos(theta))
								
			);
	camera_orientation *= rotation_matrix;
}

void pan(float theta) {
	glm::mat3 rotation_matrix(
				glm::vec3(cos(theta), 0, -sin(theta)),
				glm::vec3(0, 1, 0),
				glm::vec3(sin(theta), 0, cos(theta))
								
			);
	camera_orientation *= rotation_matrix;
}

void reset_cam() {
	camera = glm::vec3(0.0, 0.0, 4.0);
	light = glm::vec3(0.0, 0.9, 0.0);
	camera_orientation = glm::mat3(glm::vec3(1.0, 0.0, 0.0),
								   glm::vec3(0.0, 1.0, 0.0),
								   glm::vec3(0.0, 0.0, 1.0));
}

void look_at() {
	glm::vec3 forward = glm::normalize(camera - glm::vec3(0.0, 0.0, 0.0));
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0.0, 1.0, 0.0), forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));
	camera_orientation = glm::mat3(right,
								   up,
								   forward);
}

void orbit(bool orbiting) {
	if (orbiting) {
		float theta = PI / 180;
		glm::mat3 rotation_matrix = glm::mat3(
					cos(theta), 0, -sin(theta),
					0, 1, 0,
					sin(theta), 0, cos(theta)
				);
		camera = camera * rotation_matrix;
		look_at();
	}
}

std::vector<ModelTriangle> calc_vertex_normals(std::vector<ModelTriangle> triangles) {

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle tri = triangles[i];
		std::vector<glm::vec3> vertex_normals;
		for(int v = 0; v < tri.vertices.size(); v++) {
			glm::vec3 vertex = tri.normal;
			int count = 1;
			for(int j = 0; j < triangles.size(); j++) {
				ModelTriangle tri_ = triangles[j];
				for(int u = 0; u < tri_.vertices.size(); u++) {
					if(i != j && tri.vertices[v].x == tri_.vertices[u].x && tri.vertices[v].y == tri_.vertices[u].y && tri.vertices[v].z == tri_.vertices[u].z) {
						if (std::acos(glm::dot(tri.normal, tri_.normal) / (length(tri.normal) * length(tri_.normal))) < PI / 4) {
							vertex = vertex + tri_.normal;
							count = count + 1;
						}
					}
				}
			}
			vertex = vertex / float(count);
			triangles[i].vertex_normals[v] = normalize(vertex);
		}
	}

	return triangles;
}

std::vector<ModelTriangle> read_obj(std::string file_name, float scale, std::unordered_map<std::string, Colour> colours) {
	std::ifstream obj_file;
	std::string line;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> vertex_normals;
	std::vector<TexturePoint> texture_vertices;
	std::vector<ModelTriangle> output;
	std::string colour;

	obj_file.open(file_name);

	if(!obj_file) std::cout << "unable to read obj file" << std::endl;
	else std::cout << "obj file read successfully" << std::endl;

	while(std::getline(obj_file, line)) {
		if(line == "") continue;

		std::vector<std::string> token = split(line, ' ');

		if(token[0] == "v") {
			glm::vec3 vertex = glm::vec3(stof(token[1]) * scale, stof(token[2]) * scale, stof(token[3]) * scale);
			vertices.push_back(vertex);
		}
		else if(token[0] == "vn") {
			glm::vec3 vertex_normal = glm::vec3(stof(token[1]), stof(token[2]), stof(token[3]));
			vertex_normals.push_back(vertex_normal);
		} 
		else if(token[0] == "vt") {
			texture_vertices.push_back(TexturePoint(stof(token[1]), stof(token[2])));
		}
		else if(token[0] == "f") {
			std::string face1 = token[1];
			std::string face2 = token[2];
			std::string face3 = token[3];

			std::vector<std::string> new_token = split(token[1], '/');
			std::vector<std::string> new_token_2 = split(token[2], '/');
			std::vector<std::string> new_token_3 = split(token[3], '/');

			if(new_token[1] == "") {
				ModelTriangle triangle(vertices[stoi(face1) - 1], vertices[stoi(face2) - 1], vertices[stoi(face3) - 1], colours[colour]);
				if(!vertex_normals.empty()) {
					triangle.vertex_normals[0] = vertex_normals[stoi(face1) - 1]; 
					triangle.vertex_normals[1] = vertex_normals[stoi(face2) - 1]; 
					triangle.vertex_normals[2] = vertex_normals[stoi(face3) - 1]; 
				}
				triangle.normal = glm::normalize(glm::cross(glm::vec3(triangle.vertices[1] - triangle.vertices[0]), glm::vec3(triangle.vertices[2] - triangle.vertices[0])));
				output.push_back(triangle);
			}
			else {
				ModelTriangle triangle(vertices[stoi(face1) - 1], vertices[stoi(face2) - 1], vertices[stoi(face3) - 1], colours[colour]);
				triangle.texturePoints[0] = texture_vertices[stof(new_token[1]) - 1];
				triangle.texturePoints[1] = texture_vertices[stof(new_token_2[1]) - 1];
				triangle.texturePoints[2] = texture_vertices[stof(new_token_3[1]) - 1];
				output.push_back(triangle);
			}
		}
		else if(token[0] == "usemtl") {
			colour = token[1];
		}
	}

	if(vertex_normals.empty()) {
		std::cout << ".obj DID NOT contain any vertex normals" << std::endl;
		output = calc_vertex_normals(output);
	} 
	else {
		std::cout << ".obj DID contain vertex normals" << std::endl;
	}

	obj_file.close();

	return output;
}

std::unordered_map<std::string, Colour> read_mtl(std::string file_name) {
	std::ifstream mtl_file;
	std::string line;
	std::string colour_name;
	std::unordered_map<std::string, Colour> output;

	mtl_file.open(file_name);

	if(!mtl_file) std::cout << "unable to read mtl file" << std::endl;
	else std::cout << "mtl file read successfully" << std::endl;

	while(std::getline(mtl_file, line)) {
		if(line == "") continue;

		std::vector<std::string> token = split(line, ' ');

		if(token[0] == "newmtl") {
			colour_name = token[1];
		} 
		else if(token[0] == "Kd") {
			std::string r = token[1];
			std::string g = token[2];
			std::string b = token[3];

			Colour colour = Colour((stof(r) * 255), int(stof(g) * 255), int(stof(b) * 255));
			output.insert({colour_name, colour});
		}
		else if(token[0] == "map_Kd") {
			Colour colour = output[colour_name];
			colour.name = token[1];
			output[colour_name] = colour;
		}
	}

	mtl_file.close();

	return output;
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
		}
	}

}

void update(DrawingWindow &window) {
	
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		//move camera left
		if (event.key.keysym.sym == SDLK_LEFT) light.x -= 0.1;
		//move camera right
		else if (event.key.keysym.sym == SDLK_RIGHT) light.x += 0.1;
		//move camera up
		else if (event.key.keysym.sym == SDLK_UP) light.y += 0.1;
		//move camera down
		else if (event.key.keysym.sym == SDLK_DOWN) light.y -= 0.1;
		//move camera forwards
		else if (event.key.keysym.sym == SDLK_w) light.z -= 0.1;
		//move camera backwards
		else if (event.key.keysym.sym == SDLK_s) light.z += 0.1;
		//rotate camera about x axis
		else if (event.key.keysym.sym == SDLK_a) rotate_x(PI / 180);
		//rotate camera about x axis
		else if (event.key.keysym.sym == SDLK_d) rotate_x(-PI / 180);
		//rotate camera about y axis
		else if (event.key.keysym.sym == SDLK_q) rotate_y(PI / 180);
		//rotate camera about y axis
		else if (event.key.keysym.sym == SDLK_e) rotate_y(-PI / 180);
		//tilt camera up
		else if (event.key.keysym.sym == SDLK_u) tilt(PI / 180);
		//tilt camera down
		else if (event.key.keysym.sym == SDLK_j) tilt(-PI / 180);
		//pan camera left
		else if (event.key.keysym.sym == SDLK_h) pan(PI / 180);
		//pan camera right
		else if (event.key.keysym.sym == SDLK_k) pan(-PI / 180);
		//reset camera + light
		else if (event.key.keysym.sym == SDLK_r) reset_cam();
		//orbit + look at
		else if (event.key.keysym.sym == SDLK_o) orbiting = (orbiting) ? false : true;
		//draw raytrace
		else if (event.key.keysym.sym == SDLK_1) { drawing_method = "raytrace"; std::cout << "Drawing using raytrace" << std::endl; }
		//draw raytrace
		else if (event.key.keysym.sym == SDLK_2) { drawing_method = "rasterise"; std::cout << "Drawing using rasterise" << std::endl; }
		//draw wireframe
		else if (event.key.keysym.sym == SDLK_3) { drawing_method = "wireframe"; std::cout << "Drawing using wireframe" << std::endl; }
		//turn on gouraurd
		else if (event.key.keysym.sym == SDLK_g) { shading_method = "gouraurd"; std::cout << "Drawing using gouraurd" << std::endl; }
		//turn on phong
		else if (event.key.keysym.sym == SDLK_p) { shading_method = "phong"; std::cout << "Drawing using phong" << std::endl; }
    }
	else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	std::unordered_map<std::string, Colour> colours;
	std::vector<ModelTriangle> triangles;
	colours = read_mtl("cornell-box.mtl");

	triangles = read_obj("cornell-box.obj", 0.4, colours);

	//triangles = read_obj("sphere.obj", 0.4, colours);

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);

		orbit(orbiting);

		if(drawing_method == "rasterise") rasterise_draw_cornell(triangles, window);
		else if(drawing_method == "raytrace") raytrace_draw_cornell(triangles, window);
		else if(drawing_method == "wireframe") wireframe_draw_cornell(triangles, window);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
