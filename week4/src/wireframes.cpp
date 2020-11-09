#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <Utils.h>
#include <glm/glm.hpp> 
#include <glm/gtx/string_cast.hpp>

#define WIDTH 400
#define HEIGHT 400

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

int find_mid_x_value(CanvasPoint mid_point, CanvasPoint from_point, CanvasPoint to_point) {
	float x = from_point.x + ((mid_point.y - from_point.y) / (to_point.y - from_point.y) * (to_point.x - from_point.x));
	return x;
}

int find_mid_depth_value(CanvasPoint mid_point, CanvasPoint from_point, CanvasPoint to_point) {
	float depth = from_point.depth + ((mid_point.depth - from_point.depth) / (to_point.depth - from_point.depth) * (to_point.depth - from_point.depth));
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

	float stepX = (end.x - start.x)/(steps);
	float stepY = (end.y - start.y)/(steps);
	float stepDepth = (end.depth - start.depth)/(steps);
	
	CanvasPoint temp = start;
	result.push_back(temp);

	for (int i = 0; i < steps; i++) {
		temp.x = temp.x + stepX;
		temp.y = temp.y + stepY;
		temp.depth = temp.depth + stepDepth;

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

void fill_texture_top(CanvasTriangle triangle, TextureMap texture, DrawingWindow &window) {
	CanvasPoint top = triangle.v0();
	CanvasPoint mid = triangle.v1();
	CanvasPoint bot = triangle.v2();

	std::vector<CanvasPoint> left = interpolate_points(top, mid, mid.y-top.y+1);
	std::vector<TexturePoint> left_texture = interpolate_points(top.texturePoint, mid.texturePoint, mid.y-top.y+1);

	std::vector<CanvasPoint> right = interpolate_points(top, bot, mid.y-top.y+1);
	std::vector<TexturePoint> right_texture = interpolate_points(top.texturePoint, bot.texturePoint, mid.y-top.y+1);

	for (float i = 0.0; i < left.size(); i++) {
		int steps = (int) abs(left[i].x - right[i].x);
				
		std::vector<CanvasPoint> points = interpolate_points(left[i], right[i], steps+1);
		std::vector<TexturePoint> texturePoints = interpolate_points(left_texture[i], right_texture[i], steps+1);

		for (int c = 0; c < steps + 1; c++) {
			int x_coord = int(texturePoints.at(c).x);
			int y_coord = int(texturePoints.at(c).y);
			uint32_t colour = texture.pixels.at((y_coord*texture.width) + x_coord);

			window.setPixelColour(round(points[c].x), round(points[i].y), colour);
		}

	}
}

void fill_texture_bottom(CanvasTriangle triangle, TextureMap texture, DrawingWindow &window) {
	CanvasPoint mid = triangle.v0();
	CanvasPoint split = triangle.v1();
	CanvasPoint bot = triangle.v2();

	std::vector<CanvasPoint> left = interpolate_points(bot, mid, bot.y-mid.y+1);
	std::vector<TexturePoint> left_texture = interpolate_points(bot.texturePoint, mid.texturePoint, bot.y-mid.y+1);

	std::vector<CanvasPoint> right = interpolate_points(bot, split, bot.y-mid.y+1);
	std::vector<TexturePoint> right_texture = interpolate_points(bot.texturePoint, split.texturePoint, bot.y-mid.y+1);

	for (float i = 0.0; i < left.size(); i++) {
		int steps = (int) abs(left[i].x - right[i].x);
				
		std::vector<CanvasPoint> points = interpolate_points(left[i], right[i], steps+1);
		std::vector<TexturePoint> texturePoints = interpolate_points(left_texture[i], right_texture[i], steps+1);

		for (int c = 0; c < steps + 1; c++) {
			int x_coord = int(texturePoints.at(c).x);
			int y_coord = int(texturePoints.at(c).y);
			uint32_t colour = texture.pixels.at((y_coord*texture.width) + x_coord);

			window.setPixelColour(round(points[c].x), round(points[i].y), colour);
		}

	}
}

void map_texture(CanvasTriangle triangle, TextureMap texture, DrawingWindow &window) {
	triangle = sort_vertices(triangle);

	CanvasPoint top = triangle.v0();
	CanvasPoint mid = triangle.v1();
	CanvasPoint bot = triangle.v2();

	CanvasPoint split;
	split.y = mid.y;

	split.x = find_mid_x_value(mid, top, bot);

	float scale = (mid.y - top.y)/(bot.y-top.y);

	split.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	split.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

	CanvasTriangle top_triangle = CanvasTriangle(top, mid, split);
	fill_texture_top(top_triangle, texture, window);

	CanvasTriangle bottom_triangle = CanvasTriangle(mid, split, bot);
	fill_texture_bottom(bottom_triangle, texture, window);

	draw_unfilled_triangle(triangle, Colour(255,255,255), window);
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
			if(-1/points[c].depth > depths[round(points[c].x)][round(points[c].y)]) {
				depths[round(points[c].x)][round(points[c].y)] = -1/points[c].depth;
				uint32_t tri_colour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
				window.setPixelColour(round(points[c].x), round(points[c].y), tri_colour);
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
			if(-1/points[c].depth > depths[round(points[c].x)][round(points[c].y)]) {
				depths[round(points[c].x)][round(points[c].y)] = -1/points[c].depth;
				uint32_t tri_colour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
				window.setPixelColour(round(points[c].x), round(points[c].y), tri_colour);
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

	split.x = find_mid_x_value(mid, top, bot);

	split.depth = find_mid_depth_value(mid, top, bot);

	CanvasTriangle top_triangle = CanvasTriangle(top, mid, split);
	fill_cornell_top(top_triangle, colour, depths, window);

	CanvasTriangle bottom_triangle = CanvasTriangle(mid, split, bot);
	fill_cornell_bottom(bottom_triangle, colour, depths, window);

}

void draw_cornell(std::vector<ModelTriangle> triangles, DrawingWindow &window) {
	glm::vec3 camera(0.0, 0.0, 4.0);
	float focal_length = 500;
	std::vector<std::vector<float>> depths(window.width, std::vector<float> (window.height));

	for(int i = 0; i < triangles.size(); i++) {
		CanvasTriangle triangle;
		for(int j = 0; j < 3; j++) {
			int u = -(focal_length * triangles[i].vertices[j].x / (triangles[i].vertices[j].z - camera.z)) + (window.width / 2);
			int v = (focal_length * triangles[i].vertices[j].y / (triangles[i].vertices[j].z - camera.z)) + (window.height / 2);

			triangle[j] = CanvasPoint(u, v, triangles[i].vertices[j].z - camera.z);
		}
		fill_cornell(triangle, triangles[i].colour, depths, window);
	}
}

std::vector<ModelTriangle> read_obj(std::string file_name, float scale, std::vector<Colour> colours) {
	std::ifstream obj_file;
	std::string line;
	std::vector<glm::vec3> vertices;
	std::vector<ModelTriangle> output;
	Colour colour;

	obj_file.open(file_name);

	if(!obj_file) std::cout << "unable to read obj file" << std::endl;
	else std::cout << "obj file read successfully" << std::endl;

	while(std::getline(obj_file, line)) {
		if(line == "") continue;

		std::vector<std::string> token = split(line, ' ');

		if(token[0] == "v") {
			glm::vec3 vertex = glm::vec3(stof(token[1]) * scale, stof(token[2]) * scale, stof(token[3]) * scale);
			vertices.push_back(vertex);
		} else if(token[0] == "f") {
			std::string face1 = token[1];
			std::string face2 = token[2];
			std::string face3 = token[3];

			ModelTriangle triangle(vertices[stoi(face1) - 1], vertices[stoi(face2) - 1], vertices[stoi(face3) - 1], colour);
			output.push_back(triangle);
		} else if(token[0] == "usemtl") {
			for(int i = 0; i < colours.size(); i++) {
				if(colours[i].name == token[1]) {
					colour = Colour(colours[i]);
				}
			}
		}

	}

	obj_file.close();

	return output;
}

std::vector<Colour> read_mtl(std::string file_name) {
	std::ifstream mtl_file;
	std::string line;
	std::string colour_name;
	std::vector<Colour> output;

	mtl_file.open(file_name);

	if(!mtl_file) std::cout << "unable to read mtl file" << std::endl;
	else std::cout << "mtl file read successfully" << std::endl;

	while(std::getline(mtl_file, line)) {
		if(line == "") continue;

		std::vector<std::string> token = split(line, ' ');

		if(token[0] == "newmtl") {
			colour_name = token[1];
		} else if(token[0] == "Kd") {
			std::string r = token[1];
			std::string g = token[2];
			std::string b = token[3];

			Colour colour = Colour(colour_name, (stof(r) * 255), int(stof(g) * 255), int(stof(b) * 255));
			output.push_back(colour);
		}
	}

	mtl_file.close();

	return output;
}

void draw(DrawingWindow &window) {
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
	// Function for performing animation (shifting artifacts or moving the camera)
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_r) {
			std::vector<Colour> colours;
			std::vector<ModelTriangle> triangles;
			colours = read_mtl("cornell-box.mtl");
			triangles = read_obj("cornell-box.obj", 0.4, colours);
			draw_cornell(triangles, window);
		}
    }
	else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
