#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>

#define WIDTH 320
#define HEIGHT 240

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

int find_mid_value(CanvasPoint mid_point, CanvasPoint from_point, CanvasPoint to_point) {
	float x = from_point.x + ((mid_point.y - from_point.y) / (to_point.y - from_point.y) * (to_point.x - from_point.x));
	return x;
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

	//this value isn't used anymore

	//int x = find_mid_value(triangle.v1(), triangle.v0(), triangle.v2());

	//we don't need to do this because for some reason 'fill_triangle' works when you pass in the entire triangle
	
	/*CanvasTriangle top_triangle = CanvasTriangle(triangle.v0(), CanvasPoint(x,triangle.v1().y), triangle.v1());
	CanvasTriangle bottom_triangle = CanvasTriangle(CanvasPoint(x,triangle.v1().y), triangle.v1(), triangle.v2());
	fill_half_triangle(top_triangle, triangle_colour, window);
	fill_half_triangle(bottom_triangle, triangle_colour, window);*/

	fill_triangle(triangle, triangle_colour, window);

	draw_unfilled_triangle(triangle, Colour(255,255,255), window);
}

std::vector<CanvasPoint> interpolate_points(CanvasPoint start, CanvasPoint end, int steps) {
	std::vector<CanvasPoint> result;
	
	for (size_t i = 0; i < steps + 1; i++) {
		CanvasPoint p;
		if (steps == 1) steps = 2;
		p.x = start.x + ((end.x - start.x) * i / (steps));
		p.y = start.y + ((end.y - start.y) * i / (steps));
		result.push_back(p);
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

	split.x = find_mid_value(mid, top, bot);

	float scale = (mid.y - top.y)/(bot.y-top.y);

	split.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	split.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

	CanvasTriangle top_triangle = CanvasTriangle(top, mid, split);
	fill_texture_top(top_triangle, texture, window);

	CanvasTriangle bottom_triangle = CanvasTriangle(mid, split, bot);
	fill_texture_bottom(bottom_triangle, texture, window);

	draw_unfilled_triangle(triangle, Colour(255,255,255), window);
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
		else if (event.key.keysym.sym == SDLK_u	) {
			CanvasTriangle triangle = CanvasTriangle(CanvasPoint(rand() % window.width, rand() % window.height),
																							 CanvasPoint(rand() % window.width, rand() % window.height),
																							 CanvasPoint(rand() % window.width, rand() % window.height));
			Colour triangle_colour = Colour(rand() % 255, rand() % 255, rand() % 255);
			draw_unfilled_triangle(triangle, triangle_colour, window);
		}
		else if (event.key.keysym.sym == SDLK_f	) {
			CanvasTriangle triangle = CanvasTriangle(CanvasPoint(rand() % window.width, rand() % window.height),
																							 CanvasPoint(rand() % window.width, rand() % window.height),
																							 CanvasPoint(rand() % window.width, rand() % window.height));
			Colour triangle_colour = Colour(rand() % 255, rand() % 255, rand() % 255);
			draw_filled_triangle(triangle, triangle_colour, window);
		}
		else if (event.key.keysym.sym == SDLK_t	) {
			TextureMap texture = TextureMap("texture.ppm");

			CanvasPoint p1 = CanvasPoint(160, 10);
      		CanvasPoint p2 = CanvasPoint(300,230);
      		CanvasPoint p3 = CanvasPoint(10, 150);

      		p1.texturePoint = TexturePoint(195, 5);
      		p2.texturePoint = TexturePoint(395, 380);
      		p3.texturePoint = TexturePoint(65, 330);

			CanvasTriangle triangle = CanvasTriangle(p1, p2, p3);

			map_texture(triangle, texture, window);
    }
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
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
