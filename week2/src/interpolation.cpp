#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#define WIDTH 320
#define HEIGHT 240

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	std::vector<float> result;
	float step = ((to - from) / (numberOfValues - 1));
	for (size_t i = 0; i < numberOfValues; i++) {
		result.push_back(from + (step * i));
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> result;
	float first_step = ((to[0] - from[0]) / (numberOfValues - 1));
	float second_step = ((to[1] - from[1]) / (numberOfValues - 1));
	float third_step = ((to[2] - from[2]) / (numberOfValues - 1));
	for (size_t i = 0; i < numberOfValues; i++) {
		glm::vec3 temp(from[0] + (first_step * i), from[1] + (second_step * i), from[2] + (third_step * i));
		result.push_back(temp);
	}
	return result;
}

void draw(DrawingWindow &window) {
	/*window.clearPixels();
	std::vector<float> interpolation = interpolateSingleFloats(255.0, 0.0, window.width);
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = interpolation[x];
			float green = interpolation[x];
			float blue = interpolation[x];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}*/

	window.clearPixels();
	glm::vec3 top_left(255, 0, 0); //red
	glm::vec3 top_right(0, 0, 255); //blue
	glm::vec3 bottom_right(0, 255, 0); //green
	glm::vec3 bottom_left(255, 255, 0); //yellow

	std::vector<glm::vec3> first_column = interpolateThreeElementValues(top_left, bottom_left, window.height);
	std::vector<glm::vec3> last_column = interpolateThreeElementValues(top_right, bottom_right, window.height);

	for (size_t y = 0; y < window.height; y++) {
		std::vector<glm::vec3> row = interpolateThreeElementValues(first_column[y], last_column[y], window.width);
		for (size_t x = 0; x < window.width; x++) {
			float red = row[x].r;
			float green = row[x].g;
			float blue = row[x].b;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
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

	std::vector<glm::vec3> result;
	result = interpolateThreeElementValues(glm::vec3(1.0, 4.0, 9.2), glm::vec3(4.0, 1.0, 9.8), 4);
	for (size_t i = 0; i < result.size(); i++) std::cout << glm::to_string(result[i]) << " ";
	std::cout << std::endl;
}
