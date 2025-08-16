#include"common.h"
#undef NOTPREV 1
#undef NORANGE 1

std::mutex image_mutex;

bool inrange(uint x, uint y){
    return ((x-5000)*(x-5000) + (y-5000)*(y-5000)) > 6250000;
}

void compute(uint width, uint height, uint num, uint preiter, png::image<png::rgb_pixel> &image, std::vector<std::pair<uint, uint>> points, double weight){
    // seed the machine
    unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> x(0, width - 1);
    std::uniform_int_distribution<int> y(0, height - 1);
    std::uniform_int_distribution<int> index(0, points.size() - 1);
    // get start values
    uint wx = x(gen); uint wy = y(gen);
    uint tx, ty, px, py;
    // pre iterate for quality
    #ifdef NOTPREV
    std::pair<uint, uint> prev={999999999,999999998};
    #endif
    std::pair<uint, uint> chosen_point = {999999999,999999999};
    for(uint i = 0; i < preiter; ++i){
        tx = wx; ty = wy;
        rechoosepre:
        // choose the point
        chosen_point = points[index(gen)];
        #ifdef NOTPREV
        while (prev == chosen_point){
        #endif
            chosen_point = points[index(gen)];
        #ifdef NOTPREV
        }
        prev = chosen_point;
        #endif
        px = chosen_point.first; py = chosen_point.second;
        wx = px * weight + tx * (1.0 - weight);
        wy = py * weight + ty * (1.0 - weight);
        #ifdef NORANGE
        if (!inrange(wx, wy)) goto rechoosepre;
        #endif
    }
    for(uint i = 0; i < num; ++i){
        tx = wx; ty = wy;
        rechoose:
        // choose the point
        chosen_point = points[index(gen)];
        #ifdef NOTPREV
        while (prev == chosen_point){
        #endif
            chosen_point = points[index(gen)];
        #ifdef NOTPREV
        }
        prev = chosen_point;
        #endif
        px = chosen_point.first; py = chosen_point.second;
        wx = px * weight + tx * (1.0 - weight);
        wy = py * weight + ty * (1.0 - weight);
        #ifdef NORANGE
        if (!inrange(wx, wy)) goto rechoose;
        #endif
        std::lock_guard<std::mutex> lock(image_mutex);
        image[wy][wx] = png::rgb_pixel(255, 255, 255);
    }
}

int main(int argc, char ** argv){
    uint width = 6561;
    uint height = 6561;
    uint tries = 10000000;
    uint preiter = 100;
    std::vector<std::pair<uint, uint>> points = {{0,0}, {8192, 0}, {8192, 8192}, {0, 8192}, {4096, 0}, {0, 4096}, {8192, 4096}, {4096, 8192}};
    double weight = 2.0 / 3;
    png::image<png::rgb_pixel> image(width, height);
    // threads work
    uint num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    for (uint i = 0; i < num_threads; ++i) {
        threads.push_back(std::thread(compute, width, height, tries, preiter, std::ref(image), points, weight));
    }
    // Join all threads
    std::for_each(threads.begin(), threads.end(), [](std::thread &t) { t.join(); });
    try {
        image.write("img.png");
        std::cout << "Image 'img.png' created successfully." << std::endl;
    } catch (const png::error& e) {
        std::cerr << "Error writing PNG: " << e.what() << std::endl;
        return 1;
    }
}