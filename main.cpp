#include"mainmandelbrot.h"
#include"mainburningship.h"

int main(){
    mandelbrot(800, 800, -0.5, 0, 10000, 1, 200, 100);
    system("mv img.png img0.png");
    burningship(800, 800, -0.5, 0, 10000, 1, 200, 100);
}
