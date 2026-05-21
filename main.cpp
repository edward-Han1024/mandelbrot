#include"mainmandelbrot.h"
int main(){
    mandelbrot(8000, 8000, -0.5, 0, 10000, 1, 200, 100);
    system("mv img.png img0.png");
}