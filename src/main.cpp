#include "mainburningship.h"
#include "mainmandelbrot.h"

int
main ()
{
  // mandelbrot (800, 1200, -1.76918800235998352477, 0.00355923323449534867, 10000, 150000000, 4600, 2000);
  // system ("mv img.png img0.png");
  // burningship (800, 800, -0.5, 0, 10000, 1, 200, 100);
  // system ("mv img.png img-1.png");
  // mandelbrot(9000, 16000, -0.75, 0, 10000, 1, 1000, 200);
  mandelbrot(2234, 3456, -0.75, 0, 10000, 1, 1000, 20);
  return 0;
}
