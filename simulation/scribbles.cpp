#include <iostream>
#include <numeric>
#include <vector>
#include "Eigen/Dense"
#include "sph_sim.h"
#include "gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communikation with Gnuplot

using namespace Eigen;
using namespace std;

void wait_for_key(); 	// Programm halts until keypress
#define NPOINTS    50 	// length of array

double myfunction (double a, double b) {
	return a + b;
}
int main()
{
	Gnuplot g1("lines");
	g1.showonscreen(); // window output
	//
	// User defined 1d, 2d and 3d point sets
	//
	std::vector<double> x, y, y2;

	for (int i = 0; i < NPOINTS; i++)  // fill double arrays x, y, z
	{
		x.push_back((double)i);             // x[i] = i
		y.push_back((double)i * (double)i); // y[i] = i^2
	}

	g1.reset_plot();
	cout << endl << endl << "*** user-defined lists of points (x,y)" << endl;
	g1.set_grid();
	g1.set_style("points").plot_xy(x,y,"user-defined points 2d");

	wait_for_key();



}
void wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}