#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <QApplication>
#include <filesystem>
#include "./TracerVisualiserWindow.h"
#include "./Timer.h"
#include "./ReebSpace.h"

#include "CGALTypedefs.h"

#include <CGAL/Union_find.h>

using namespace std;
namespace fs = std::filesystem;


int main(int argc, char* argv[])
{

// Initialize the union-find structure
    //CGAL::Union_find<set<int>> uf;

    //// Create elements
    //CGAL::Union_find<set<int>>::handle a = uf.push_back({1, 2});  // Adds 11 and returns the index
    //CGAL::Union_find<set<int>>::handle b = uf.push_back({2, 3});  // Adds 12 and returns the index

    //cout << "Number of sets before unification: " << uf.number_of_sets() << endl;

    //uf.unify_sets(a, b);  // Use the indices returned by push_back

    //cout << "Number of sets after unification: " << uf.number_of_sets() << endl;

    //return 0;

    if (argc != 2)
    {
        cout << "The usage is ./tv9k <input_file>." << endl;
        return 0;
    }

    std::string filename = argv[1];

    fs::path filePath(argv[1]);
    
    if (!fs::exists(filePath)) 
    {
        std::cerr << "Error: File does not exist: " << filename << std::endl;
        return 0;
    }

    // Set up the data class
    Data* data = new Data();

    std::string extension = filePath.extension().string();
    if (extension == ".vtu") 
    {
        data->readDataVTK(filename);
    } 
    else if (extension == ".txt") 
    {
        data->readData(filename);
    } 
    else 
    {
        std::cerr << "Error: Unsupported file type: " << extension << std::endl;
    }

    // Compute the 2D arrangement
    ReebSpace::computeArrangement(data);


    // Compute the upper and lower links of each edge
    Timer::start();
    ReebSpace::computeUpperLowerLink(data);
    Timer::stop("Computed upper and lower link          :");


    // Compute the adjacency of triangles in the mesh
    Timer::start();
    ReebSpace::computeTriangleAdjacency(data);
    Timer::stop("Computed triangle adjacency            :");

    Timer::start();
    ReebSpace::testTraverseArrangement(data);
    Timer::stop("Computed empty traversal               :");

    // Compute the preimageGraphs of each face in the arrangement
    Timer::start();
    ReebSpace::computePreimageGraphs(data);
    Timer::stop("Computed preimage graph                :");


    // Compute the Reeb space (connected components of preimage graph components)
    Timer::start();
    ReebSpace::computeReebSpace(data);
    Timer::stop("Computed H and RS                      :");

    //data->pl = Point_location(data->arr);
    Timer::start();
    data->pl = std::make_unique<Point_location>(data->arr);
    Timer::stop("Arrangement search structure           :");

    // Set up QT Application
    QApplication app(argc, argv);
    glutInit(&argc, argv);

    // Create the widget
    TracerVisualiserWindow* window = new TracerVisualiserWindow(NULL, data);
    window->setWindowTitle("Fiber Visualiser");

    // Make the window full screen by default
    window->showMaximized();

    // Show the label
    window->show();

    // start it running
    app.exec();

    // clean up
    delete window;
    delete data;

    // return to caller
    return 0;
} // main()
