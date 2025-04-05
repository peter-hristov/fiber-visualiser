#include "src/DisjointSet.h"
#include <unordered_map>
#include <set>
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


#include "./utility/CLI11.hpp"


using namespace std;
namespace fs = std::filesystem;


int main(int argc, char* argv[])
{

    //std::unordered_map<set<int>, int, MyHash<set<int>> test;

    //std::unordered_map<set<int>, int, MyHash<set<int>>> testDta;

    //set<int> a = {1,2,3};
    //set<int> b = {4,2,3};

    //testDta[a] = 2;
    //testDta[b] = 3;

    //assert(testDta[a] == testDta[b]);

    //cout << testDta[a] << " " << testDta[b];


    //return 0;

    // Initialize the union-find structure
    //CGAL::Union_find<set<int>> uf;

    //// Create elements
    //CGAL::Union_find<set<int>>::handle a = uf.push_back({1, 2});  // Adds 11 and returns the index
    //CGAL::Union_find<set<int>>::handle b = uf.push_back({2, 3});  // Adds 12 and returns the index

    //cout << "Number of sets before unification: " << uf.number_of_sets() << endl;

    //uf.unify_sets(a, b);  // Use the indices returned by push_back

    //cout << "Number of sets after unification: " << uf.number_of_sets() << endl;

    //return 0;


    // Parse the command line arguments
    CLI::App cliApp("Fiber Visualiser");

    string filename = "a";
    cliApp.add_option("--file, -f", filename, "Input data filename. Has to be either .txt of .vti.")->required();

    bool performanceRun = false;
    cliApp.add_flag("--performanceRun, -p", performanceRun, "Only compute the Reeb space, no graphics..");

    bool discardFiberSeeds = false;
    cliApp.add_flag("--discardPreimageGraphs, -d", discardFiberSeeds, "Only compute the Reeb space, no graphics..");

    CLI11_PARSE(cliApp, argc, argv);

    // For convenience
    if (performanceRun == true)
    {
        discardFiberSeeds = true;
    }




    fs::path filePath(filename);
    
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

    //std::cout << "Press Enter to continue...";
    //std::cin.get();  // waits for Enter key

    Timer::start();
    ReebSpace::computeUpperLowerLink(data);
    Timer::stop("Computed upper and lower link          :");

    Timer::start();
    ReebSpace::computeTriangleAdjacency(data);
    Timer::stop("Computed triangle adjacency            :");


    //cout << "Triangles to Index " << endl;
    //for (const auto &[triangle, triangleId] : data->triangleToIndex)
    //{
        //cout << "Triangle = ";
        //for (const auto v : triangle)
        //{
            //cout << v << " ";
        //}
        //cout << "  ID = " << triangleId << endl;
    //}

    //cout << "\n\nIndex to triangle" << endl;
    //for (int i = 0 ; i < data->indexToTriangle.size() ; i++)
    //{
        //cout << "  ID = " << i << " ";

        //cout << "Triangle = ";
        //for (const auto v : data->indexToTriangle[i])
        //{
            //cout << v << " ";
        //}

        //cout << endl;
    //}


    Timer::start();
    ReebSpace::testTraverseArrangement(data);
    Timer::stop("Computed empty traversal               :");

    Timer::start();
    ReebSpace::computePreimageGraphs(data, discardFiberSeeds);
    Timer::stop("Computed {G_F} and H                   :");

    //Timer::start();
    //ReebSpace::computeCorrespondenceGraph(data);
    //Timer::stop("Computed H                             :");

    Timer::start();
    ReebSpace::computeReebSpacePostprocess(data);
    Timer::stop("Computed RS(f) Postprocess             :");
    //Timer::stop("Computed RS(f)                         :");

    //std::cout << "Press Enter to continue...";
    //std::cin.get();  // waits for Enter key

    if (performanceRun == true)
    {
        return 0;
    }

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
