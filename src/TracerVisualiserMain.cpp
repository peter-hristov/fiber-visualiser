#include <cstddef>
#include <filesystem>

#include <GL/glut.h>
#include <QApplication>

#include "./io.h"
#include "./Timer.h"
#include "./ReebSpace.h"
#include "./CGALTypedefs.h"
#include "./Arrangement.h"
#include "./utility/CLI11.hpp"
#include "./TracerVisualiserWindow.h"

using namespace std;

int main(int argc, char* argv[])
{
    CLI::App cliApp("Reeb Space Fiber Visualiser");

    string filename;
    cliApp.add_option("--file, -f", filename, "Input data filename. Has to be either .txt of .vti.")->required();

    bool performanceRun = false;
    cliApp.add_flag("--performanceRun, -p", performanceRun, "Only compute the Reeb space, no graphics..");

    bool discardFiberSeeds = false;
    cliApp.add_flag("--discardPreimageGraphs, -d", discardFiberSeeds, "Discard the seeds for generating fibers based on sheets, discard to save a bit of memory (not too much).");

    float perturbationEpsilon = 1e-2;
    cliApp.add_option("--epsilon, -e", perturbationEpsilon, "Strength of the numerial perturbation in the range [-e, e].");

    string outputSheetPolygonsFilename;
    cliApp.add_option("--outputSheetPolygons, -o", outputSheetPolygonsFilename, "Filename where to output the coordinates of the polygons that represent each sheet.");

    int fiberSampling = 1;
    cliApp.add_option("--fiberSampling, -s", fiberSampling, "When saving fibers per component, how many do we save. Default is to save the centroid, otherwise sample along the boundary.");

    int sheetOutputCount = 10;
    cliApp.add_option("--sheetOutputCount", sheetOutputCount, "How many sheets to sample for automatic feature extraction.");

    string outputSheetFibersFolder;
    cliApp.add_option("--outputSheetFibersFolder", outputSheetFibersFolder, "Folder in which to ouput fiber for each sheet.");

    //string outputFibersFilename = "./fibers.vtp";
    //cliApp.add_option("--outputFibers", outputSheetPolygonsFilename, "Filename where to save the visible fiber components. Must be .vtp");

    CLI11_PARSE(cliApp, argc, argv);

    // For convenience
    if (performanceRun == true)
    {
        discardFiberSeeds = true;
    }

    // Read, perturb and sort the indices of the vertices lexicographically (by their range position).
    TetMesh tetMesh;
    try
    {
        tetMesh = io::readData(filename);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    Timer::start();
    tetMesh.perturbRangeValues(perturbationEpsilon);
    tetMesh.sortVertices();
    tetMesh.computeBoundingBoxes();
    tetMesh.computeCombinatorialStructure();
    tetMesh.computeUpperLowerLinkAndStar();
    tetMesh.computeSingularEdgeTypes();
    Timer::stop("Input mesh postprocessing              :");

    Timer::start();
    Arrangement arrangement;
    arrangement.computeArrangement(tetMesh);
    Timer::stop("Arrangement                            :");

    Timer::start();
    arrangement.computePointLocationDataStructure();
    Timer::stop("Arrangement search structure           :");

    Timer::start();
    ReebSpace reebSpace;
    reebSpace.computeTraversal(tetMesh, arrangement, discardFiberSeeds);
    Timer::stop("Computed {G_F} and H                   :");

    std::cout << "Postprocessing..." << std::endl;
    Timer::start();
    reebSpace.computeSheetGeometry(tetMesh, arrangement);
    reebSpace.computeSheetArea(tetMesh, arrangement);
    reebSpace.printTopSheets(tetMesh, arrangement, 20);
    Timer::stop("Computed RS(f) Postprocess             :");


    // Save all the polygons
    if (false == outputSheetPolygonsFilename.empty())
    {
        try
        {
            io::saveSheets(tetMesh, arrangement, reebSpace, outputSheetPolygonsFilename);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        }
    }

    if (performanceRun == true)
    {
        return 0;
    }


    // Package all my data for visualisation
    Data data(tetMesh, arrangement, reebSpace);

    if (false == outputSheetFibersFolder.empty())
    {
        data.generatefFaceFibersForSheets(sheetOutputCount, fiberSampling, outputSheetFibersFolder);
    }

    //ReebSpace::countIntersectionsTypes(data);


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

    // return to caller
    return 0;
} // main()
