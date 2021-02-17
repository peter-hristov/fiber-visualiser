#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <QApplication>
#include <QVBoxLayout>
#include <memory>

#include "./ScatterPlot.h"
#include "./TracerVisualiserWindow.h"
#include "./external/CLI11.hpp"
#include "./utility/MergeTree.h"
#include "./utility/Utility.h"

#include "./GlobalConfig.h"

using namespace std;

const string
checkFunction(string itype)
{
    if (itype == "none" || itype == "bilinear" || itype == "nearest") {
        return "";
    } else {
        return "The interpolation type must be one of none, bilinear or nearest.";
    }
}

int
main(int argc, char* argv[])
{
    tv9k::InputInformation input;

    CLI::App cliApp("K3Pi goofit fitter");
    cliApp.add_option("-f,--file,file", input.filename, "File name")->required();
    cliApp.add_option("-i,--isovalue-property", input.attributeNames[0], "Isovalue Property Name")->required();
    cliApp
      .add_option("-u,--fiber-surface-property-u", input.attributeNames[1], "Fiber Surface Property U (horizontal)")
      ->required();
    cliApp.add_option("-v,--fiber-surface-property-v", input.attributeNames[2], "Fiber Surface Property V (vertical)")
      ->required();
    cliApp.add_option("-x,--x-dimension-name", input.xName, "X Dimension Name");
    cliApp.add_option("-y,--y-dimension-name", input.yName, "Y Dimension Name");
    cliApp.add_option("-z,--z-dimension-name", input.zName, "Z Dimension Name");
    cliApp.add_option("-t,--t-dimension-name", input.tName, "T Dimension Name");
    cliApp.add_option("--xmin", input.xMin, "The lower bound to clip in the X Dimension.");
    cliApp.add_option("--xmax", input.xMax, "The upper bound to clip in the X Dimension.");
    cliApp.add_option("--ymin", input.yMin, "The lower bound to clip in the Y Dimension.");
    cliApp.add_option("--ymax", input.yMax, "The upper bound to clip in the Y Dimension.");
    cliApp.add_option("--zmin", input.zMin, "The lower bound to clip in the Z Dimension.");
    cliApp.add_option("--zmax", input.zMax, "The upper bound to clip in the Z Dimension.");
    cliApp.add_option("--tmin", input.tMin, "The lower bound to clip in the T Dimension.");
    cliApp.add_option("--tmax", input.tMax, "The upper bound to clip in the T Dimension.");
    cliApp.add_option("--downsampleX",
                      input.downsampleX,
                      "Each dimension is downsampled by this factor. Thus data "
                      "size effectively downsampled d^3, where d is he "
                      "downsample value. Default is 1, or no downsampling.");
    cliApp.add_option("--downsampleY",
                      input.downsampleY,
                      "Each dimension is downsampled by this factor. Thus data "
                      "size effectively downsampled d^3, where d is he "
                      "downsample value. Default is 1, or no downsampling.");
    cliApp.add_option("--downsampleZ",
                      input.downsampleZ,
                      "Each dimension is downsampled by this factor. Thus data "
                      "size effectively downsampled d^3, where d is he "
                      "downsample value. Default is 1, or no downsampling.");

    cliApp.add_option("--verticalLines",
                      input.verticalLines,
                      "Vertical lines for the U Field at which to draw a constant line in the scatterplot");
    cliApp.add_option("--horizontalLines",
                      input.horizontalLines,
                      "Horizontal lines for the V Field at which to draw a constant line in the scatterplot");

    cliApp.add_option("--resolution", input.scatterplotResolution, "The resolution of the scatterplot.");

    string interpolationType = "none";
    cliApp
      .add_option("--interpolation-type",
                  interpolationType,
                  "Do not set unless you have performance issues! Type of interpolation to "
                  "use for computing the Fiber Surface distance field. Available options "
                  "are none, nearest, and bilinear. Default is none, use bilinear and none "
                  "for faster fiber surface computation with decreased accuracy and more "
                  "artifacts. Nearest is faster than bilinear.")
      ->check(checkFunction);

    Data* data = new Data();

    cliApp.add_flag("--cacheJoinTree", data->cacheJoinTree,
                  "If provided jointree will be cached (and read from cache) for the iso field "
                  "(and all other fields if combined with the `-m` flag). The join tree data "
                  "will be written to `<input_file>.tree`.");

    cliApp.add_flag(
      "--dynamicPolygon", data->dynamicPolygon, "Dynamically render fiber surface polygon as you move points.");
    cliApp.add_flag("--render3DLabels, -r", data->render3DLabels, "Render 3D labels with distance in the plot.");
    cliApp.add_flag("--continuousScatterPlot, -c", data->continuousScatterPlot, "Render Continuous Scatterplot.");
    cliApp.add_flag("--flatNormals, -n", data->flatNormals, "Use flat normals for debugging.");

    cliApp.add_flag("--frontloadComputation, -m",
                    data->precomputeMergeTrees,
                    "Frontload the merge computation so it's smooth when you "
                    "switch isofields.");

    CLI11_PARSE(cliApp, argc, argv);

    if (input.filename.find(".vti") != std::string::npos) {
        try {

#ifdef VTK_DIR
            // data->readVtkData(filename, attributeNames[0], attributeNames[1],
            // attributeNames[2], xName, yName, zName, tName, xMin, xMax, yMin, yMax,
            // zMin, zMax, tMin, tMax, downsampleX, downsampleY, downsampleZ);
#endif
        } catch (string message) {
            cerr << message;
            return 1;
        }
    } else if (input.filename.find(".nc") != std::string::npos) {
        data->readNcData(input);
    } else {
        printf("File must be be either vtk or nc. Not other format are supported.");
    }

    input.explodeLinesToVector();

    // data->computeMinMax();

    // Compute visited array and merge tree

    // Utility::startTimer();
    // printf("\n\n------------------------------------ Computing Continuous "
    //"Scatterplot. \n");
    // data->plot->computeDensityDiscrete(*data->uField, *data->vField, 0);
    // data->plot->populateImage(0);

    // printf("------------------------------------ Done - %3ld.%06ld seconds.\n",
    // Utility::endTimer().first,
    // Utility::endTimer().second);

    QApplication app(argc, argv);
    glutInit(&argc, argv);

    // create model (polygon) as a triangle
    //	GLPolygon *polygon = new GLPolygon();

    // create a master widget
    TracerVisualiserWindow* window = new TracerVisualiserWindow(NULL, data, input, interpolationType);
    window->setWindowTitle("Tracer Visualiser 9000");

    // create a controller to hook things up
    //	GLPolygonController *controller = new GLPolygonController(window,
    // polygon);

    // make the window full screen by default
    window->showMaximized();

    // show the label
    window->show();

    // start it running
    app.exec();

    // clean up
    //	delete controller;
    delete window;
    delete data;

    // return to caller
    return 0;
} // main()
