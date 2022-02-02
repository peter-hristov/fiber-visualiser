#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <QApplication>

#include "./TracerVisualiserWindow.h"
#include "./external/CLI11.hpp"
#include "./GlobalConfig.h"

using namespace std;


int main(int argc, char* argv[])
{
    // Read in the input file
    tv9k::InputInformation input;
    input.filename = argv[1];

    // Read in data file
    Data* data = new Data();
    data->readNcData(input);


    // Parse the horizontal and vertical lines
    input.explodeLinesToVector();

    // Set up QT Application
    QApplication app(argc, argv);
    glutInit(&argc, argv);

    // Create the widget
    TracerVisualiserWindow* window = new TracerVisualiserWindow(NULL, data, input);
    window->setWindowTitle("Tracer Visualiser 9000");

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
