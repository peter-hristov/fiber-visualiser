#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <QApplication>
#include "./TracerVisualiserWindow.h"

#include "./ReebSpace.h"

#include "CGALTypedefs.h"


using namespace std;


int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "The usage is ./tv9k <input_file>." << endl;
        return 0;
    }

    std::string filename = argv[1];

    // Set up the data class
    Data* data = new Data();

    // Read in data file
    data->readData(filename);

    // Compute domain/range data bounds
    data->computeMinMaxRangeDomainCoordinates();

    // Compute the 2D arrangement
    data->arr = ReebSpace::computeArrangement(data);

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
