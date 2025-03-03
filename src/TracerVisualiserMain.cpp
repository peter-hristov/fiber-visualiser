#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <QApplication>
#include <filesystem>
#include "./TracerVisualiserWindow.h"

#include "./ReebSpace.h"

#include "CGALTypedefs.h"


using namespace std;
namespace fs = std::filesystem;


int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "The usage is ./tv9k <input_file>." << endl;
        return 0;
    }

    std::string filename = argv[1];



    fs::path filePath(filename);
    
    if (!fs::exists(filePath)) {
        std::cerr << "Error: File does not exist: " << filename << std::endl;
    }

    // Set up the data class
    Data* data = new Data();

    std::string extension = filePath.extension().string();
    if (extension == ".vtu") 
    {
        data->readDataVTK(filename);
    } else if (extension == ".txt") 
    {
        data->readData(filename);
    } else 
    {
        std::cerr << "Error: Unsupported file type: " << extension << std::endl;
    }












    // Read in data file

    // Compute domain/range data bounds
    data->computeMinMaxRangeDomainCoordinates();

    // Compute the 2D arrangement
    ReebSpace::computeArrangement(data);
    ReebSpace::computeUpperLowerLink(data);

    ReebSpace::BFS(data);

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
