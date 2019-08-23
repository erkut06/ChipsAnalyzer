#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLWidget.h>

#include "ChipsAnalyzer.h"
#include "vtkFileOutputWindow.h"

int main(int argc, char** argv)
{

	vtkFileOutputWindow* w = vtkFileOutputWindow::New();
	std::string str("");
	str += "vtkOutput.txt";
	w->SetFileName(str.c_str());
	vtkOutputWindow::SetInstance(w);
	// needed to ensure appropriate OpenGL context is created for VTK rendering.
	QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());

	// QT Stuff
	QApplication app(argc, argv);

	ChipsAnalyzer chipsAnalyzer;
	chipsAnalyzer.show();
	w->Delete();
	return app.exec();

}
