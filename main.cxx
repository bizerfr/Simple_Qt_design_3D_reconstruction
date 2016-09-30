#include "QVTKApplication.h"

#include "GUI.h"

int main(int argc, char** argv)
{
  QVTKApplication app(argc, argv);
  GUI widget;

  widget.show();

  return app.exec();
}
