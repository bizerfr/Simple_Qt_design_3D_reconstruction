#ifndef _GUI_h
#define _GUI_h

#include <QMainWindow>
#include "ui_GUI.h"
#include <vtkSmartPointer.h>


class vtkRenderer;
class vtkEventQtSlotConnect;
class vtkObject;
class vtkCommand;

class vtkVolume;
class vtkActor;

class GUI : public QMainWindow, public Ui::GUI
{
  Q_OBJECT
public:
  GUI();
  ~GUI();

public slots:
  void updateCoords(vtkObject*);
  void popup(vtkObject * obj, unsigned long,
             void * client_data, void *,
             vtkCommand * command);
  void color1(QAction*);
  void color2(QAction*);

  void on_actionInput_image_size_triggered();
  void on_actionInput_the_number_of_image_series_triggered();
  void on_actionInput_data_space_triggered();

  void on_actionRayCasting_triggered();
  void on_actionTextureMapper2D_triggered();
  void on_actionTextureMapper3D_triggered();
  void on_actionMarchingCubes_triggered();
  void on_actionNewMarchingCubes_triggered();


protected:
  vtkRenderer* Ren1;
  vtkRenderer* Ren2;
  vtkEventQtSlotConnect* Connections;

  vtkSmartPointer<vtkVolume> volume;
  vtkSmartPointer<vtkActor> actor;
  int ImageLength;
  int ImageWidth;
  int NumberOfImages;

  int XDataSpace;
  int YDataSpace;
  int ZDataSpace;
};

#endif // _GUI_h

