
#ifndef vtkMarchingCubesSpaceCases_h
#define vtkMarchingCubesSpaceCases_h
//
// marching cubes space case table for generating isosurfaces
//
#include "vtkCommonDataModelModule.h" // For export macro
#include "vtkSystemIncludes.h"



typedef int EDGE_LIST;
//struct VTKCOMMONDATAMODEL_EXPORT vtkMarchingCubesSpaceSpaceSpaceTriangleCases  define VTKFILTERSCORE_EXPORT __declspec(dllimport),与dll调用有关
struct  vtkMarchingCubesSpaceTriangleCases
{
  EDGE_LIST edges[16];
  static vtkMarchingCubesSpaceTriangleCases* GetCases();
};


struct  vtkMarchingCubesSpaceNeighbourCases
{
  bool neighbours[6];
  static vtkMarchingCubesSpaceNeighbourCases* GetCases();
};




#endif
// VTK-HeaderTest-Exclude: vtkMarchingCubesSpaceTriangleCases.h
