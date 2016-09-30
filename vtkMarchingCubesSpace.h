

#ifndef vtkMarchingCubesSpace_h
#define vtkMarchingCubesSpace_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

#include "vtkContourValues.h" // Needed for direct access to ContourValues

class vtkIncrementalPointLocator;

//class VTKFILTERSCORE_EXPORT vtkMarchingCubesSpace : public vtkPolyDataAlgorithm  //#      define VTKFILTERSCORE_EXPORT __declspec(dllimport),与dll调用有关
class  vtkMarchingCubesSpace : public vtkPolyDataAlgorithm
{
public:
  static vtkMarchingCubesSpace *New();
  vtkTypeMacro(vtkMarchingCubesSpace,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Methods to set contour values
  void SetValue(int i, double value);
  double GetValue(int i);
  double *GetValues();
  void GetValues(double *contourValues);
  void SetNumberOfContours(int number);
  int GetNumberOfContours();
  void GenerateValues(int numContours, double range[2]);
  void GenerateValues(int numContours, double rangeStart, double rangeEnd);

  // Because we delegate to vtkContourValues
  unsigned long int GetMTime();

  // Description:
  // Set/Get the computation of normals. Normal computation is fairly
  // expensive in both time and storage. If the output data will be
  // processed by filters that modify topology or geometry, it may be
  // wise to turn Normals and Gradients off.
  vtkSetMacro(ComputeNormals,int);
  vtkGetMacro(ComputeNormals,int);
  vtkBooleanMacro(ComputeNormals,int);

  // Description:
  // Set/Get the computation of gradients. Gradient computation is
  // fairly expensive in both time and storage. Note that if
  // ComputeNormals is on, gradients will have to be calculated, but
  // will not be stored in the output dataset.  If the output data
  // will be processed by filters that modify topology or geometry, it
  // may be wise to turn Normals and Gradients off.
  vtkSetMacro(ComputeGradients,int);
  vtkGetMacro(ComputeGradients,int);
  vtkBooleanMacro(ComputeGradients,int);

  // Description:
  // Set/Get the computation of scalars.
  vtkSetMacro(ComputeScalars,int);
  vtkGetMacro(ComputeScalars,int);
  vtkBooleanMacro(ComputeScalars,int);

  // Description:
  // Overide the default locator.  Useful for changing the number of
  // bins for performance or specifying a more aggressive locator.
  void SetLocator(vtkIncrementalPointLocator *locator);
  vtkGetObjectMacro(Locator,vtkIncrementalPointLocator);

  // Description:
  // Create default locator. Used to create one when none is
  // specified. The locator is used to merge coincident points.
  void CreateDefaultLocator();

protected:
  vtkMarchingCubesSpace();
  ~vtkMarchingCubesSpace();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int port, vtkInformation *info);

  vtkContourValues *ContourValues;
  int ComputeNormals;
  int ComputeGradients;
  int ComputeScalars;
  vtkIncrementalPointLocator *Locator;
private:
  vtkMarchingCubesSpace(const vtkMarchingCubesSpace&);  // Not implemented.
  void operator=(const vtkMarchingCubesSpace&);  // Not implemented.
};

// Description:
// Set a particular contour value at contour number i. The index i ranges
// between 0<=i<NumberOfContours.
inline void vtkMarchingCubesSpace::SetValue(int i, double value)
{this->ContourValues->SetValue(i,value);}

// Description:
// Get the ith contour value.
inline double vtkMarchingCubesSpace::GetValue(int i)
{return this->ContourValues->GetValue(i);}

// Description:
// Get a pointer to an array of contour values. There will be
// GetNumberOfContours() values in the list.
inline double *vtkMarchingCubesSpace::GetValues()
{return this->ContourValues->GetValues();}

// Description:
// Fill a supplied list with contour values. There will be
// GetNumberOfContours() values in the list. Make sure you allocate
// enough memory to hold the list.
inline void vtkMarchingCubesSpace::GetValues(double *contourValues)
{this->ContourValues->GetValues(contourValues);}

// Description:
// Set the number of contours to place into the list. You only really
// need to use this method to reduce list size. The method SetValue()
// will automatically increase list size as needed.
inline void vtkMarchingCubesSpace::SetNumberOfContours(int number)
{this->ContourValues->SetNumberOfContours(number);}

// Description:
// Get the number of contours in the list of contour values.
inline int vtkMarchingCubesSpace::GetNumberOfContours()
{return this->ContourValues->GetNumberOfContours();}

// Description:
// Generate numContours equally spaced contour values between specified
// range. Contour values will include min/max range values.
inline void vtkMarchingCubesSpace::GenerateValues(int numContours, double range[2])
{this->ContourValues->GenerateValues(numContours, range);}

// Description:
// Generate numContours equally spaced contour values between specified
// range. Contour values will include min/max range values.
inline void vtkMarchingCubesSpace::GenerateValues(int numContours, double
                                             rangeStart, double rangeEnd)
{this->ContourValues->GenerateValues(numContours, rangeStart, rangeEnd);}



struct Cube{
	int oi;
	int oj;
	int ok;
	int index;
	bool Flag;
	Cube(int i=0,int j=0,int k=0,int id=0,bool f=0):oi(i),oj(j),ok(k),index(id),Flag(f){}
};



#endif
