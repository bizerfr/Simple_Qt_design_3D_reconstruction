
#include "vtkMarchingCubesSpace.h"
#include "vtkCellArray.h"
#include "vtkCharArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkMarchingCubesSpaceCases.h"
#include "vtkMath.h"
#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkShortArray.h"
#include "vtkStructuredPoints.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkIncrementalPointLocator.h"

#include<queue>
#include<vector>

vtkStandardNewMacro(vtkMarchingCubesSpace);

// Description:
// Construct object with initial range (0,1) and single contour value
// of 0.0. ComputeNormal is on, ComputeGradients is off and ComputeScalars is on.
vtkMarchingCubesSpace::vtkMarchingCubesSpace()
{
  this->ContourValues = vtkContourValues::New();
  this->ComputeNormals = 1;
  this->ComputeGradients = 0;
  this->ComputeScalars = 1;
  this->Locator = NULL;
  this->SetNumberOfInputPorts(2);
}

vtkMarchingCubesSpace::~vtkMarchingCubesSpace()
{
  this->ContourValues->Delete();
  if ( this->Locator )
    {
    this->Locator->UnRegister(this);
    this->Locator = NULL;
    }
}

// Description:
// Overload standard modified time function. If contour values are modified,
// then this object is modified as well.
unsigned long vtkMarchingCubesSpace::GetMTime()
{
  unsigned long mTime=this->Superclass::GetMTime();
  unsigned long mTime2=this->ContourValues->GetMTime();

  mTime = ( mTime2 > mTime ? mTime2 : mTime );
  if (this->Locator)
    {
    mTime2=this->Locator->GetMTime();
    mTime = ( mTime2 > mTime ? mTime2 : mTime );
    }

  return mTime;
}

// Calculate the gradient using central difference.
// NOTE: We calculate the negative of the gradient for efficiency
template <class T>
void vtkMarchingCubesSpaceComputePointGradient(int i, int j, int k, T *s, int dims[3],
                                          vtkIdType sliceSize, double spacing[3], double n[3])
{
  double sp, sm;

  // x-direction
  if ( i == 0 )
    {
    sp = s[i+1 + j*dims[0] + k*sliceSize];    // dims[0] stores xmax,dims[1] stores ymax,dims[2] stores zmax
    sm = s[i + j*dims[0] + k*sliceSize];
    n[0] = (sm - sp) / spacing[0];  //spacing[0] stores delta x;spacing[2] stores delta y,spacing[3] stores delta z
    }
  else if ( i == (dims[0]-1) )
    {
    sp = s[i + j*dims[0] + k*sliceSize];
    sm = s[i-1 + j*dims[0] + k*sliceSize];
    n[0] = (sm - sp) / spacing[0];
    }
  else
    {
    sp = s[i+1 + j*dims[0] + k*sliceSize];
    sm = s[i-1 + j*dims[0] + k*sliceSize];
    n[0] = 0.5 * (sm - sp) / spacing[0];
    }

  // y-direction
  if ( j == 0 )
    {
    sp = s[i + (j+1)*dims[0] + k*sliceSize];
    sm = s[i + j*dims[0] + k*sliceSize];
    n[1] = (sm - sp) / spacing[1];
    }
  else if ( j == (dims[1]-1) )
    {
    sp = s[i + j*dims[0] + k*sliceSize];
    sm = s[i + (j-1)*dims[0] + k*sliceSize];
    n[1] = (sm - sp) / spacing[1];
    }
  else
    {
    sp = s[i + (j+1)*dims[0] + k*sliceSize];
    sm = s[i + (j-1)*dims[0] + k*sliceSize];
    n[1] = 0.5 * (sm - sp) / spacing[1];
    }

  // z-direction
  if ( k == 0 )
    {
    sp = s[i + j*dims[0] + (k+1)*sliceSize];
    sm = s[i + j*dims[0] + k*sliceSize];
    n[2] = (sm - sp) / spacing[2];
    }
  else if ( k == (dims[2]-1) )
    {
    sp = s[i + j*dims[0] + k*sliceSize];
    sm = s[i + j*dims[0] + (k-1)*sliceSize];
    n[2] = (sm - sp) / spacing[2];
    }
  else
    {
    sp = s[i + j*dims[0] + (k+1)*sliceSize];
    sm = s[i + j*dims[0] + (k-1)*sliceSize];
    n[2] = 0.5 * (sm - sp) / spacing[2];
    }
}






/*//
// Contouring filter specialized for volumes and "short int" data values.
//
template <class T>
void vtkMarchingCubesSpaceComputeGradient(vtkMarchingCubesSpace *self,T *scalars, int dims[3],
                                     double origin[3], double spacing[3],
                                     vtkIncrementalPointLocator *locator,
                                     vtkDataArray *newScalars,
                                     vtkDataArray *newGradients,
                                     vtkDataArray *newNormals,
                                     vtkCellArray *newPolys, double *values,
                                     int numValues)
{
  double s[8], value;
  int i, j, k;
  vtkIdType sliceSize;
  static int CASE_MASK[8] = {1,2,4,8,16,32,64,128};  //即0000 0001,000 0010等等，与index或运算，index即cube八个顶点状态
  vtkMarchingCubesSpaceTriangleCases *triCase, *triCases;  //结构体vtkMarchingCubesSpaceTriangleCases含有edges[16](cube边所对应的triangle)和指向结构体指针函数
  EDGE_LIST  *edge;
  int contNum, jOffset, ii, index, *vert;
  vtkIdType kOffset, idx;
  vtkIdType ptIds[3];
  int ComputeNormals = newNormals != NULL;
  int ComputeGradients = newGradients != NULL;
  int ComputeScalars = newScalars != NULL;
  int NeedGradients;
  int extent[6];
  double t, *x1, *x2, x[3], *n1, *n2, n[3], min, max;
  double pts[8][3], gradients[8][3], xp, yp, zp;      //pts[8][3]中，8是cube八个顶点，3是对应x,y,z坐标
  static int edges[12][2] = { {0,1}, {1,2}, {3,2}, {0,3},
                              {4,5}, {5,6}, {7,6}, {4,7},
                              {0,4}, {1,5}, {3,7}, {2,6}};  // 0-7是cube顶点编号，{0,1}表示由0号顶点和1号顶点之间的棱

  vtkInformation *inInfo = self->GetExecutive()->GetInputInformation(0, 0); //Get this algorithm's executive,Get the pipeline information for the given input connection.
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent);

  triCases =  vtkMarchingCubesSpaceTriangleCases::GetCases();  //得到vtkMarchingCubesSpaceTriangleCases.cxx中静态成员“构型—三角剖分”查找表的头指针

//
// Get min/max contour values
//
  if ( numValues < 1 )
    {
    return;
    }
  for ( min=max=values[0], i=1; i < numValues; i++)
    {
    if ( values[i] < min )
      {
      min = values[i];
      }
    if ( values[i] > max )
      {
      max = values[i];
      }
    }
//
// Traverse all voxel cells, generating triangles and point gradients
// using marching cubes algorithm.
//
  sliceSize = dims[0] * dims[1];
  for ( k=0; k < (dims[2]-1); k++)
    {
    self->UpdateProgress (k / static_cast<double>(dims[2] - 1));
    if (self->GetAbortExecute())
      {
      break;
      }
    kOffset = k*sliceSize;
    pts[0][2] = origin[2] + (k+extent[4]) * spacing[2];
    zp = pts[0][2] + spacing[2];
    for ( j=0; j < (dims[1]-1); j++)
      {
      jOffset = j*dims[0];
      pts[0][1] = origin[1] + (j+extent[2]) * spacing[1];
      yp = pts[0][1] + spacing[1];
      for ( i=0; i < (dims[0]-1); i++)
        {
        //get scalar values
        idx = i + jOffset + kOffset;
        s[0] = scalars[idx];                          //point(i,j,k)对应标量
        s[1] = scalars[idx+1];                       //point(i+1,j,k）对应标量
        s[2] = scalars[idx+1 + dims[0]];             //point(i+1,j+1,k)对应标量
        s[3] = scalars[idx + dims[0]];               //point(i,j+1,k)对应标量
        s[4] = scalars[idx + sliceSize];             //point(i,j,k+1)对应标量
        s[5] = scalars[idx+1 + sliceSize];           //point(i+1,j,k+1)对应标量
        s[6] = scalars[idx+1 + dims[0] + sliceSize]; //point(i+1,j+1,k+1)对应标量
        s[7] = scalars[idx + dims[0] + sliceSize];   //point(i,j+1,k+1)对应标量

        if ( (s[0] < min && s[1] < min && s[2] < min && s[3] < min &&
        s[4] < min && s[5] < min && s[6] < min && s[7] < min) ||
        (s[0] > max && s[1] > max && s[2] > max && s[3] > max &&
        s[4] > max && s[5] > max && s[6] > max && s[7] > max) )
          {
          continue; // no contours possible
          }

        //create voxel points
        pts[0][0] = origin[0] + (i+extent[0]) * spacing[0];
        xp = pts[0][0] + spacing[0];

        pts[1][0] = xp;
        pts[1][1] = pts[0][1];
        pts[1][2] = pts[0][2];

        pts[2][0] = xp;
        pts[2][1] = yp;
        pts[2][2] = pts[0][2];

        pts[3][0] = pts[0][0];
        pts[3][1] = yp;
        pts[3][2] = pts[0][2];

        pts[4][0] = pts[0][0];
        pts[4][1] = pts[0][1];
        pts[4][2] = zp;

        pts[5][0] = xp;
        pts[5][1] = pts[0][1];
        pts[5][2] = zp;

        pts[6][0] = xp;
        pts[6][1] = yp;
        pts[6][2] = zp;

        pts[7][0] = pts[0][0];
        pts[7][1] = yp;
        pts[7][2] = zp;           //cube八个顶点坐标

        NeedGradients = ComputeGradients || ComputeNormals;

        //create gradients if needed
        if (NeedGradients)
          {
          vtkMarchingCubesSpaceComputePointGradient(i,j,k, scalars, dims, sliceSize, spacing, gradients[0]);  //前面定义了double gradients[8][3]
          vtkMarchingCubesSpaceComputePointGradient(i+1,j,k, scalars, dims, sliceSize, spacing, gradients[1]);
          vtkMarchingCubesSpaceComputePointGradient(i+1,j+1,k, scalars, dims, sliceSize, spacing, gradients[2]);
          vtkMarchingCubesSpaceComputePointGradient(i,j+1,k, scalars, dims, sliceSize, spacing, gradients[3]);
          vtkMarchingCubesSpaceComputePointGradient(i,j,k+1, scalars, dims, sliceSize, spacing, gradients[4]);
          vtkMarchingCubesSpaceComputePointGradient(i+1,j,k+1, scalars, dims, sliceSize, spacing, gradients[5]);
          vtkMarchingCubesSpaceComputePointGradient(i+1,j+1,k+1, scalars, dims, sliceSize, spacing, gradients[6]);
          vtkMarchingCubesSpaceComputePointGradient(i,j+1,k+1, scalars, dims, sliceSize, spacing, gradients[7]);
          }
        for (contNum=0; contNum < numValues; contNum++)
          {
          value = values[contNum];
          // Build the case table
          for ( ii=0, index = 0; ii < 8; ii++)
            {
            if ( s[ii] >= value )
              {
              index |= CASE_MASK[ii];  //index是cube的索引号（0-255，表示vertex八个状态）
              }
            }
          if ( index == 0 || index == 255 ) //no surface
            {
            continue;
            }

          triCase = triCases+ index; //triCases是vtkMarchingCubesSpaceTriangleCases.cxx中“构型—三角剖分”查找表的头指针,
          edge = triCase->edges; //edges是vtkMarchingCubesSpaceTriangleCases结构体成员EDGE_LIST edges[16];edges[16]每个元素对应cube边的序号;edge即edges[16]中头指针

          for ( ; edge[0] > -1; edge += 3 )  //每三个边就是一个triangle，所以edge+=3
            {
            for (ii=0; ii<3; ii++) //insert triangle
			 {
              vert = edges[edge[ii]];   //函数开始定义edges[12][2] = {{0,1},{1,2},{3,2},{0,3},{4,5},{5,6},{7,6},{4,7},{0,4},{1,5},{3,7},{2,6}};vert是每条棱对应的头指针
              //t = (value - s[vert[0]]) / (s[vert[1]] - s[vert[0]]);
              x1 = pts[vert[0]];   //函数开始定义pts[8][3]
              x2 = pts[vert[1]];
			  x[0] =0.5*(x1[0] + x2[0]);
              x[1] =0.5*(x1[1] + x2[1]);
              x[2] =0.5*(x1[2] + x2[2]);       //线性插值改为取中点

              // check for a new point
              if ( locator->InsertUniquePoint(x, ptIds[ii]) ) //Determine whether point given by x[3] has been inserted into points list
                  {                                           //Return 0 if point was already in the list, otherwise return 1.-- 
                  if (NeedGradients)                          //---If the point was not in the list, it will be ADDED
                    {
                    n1 = gradients[vert[0]];
                    n2 = gradients[vert[1]];
                    n[0] = 0.5*(n1[0] + n2[0]);
                    n[1] = 0.5*(n1[1] + n2[1]);
                    n[2] = 0.5*(n1[2] + n2[2]);    //梯度计算改为取平均
                    }
                  if (ComputeScalars)
                    {
                    newScalars->InsertTuple(ptIds[ii],&value);
                    }
                  if (ComputeGradients)
                    {
                    newGradients->InsertTuple(ptIds[ii],n);
                    }
                  if (ComputeNormals)
                    {
                    vtkMath::Normalize(n);
                    newNormals->InsertTuple(ptIds[ii],n);
                    }
                  }
			}
            // check for degenerate triangle
            if ( ptIds[0] != ptIds[1] &&
                 ptIds[0] != ptIds[2] &&
                 ptIds[1] != ptIds[2] )
                {
                newPolys->InsertNextCell(3,ptIds);
                }
            }//for each triangle
          }//for all contours
        }//for i
      }//for j
    }//for k
}*/



//
// My vtkMarchingCubesSpaceComputeGradient
//
template <class T>
void vtkMarchingCubesSpaceComputeGradient(vtkMarchingCubesSpace *self,T *scalars, T *scalars1,int dims[3],   //scalars1作为原始图像的标量输入
										  double origin[3], double spacing[3],
										  vtkIncrementalPointLocator *locator,
										  vtkDataArray *newScalars,
										  vtkDataArray *newGradients,
										  vtkDataArray *newNormals,
										  vtkCellArray *newPolys, double *values,
										  int numValues)
{
	double s[8], value;
	int i, j, k;
	vtkIdType sliceSize;
	static int CASE_MASK[8] = {1,2,4,8,16,32,64,128};  //即0000 0001,000 0010等等，与index或运算，index即cube八个顶点状态
	vtkMarchingCubesSpaceTriangleCases *triCase, *triCases;  //结构体vtkMarchingCubesSpaceTriangleCases含有edges[16](cube边所对应的triangle)和指向结构体指针函数
	EDGE_LIST  *edge;
	int contNum, jOffset, ii, index, *vert;
	vtkIdType kOffset, idx;
	vtkIdType ptIds[3];
	int ComputeNormals = newNormals != NULL;
	int ComputeGradients = newGradients != NULL;
	int ComputeScalars = newScalars != NULL;
	int NeedGradients;
	int extent[6];
	double t, *x1, *x2, x[3], *n1, *n2, n[3], min, max;
	double pts[8][3], gradients[8][3], xp, yp, zp;      //pts[8][3]中，8是cube八个顶点，3是对应x,y,z坐标
	static int edges[12][2] = { {0,1}, {1,2}, {3,2}, {0,3},
	{4,5}, {5,6}, {7,6}, {4,7},
	{0,4}, {1,5}, {3,7}, {2,6}};  // 0-7是cube顶点编号，{0,1}表示由0号顶点和1号顶点之间的棱

	vtkMarchingCubesSpaceNeighbourCases *neiborCase, *neiborCases;  //结构体vtkMarchingCubesSpaceNeighbourCases含有neighbours[6](cube的邻居)和指向结构体指针函数
	bool  *neighbour;
	std::queue<Cube> CubeQueue;
	std::vector<Cube> List(dims[0]*dims[1]*dims[2]);	


	vtkInformation *inInfo = self->GetExecutive()->GetInputInformation(0, 0); //Get this algorithm's executive,Get the pipeline information for the given input connection.
	inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent);

	triCases =  vtkMarchingCubesSpaceTriangleCases::GetCases();  //得到vtkMarchingCubesSpaceCases.cxx中静态成员“构型—三角剖分”查找表的头指针
	neiborCases =  vtkMarchingCubesSpaceNeighbourCases::GetCases();  //得到vtkMarchingCubesSpaceCases.cxx中静态成员邻接查找表的头指针

	//
	// Get min/max contour values
	//
	if ( numValues < 1 )
	{
		return;
	}
	for ( min=max=values[0], i=1; i < numValues; i++)
	{
		if ( values[i] < min )
		{
			min = values[i];
		}
		if ( values[i] > max )
		{
			max = values[i];
		}
	}


	//
	// 发现index不为0或255的立方体
	//
	sliceSize = dims[0] * dims[1];
	for ( k=0; k < (dims[2]-1); k++)
	{
		self->UpdateProgress (k / static_cast<double>(dims[2] - 1));
		if (self->GetAbortExecute())
		{
			break;
		}
		kOffset = k*sliceSize;
		pts[0][2] = origin[2] + (k+extent[4]) * spacing[2];
		zp = pts[0][2] + spacing[2];
		for ( j=0; j < (dims[1]-1); j++)
		{
			jOffset = j*dims[0];
			pts[0][1] = origin[1] + (j+extent[2]) * spacing[1];
			yp = pts[0][1] + spacing[1];
			for ( i=0; i < (dims[0]-1); i++)
			{
				//get scalar values
				idx = i + jOffset + kOffset;
				s[0] = scalars[idx];                          //point(i,j,k)对应标量
				s[1] = scalars[idx+1];                       //point(i+1,j,k）对应标量
				s[2] = scalars[idx+1 + dims[0]];             //point(i+1,j+1,k)对应标量
				s[3] = scalars[idx + dims[0]];               //point(i,j+1,k)对应标量
				s[4] = scalars[idx + sliceSize];             //point(i,j,k+1)对应标量
				s[5] = scalars[idx+1 + sliceSize];           //point(i+1,j,k+1)对应标量
				s[6] = scalars[idx+1 + dims[0] + sliceSize]; //point(i+1,j+1,k+1)对应标量
				s[7] = scalars[idx + dims[0] + sliceSize];   //point(i,j+1,k+1)对应标量

				if ( (s[0] < min && s[1] < min && s[2] < min && s[3] < min &&
					s[4] < min && s[5] < min && s[6] < min && s[7] < min) ||
					(s[0] > max && s[1] > max && s[2] > max && s[3] > max &&
					s[4] > max && s[5] > max && s[6] > max && s[7] > max) )
				{
					continue; // no contours possible
				} 

				for (contNum=0; contNum < numValues; contNum++)
				{
					value = values[contNum];
					// Build the case table
					for ( ii=0, index = 0; ii < 8; ii++)
					{
						if ( s[ii] >= value )
						{
							index |= CASE_MASK[ii];  //index是cube的索引号（0-255，表示vertex八个状态）
						}
					}
					if ( index == 0 || index == 255 ) //no surface
					{
						continue;
					}
					goto FindSeedPoint;
				}//for all contours
			}//for i
		}//for j
	}//for k



FindSeedPoint:
	bool flag=0;Cube CubeTemp;
	CubeQueue.push(Cube(i,j,k,index,flag));

	while(!CubeQueue.empty()){
		CubeTemp=CubeQueue.front();  	 
		i=CubeTemp.oi;j=CubeTemp.oj;k=CubeTemp.ok;index=CubeTemp.index;
		self->UpdateProgress (k / static_cast<double>(dims[2] - 1));
		if (self->GetAbortExecute())
		{
			break;
		}
		kOffset = k*sliceSize;
		pts[0][2] = origin[2] + (k+extent[4]) * spacing[2];
		zp = pts[0][2] + spacing[2];

		jOffset = j*dims[0];
		pts[0][1] = origin[1] + (j+extent[2]) * spacing[1];
		yp = pts[0][1] + spacing[1];


		idx = i + jOffset + kOffset;
		CubeQueue.pop();
		if(List[idx].Flag==1)
			continue;
		List[idx].Flag=1;


		s[0] = scalars[idx];                          //point(i,j,k)对应标量
		s[1] = scalars[idx+1];                       //point(i+1,j,k）对应标量
		s[2] = scalars[idx+1 + dims[0]];             //point(i+1,j+1,k)对应标量
		s[3] = scalars[idx + dims[0]];               //point(i,j+1,k)对应标量
		s[4] = scalars[idx + sliceSize];             //point(i,j,k+1)对应标量
		s[5] = scalars[idx+1 + sliceSize];           //point(i+1,j,k+1)对应标量
		s[6] = scalars[idx+1 + dims[0] + sliceSize]; //point(i+1,j+1,k+1)对应标量
		s[7] = scalars[idx + dims[0] + sliceSize];   //point(i,j+1,k+1)对应标量

		if ( (s[0] < min && s[1] < min && s[2] < min && s[3] < min &&
			s[4] < min && s[5] < min && s[6] < min && s[7] < min) ||
			(s[0] > max && s[1] > max && s[2] > max && s[3] > max &&
			s[4] > max && s[5] > max && s[6] > max && s[7] > max) )
		{
			continue; // no contours possible
		}


		//create voxel points
		pts[0][0] = origin[0] + (i+extent[0]) * spacing[0];
		xp = pts[0][0] + spacing[0];

		pts[1][0] = xp;
		pts[1][1] = pts[0][1];
		pts[1][2] = pts[0][2];

		pts[2][0] = xp;
		pts[2][1] = yp;
		pts[2][2] = pts[0][2];

		pts[3][0] = pts[0][0];
		pts[3][1] = yp;
		pts[3][2] = pts[0][2];

		pts[4][0] = pts[0][0];
		pts[4][1] = pts[0][1];
		pts[4][2] = zp;

		pts[5][0] = xp;
		pts[5][1] = pts[0][1];
		pts[5][2] = zp;

		pts[6][0] = xp;
		pts[6][1] = yp;
		pts[6][2] = zp;

		pts[7][0] = pts[0][0];
		pts[7][1] = yp;
		pts[7][2] = zp;

		NeedGradients = ComputeGradients || ComputeNormals;

		//create gradients if needed
		if (NeedGradients)
		{
			vtkMarchingCubesSpaceComputePointGradient(i,j,k, scalars1, dims, sliceSize, spacing, gradients[0]);   //前面定义了double gradients[8][3]
			vtkMarchingCubesSpaceComputePointGradient(i+1,j,k, scalars1, dims, sliceSize, spacing, gradients[1]);
			vtkMarchingCubesSpaceComputePointGradient(i+1,j+1,k, scalars1, dims, sliceSize, spacing, gradients[2]);
			vtkMarchingCubesSpaceComputePointGradient(i,j+1,k, scalars1, dims, sliceSize, spacing, gradients[3]);
			vtkMarchingCubesSpaceComputePointGradient(i,j,k+1, scalars1, dims, sliceSize, spacing, gradients[4]);
			vtkMarchingCubesSpaceComputePointGradient(i+1,j,k+1, scalars1, dims, sliceSize, spacing, gradients[5]);
			vtkMarchingCubesSpaceComputePointGradient(i+1,j+1,k+1, scalars1, dims, sliceSize, spacing, gradients[6]);
			vtkMarchingCubesSpaceComputePointGradient(i,j+1,k+1, scalars1, dims, sliceSize, spacing, gradients[7]);
		}

		for (contNum=0; contNum < numValues; contNum++)
		{
			value = values[contNum];
			// Build the case table
			for ( ii=0, index = 0; ii < 8; ii++)
			{
				if ( s[ii] >= value )
				{
					index |= CASE_MASK[ii];  //index是cube的索引号（0-255，表示vertex八个状态）
				}
			}
			if ( index == 0 || index == 255 ) //no surface
			{
				continue;
			}

			triCase = triCases+ index; //triCases是vtkMarchingCubesSpaceCases.cxx中“构型—三角剖分”查找表的头指针,
			edge = triCase->edges;    //edges是vtkMarchingCubesSpaceCases结构体成员EDGE_LIST edges[16];edges[16]每个元素对应cube边的序号;edge即edges[16]中头指针
			neiborCase = neiborCases+ index; //neiborCases是vtkMarchingCubesSpaceCases.cxx中邻接查找表的头指针,
			neighbour = neiborCase->neighbours;    //neighbours是vtkMarchingCubesSpaceNeighbourCases结构体成员bool neighbours[6];neighbour即neighbours[6]中头指针

			//cout<<index<<endl;
			for ( ; edge[0] > -1; edge += 3 )  //每三个边就是一个triangle，所以edge+=3
			{
				for (ii=0; ii<3; ii++) //insert triangle
				{
					vert = edges[edge[ii]];   //函数开始定义edges[12][2] = {{0,1},{1,2},{3,2},{0,3},{4,5},{5,6},{7,6},{4,7},{0,4},{1,5},{3,7},{2,6}};vert是每条棱对应的头指针
					//t = (value - s[vert[0]]) / (s[vert[1]] - s[vert[0]]);
					x1 = pts[vert[0]];   //函数开始定义pts[8][3]
					x2 = pts[vert[1]];
					x[0] =0.5*(x1[0] + x2[0]);
					x[1] =0.5*(x1[1] + x2[1]);
					x[2] =0.5*(x1[2] + x2[2]);       //线性插值改为取中点

					// check for a new point
					if ( locator->InsertUniquePoint(x, ptIds[ii]) )
					{
						if (NeedGradients)
						{
							n1 = gradients[vert[0]];
							n2 = gradients[vert[1]];
							n[0] = 0.5*(n1[0] + n2[0]);
							n[1] = 0.5*(n1[1] + n2[1]);
							n[2] = 0.5*(n1[2] + n2[2]);    //梯度计算改为去平均
						}
						if (ComputeScalars)
						{
							newScalars->InsertTuple(ptIds[ii],&value);
						}
						if (ComputeGradients)
						{
							newGradients->InsertTuple(ptIds[ii],n);
						}
						if (ComputeNormals)
						{
							vtkMath::Normalize(n);
							newNormals->InsertTuple(ptIds[ii],n);
						}
					}
				}
				// check for degenerate triangle
				if ( ptIds[0] != ptIds[1] &&
					ptIds[0] != ptIds[2] &&
					ptIds[1] != ptIds[2] )
				{
					newPolys->InsertNextCell(3,ptIds);
				}
			}//for each triangle
		}//for all contours


		if(neighbour[0]){     //front direction
			if(i<(dims[0]-1)){        
				List[idx+1].oi=i+1;List[idx+1].oj=j;List[idx+1].ok=k;List[idx+1].index=index;
				CubeQueue.push(List[idx+1]);
			}
		}
		if(neighbour[1]){   //back direction
			if(i>0){          
				List[idx-1].oi=i-1;List[idx-1].oj=j;List[idx-1].ok=k;List[idx-1].index=index;
				CubeQueue.push(List[idx-1]);
			}
		}
		if(neighbour[2]){   //right direction
			if(j<(dims[1]-1)){         
				List[idx + dims[0]].oi=i;List[idx + dims[0]].oj=j+1;List[idx + dims[0]].ok=k;List[idx + dims[0]].index=index;
				CubeQueue.push(List[idx + dims[0]]);
			}
		}
		if(neighbour[3]){   //left direction
			if(j>0){          
				List[idx - dims[0]].oi=i;List[idx - dims[0]].oj=j-1;List[idx - dims[0]].ok=k;List[idx - dims[0]].index=index;
				CubeQueue.push(List[idx - dims[0]]);
			}
		}
		if(neighbour[4]){       //up direction 
			if(k<(dims[2]-2)){        // 注意k<dims[2]-1的话，正方体相对原点到顶层会超界    
				List[idx + sliceSize].oi=i;List[idx + sliceSize].oj=j;List[idx + sliceSize].ok=k+1;List[idx + sliceSize].index=index;		
				CubeQueue.push(List[idx + sliceSize]);
			}
		}
		if(neighbour[5]){  //left direction
			if(k>0){         
				List[idx - sliceSize].oi=i;List[idx - sliceSize].oj=j;List[idx - sliceSize].ok=k-1;List[idx - sliceSize].index=index;
				CubeQueue.push(List[idx - sliceSize]);
			}	 
		}
	}//exit the queue


}






//
// Contouring filter specialized for volumes and "short int" data values.
//
int vtkMarchingCubesSpace::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkImageData *input = vtkImageData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkImageData *input1 = vtkImageData::SafeDownCast(
	inInfo1->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPoints *newPts;
  vtkCellArray *newPolys;
  vtkFloatArray *newScalars;
  vtkFloatArray *newNormals;
  vtkFloatArray *newGradients;
  vtkPointData *pd; vtkPointData *pd1; 
  vtkDataArray *inScalars;vtkDataArray *inScalars1;
  int dims[3], extent[6];
  vtkIdType estimatedSize;
  double spacing[3], origin[3];
  double bounds[6];
  int numContours=this->ContourValues->GetNumberOfContours();
  double *values=this->ContourValues->GetValues();

  vtkDebugMacro(<< "Executing marching cubes");

//
// Initialize and check input
//
  pd=input->GetPointData();pd1=input1->GetPointData(); 
  if (pd ==NULL)
    {
    vtkErrorMacro(<<"PointData is NULL");
    return 1;
    }
  inScalars=pd->GetScalars();inScalars1=pd1->GetScalars();
  if ( inScalars == NULL )
    {
    vtkErrorMacro(<<"Scalars must be defined for contouring");
    return 1;
    }

  if ( input->GetDataDimension() != 3 )
    {
    vtkErrorMacro(<<"Cannot contour data of dimension != 3");
    return 1;
    }
  input->GetDimensions(dims);
  input->GetOrigin(origin);
  input->GetSpacing(spacing);

  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);

  // estimate the number of points from the volume dimensions
  estimatedSize = static_cast<vtkIdType>(
    pow(1.0*dims[0]*dims[1]*dims[2], 0.75));
  estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
  if (estimatedSize < 1024)
    {
    estimatedSize = 1024;
    }
  vtkDebugMacro(<< "Estimated allocation size is " << estimatedSize);
  newPts = vtkPoints::New(); newPts->Allocate(estimatedSize,estimatedSize/2);
  // compute bounds for merging points
  for ( int i=0; i<3; i++)
    {
    bounds[2*i] = origin[i] + extent[2*i] * spacing[i];
    bounds[2*i+1] = origin[i] + extent[2*i+1] * spacing[i];
    }
  if ( this->Locator == NULL )
    {
    this->CreateDefaultLocator();
    }
  this->Locator->InitPointInsertion (newPts, bounds, estimatedSize);

  if (this->ComputeNormals)
    {
    newNormals = vtkFloatArray::New();
    newNormals->SetNumberOfComponents(3);
    newNormals->Allocate(3*estimatedSize,3*estimatedSize/2);
    }
  else
    {
    newNormals = NULL;
    }

  if (this->ComputeGradients)
    {
    newGradients = vtkFloatArray::New();
    newGradients->SetNumberOfComponents(3);
    newGradients->Allocate(3*estimatedSize,3*estimatedSize/2);
    }
  else
    {
    newGradients = NULL;
    }

  newPolys = vtkCellArray::New();
  newPolys->Allocate(newPolys->EstimateSize(estimatedSize,3));

  if (this->ComputeScalars)
    {
    newScalars = vtkFloatArray::New();
    newScalars->Allocate(estimatedSize,estimatedSize/2);
    }
  else
    {
    newScalars = NULL;
    }

  if (inScalars->GetNumberOfComponents() == 1 )
    {
    void* scalars = inScalars->GetVoidPointer(0);void* scalars1 = inScalars1->GetVoidPointer(0);
    switch (inScalars->GetDataType())
      {
      vtkTemplateMacro(
        vtkMarchingCubesSpaceComputeGradient(this, static_cast<VTK_TT*>(scalars), static_cast<VTK_TT*>(scalars1),
                                        dims,origin,spacing,this->Locator,
                                        newScalars,newGradients,
                                        newNormals,newPolys,values,
                                        numContours)
        );
      } //switch
    }

  else //multiple components - have to convert
    {
    vtkIdType dataSize = static_cast<vtkIdType>(dims[0]) * dims[1] * dims[2];
    vtkDoubleArray *image=vtkDoubleArray::New();                          vtkDoubleArray *image1=vtkDoubleArray::New();
    image->SetNumberOfComponents(inScalars->GetNumberOfComponents());     image1->SetNumberOfComponents(inScalars1->GetNumberOfComponents());
    image->SetNumberOfTuples(image->GetNumberOfComponents()*dataSize);    image1->SetNumberOfTuples(image1->GetNumberOfComponents()*dataSize);
    inScalars->GetTuples(0,dataSize,image);                               inScalars1->GetTuples(0,dataSize,image);

    double *scalars = image->GetPointer(0);                               double *scalars1 = image1->GetPointer(0);
    vtkMarchingCubesSpaceComputeGradient(this,scalars,scalars1,dims,origin,spacing,this->Locator,
                  newScalars,newGradients,
                  newNormals,newPolys,values,numContours);
    image->Delete();
    }

  vtkDebugMacro(<<"Created: "
               << newPts->GetNumberOfPoints() << " points, "
               << newPolys->GetNumberOfCells() << " triangles");
  //
  // Update ourselves.  Because we don't know up front how many triangles
  // we've created, take care to reclaim memory.
  //
  output->SetPoints(newPts);
  newPts->Delete();

  output->SetPolys(newPolys);
  newPolys->Delete();

  if (newScalars)
    {
    int idx = output->GetPointData()->AddArray(newScalars);
    output->GetPointData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    newScalars->Delete();
    }
  if (newGradients)
    {
    output->GetPointData()->SetVectors(newGradients);
    newGradients->Delete();
    }
  if (newNormals)
    {
    output->GetPointData()->SetNormals(newNormals);
    newNormals->Delete();
    }
  output->Squeeze();
  if (this->Locator)
    {
    this->Locator->Initialize(); //free storage
    }

  return 1;
}

// Description:
// Specify a spatial locator for merging points. By default,
// an instance of vtkMergePoints is used.
void vtkMarchingCubesSpace::SetLocator(vtkIncrementalPointLocator *locator)
{
  if ( this->Locator == locator )
    {
    return;
    }

  if ( this->Locator )
    {
    this->Locator->UnRegister(this);
    this->Locator = NULL;
    }

  if (locator)
    {
    locator->Register(this);
    }

  this->Locator = locator;
  this->Modified();
}

void vtkMarchingCubesSpace::CreateDefaultLocator()
{
  if ( this->Locator == NULL)
    {
    this->Locator = vtkMergePoints::New();
    }
}

// 1 input
/*int vtkMarchingCubesSpace::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}*/

// 2 inputs
int vtkMarchingCubesSpace::FillInputPortInformation(int port, vtkInformation *info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
    }
  else if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
    }
  return 0;
}

void vtkMarchingCubesSpace::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  this->ContourValues->PrintSelf(os,indent.GetNextIndent());

  os << indent << "Compute Normals: " << (this->ComputeNormals ? "On\n" : "Off\n");
  os << indent << "Compute Gradients: " << (this->ComputeGradients ? "On\n" : "Off\n");
  os << indent << "Compute Scalars: " << (this->ComputeScalars ? "On\n" : "Off\n");

  if ( this->Locator )
    {
    os << indent << "Locator:" << this->Locator << "\n";
    this->Locator->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Locator: (none)\n";
    }
}


