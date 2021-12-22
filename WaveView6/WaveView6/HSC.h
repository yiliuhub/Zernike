#pragma once
#include <string>

#include "CImg.h"
using namespace cimg_library;

typedef unsigned int  UInt;

const float PIXEL_SIZE_MM	= 4.65E-3;
const float FOCUS_LENGTH_MM	= 5.2;
const float NORMAL_RADIUS_MM	= 1.3; 

#define INVALID_CENTROID 100;
#define MIISED_CENTROID 101;
#define CENTROID_FLAG	102;
#define CENTROID_REGION 103;
#define HS_STATIONARY_CENTROID 104;

class HSC
{
public:
  HSC(void);
  ~HSC(void);

  template<typename T>
  bool  SetAndCalcCentroidWithConfig(
							T* data, int width, int height,		// image data
							float* pixelX, float* pixelY,			// reference data
							int* gridX, int* gridY, int* flag, int refSize,		// reference data
							float* locX, float* locY, float* dx, float* y, int& dataSize,	// output data size
							//float Pixel_Size_MM = 4.65E-3, float Focus_Length_MM = 5.2, float Normal_Radius_MM	= 1.3	// hardware parameters
							float Pixel_Size_MM, float Focus_Length_MM, float Normal_Radius_MM,
							int brightSpotRadius
						  );
  
  template<typename T>
  bool  SetAndCalcCentroid(
							T* data, int width, int height,		// image data
							float* pixelX, float* pixelY,			// reference data
							int* gridX, int* gridY, int* flag, int refSize,		// reference data
							float* locX, float* locY, float* dx, float* y, int& dataSize	// output data size
						  );
 
  bool  CalcKeratometerRingSpots( float* data, int width, int height,
									float* locX, float* locY, int& dataSize,	// output data size
									int refRadius=0, int refLocX=0, int refLocY=0
								  );

  bool  CalcKeratometerRingSpotsWithReference(float* data, int width, int height,
											  int refRadius, int refLocX, int refLocY,
											  float* locX, float* locY, int& dataSize	// output data size
											  );

  float CalcContrast( float* imageData, int imageWidth, int imagHeight, int refRadius, int refLocX, int refLocY );

  template<typename T>
  bool  SetStreamData( T* data, int width, int height );
  void  SetHardwareConfig(float Pixel_Size_MM, float Focus_Length_MM, float Normal_Radius_MM);
  void  SetBrightSpotRadius( int Bright_Spot_Radius );
  void  SetReferenceData( float* pixelX, float* pixelY, 
						  int* gridX, int* gridY, int* flag, int size );
  void  GetCentroidResult( float* locX, float* locY, float* dx, float* dy, int& dataSize );

  bool  OnExcute();
  bool  OpenImageFromFile(std::string filename);
  bool  GoHSCentroid();
  void  ResultDisplay();
  void  Clear();
  void  OutputResult( std::string outputName );
  void  OpenReferenceCentroid();
  void  SaveAsReference();
  void  Test();

private:
  bool  FindGridPicth();
  bool  FindGridPicthNew();
  bool  LocateMissingCentroid();
  float FindLocalMaxPixelPos(UInt x0, UInt y0, UInt x1, UInt y1, UInt& maxX, UInt& maxY );
  float FindLocalMaxPixelPos1(UInt x0, UInt y0, UInt x1, UInt y1, UInt& maxX, UInt& maxY );
  template<typename T>
  float FindLocalMaxPixelPos1(CImg<T> img, int x0, int y0, int x1, int y1, int& maxX, int& maxY );
  float FindLocalMaxPixelPos2(UInt& maxX, UInt& maxY );
  void  SearchInDirection(const UInt startX, const UInt startY, const UInt startIndexX, const UInt startIndexY, const int dirX, const int dirY);
  void  RadiantSearch();
  void  OutputImage();
  void  UpdateResultImage();
  void  CalcSubpixel();
  template<typename T>
  void  CalcSubpixel(CImg<T> img, int locX, int locY, float& subLocX, float& subLocY);
  void  CentroidValidation();
  void  LocateSearchCenter();
  bool  CalcStepSizeAndStartIndex();
  void  DrawCentroidDisplacement();
  bool  ReferenceAssociation();
  bool  ReferenceAssociationNew();
  template<typename T>
  bool  ReferenceAssociationNew(CImg<T> ref4x4, CImg<T> targetImg, int& associationX, int& associationY);
  void  SmoothAndDenoise();
  void  ConvertResultToMM();
  
  void CalcKeratometerRingCenter( 
		 float* data,int width, int height, 
		 int refRadius, int refLocX, int refLocY , 
		 int& centerX, int& centerY 
	   );
  
  void CalcKeratometerRingMirrorCenter( 
	     float* data,int width, int height, 
		 int refRadius, int refLocX, int refLocY , 
		 int& centerX, int& centerY 
	   );

  void CalcKeratometerRingRadius( 
		 float* data,int width, int height, 
		 int centerX, int centerY, float& radius 
	   );

private:

  CImg<float> m_imgInput;
  CImg<float> m_image;  
  CImg<float> m_resultImage;
  CImg<float> m_resultMatrix;
  CImg<float> m_referenceMatrix;
  CImg<float> m_resultInMM;


  UInt m_stepX;
  UInt m_stepY;
  UInt m_centerX;
  UInt m_centerY;
  UInt m_refCenterX;
  UInt m_refCenterY;
  double m_discCenterX;
  double m_discCenterY;
  UInt m_LSRadius;                  // local search radius
  UInt m_startIndexX;
  UInt m_startIndexY;
  float m_CENTROID_THRESHOLD;        // centroid threshold
  float m_minGridSize;              // minimum grid pitch size
  int  m_referenceAssocX;
  int  m_referenceAssocY;
  int  m_stationaryX;
  int  m_stationaryY;

  // hardware parameters
  float m_PIXEL_SIZE_MM;
  float m_FOCUS_LENGTH_MM;
  float m_NORMAL_RADIUS_MM; 

  // bright spot removing parameters
  int m_brightSpotRadius;
  int m_removeCenterX;
  int m_removeCenterY;
  
};
