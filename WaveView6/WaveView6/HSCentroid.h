#pragma once
#include <string>
#include "stdlib.h"
#include "BasicDefine.h"

using namespace LLW;

#include "CImg.h"
using namespace cimg_library;

typedef unsigned int  UInt;

#define PIXLE_SIZE_MM		4.65E-3;
#define NORMAL_RADIUS_MM	1.3; 
#define FOCUS_LENGTH_MM		5.2; //

enum Centroid_Enum
{
	INVALID_CENTROID = 100,
	CENTROID_FLAG,
	STATIONARY_CENTROID
};


namespace LLW 
{
	namespace Centroid
	{
		class HSCentroid
		{
		public:
		  HSCentroid(void);
		  ~HSCentroid(void);

		  
		  template<typename T>
		  void  SetAndCalcCentroid(
									T* data, int width, int height,		// image data
									float* pixelX, float* pixelY,			// reference data
									int* gridX, int* gridY, int* flag, int refSize,		// reference data
									float* locX, float* locY, float* dx, float* y, int& dataSize	// output data size
								  );
		  
		  template<typename T>
		  void  SetStreamData( T* data, int width, int height );

		  void  SetReferenceData( float* pixelX, float* pixelY, 
								  int* gridX, int* gridY, int* flag, int size );

		  void  GetCentroidResult( float* locX, float* locY, float* dx, float* y, int& dataSize );
		  bool  OnExcute();
		  bool  OpenImageFromFile(const std::string filename);
		  bool  GoHSCentroid();
		  void  ResultDisplay();
		  void  Clear();
		  void  OutputResult( std::string outputName );
		  void  OpenReferenceCentroid();
		  void  SaveAsReference();
		  void  Test();

		private:
		  bool  FindGridPicth();
		  bool 	FindGridPicthNew();
		  bool  LocateMissingCentroid();
		  float FindLocalMaxPixelPos(int x0, int y0, int x1, int y1, int& maxX, int& maxY );
		  float FindLocalMaxPixelPos1(int x0, int y0, int x1, int y1, int& maxX, int& maxY );
		  template<typename T>
		  float FindLocalMaxPixelPos1(CImg<T> img, int x0, int y0, int x1, int y1, int& maxX, int& maxY );

		  float FindLocalMaxPixelPos2(UInt& maxX, UInt& maxY );
		  void  SearchInDirection(const UInt startX, const UInt startY, const UInt startIndexX, const UInt startIndexY, const int dirX, const int dirY);
		  void  RadiantSearch();
		  void  OutputImage();
		  void  UpdateResultImage();
		  void  CalcSubpixel();
		  void  CentroidValidation();
		  void  LocateSearchCenter();
		  bool  CalcStepSizeAndStartIndex();
		  void  DrawCentroidDisplacement();
		  bool  ReferenceAssociation();
		  void  SmoothAndDenoise();
		  void  ConvertResultToMM();
		  

		private:
		  CImg<float> m_image;  
		  CImg<float> m_resultImage;
		  CImg<float> m_centroidMatrix;
		  CImg<float> m_referenceMatrix;
		  CImg<float> m_resultInMM;
		  CImg<float> m_resultMatrix;


		  LLW::Point2d m_step;
		  LLW::Point2n m_center;
		  LLW::Point2n m_startIndex;


		  //UInt m_stepX;
		  //UInt m_stepY;
		  //UInt m_centerX;
		  //UInt m_centerY;
		  //UInt m_refCenterX;
		  //UInt m_refCenterY;
		  double m_discCenterX;
		  double m_discCenterY;
		  UInt m_LSRadius;                  // local search radius
		  UInt m_startIndexX;
		  UInt m_startIndexY;
		  float m_CENTROID_THRESHOLD;       // centroid threshold
		  float m_minGridSize;              // minimum grid pitch size
		  int  m_referenceAssocX;
		  int  m_referenceAssocY;
		  int  m_stationaryX;
		  int  m_stationaryY;
		};
	}
}