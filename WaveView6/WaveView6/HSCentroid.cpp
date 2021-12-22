#include "StdAfx.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "HSCentroid.h"

using namespace std;
#include "CImg.h"
using namespace cimg_library;

namespace LLW {
 namespace Centroid	{

	HSCentroid::HSCentroid(void)
	{
		m_image.assign();  
		m_resultImage.assign();
		m_centroidMatrix.assign();
		m_referenceMatrix.assign();
		m_resultInMM.assign();

	}

	HSCentroid::~HSCentroid(void)
	{
		m_image.assign();
		m_resultImage.assign();
		m_centroidMatrix.assign();
		m_referenceMatrix.assign();
		m_resultInMM.assign();
	}

	bool 
	HSCentroid::OpenImageFromFile(const std::string filename)
	{
		CImg<float> tempImage = CImg<float>(filename.c_str()).get_channel(0);
		SetStreamData<float>(tempImage.data(), tempImage.width(), tempImage.height() );

		SmoothAndDenoise();

		m_image.display("fjdksl", true);

		FindGridPicth();
		FindGridPicthNew();
	
		return true;
	}

	template<typename T>
	void HSCentroid::SetStreamData( T* data, int width, int height )
	{
		m_image.assign(data, width, height);
	}

	void HSCentroid::SmoothAndDenoise()
	{
		float diffsion = 2.0;
		float threshold = 100;
		int erodeVal = 10;
		int dilateVal = 10;

		CImg<float> res = m_image.get_blur(diffsion);

		int x[2] = { (int)(0.25*m_image.width()), (int)(0.75*m_image.width()) };
		int y[2] = { (int)(0.25*m_image.height()), (int)(0.75*m_image.height()) };

		float variance = m_image.get_crop(x[0], y[0], x[1], y[1]).variance();
		threshold =  sqrt(variance) + 2*m_image.get_crop(x[0], y[0], x[1], y[1]).mean();
		res.threshold(threshold, false, false);

		res.erode(erodeVal);
		res.dilate(dilateVal);

		cimg_forXY(m_image, x,y)
		{ 
		  if(res(x,y))
			  m_image(x,y) = 0;
		}

		m_image.blur(diffsion);
	}

	float HSCentroid::FindLocalMaxPixelPos1(int x0, int y0, int x1, int y1, int& maxX, int& maxY )
	{
	  CImg<float> cropImg = m_image.get_crop(x0, y0, x1, y1);
	  float *ptr = cropImg._data;
	  float max_value = -9999999;
	  long count = 0;
	  long maxIndex = 0;
	  long size = cropImg.size();

	  for(int i = 0; i<size;i++)
	  {
		const float val = *ptr++;
		if (val>max_value) {max_value = val; maxIndex = count; }
		count++;
	  }

	  maxX = x0 + maxIndex % cropImg.width();
	  maxY = y0 + (int)(maxIndex / cropImg.width());

	  return max_value;
	}

	template<typename T>
	float HSCentroid::FindLocalMaxPixelPos1(CImg<T> img, int x0, int y0, int x1, int y1, int& maxX, int& maxY )
	{
	  CImg<float> cropImg = img.get_crop(x0, y0, x1, y1);
	  float *ptr = cropImg._data;
	  float max_value = -9999999;
	  long count = 0;
	  long maxIndex = 0;
	  long size = cropImg.size();

	  for(int i = 0; i<size;i++)
	  {
		const float val = *ptr++;
		if (val>max_value) {max_value = val; maxIndex = count; }
		count++;
	  }

	  maxX = x0 + maxIndex % cropImg.width();
	  maxY = y0 + (int)(maxIndex / cropImg.width());

	  return max_value;
	}

	float HSCentroid::FindLocalMaxPixelPos(int x0, int y0, int x1, int y1, int& maxX, int& maxY )
	{
	  CImg<float> cropImg = m_image.get_crop(x0, y0, x1, y1);
	  float max_value = -9999999;
	  long count = 0;
	  long maxIndex = 0;

	  cimg_for(cropImg,ptrs,float) {
		count++;
		const float val = *ptrs;
		if (val>max_value) {max_value = val; maxIndex = count; }
	  }

	  long reversedNum = cropImg.size() - maxIndex;
	  maxX = x0 + reversedNum % cropImg.width();
	  maxY = y0 + (int)(reversedNum / cropImg.width());

	  return max_value;
	}

	//////////////////////////////////////////////////////////////////////////////////
	// this function does two things: find the grid picthes of centroid spots 
	// in x and y direction, and figure out the result matrix dimmension
	//////////////////////////////////////////////////////////////////////////////////
	bool 
	HSCentroid::FindGridPicth()
	{
	  bool bOk = true;
	  
	  int x[2] = {(int)( 0.4*m_image.width() ),  (int)( 0.6*m_image.width() ) };
	  int y[2] = {(int)( 0.6*m_image.height() ), (int)( 0.8*m_image.height() ) };

	  int count = 0;
	  Point2n localMax;

	  float maxVal = FindLocalMaxPixelPos1(x[0], y[0], x[1], y[1], m_center.x, m_center.y);

	  // stepX
	  int deltaX=5,	deltaY = 5;
	  x[0] = m_center.x;			y[0] = m_center.y - deltaY;
	  x[1] = x[0] + 4*deltaX;		y[1] = m_center.y + deltaY;

	  bOk = false;
	  while( !bOk && x[0] < this->m_image.width() - 4*(int)deltaX)
	  {
		x[0] += deltaX;
		x[1] = x[0] + 4*deltaX;

		maxVal = FindLocalMaxPixelPos1( x[0], y[0], x[1], y[1], localMax.x, localMax.y );

		if(		(maxVal < m_CENTROID_THRESHOLD)
			||	(int)localMax.x - x[0] < deltaX
			||	(x[1] - (int)localMax.x < deltaX) )
		{
		  bOk = false;
		}
		else
		  bOk = true;
	  }

	  if(bOk > 0)
		  m_step.x = localMax.x - m_center.x;

	  // reset bOk
	  bOk = false;

	  // stepY
	  x[0] = m_center.x - deltaX;      y[0] = m_center.y;
	  x[1] = m_center.x + deltaX;      y[1] = m_center.y + 4*deltaY;

	  while(!bOk && y[0] < m_image.height() - 4*deltaY)
	  {
		y[0] += deltaX;
		y[1] = y[0] + 4*deltaX;

		maxVal = FindLocalMaxPixelPos1(x[0], y[0], x[1], y[1], localMax.x, localMax.y);

		if( maxVal < m_CENTROID_THRESHOLD
			||(int)localMax.y - y[0] < deltaY
			||(y[1] - (int)localMax.y < deltaY) )
		{
		  bOk = false;
		}
		else
		  bOk = true;
	  }

	  if(bOk > 0)
	  {
		  m_step.y = localMax.y - m_center.y;
	  }

	  return bOk;
	}

	bool 
	HSCentroid::FindGridPicthNew()
	{
		bool bOk = true;

		int x[2] = {(int)( 0.4*m_image.width() ),  (int)( 0.6*m_image.width() ) };
		int y[2] = {(int)( 0.4*m_image.height() ), (int)( 0.6*m_image.height() ) };

		CImg<float> sampleImage = m_image.get_crop(x[0], y[0], x[1], y[1]);

		// vector saving pixel values along y axis
		CImg<float> VSum(x[1]-x[0]+1, 1);
		VSum.fill(0);

		// vector saving pixel values along x axis
		CImg<float> HSum(y[1]-y[0]+1, 1);
		HSum.fill(0);

		// vector saving local maximum flags
		CImg<unsigned char> VFlag(x[1]-x[0]+1, 1);
		CImg<unsigned char> HFlag(y[1]-y[0]+1, 1);
		VFlag.fill(0);
		HFlag.fill(0);

		cimg_forX(sampleImage, x)
		{
			VSum(x, 0) = sampleImage.get_crop(x, 0, x, sampleImage.height()-1).sum();
		}

		cimg_forY(sampleImage, y)
		{
			HSum(y, 0) = sampleImage.get_crop(0, y, sampleImage.width()-1, y).sum();
		}

		int radius = 10; 
		int step = 5;
		int localMaxX = 0;
		int localMaxY = 0;
		int localMaxCount = 0;

		m_step.x = 0;
		float startPos = 0;
		float endPos = 0;

		// find local maximum of intensity sum
		for(int i=radius; i<x[1]-x[0]-radius; i += step)
		{
			float maxValue = FindLocalMaxPixelPos1<float>(VSum, i-radius, 0, i+radius, 0, localMaxX, localMaxY);
			if(abs(i - localMaxX) < step)
			{
				VFlag(localMaxX, 0) = 1;
				endPos = localMaxX;
				++localMaxCount;
				i += step;

				if(localMaxCount == 1)
				{
					startPos = localMaxX;
				}
			}
		}

		assert( localMaxCount > 1 );

		m_step.x = (endPos - startPos)/(localMaxCount-1);

		m_step.y = 0;
		startPos = 0;
		endPos = 0;
		localMaxCount = 0;

		// find local maximum of intensity sum
		for(int i=radius; i<y[1]-y[0]-radius; i += step)
		{
			float maxValue = FindLocalMaxPixelPos1(HSum, i-radius, 0, i+radius, 0, localMaxX, localMaxY);
			if(abs(i - localMaxX) < step)
			{
				HFlag(localMaxX, 0) = 1;
				endPos = localMaxX;
				++localMaxCount;
				i += step;

				if(localMaxCount == 1)
				{
					startPos = localMaxX;
				}
			}
		}

		assert( localMaxCount > 1 );

		m_step.y = (endPos - startPos)/(localMaxCount-1);

		return bOk;
	}

	void HSCentroid::RadiantSearch()
	{
		SearchInDirection(m_center.x, m_center.y, m_startIndex.x, m_startIndex.y, 0,  1); // up
		SearchInDirection(m_center.x, m_center.y, m_startIndex.x, m_startIndex.y, 0, -1); // down

		for( int i = 0; i < m_resultMatrix.height(); i++)
		{
			if( m_resultMatrix( m_startIndexX, i ) > 0 )
			{
				int cx = (int)( m_resultMatrix( m_startIndex.x, i, 0, 1 ) + 0.5 );
				int cy = (int)( m_resultMatrix( m_startIndex.x, i, 0, 2 ) + 0.5 );

				this->SearchInDirection( cx, cy, m_startIndex.x, i,  1, 0 ); // to the right
				this->SearchInDirection( cx, cy, m_startIndex.x, i, -1, 0 ); // to the left
			}
		}
	}

	void  HSCentroid::SearchInDirection(const UInt startX, const UInt startY, const UInt startIndexX, const UInt startIndexY, const int dirX, const int dirY)
	{
		int Minimum_Grid = 10;
		int dx = cimg::sign(dirX);
		int dy = cimg::sign(dirY);

		if( (dirX == 0 && dirY == 0) || 
			(abs(dirX) > 1) || (abs(dirY) > 1) ||
			(m_step.x < m_minGridSize) || (m_step.y < m_minGridSize)  )
		{
			return;
		}

		int stepX = dirX*m_step.x;
		int stepY = dirY*m_step.y;

		int centerX = startX + stepX;
		int centerY = startY + stepY;

		int resultIndexX = startIndexX + dirX;
		int resultIndexY = startIndexY + dirY;

		int localMaxX, localMaxY;

		int max_Off_Distance = (int)(__min(m_step.x, m_step.y) / 5.0 );

		while(centerX > m_LSRadius && centerX < m_image.width()-m_LSRadius &&
			  centerY > m_LSRadius && centerY < m_image.height()-m_LSRadius)
		{
			if( resultIndexX >= 0 && resultIndexY >= 0 &&
				resultIndexX < m_resultMatrix.width() &&
				resultIndexY < m_resultMatrix.height() )
			{
				float maxVal = FindLocalMaxPixelPos1( centerX-m_LSRadius, centerY-m_LSRadius, 
													  centerX+m_LSRadius, centerY+m_LSRadius,
													  localMaxX, localMaxY );
				if( maxVal > m_CENTROID_THRESHOLD && 
				    abs((int)localMaxX - centerX) <= max_Off_Distance &&
				    abs((int)localMaxY - centerY) <= max_Off_Distance   )
				{
					int x = localMaxX;
					int y = localMaxY;

					float moduledValue = 1000*((resultIndexX+resultIndexY)%2);
					m_resultMatrix( resultIndexX, resultIndexY, 0, 0 ) = 3000 + moduledValue;
					m_resultMatrix( resultIndexX, resultIndexY, 0, 1 ) = localMaxX;
					m_resultMatrix( resultIndexX, resultIndexY, 0, 2 ) = localMaxY;

					//centerX = (int)(0.5*(centerX + localMaxX)) + stepX;
					//centerY = (int)(0.5*(centerY + localMaxY)) + stepY;

					centerX = localMaxX + stepX;
					centerY = localMaxY + stepY;
				  
				}
				else
				{
					m_resultMatrix( resultIndexX, resultIndexY, 0, 0 ) = INVALID_CENTROID;
					m_resultMatrix( resultIndexX, resultIndexY, 0, 1 ) = centerX;
					m_resultMatrix( resultIndexX, resultIndexY, 0, 2 ) = centerY;

					centerX += stepX;
					centerY += stepY;
				}

				// update centroid index in result matrix along searching direction
				resultIndexX += dirX;
				resultIndexY += dirY;
			}
			else
			{
				break;
			}
	        
		}
	}

}}


/*
void HSCentroid::Clear()
{
  m_stepX = 0;
  m_stepY = 0;
  m_LSRadius = 5;
  m_minGridSize = 10;
  m_CENTROID_THRESHOLD = 8.0f;

  m_image.assign();
  m_imgInput.assign();
  m_resultImage.assign();

  m_resultMatrix.assign(45, 45, 1, 3);
  m_resultMatrix.fill(0);

  m_referenceMatrix.assign();
}

template<typename T>
void HSCentroid::SetStreamData( T* data, int width, int height )
{
	m_imgInput.assign(data, width, height, 1, 1);
	m_image = m_imgInput.get_channel(0);
}

	
bool HSCentroid::OnExcute()
{
	bool bOk = true;
	
	SmoothAndDenoise();

	m_CENTROID_THRESHOLD = 2*m_image.mean();

	// color image
	m_resultImage = CImg<float>( m_image.width(), m_image.height(), 1, 3 );
	
	m_resultImage.fill(0);
	m_resultMatrix.fill(0);

	bOk = FindGridPicth();
	bOk = LocateMissingCentroid();

	RadiantSearch();

	CentroidValidation();

	CalcSubpixel();

	bOk &= ReferenceAssociation();

	ConvertResultToMM();

  
  //int found = filename.find_last_of(".");
  //std::string outputName = filename.substr(0, found);
  //outputName = outputName + ".txt";
  //OutputResult( outputName );

  return true;
}

void HSCentroid::SmoothAndDenoise()
{
	float diffsion = 2.0;
	float threshold = 100;
	int erodeVal = 10;
	int dilateVal = 10;

	CImg<float> res = m_image.get_blur(diffsion);
	res.threshold(threshold, false, false);

	res.erode(erodeVal);
	res.dilate(dilateVal);

	cimg_forXY(m_image, x,y)
	{ 
	  if(res(x,y))
		  m_image(x,y) = 0;
	}

	m_image.blur(diffsion);
}


bool 
HSCentroid::OpenImageFromFile(std::string filename)
{
	OpenReferenceCentroid();
	m_resultMatrix.fill(0);

	CImg<float> tempImg = CImg<float>(filename.c_str());
	this->SetStreamData(tempImg.data(), tempImg.width(), tempImg.height());

	this->SmoothAndDenoise();

	this->OnExcute();

	DrawCentroidDisplacement();

	size_t found = filename.find_last_of(".");
	std::string outputName = filename.substr(0, found);
	outputName = outputName + ".txt";

	OutputResult( outputName );

	return true;
}

//////////////////////////////////////////////////////////////////////////////////
// this function does two things: find the grid picthes of centroid spots 
// in x and y direction, and figure out the result matrix dimmension
//////////////////////////////////////////////////////////////////////////////////
bool 
HSCentroid::FindGridPicth()
{
  bool bOk = true;
  
  int x[2] = {(int)( 0.4*m_image.width() ),  (int)( 0.6*m_image.width() ) };
  int y[2] = {(int)( 0.6*m_image.height() ), (int)( 0.8*m_image.height() ) };

  int count = 0;
  UInt localMaxX, localMaxY;

  float maxVal = FindLocalMaxPixelPos1(x[0], y[0], x[1], y[1], m_centerX, m_centerY);

  // stepX
  int deltaX=5, deltaY = 5;
  x[0] = m_centerX;           y[0] = m_centerY - deltaY;
  x[1] = x[0] + 4*deltaX;     y[1] = m_centerY + deltaY;

  bOk = false;
  while( !bOk && x[0] < this->m_image.width() - 4*(int)deltaX)
  {
    x[0] += deltaX;
    x[1] = x[0] + 4*deltaX;

    maxVal = FindLocalMaxPixelPos(x[0], y[0], x[1], y[1], localMaxX, localMaxY);

    if( (maxVal < m_CENTROID_THRESHOLD) || ((int)localMaxX - x[0] < deltaX) || (x[1] - (int)localMaxX < deltaX) )
    {
      bOk = false;
    }
    else
      bOk = true;
  }

  if(bOk > 0)
    m_stepX = localMaxX - m_centerX;

  // reset bOk
  bOk = false;

  // stepY
  x[0] = m_centerX - deltaX;      y[0] = m_centerY;
  x[1] = m_centerX + deltaX;      y[1] = m_centerY + 4*deltaY;

  while(!bOk && y[0] < m_image.height() - 4*deltaY)
  {
    y[0] += deltaX;
    y[1] = y[0] + 4*deltaX;

    maxVal = FindLocalMaxPixelPos(x[0], y[0], x[1], y[1], localMaxX, localMaxY);

    if( (maxVal < m_CENTROID_THRESHOLD) ||((int)localMaxY - y[0] < deltaY) || (y[1] - (int)localMaxY < deltaY) )
    {
      bOk = false;
    }
    else
      bOk = true;
  }

  if(bOk > 0)
  {
    m_stepY = localMaxY - m_centerY;
  }

  return bOk;
}

float HSCentroid::FindLocalMaxPixelPos(UInt x0, UInt y0, UInt x1, UINT y1, UInt& maxX, UInt& maxY )
{
  CImg<float> cropImg = m_image.get_crop(x0, y0, x1, y1);
  float max_value = -9999999;
  long count = 0;
  long maxIndex = 0;

  cimg_for(cropImg,ptrs,float) {
    count++;
    const float val = *ptrs;
    if (val>max_value) {max_value = val; maxIndex = count; }
  }

  long reversedNum = cropImg.size() - maxIndex;
  maxX = x0 + reversedNum % cropImg.width();
  maxY = y0 + (int)(reversedNum / cropImg.width());

  return max_value;
}

float HSCentroid::FindLocalMaxPixelPos1(UInt x0, UInt y0, UInt x1, UINT y1, UInt& maxX, UInt& maxY )
{
  CImg<float> cropImg = m_image.get_crop(x0, y0, x1, y1);
  float *ptr = cropImg._data;
  float max_value = -9999999;
  long count = 0;
  long maxIndex = 0;
  long size = cropImg.size();

  for(int i = 0; i<size;i++)
  {
    const float val = *ptr++;
    if (val>max_value) {max_value = val; maxIndex = count; }
    count++;
  }

  maxX = x0 + maxIndex % cropImg.width();
  maxY = y0 + (int)(maxIndex / cropImg.width());

  return max_value;
}


float HSCentroid::FindLocalMaxPixelPos2(UInt& maxX, UInt& maxY )
{
  float max_value = -9999999;
  long count = 0;
  long maxIndex = 0;

  cimg_for(m_image,ptrs,float) {
    count++;
    const float val = *ptrs;
    if (val>max_value) {max_value = val; maxIndex = count; }
  }

  long reversedNum = m_image.size() - maxIndex;
  maxX = reversedNum % m_image.width();
  maxY = (int)(reversedNum / m_image.width());

  return max_value;
}

void  HSCentroid::SearchInDirection(const UInt startX, const UInt startY, const UInt startIndexX, const UInt startIndexY, const int dirX, const int dirY)
{
    int Minimum_Grid = 10;
    int dx = cimg::sign(dirX);
    int dy = cimg::sign(dirY);

    if( (dirX == 0 && dirY == 0) || 
		(abs(dirX) > 1) || (abs(dirY) > 1) ||
		(m_stepX < m_minGridSize) || (m_stepY < m_minGridSize)  )
    {
        return;
    }

    int stepX = dirX*m_stepX;
    int stepY = dirY*m_stepY;

    int centerX = startX + stepX;
    int centerY = startY + stepY;

    int resultIndexX = startIndexX + dirX;
    int resultIndexY = startIndexY + dirY;

    UInt localMaxX, localMaxY;

	int max_Off_Distance = (int)(__min(m_stepX, m_stepY) / 5.0 );

    while(centerX > m_LSRadius && centerX < m_image.width()-m_LSRadius &&
          centerY > m_LSRadius && centerY < m_image.height()-m_LSRadius)
    {
		if( resultIndexX >= 0 && resultIndexY >= 0 &&
			resultIndexX < m_resultMatrix.width() &&
			resultIndexY < m_resultMatrix.height() )
        {
			float maxVal = FindLocalMaxPixelPos1( centerX-m_LSRadius, centerY-m_LSRadius, 
                                              centerX+m_LSRadius, centerY+m_LSRadius,
                                              localMaxX, localMaxY );
			if( maxVal > m_CENTROID_THRESHOLD && 
					  abs((int)localMaxX - centerX) <= max_Off_Distance &&
					  abs((int)localMaxY - centerY) <= max_Off_Distance   )
			{
			  int x = localMaxX;
			  int y = localMaxY;
	    
			  float moduledValue = 1000*((resultIndexX+resultIndexY)%2);
			  m_resultMatrix( resultIndexX, resultIndexY, 0, 0 ) = 3000 + moduledValue;
			  m_resultMatrix( resultIndexX, resultIndexY, 0, 1 ) = localMaxX;
			  m_resultMatrix( resultIndexX, resultIndexY, 0, 2 ) = localMaxY;
	      
			  //centerX = (int)(0.5*(centerX + localMaxX)) + stepX;
			  //centerY = (int)(0.5*(centerY + localMaxY)) + stepY;

			  centerX = localMaxX + stepX;
			  centerY = localMaxY + stepY;
			  
			}
			else
			{
				m_resultMatrix( resultIndexX, resultIndexY, 0, 0 ) = INVALID_CENTROID;
				m_resultMatrix( resultIndexX, resultIndexY, 0, 1 ) = centerX;
				m_resultMatrix( resultIndexX, resultIndexY, 0, 2 ) = centerY;

				centerX += stepX;
				centerY += stepY;
			}

			// update centroid index in result matrix along searching direction
			resultIndexX += dirX;
			resultIndexY += dirY;
		}
		else
		{
			break;
		}
        
    }
}

void HSCentroid::RadiantSearch()
{
  	SearchInDirection(m_centerX, m_centerY, m_startIndexX, m_startIndexY, 0,  1); // up
    SearchInDirection(m_centerX, m_centerY, m_startIndexX, m_startIndexY, 0, -1); // down

	for( int i = 0; i < m_resultMatrix.height(); i++)
	{
		if( m_resultMatrix( m_startIndexX, i ) > 0 )
		{
			int cx = (int)( m_resultMatrix( m_startIndexX, i, 0, 1 ) + 0.5 );
			int cy = (int)( m_resultMatrix( m_startIndexX, i, 0, 2 ) + 0.5 );

			this->SearchInDirection( cx, cy, m_startIndexX, i,  1, 0 ); // to the right
			this->SearchInDirection( cx, cy, m_startIndexX, i, -1, 0 ); // to the left
		}
	}
}

bool HSCentroid::LocateMissingCentroid()
{
  bool bOk = false;
  int stationaryX = 0,
      stationaryY = 0;
  
  int imgCenterX = 0.5 * m_image.width();
  int imgCenterY = 0.5 * m_image.height();

  int x[2] = { imgCenterX - 4*m_stepX, imgCenterX - 2*m_stepX };
  int y[2] = { imgCenterY - 4*m_stepY, imgCenterY - 2*m_stepY };

  float stationaryVal = HS_STATIONARY_CENTROID;

  if(m_referenceMatrix.width() && m_referenceMatrix.height() )
  {
    for(int i=0; !bOk && i < m_referenceMatrix.height(); i++)   
	{
	    for(int j=0; !bOk && j<this->m_referenceMatrix.width(); j++)  
		{
		    float rVal = m_referenceMatrix(j, i, 0, 0);
			if( rVal == stationaryVal ) 
			{ //HS_STATIONARY_CENTROID  
				stationaryX = j;
			    stationaryY = i;
			    bOk = true;
		    }
	    }
    }

	if(bOk)
	{
		if( m_referenceMatrix( stationaryX+1, stationaryY+1, 0 ) == stationaryVal )
		{
			x[0] = (int)m_referenceMatrix(stationaryX, stationaryY, 1);
			y[0] = (int)m_referenceMatrix(stationaryX, stationaryY, 2);
			x[1] = (int)m_referenceMatrix(stationaryX+1, stationaryY+1, 1);
			y[1] = (int)m_referenceMatrix(stationaryX+1, stationaryY+1, 2);

			int locDx = (int)(0.5*(x[1] - x[0]));
			int locDy = (int)(0.5*(x[1] - x[0]));

			x[0] += locDx;
			y[0] += locDy;

			x[1] += locDx;
			y[1] += locDy;

			UInt localMaxX, localMaxY;
			float local_max = FindLocalMaxPixelPos1( x[0], y[0], x[1], y[1], localMaxX, localMaxY );

			if( local_max > m_CENTROID_THRESHOLD )
			{
			  this->m_centerX = localMaxX;
			  this->m_centerY = localMaxY;

			  // calc the index of centroid @ (m_centerX, m_centerY) in m_resultMatrix
			  m_startIndexX = stationaryX + 1;
			  m_startIndexY = stationaryY + 1;

			  float moduledValue = 1000*( ( m_startIndexX+m_startIndexY ) % 2 );

			  m_resultMatrix( stationaryX+1, stationaryY+1, 0, 0 ) = 3000 + moduledValue;
			  m_resultMatrix( stationaryX+1, stationaryY+1, 0, 1 ) = localMaxX;
			  m_resultMatrix( stationaryX+1, stationaryY+1, 0, 2 ) = localMaxY;

			  bOk = true;
			}
			else
			{
			  bOk = false;
			}

      }
    }
  }


  if(bOk)
  {
    return bOk;
  }


  UInt localMaxX, localMaxY;
  float REGION_MAX = FindLocalMaxPixelPos1( x[0], y[0], x[1], y[1], localMaxX, localMaxY );

  float startX = localMaxX;
  float startY = localMaxY;
  
  float prevX, prevY;
  UINT  closestToCenterX = localMaxX,
        closestToCenterY = localMaxY;

  float minMaxDist = abs( (int)closestToCenterX - imgCenterX ) 
                    + abs( (int)closestToCenterY - imgCenterY );

  for( y[0] = startY - m_LSRadius; !bOk && y[0] < imgCenterY + 2*m_stepY; y[0] += m_stepY )
  {
    y[1] = y[0] + 2*m_LSRadius;

    for( x[0] = startX - m_LSRadius; x[0] < imgCenterX + 2*m_stepX; x[0] += m_stepX )
    {
      x[1] = x[0] + 2*m_LSRadius;
      float local_max = FindLocalMaxPixelPos1( x[0], y[0], x[1], y[1], localMaxX, localMaxY );
	  
	    float curMaxDist = abs((int)localMaxX - imgCenterX) + abs((int)localMaxY - imgCenterY);
	    minMaxDist = __max(minMaxDist, curMaxDist);

      if( local_max < 0.1*REGION_MAX )
      {
        bOk = true;
        m_centerX = prevX + m_stepX;
        m_centerY = prevY;
        break;
      }
      else
      {
        prevX = localMaxX;
        prevY = localMaxY;

		    if( curMaxDist < minMaxDist )
        {
          closestToCenterX = localMaxX;
          closestToCenterY = localMaxY;
        }
      }

    }
  }

  if(bOk)
  {
    UINT locMaxX[4], locMaxY[4];
    float 
    local_max = FindLocalMaxPixelPos1( m_centerX-m_stepX-m_LSRadius, m_centerY-m_LSRadius, 
                                       m_centerX-m_stepX+m_LSRadius, m_centerY+m_LSRadius, locMaxX[0], locMaxY[0] );

    local_max = FindLocalMaxPixelPos1( m_centerX+m_stepX-m_LSRadius, m_centerY-m_LSRadius, 
                                       m_centerX+m_stepX+m_LSRadius, m_centerY+m_LSRadius, locMaxX[1], locMaxY[1] );

    local_max = FindLocalMaxPixelPos1( m_centerX-m_LSRadius, m_centerY-m_stepY-m_LSRadius, 
                                       m_centerX+m_LSRadius, m_centerY-m_stepY+m_LSRadius, locMaxX[2], locMaxY[2] );

    local_max = FindLocalMaxPixelPos1( m_centerX-m_LSRadius, m_centerY+m_stepY-m_LSRadius,
                                       m_centerX+m_LSRadius, m_centerY+m_stepY+m_LSRadius, locMaxX[3], locMaxY[3] );

    m_centerX = UINT( 0.25*( locMaxX[0] + locMaxX[1] + locMaxX[2] + locMaxX[3] ) + 0.5 );
    m_centerY = UINT( 0.25*( locMaxY[0] + locMaxY[1] + locMaxY[2] + locMaxY[3] ) + 0.5 );
  }
  else
  {
    m_centerX = closestToCenterX;
    m_centerY = closestToCenterY;
  }


  // calc the index of centroid @ (m_centerX, m_centerY) in m_resultMatrix
  m_startIndexX = cimg::max<UInt, int>( m_centerX, m_image.width()-m_centerX ) / m_stepX + 1;
  m_startIndexY = cimg::max<UInt, int>( m_centerY, m_image.height()-m_centerY ) / m_stepY + 1;

  m_resultMatrix = CImg<float>(2*m_startIndexX +1, 2*m_startIndexY +1, 1, 3);
  m_resultMatrix.fill(0);

  m_resultMatrix( m_startIndexX, m_startIndexY ) = 1;
  m_resultMatrix( m_startIndexX, m_startIndexY, 1 ) = m_centerX;
  m_resultMatrix( m_startIndexX, m_startIndexY, 2 ) = m_centerY;

  for(int i = -2; i <= 2; i++)
  {
    m_resultImage( m_centerX+i, m_centerY, 0 ) = m_resultImage( m_centerX, m_centerY+i, 0 ) = CENTROID_FLAG;
    m_resultImage( m_centerX+i, m_centerY, 1 ) = m_resultImage( m_centerX, m_centerY+i, 1 ) = CENTROID_FLAG;
    m_resultImage( m_centerX+i, m_centerY, 2 ) = m_resultImage( m_centerX, m_centerY+i, 2 ) = 55;
  }

  if(!bOk)  // no artifact (missing centroid)
  {
    float moduledValue = 1000*( ( m_startIndexX+m_startIndexY )%2 );
    m_resultMatrix( m_startIndexX, m_startIndexY, 0, 0 ) = 3000 + moduledValue;
    m_resultMatrix( m_startIndexX, m_startIndexY, 0, 1 ) = m_centerX;
    m_resultMatrix( m_startIndexX, m_startIndexY, 0, 2 ) = m_centerY;

    for(int i = -2; i <= 2; i++)
    {
      m_resultImage( m_centerX+i, m_centerY, 0 ) = m_resultImage( m_centerX, m_centerY+i, 0 ) = CENTROID_FLAG;
      m_resultImage( m_centerX+i, m_centerY, 1 ) = m_resultImage( m_centerX, m_centerY+i, 1 ) = m_centerX;
      m_resultImage( m_centerX+i, m_centerY, 2 ) = m_resultImage( m_centerX, m_centerY+i, 2 ) = m_centerY;
    }
  }

  return bOk;
}
void HSCentroid::CentroidValidation()
{
	float centroidVal = CENTROID_FLAG;
	for( int y = 0; y < m_resultMatrix.height(); y++ )
	{
		for ( int x = 0; x < m_resultMatrix.width(); x++ )
		{
			if( m_resultMatrix( x, y, 0) >= centroidVal )
			{
				int i = (int)( m_resultMatrix( x, y, 1) + 0.5 );
				int j = (int)( m_resultMatrix( x, y, 2) + 0.5 );

				if( m_image(i, j) <= m_image(i-1, j-1) ||
					m_image(i, j) <= m_image(i  , j-1) ||
					m_image(i, j) <= m_image(i+1, j-1) ||
					m_image(i, j) <= m_image(i-1, j  ) ||
					m_image(i, j) <= m_image(i+1, j  ) ||
					m_image(i, j) <= m_image(i-1, j+1) ||
					m_image(i, j) <= m_image(i  , j+1) ||
					m_image(i, j) <= m_image(i+1, j+1) 
				  )
				{
					m_resultMatrix( x, y, 0) = INVALID_CENTROID;
				}
			}
		}
	}
}
void HSCentroid::CalcSubpixel()
{

	for(int y = 0; y < m_resultMatrix.height(); y++)
	{
		for(int x = 0; x < m_resultMatrix.width(); x++)
		{
			float Valid_Centroid = CENTROID_FLAG;

			if( m_resultMatrix(x, y, 0, 0) >= Valid_Centroid )
			{
				int lx = (int)(m_resultMatrix( x, y, 0, 1 ) + 0.5);
				int ly = (int)(m_resultMatrix( x, y, 0, 2 ) + 0.5);

				if( lx > 1 && ly > 1 &&
					lx < m_image.width() -1 &&
					ly < m_image.height() -1 
				)
				{
					float denomX =  2*m_image( lx, ly ) - ( m_image( lx+1, ly ) + m_image( lx-1, ly ) );
					float denomY =  2*m_image( lx, ly ) - ( m_image( lx, ly+1 ) + m_image( lx, ly-1 ) );
					float denom45 =  2*m_image( lx, ly ) - ( m_image( lx+1, ly+1 ) + m_image( lx-1, ly-1 ) );
					float denom135 =  2*m_image( lx, ly ) - ( m_image( lx-1, ly+1 ) + m_image( lx+1, ly-1 ) );

					float dx0 = m_image( lx+1, ly ) - m_image( lx-1, ly );
					float dx45 = m_image( lx+1, ly+1 ) - m_image( lx-1, ly-1 );
					float dx135 = m_image( lx+1, ly-1 ) - m_image( lx-1, ly+1 );

					float dy0 = m_image( lx, ly+1) - m_image( lx, ly-1 );
					float dy45 = m_image( lx+1, ly+1 ) - m_image( lx-1, ly-1 );
					float dy135 = m_image( lx-1, ly+1 ) - m_image( lx+1, ly-1 );

					float subpixelX = lx;
					float subpixelY = ly;

					if( denomX > 0)
						subpixelX += 0.5*dx0/denomX;

					if( denomY > 0)
						subpixelY += 0.5*dy0/denomY;

					if(denom45 > 0)
					{
						subpixelX += 0.25 * dx45 / denom45;
						subpixelY += 0.25 * dy45 / denom45;
					}

					if(denom135 > 0)
					{
						subpixelX += 0.25 * dx135 / denom135;
						subpixelY += 0.25 * dy135 / denom135;
					}

					m_resultMatrix(x, y, 0, 1) = subpixelX;
					m_resultMatrix(x, y, 0, 2) = subpixelY;
				}
			}
		}
	}

}

void HSCentroid::ResultDisplay()
{
	if(this->m_imgInput.data() == NULL)
		return;

    m_resultImage.display("Result Image", true);
    m_resultMatrix.display("Result Matrix", true);

	return;

	std::string locStr("");
	char strX[10], strY[10];
  
  CImgDisplay disp(m_resultImage.width(), m_resultImage.height(), "Result Image");//disp(m_resultImage.width(), m_resultImage.height(), "My Display");
  while(!disp.is_closed()){
    if(disp.button()&1) {
      int x = disp.mouse_x(), y = disp.mouse_y();
      unsigned char foreground_color[] = {0, 0, 255};
      unsigned char background_color[] = {200, 200, 0};
      //m_imgInput.draw_text(100, 500, "+", purple);
      float ratioX = (float)m_imgInput.width()/disp.width();
      float ratioY = (float)m_imgInput.height()/disp.height();
      x = (int)(x*ratioX+0.5); 
      y = (int)(y*ratioY+0.5); 
	  
	  int fontSize = 24;(int)(13*__max(ratioX, ratioY));

	  locStr.clear();
	  sprintf(strX, "%.3f", (float)m_resultImage(x, y, 0, 1));
	  sprintf(strY, "%.3f", (float)m_resultImage(x, y, 0, 2));
	  locStr = locStr + "(" + strX + ", " + strY + ")";

	  CImg<float> img(m_imgInput);

	  if( x > 0.5*disp.width())
		x = x - locStr.size();
	  else
		x = x+5;
	  if( y < 0.5*disp.height() )
		y = y + 10;
	  else
		y = y-25;

	  img.draw_text(x,y,locStr.c_str(),foreground_color,background_color,1,fontSize);
      img.display(disp);
    }
    if( disp.is_resized()) disp.resize();
    disp.wait(25);
  }
  //delete[] str;

  return;
}

void HSCentroid::UpdateResultImage()
{
	for(int y = 0; y < m_resultMatrix.height(); y++)
	{
		for(int x = 0; x < m_resultMatrix.width(); x++)
		{
			float Valid_Centroid = CENTROID_FLAG;
			if(m_resultMatrix(x, y, 0, 0) >= Valid_Centroid )
			{
				int lx = (int)(m_resultMatrix(x, y, 0, 1) + 0.5);
				int ly = (int)(m_resultMatrix(x, y, 0, 2) + 0.5);

				if(lx > 3 && lx < m_resultImage.width() -3 &&
					ly > 3 && ly < m_resultImage.height() -3)
				{
					for(int i = -3; i <= 3; i++)
						for(int j = -3; j <=3; j++)
						{
							m_resultImage(lx+j, ly+i, 0, 0) = CENTROID_REGION;
							m_resultImage(lx+j, ly+i, 0, 1) = m_resultMatrix(x, y, 0, 1);
							m_resultImage(lx+j, ly+i, 0, 2) = m_resultMatrix(x, y, 0, 2);
						}
				}
			}
		}
	}
}
void HSCentroid::OpenReferenceCentroid()
{
	FILE* pFile = fopen("ReferenceCentroidData.txt", "r");
	
	if(!pFile) 
	{
		m_referenceMatrix.assign();
		return;
	}

	int width = 45, height = 45;
	m_referenceMatrix.assign(width, height, 1, 3);
	m_referenceMatrix.display("Reference Centroid Info", true);
	m_referenceMatrix.fill(0);

	float r, g, b;
	int row, col;

	int maxWidth = 0, maxHeight = 0;
	int size = 0;

	do
	{
		size = fscanf (pFile, "%f %f %f %d %d", &r, &g, &b, &row, &col);

		if(size == 5) 
		{
			if( col >= 0 && col < width && 
			  row >= 0 && row < height ) 
			{
				m_referenceMatrix(col, row, 0) = r;
				m_referenceMatrix(col, row, 1) = g;
				m_referenceMatrix(col, row, 2) = b;
			}

			maxWidth = __max(maxWidth, col);
			maxHeight = __max(maxHeight, row);

		}
		else
		{
			break;
		}

	}while(size == 5);

	int imageW = __min(width, maxWidth);
	int imageH = __min(height, maxHeight);

	if(imageW > 0 && imageH > 0)
	{
		m_referenceMatrix = m_referenceMatrix.get_crop(0, 0, imageW, imageH);
		m_referenceMatrix.display("Reference Centroid Info", true);
	}
	else
	{
		m_referenceMatrix.assign();
	}


  bool bFound = false;
  float stationaryVal = HS_STATIONARY_CENTROID;

  if(m_referenceMatrix.width() && m_referenceMatrix.height() )
  {
    for(int i=0; !bFound && i < m_referenceMatrix.height(); i++)   {
	    for(int j=0; !bFound && j<this->m_referenceMatrix.width(); j++)  {

		    float rVal = m_referenceMatrix(j, i, 0, 0);
        if( rVal == stationaryVal ) { //HS_STATIONARY_CENTROID  
			    m_stationaryX = j;
			    m_stationaryY = i;
			    bFound = true;
		    }
	    }
    }

	  float centroidVal = CENTROID_FLAG;

    // make sure the stationary centroids are 4x4 matrix
    if(bFound) {
      for(int i = 0; bFound && i<4; i++)  {
	      for(int j = 0; bFound && j<4; j++)  {
		      if(m_referenceMatrix(m_stationaryX+j, m_stationaryY+i, 0, 0) >= centroidVal )
			      if(m_referenceMatrix(m_stationaryX+j, m_stationaryY+i, 0, 0) != stationaryVal)
			      {
				      bFound = false;
			      }
	      }
      }
    } // if(bFound)
  }

	if(!bFound)
	{
		m_stationaryX = -1;
		m_stationaryY = -1;
	}

	fclose(pFile);
}

void HSCentroid::SaveAsReference()
{
  m_referenceMatrix.assign( m_resultMatrix );

	FILE* pFile = fopen("ReferenceCentroidData.txt", "w+");

	FILE* pRefDll = fopen("ReferenceCentroidForDll.txt", "w+");

	int stationary_x = 20;
	int stationary_y = 14;
	float CentroidVal = CENTROID_FLAG;
	float StationaryVal = HS_STATIONARY_CENTROID;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			m_resultMatrix(stationary_x+i, stationary_y+j, 0) = HS_STATIONARY_CENTROID;
		}
	}

	int width = m_resultMatrix.width();
	int height = m_resultMatrix.height();

	for(int j = 0; j < height; j++)
	{
		for(int i = 0; i < width; i++)
		{
			fprintf(pFile, "%6.3f %6.3f %6.3f %d %d\n",  
							m_resultMatrix(i, j, 0),
							m_resultMatrix(i, j, 1),
							m_resultMatrix(i, j, 2),
							j, i);
			
			if( m_resultMatrix(i, j, 0) >= CentroidVal )
			{
				int stationary_flag = 0;
				if( m_resultMatrix(i, j, 0) == StationaryVal )
				{
					stationary_flag = 1;
				}

				fprintf(	pRefDll, "%f,%f,%d,%d,%d\n",  
							m_resultMatrix(i, j, 1),
							m_resultMatrix(i, j, 2),
							i, 
							j,
							stationary_flag
						);
			}
		}
	}

	fclose(pFile);
	fclose(pRefDll);
}

void HSCentroid::LocateSearchCenter()
{
}

bool HSCentroid::CalcStepSizeAndStartIndex()
{
	bool bOk = true;
	bool bFound = false;
	float stationaryVal = HS_STATIONARY_CENTROID;
	int stationaryX = -1, 
		stationaryY = -1;

	if(m_referenceMatrix.width() && m_referenceMatrix.height() )
	{
		for(int i=0; !bFound && i < m_referenceMatrix.height(); i++)   
		{
			for(int j=0; !bFound && j<this->m_referenceMatrix.width(); j++)  
			{
				float rVal = m_referenceMatrix(j, i, 0, 0);
				if( rVal == stationaryVal ) 
				{ //HS_STATIONARY_CENTROID  
					stationaryX = j;
					stationaryY = i;
					bFound = true;
				}
			}
		}

		float centroidVal = CENTROID_FLAG;

		if(bFound) 
		{
			for(int i = 0; bOk && i<4; i++)  
			{
				for(int j = 0; bOk && j<4; j++)  
				{
					if(m_referenceMatrix(stationaryX+j, stationaryY+i, 0, 0) >= centroidVal )
						if(m_referenceMatrix(stationaryX+j, stationaryY+i, 0, 0) != stationaryVal)
						{
							bOk = false;
						}
				}
			}
		} // if(bFound)
	}

	// initialize x[] and y[]
	int x[2] = {(int)( 0.4*m_image.width() ),  (int)( 0.6*m_image.width() ) };
	int y[2] = {(int)( 0.6*m_image.height() ), (int)( 0.8*m_image.height() ) };

	if(bOk && bFound) // if found stationary centroids in reference matrix, use them
	{
		x[0] = (int)m_referenceMatrix(stationaryX, stationaryY, 1);
		y[0] = (int)m_referenceMatrix(stationaryX, stationaryY, 2);
		x[1] = (int)m_referenceMatrix(stationaryX+1, stationaryY+1, 1);
		y[1] = (int)m_referenceMatrix(stationaryX+1, stationaryY+1, 2);

		int locDx = (int)(0.5*(x[1] - x[0]));
		int locDy = (int)(0.5*(y[1] - y[0]));

		x[0] -= locDx;
		y[0] -= locDy;

		x[1] += locDx;
		y[1] += locDy;
	}
	
  int count = 0;
  UInt localMaxX, localMaxY;
  float subPixelX, subPixelY;

  float maxVal = FindLocalMaxPixelPos1(x[0], y[0], x[1], y[1], m_centerX, m_centerY);

  m_startIndexX = 20;
  m_startIndexY = 20;

  if(maxVal > m_CENTROID_THRESHOLD )
  {
	  m_resultMatrix( m_startIndexX, m_startIndexY, 0) = CENTROID_FLAG;
	  m_resultMatrix( m_startIndexX, m_startIndexY, 1) = m_centerX;
	  m_resultMatrix( m_startIndexX, m_startIndexY, 2) = m_centerY;
  }
  else
  {
	  m_resultMatrix( m_startIndexX, m_startIndexY, 0) = 1;
	  m_resultMatrix( m_startIndexX, m_startIndexY, 1) = m_centerX;
	  m_resultMatrix( m_startIndexX, m_startIndexY, 2) = m_centerY;
  }

  // stepX
  int deltaX=5, deltaY = 5;
  x[0] = m_centerX;           y[0] = m_centerY - deltaY;
  x[1] = x[0] + 4*deltaX;     y[1] = m_centerY + deltaY;

  bOk = false;
  while( !bOk && x[0] < this->m_image.width() - 4*(int)deltaX)
  {
    x[0] += deltaX;
    x[1] = x[0] + 4*deltaX;

    maxVal = FindLocalMaxPixelPos(x[0], y[0], x[1], y[1], localMaxX, localMaxY);

    if( (maxVal > m_CENTROID_THRESHOLD) &&
        ((int)localMaxX - x[0] > deltaX) && 
        (x[1] - (int)localMaxX > deltaX) )
    {
      bOk = true;
    }
  }

  if(bOk) {
    m_stepX = localMaxX - m_centerX;
  }

  if(bOk)
  {
    // reset bOk
    bOk = false;

    // stepY
    x[0] = m_centerX - deltaX;      y[0] = m_centerY;
    x[1] = m_centerX + deltaX;      y[1] = m_centerY + 4*deltaY;

    while(!bOk && y[0] < m_image.height() - 4*deltaY)
    {
      y[0] += deltaX;
      y[1] = y[0] + 4*deltaX;

      maxVal = FindLocalMaxPixelPos(x[0], y[0], x[1], y[1], localMaxX, localMaxY);

      if( (maxVal > m_CENTROID_THRESHOLD) &&
          ((int)localMaxY - y[0] > deltaY) && 
          (y[1] - (int)localMaxY > deltaY) )  {
        bOk = true;
      }
    }
  }

  if(bOk) {
    m_stepY = localMaxY - m_centerY;
  }

  if(bOk)
  {
	  this->m_LSRadius = __max( m_LSRadius, (int)(__min(m_stepY, m_stepX)/4.0 ) );
  }

  return bOk;
}

bool HSCentroid::ReferenceAssociation()
{
  bool bOk = true;

  if( !m_referenceMatrix.data() ||  !m_resultMatrix.data() )
  {
    bOk = false;
  }

  bool bFound = false;
  int Size = 4;
  int stationaryX = -1;
  int stationaryY = -1;
  float Centroid_Value = CENTROID_FLAG;

  if(m_referenceMatrix.width() && m_referenceMatrix.height() )
  {
    float stationaryVal = HS_STATIONARY_CENTROID;

    for(int i=0; !bFound && i < m_referenceMatrix.height(); i++)   {
	    for(int j=0; !bFound && j<this->m_referenceMatrix.width(); j++)  {

		    float rVal = m_referenceMatrix(j, i, 0, 0);
        if( rVal == stationaryVal ) { //HS_STATIONARY_CENTROID  
			    stationaryX = j;
			    stationaryY = i;
			    bFound = true;
		    }
	    }
    }

	double sumX = 0.0, sumY = 0.0;
    if(bFound) {
      for( int i = 0; bOk && i < Size; i++ )  {
	      for( int j = 0; bOk && j < Size; j++ )  {
		      if(m_referenceMatrix(stationaryX+j, stationaryY+i, 0, 0) >= Centroid_Value )
			  {
				  if(m_referenceMatrix(stationaryX+j, stationaryY+i, 0, 0) != stationaryVal)
			      {
				      bOk = false;
			      }
				  else
				  {
					  sumX += m_referenceMatrix(stationaryX+j, stationaryY+i, 0, 1);
					  sumY += m_referenceMatrix(stationaryX+j, stationaryY+i, 0, 2);
				  }
			  }
	      }
      }
    } // if(bFound)

	if( bOk )
	{
		m_discCenterX = sumX / ( Size * Size );
		m_discCenterY = sumY / ( Size * Size );
	}

  }


  int boarder = 3;
  boarder = __min( boarder, stationaryX);
  boarder = __min( boarder, m_referenceMatrix.width() - stationaryX );

  boarder = __min( boarder, stationaryY);
  boarder = __min( boarder, m_referenceMatrix.height() - stationaryY );

  int matchingX;
  int matchingY;
  float minDiff = 1.0e9;

  if(bFound)
  {
    CImg<float> refer4x4 = m_referenceMatrix.get_crop(stationaryX, stationaryY, stationaryX+Size -1, stationaryY+Size-1);

    for( int y = stationaryY - boarder; y < stationaryY + boarder; y++ )
    {
      for( int x = stationaryX -boarder; x < stationaryX + boarder; x++ )
      {
        float curDiff = 0.0;

        CImg<float> target4x4 = m_resultMatrix.get_crop( x, y, x + Size - 1, y + Size -1);
        CImg<float> diff4x4 = ( target4x4 - refer4x4 ).channels( 1, 2 );
        
        cimg_forXY(target4x4, i, j)
        {
          if( target4x4( i, j, 0 ) < Centroid_Value )
          {
            diff4x4( i, j, 0 ) = diff4x4( i, j, 1 ) = 0.0;
          }
        }
        curDiff = diff4x4.pow( 2 ).sum();

        if( curDiff < minDiff )
        {
          minDiff = curDiff;
          matchingX = x;
          matchingY = y;
        }
      }
     }
   }

  if(bOk)
  {
    m_referenceAssocX = matchingX;
    m_referenceAssocY = matchingY;
  }
  else
  {
    m_referenceAssocX = -1;
    m_referenceAssocY = -1;
  }

  return bOk;
}


void HSCentroid::DrawCentroidDisplacement()
{
	m_resultImage.assign();
	m_resultImage.assign( m_imgInput );


	int red[3] = {255, 100, 50};
	int yellow[3] = {255, 255, 50};
	int black[3] = {0, 0, 0};


	if( m_stationaryX < 0 || m_stationaryY < 0 )
	{
	  m_resultImage.draw_text(100, 300, "No Reference Data!", red, black, 1.0, 48 );
	  return;
	}

	int DX = m_stationaryX - m_referenceAssocX;
	int DY = m_stationaryY - m_referenceAssocY;

	int width = __min(m_referenceMatrix.width(), m_resultMatrix.width());
	int height = __min(m_referenceMatrix.height(), m_resultMatrix.height());

	for(int y = 0; y < height; y++)
	{
		for(int x = 0; x < width; x++)
		{
			if( m_resultMatrix( x, y, 0 ) >= 3000 && 
			m_referenceMatrix( x+DX, y+DY, 0) >= 102 )
			{
				if( x+DX >= 0 && x+DX < width &&
					y+DY >= 0 && y+DY < height )
				{
					int w[2] = {(int)m_referenceMatrix(x+DX, y+DY, 1), (int)m_resultMatrix(x, y, 1)};
					int h[2] = {(int)m_referenceMatrix(x+DX, y+DY, 2), (int)m_resultMatrix(x, y, 2)};

					m_resultImage.draw_circle(w[0], h[0], 2, yellow, 0.2);
					m_resultImage.draw_line(w[0], h[0], w[1], h[1], red, 0.5);
				}
			}
		}
	}

	m_resultImage.display("result image", true);
}


void HSCentroid::OutputResult( std::string outputName )
{
	FILE* pf = fopen(outputName.c_str(), "w");

	if( m_stationaryX < 0 || m_stationaryY < 0 ||
		m_discCenterX < 0 || m_discCenterY < 0 )
	{
		fclose(pf);
		return;
	}

	double Normal_Radius_mm = NORMAL_RADIUS_MM;
	double Pixel_Size_mm = PIXLE_SIZE_MM;
	double Focus_Length_mm = FOCUS_LENGTH_MM;

	double Normal_Radius = Normal_Radius_mm / Pixel_Size_mm;
	double Focus_Length = Focus_Length_mm / Pixel_Size_mm;

	double R2 = Normal_Radius * Normal_Radius;

	int DX = m_stationaryX - m_referenceAssocX;
	int DY = m_stationaryY - m_referenceAssocY;

	int width = __min(m_referenceMatrix.width(), m_resultMatrix.width());
	int height = __min(m_referenceMatrix.height(), m_resultMatrix.height());

	for(int y = 0; y < height; y++)
	{
		for(int x = 0; x < width; x++)
		{
			if( m_resultMatrix( x, y, 0 ) >= 3000 && 
			m_referenceMatrix( x+DX, y+DY, 0) >= 102 )
			{
				if( x+DX >= 0 && x+DX < width &&
					y+DY >= 0 && y+DY < height )
				{
					int w[2] = {(int)m_referenceMatrix(x+DX, y+DY, 1), (int)m_resultMatrix(x, y, 1)};
					int h[2] = {(int)m_referenceMatrix(x+DX, y+DY, 2), (int)m_resultMatrix(x, y, 2)};

					double distToCenterX = m_referenceMatrix(x+DX, y+DY, 1) - m_discCenterX;
					double distToCenterY = m_referenceMatrix(x+DX, y+DY, 2) - m_discCenterY;

					if( R2 >= (distToCenterX*distToCenterX + distToCenterY*distToCenterY) )
					{
						fprintf(pf,"%f,%f,%f,%f\n", 
									distToCenterX / Normal_Radius, 
									distToCenterY / Normal_Radius,
									(m_referenceMatrix(x+DX, y+DY, 1) - m_resultMatrix(x, y, 1))/Focus_Length, 
									(m_referenceMatrix(x+DX, y+DY, 2) - m_resultMatrix(x, y, 2))/Focus_Length
									);
					}
				}
			}
		}
	}

	fprintf(pf, "\n");
	fclose(pf);
}

void HSCentroid::GetCentroidResult( float* locX, float* locY, float* dx, float* dy, int& dataSize )
{
	dataSize = 0;
	dataSize = m_resultInMM.height();

	for( int i = 0; i < dataSize; i++)
	{
		// reference location
		locX[i] = m_resultInMM(0, i);
		locY[i] = m_resultInMM(1, i);

		// displacement
		dx[i] = m_resultInMM(2, i);
		dy[i] = m_resultInMM(3, i);
	}
}

template<typename T>
void  HSCentroid::SetAndCalcCentroid(
							T* data, int width, int height,		// image data
							float* pixelX, float* pixelY,			// reference data
							int* gridX, int* gridY, int* flag, int refSize,		// reference data
							float* locX, float* locY, float* dx, float* dy, int& dataSize	// output data size
						  )
{
	SetStreamData( data, width, height );
	SetReferenceData( pixelX, pixelY, gridX, gridY,  flag, refSize );
	this->SmoothAndDenoise();
	this->OnExcute();

	GetCentroidResult(locX, locY, dx, dy, dataSize );
	//DrawCentroidDisplacement();

}
  

void  HSCentroid::SetReferenceData( float* pixelX, float* pixelY, 
							 int* gridX, int* gridY, int* flag, int size )
{
	m_referenceMatrix.assign();

	int widthMin = 1;
	int heightMin = 1;
	int width = 0;
	int height = 0;
	for( int i = 0; i < size; i++ )
	{
		widthMin = __min( widthMin, gridX[i] );
		heightMin = __min( heightMin, gridY[i] );

		width = __max( width, gridX[i] );
		height = __max( height, gridY[i] );
	}

	if( widthMin >= 0 && heightMin >= 0 && width > 0 && height > 0)	
	{
		width += 1;
		height += 1;
		m_referenceMatrix.assign(width, height, 1, 3 );
		m_referenceMatrix.fill( 0 );
	}
	else{
		return;
	}

	// initialize upper-left stationary centroid location in grid
	m_stationaryX = width;
	m_stationaryY = height;

	int stationaryCount = 0;
	for( int i = 0; i < size; i++ )
	{
		m_referenceMatrix( gridX[i], gridY[i], 0 ) = CENTROID_FLAG;
		m_referenceMatrix( gridX[i], gridY[i], 1 ) = pixelX[i];
		m_referenceMatrix( gridX[i], gridY[i], 2 ) = pixelY[i];

		if( flag[i] )
		{
			stationaryCount++;

			// set stationary centroid
			m_referenceMatrix( gridX[i], gridY[i], 0 ) = HS_STATIONARY_CENTROID;

			// find upper-left centroid
			if( gridX[i] < m_stationaryX || 
				gridY[i] < m_stationaryY  )
			{
				m_stationaryX = gridX[i];
				m_stationaryY = gridY[i];
			}
		}
	}


	if( stationaryCount != 16 ||
		m_stationaryX < 1 || m_stationaryX >= width ||
		m_stationaryY < 1 || m_stationaryY >= height )
	{
		m_stationaryX = -1;
		m_stationaryY = -1;

		m_referenceMatrix.assign();

		return;
	}

	// check to see if this is an valid reference by checking 
	// the distribution of the 16 stationary centroids form 4x4 matrix
	bool bValidReference = true;
	float stationaryVal = HS_STATIONARY_CENTROID;

	for(int i = 0; bValidReference && i < 4; i++)  {
		for(int j = 0; bValidReference && j < 4; j++)  {
			if(m_referenceMatrix(m_stationaryX+j, m_stationaryY+i, 0, 0) != stationaryVal)
			{
				bValidReference = false;
			}
		}
	}

	if(!bValidReference)
	{
		m_stationaryX = -1;
		m_stationaryY = -1;

		m_referenceMatrix.assign();
	}

 }
void HSCentroid::ConvertResultToMM()
{
	m_resultInMM.assign();

	if( m_stationaryX < 0 || m_stationaryY < 0 ||
		m_discCenterX < 0 || m_discCenterY < 0 )
	{
		return;
	}

	int initSize = 900;
	m_resultInMM.assign(4, initSize, 1, 1);
	m_resultInMM.fill(0);


	double Normal_Radius_mm = NORMAL_RADIUS_MM;
	double Pixel_Size_mm = PIXLE_SIZE_MM;
	double Focus_Length_mm = FOCUS_LENGTH_MM;

	double Normal_Radius = Normal_Radius_mm / Pixel_Size_mm;
	double Focus_Length = Focus_Length_mm / Pixel_Size_mm;

	double R2 = Normal_Radius * Normal_Radius;

	int DX = m_stationaryX - m_referenceAssocX;
	int DY = m_stationaryY - m_referenceAssocY;

	int width = __min(m_referenceMatrix.width(), m_resultMatrix.width());
	int height = __min(m_referenceMatrix.height(), m_resultMatrix.height());

	int count = 0;
	for(int y = 0; y < height; y++)
	{
		for(int x = 0; x < width; x++)
		{
			if( m_resultMatrix( x, y, 0 ) >= 3000 && 
			m_referenceMatrix( x+DX, y+DY, 0) >= 102 )
			{
				if( x+DX >= 0 && x+DX < width &&
					y+DY >= 0 && y+DY < height )
				{
					int w[2] = {(int)m_referenceMatrix(x+DX, y+DY, 1), (int)m_resultMatrix(x, y, 1)};
					int h[2] = {(int)m_referenceMatrix(x+DX, y+DY, 2), (int)m_resultMatrix(x, y, 2)};

					double distToCenterX = m_referenceMatrix(x+DX, y+DY, 1) - m_discCenterX;
					double distToCenterY = m_referenceMatrix(x+DX, y+DY, 2) - m_discCenterY;

					if( R2 >= (distToCenterX*distToCenterX + distToCenterY*distToCenterY) )
					{
					
						//fprintf(pf,"%f,%f,%f,%f\n", 
						//			distToCenterX / Normal_Radius, 
						//			distToCenterY / Normal_Radius,
						//			(m_referenceMatrix(x+DX, y+DY, 1) - m_resultMatrix(x, y, 1))/Focus_Length, 

						//			(m_referenceMatrix(x+DX, y+DY, 2) - m_resultMatrix(x, y, 2))/Focus_Length
						//			);
						
						m_resultInMM(0, count) = distToCenterX / Normal_Radius;
						m_resultInMM(1, count) = distToCenterY / Normal_Radius;
						m_resultInMM(2, count) = (m_referenceMatrix(x+DX, y+DY, 1) - m_resultMatrix(x, y, 1))/Focus_Length;
						m_resultInMM(3, count) = (m_referenceMatrix(x+DX, y+DY, 2) - m_resultMatrix(x, y, 2))/Focus_Length;
						++count;
					}
				}
			}
		}
	}

	m_resultInMM = m_resultInMM.get_crop(0, 0, 3, count-1);
	//m_resultInMM.display("result in MM", true);
}

void HSCentroid::Test()
 {
	 CImg<float> image(500, 500, 1, 3);

	 int red[3] = {255, 0,0};
	 int gray[3] = {55, 55,50};

	 image.fill(0);

	 image(250, 250, 0, 2) = 255;

	 for( int i = 200; i>100; --i)
	 {
		 float oppac = 0.3 + 0.7*(1- 2*abs(0.01*i - (int)(0.01*i) - 0.5));

		 image.draw_circle(256, 256, i, red, oppac, 0); 
	 }

	 image.display();
 }

 */