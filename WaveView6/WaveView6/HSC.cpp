#include "StdAfx.h"
#include "HSC.h"

#include <assert.h>
#include <float.h>

HSC::HSC(void)
:m_stepX(0)
,m_stepY(0)
,m_LSRadius(5)
,m_minGridSize(10)
,m_discCenterX(-1.0)
,m_discCenterY(-1.0)
,m_CENTROID_THRESHOLD(8.0f)
,m_PIXEL_SIZE_MM(4.65E-3)
,m_FOCUS_LENGTH_MM(5.2)
,m_NORMAL_RADIUS_MM(1.3)
,m_brightSpotRadius(75)
{
  m_imgInput.assign();
  m_image.assign();
  m_resultImage.assign();
  m_resultMatrix.assign(45, 45, 1, 3);
  m_resultMatrix.fill(0);
  m_referenceMatrix.assign();
}

HSC::~HSC(void)
{
  Clear();
}

void HSC::Clear()
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
bool HSC::SetStreamData( T* data, int width, int height )
{
	bool bOk = true;
  
	if( width < 1 || height < 1 ||
	  data == NULL )
	{
		m_imgInput.assign();
		m_image.assign();
		bOk = false;
	}
	else
	{
		m_imgInput.assign(data, width, height, 1, 1);
		m_image = m_imgInput.get_channel(0);
	}

	return bOk;
}

	
bool HSC::OnExcute()
{
	SmoothAndDenoise();

	float diffusion = 3.0;
	m_image.blur( diffusion );

	m_CENTROID_THRESHOLD = m_image.mean();

	// color image
	m_resultImage = CImg<float>( m_image.width(), m_image.height(), 1, 3 );
	
	m_resultImage.fill(0);
	m_resultMatrix.fill(0);

	bool bOk = FindGridPicthNew();

	if(bOk)
	{
		RadiantSearch();

		CentroidValidation();

		CalcSubpixel();

		bOk &= ReferenceAssociation();
		bOk &= ReferenceAssociationNew();

		ConvertResultToMM();
	}

	return bOk;
}

void HSC::SmoothAndDenoise()
{
/*
	float diffsion = 2.0;
	float threshold = 100;
	int erodeVal = 10;
	int dilateVal = 20;

	CImg<float> res = m_image.get_blur(diffsion);
	res.threshold(threshold, false, false);

	res.erode(erodeVal);
	res.display("after erode", true);
	res.dilate(dilateVal);
	res.display("after dialate", true);

	cimg_forXY(m_image, x,y)
	{ 
	  if(res(x,y))
		  m_image(x,y) = 0;
	}
*/
	int x0 = (int)(0.5*m_image.width());
	int y0 = (int)(0.5*m_image.height());
	int R = m_brightSpotRadius;
	int black[] = {0, 0, 0};
	m_image.draw_circle(x0, y0, R, black, 1);
}


bool 
HSC::OpenImageFromFile(std::string filename)
{
	CImg<float> ringImg = CImg<float>(filename.c_str());

	OpenReferenceCentroid();
	m_resultMatrix.fill(0);

	CImg<float> tempImg = CImg<float>(filename.c_str());
	SetStreamData(tempImg.data(), tempImg.width(), tempImg.height());

	bool bOk = this->OnExcute();

	if(bOk)
	{
		DrawCentroidDisplacement();

		size_t found = filename.find_last_of(".");
		std::string outputName = filename.substr(0, found);
		outputName = outputName + ".txt";

		OutputResult( outputName );
	}

  
  int  qSize = m_imgInput.width() * m_imgInput.height();
  float* data = m_imgInput.data();

  for (int i = 0; i <  qSize; i++) {
    data[i] = this->m_resultImage.data()[i];
  }

  CImg<float> tempImage( data, m_imgInput.width(), m_imgInput.height() );
  tempImage.display("Result Image", true);

	return bOk;
}

//////////////////////////////////////////////////////////////////////////////////
// this function does two things: find the grid picthes of centroid spots 
// in x and y direction, and figure out the result matrix dimmension
//////////////////////////////////////////////////////////////////////////////////
bool 
HSC::FindGridPicth()
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

bool 
HSC::FindGridPicthNew()
{
	bool bOk = true;

	int x[2] = {(int)( 0.3*m_image.width() ),  (int)( 0.7*m_image.width() ) };
    int y[2] = {(int)( 0.3*m_image.height() ), (int)( 0.7*m_image.height() ) };

    CImg<float> sampleImage = m_image.get_crop(x[0], y[0], x[1], y[1]);

	//sampleImage.display();

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

	m_stepX = 0;
	float startPos = 0;
	float endPos = 0;

	float minStep = x[1]-x[0], 
		  maxStep = 0;

	int localMaxCountX = 0;

	// find local maximum of intensity sum
	for(int i=radius; i<x[1]-x[0]-radius; i += step) {
		float maxValue = FindLocalMaxPixelPos1<float>(  VSum, 
														i-radius, 
														0, 
														i+radius, 
														0, 
														localMaxX, 
														localMaxY );

		if(abs(i - localMaxX) < step) {
			if(localMaxCountX > 0) {
				float curStep = localMaxX - endPos;
				minStep = __min(minStep, curStep);
				maxStep = __max(maxStep, curStep);
			}
			else {
				startPos = endPos = localMaxX;
			}

			VFlag(localMaxX, 0) = 1;
			endPos = localMaxX;
			++localMaxCountX;
			i += step;
		}
	}

	if( localMaxCountX < 2 || 
		(maxStep-minStep > step )
	) {
		return false;
	}

	m_stepX = (endPos - startPos)/(localMaxCountX-1);
	m_centerX = x[0] + startPos;

	m_stepY = 0;
	startPos = 0;
	endPos = 0;
	int localMaxCountY = 0;

	minStep = x[1]-x[0], maxStep = 0;

	// find local maximum of intensity sum
	for(int i=radius; i<y[1]-y[0]-radius; i += step) {
		float maxValue = FindLocalMaxPixelPos1(HSum, i-radius, 0, i+radius, 0, localMaxX, localMaxY);
		if(abs(i - localMaxX) < step) {
			if(localMaxCountY > 0) {
				float curStep = localMaxX - endPos;
				minStep = __min(minStep, curStep);
				maxStep = __max(maxStep, curStep);
			}
			else {
				startPos = endPos = localMaxX;
			}

			HFlag(localMaxX, 0) = 1;
			endPos = localMaxX;
			++localMaxCountY;
			i += step;
		}
	}

	if(localMaxCountY < 2 || (maxStep-minStep > step)) {
		return false;
	}

	m_stepY = (endPos - startPos)/(localMaxCountY-1);
	m_centerY = y[0] + startPos;

	float minX = abs( x[0] - m_referenceMatrix( m_stationaryX, m_stationaryY, 1 ) );
	float minY = abs( y[0] - m_referenceMatrix( m_stationaryX, m_stationaryY, 2 ) );
	int minIndexX = 0;
	int minIndexY = 0;

	for( int i = 0; i < VFlag.width(); i++ ) {
		if( VFlag( i, 0 ) == 1 ) {
			float curDiff = abs( x[0]+i - m_referenceMatrix( m_stationaryX, m_stationaryY, 1 ) );
			if( curDiff < minX ) {
				minX = curDiff;
				minIndexX = x[0] + i;
			}
		}
	}

	for( int i = 0; i < HFlag.width(); i++ ) {
		if( HFlag( i, 0 ) == 1 ) {
			float curDiff = abs( y[0]+i - m_referenceMatrix( m_stationaryX, m_stationaryY, 2 ) );
			if( curDiff < minY ) {
				minY = curDiff;
				minIndexY = y[0] + i;
			}
		}
	}

	m_centerX = minIndexX;
	m_centerY = minIndexY;

	m_startIndexX = m_stationaryX;
	m_startIndexY = m_stationaryY;


	// grid intersection
	CImg<float> gridImage( m_image );
	int black[] = {0,0,0};
	int r = 3;

	for(int row = 0; row < HFlag.width(); row++ ) {
		if( HFlag(row) ) {
			for( int col = 0; col < VFlag.width(); col++ ) {
				if( VFlag( col ) )
				gridImage.draw_circle(x[0] + col, y[0] + row, r, black);
			}
		}
	}

	gridImage.display(" intersection grid", true );

/*
	m_resultMatrix.assign( m_referenceMatrix );	// same dimension as reference matrix
	m_resultMatrix.fill( 0 );

	int size = 4;

	CImg<float> stationaryMatrix = m_referenceMatrix.get_crop(  m_stationaryX, 
																m_stationaryX,
																m_stationaryX+size-1,
																m_stationaryX+size-1
															);

	float minVal = FLT_MAX;
	float val = 0;
	int newCenterX = m_centerX;
	int newCenterY = m_centerY;

	for( int i = 0; i <= localMaxCountX - size; i++) 
	{
		for( int j = 0; j <= localMaxCountY - size; j++)
		{
			val = 0;
			cimg_forXY(stationaryMatrix, ix, jy)
			{
				float dx = stationaryMatrix( ix, jy, 1 ) - (m_centerX + (i+ix)*m_stepX);
				float dy = stationaryMatrix( ix, jy, 2 ) - (m_centerY + (j+jy)*m_stepY);
				val += dx*dx + dy*dy;
			}
			
			if( val < minVal ) {
				minVal = val;
				newCenterX = m_centerX + i*m_stepX;
				newCenterY = m_centerY + j*m_stepY;
			}
		}
	}

	m_centerX = newCenterX;
	m_centerY = newCenterY;

	m_startIndexX = m_stationaryX;
	m_startIndexY = m_stationaryY;

	m_referenceAssocX = m_stationaryX;
	m_referenceAssocY = m_stationaryY;
*/

	float maxValue = FindLocalMaxPixelPos1<float>( m_image, m_centerX-radius, m_centerY-radius, 
												   m_centerX+radius, m_centerY+radius, 
												   localMaxX, localMaxY );

	m_resultMatrix( m_startIndexX, m_startIndexY, 0, 1 ) = localMaxX;
	m_resultMatrix( m_startIndexX, m_startIndexY, 0, 2 ) = localMaxY;

	if( maxValue > m_CENTROID_THRESHOLD ) {
		m_resultMatrix( m_startIndexX, m_startIndexY, 0, 0 ) = 3000;
	}
	else {
		m_resultMatrix( m_startIndexX, m_startIndexY, 0, 0 ) = INVALID_CENTROID;
	}

	return true;
}


float HSC::FindLocalMaxPixelPos(UInt x0, UInt y0, UInt x1, UINT y1, UInt& maxX, UInt& maxY )
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

template<typename T>
float HSC::FindLocalMaxPixelPos1(CImg<T> img, int x0, int y0, int x1, int y1, int& maxX, int& maxY )
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



float HSC::FindLocalMaxPixelPos1(UInt x0, UInt y0, UInt x1, UINT y1, UInt& maxX, UInt& maxY )
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


float HSC::FindLocalMaxPixelPos2(UInt& maxX, UInt& maxY )
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

void  HSC::SearchInDirection(const UInt startX, const UInt startY, const UInt startIndexX, const UInt startIndexY, const int dirX, const int dirY)
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
	    
			  float moduledValue = (float)(1000*((resultIndexX+resultIndexY)%2));
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

void HSC::RadiantSearch()
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

bool HSC::LocateMissingCentroid()
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
				if( rVal == stationaryVal ) //HS_STATIONARY_CENTROID 
				{  
					stationaryX = j;
					stationaryY = i;
					bOk = true;
				}
			}
		}

	  if(bOk)
	  {
		  x[0] = (int)m_referenceMatrix(stationaryX, stationaryY, 1);
		  y[0] = (int)m_referenceMatrix(stationaryX, stationaryY, 2);
		  x[1] = (int)m_referenceMatrix(stationaryX+1, stationaryY+1, 1);
		  y[1] = (int)m_referenceMatrix(stationaryX+1, stationaryY+1, 2);

		  int locDx = (int)(0.5*(x[1] - x[0]));
		  int locDy = (int)(0.5*(x[1] - x[0]));

		  x[0] -= locDx;
		  y[0] -= locDy;

		  x[1] -= locDx;
		  y[1] -= locDy;

		  UInt localMaxX, localMaxY;
		  float local_max = FindLocalMaxPixelPos1( x[0], y[0], x[1], y[1], localMaxX, localMaxY );

		  if( local_max > m_CENTROID_THRESHOLD )
		  {
		    this->m_centerX = localMaxX;
		    this->m_centerY = localMaxY;

		    // calc the index of centroid @ (m_centerX, m_centerY) in m_resultMatrix
		    m_startIndexX = stationaryX ;
		    m_startIndexY = stationaryY ;

		    float moduledValue = 1000*( ( m_startIndexX + m_startIndexY ) % 2 );

			m_resultMatrix( m_startIndexX, m_startIndexY, 0, 0 ) = 3000 + moduledValue;
			m_resultMatrix( m_startIndexX, m_startIndexY, 0, 1 ) = localMaxX;
			m_resultMatrix( m_startIndexX, m_startIndexY, 0, 2 ) = localMaxY;

		    bOk = true;
		  }
		  else
		  {
			x[0] = (int)m_referenceMatrix(stationaryX+2, stationaryY+2, 1);
		    y[0] = (int)m_referenceMatrix(stationaryX+2, stationaryY+2, 2);
		    x[1] = (int)m_referenceMatrix(stationaryX+3, stationaryY+3, 1);
		    y[1] = (int)m_referenceMatrix(stationaryX+3, stationaryY+3, 2);

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
			  m_startIndexX = stationaryX + 3;
			  m_startIndexY = stationaryY + 3;

			  float moduledValue = 1000*( ( m_startIndexX+m_startIndexY ) % 2 );

			  m_resultMatrix( m_startIndexX, m_startIndexY, 0, 0 ) = 3000 + moduledValue;
			  m_resultMatrix( m_startIndexX, m_startIndexY, 0, 1 ) = localMaxX;
			  m_resultMatrix( m_startIndexX, m_startIndexY, 0, 2 ) = localMaxY;

			  bOk = true;
			}
		}
    }
  }

  return bOk;
}
void HSC::CentroidValidation()
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
void HSC::CalcSubpixel()
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


template<typename T>
void  HSC::CalcSubpixel(CImg<T> img, int locX, int locY, float& subLocX, float& subLocY)
{

	subLocX = locX;
	subLocY = locY;

	int lx = locX;
	int ly = locY;

	if( lx > 1 && ly > 1 &&
		lx < img.width() -1 &&
		ly < img.height() -1 
	)
	{
		float denomX =  2*img( lx, ly ) - ( img( lx+1, ly ) + img( lx-1, ly ) );
		float denomY =  2*img( lx, ly ) - ( img( lx, ly+1 ) + img( lx, ly-1 ) );
		float denom45 =  2*img( lx, ly ) - ( img( lx+1, ly+1 ) + img( lx-1, ly-1 ) );
		float denom135 =  2*img( lx, ly ) - ( img( lx-1, ly+1 ) + img( lx+1, ly-1 ) );

		float dx0 = img( lx+1, ly ) - img( lx-1, ly );
		float dx45 = img( lx+1, ly+1 ) - img( lx-1, ly-1 );
		float dx135 = img( lx+1, ly-1 ) - img( lx-1, ly+1 );

		float dy0 = img( lx, ly+1) - img( lx, ly-1 );
		float dy45 = img( lx+1, ly+1 ) - img( lx-1, ly-1 );
		float dy135 = img( lx-1, ly+1 ) - img( lx+1, ly-1 );

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

		subLocX = subpixelX;
		subLocY = subpixelY;
	}
}

void HSC::ResultDisplay()
{
	if(this->m_imgInput.data() == NULL)
		return;

	m_image.display();
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

void HSC::UpdateResultImage()
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
void HSC::OpenReferenceCentroid()
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

void HSC::SaveAsReference()
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



bool HSC::ReferenceAssociationNew()
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

  // check if reference matrix has the right number of stationary centroids
  if(m_referenceMatrix.width() && m_referenceMatrix.height() )
  {
    float stationaryVal = HS_STATIONARY_CENTROID;

    // find the upper-left stationary centroid
    for(int i=0; !bFound && i < m_referenceMatrix.height(); i++)   
    {
	    for(int j=0; !bFound && j<this->m_referenceMatrix.width(); j++)  
      {
		    float rVal = m_referenceMatrix(j, i, 0, 0);
			if( rVal == stationaryVal ) { //HS_STATIONARY_CENTROID  
			    stationaryX = j;
			    stationaryY = i;
			    bFound = true;
		    }
	    }
    }

	double sumX = 0.0, sumY = 0.0;
    if(bFound) 
    {
      for( int i = 0; bOk && i < Size; i++ )  
      {
	      for( int j = 0; bOk && j < Size; j++ )  
        {
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
	else
	{
		return false;
	}

    // we found SizexSize stationary centroids
	if( bOk )
	{
		m_discCenterX = sumX / ( Size * Size );
		m_discCenterY = sumY / ( Size * Size );
	}
    else
    {
		return bOk;
    }
  }
  else
  {
	  return false;
  }

  float refStepX = ( m_referenceMatrix(stationaryX+Size-1, stationaryY+Size-1, 1)
                      - m_referenceMatrix(stationaryX, stationaryY, 1))/(Size-1);
  float refStepY = ( m_referenceMatrix(stationaryX+Size-1, stationaryY+Size-1, 2)
                      - m_referenceMatrix(stationaryX, stationaryY, 2))/(Size-1);

  int startIndexX, startIndexY;
  bOk = false;

  for( int i = 0; !bOk && i < m_resultMatrix.height(); i++ )  
  {
    for( int j = 0; !bOk && j < m_resultMatrix.width(); j++ )  
      {
        //if(m_resultMatrix(j, i, 0) > Centroid_Value )
        {
          if( abs(m_discCenterX-m_resultMatrix(j, i, 1)) < m_stepX &&
              abs(m_discCenterY-m_resultMatrix(j, i, 2)) < m_stepY )
          {
            startIndexX = j-1;
            startIndexY = i-1;
            bOk = true;
          }
        }
      }
  }

  if(!bOk)
  {
    return bOk;
  }

  int boarder = 3;
  boarder = __min( boarder, startIndexX);
  boarder = __min( boarder, m_referenceMatrix.width() - startIndexX );

  boarder = __min( boarder, startIndexY);
  boarder = __min( boarder, m_referenceMatrix.height() - startIndexY );

  int matchingX;
  int matchingY;
  float minDiff = 1.0e9;

  int bestCount = 1;
  if(bOk)
  {
    CImg<float> refer4x4 = m_referenceMatrix.get_crop(stationaryX, stationaryY, stationaryX+Size -1, stationaryY+Size-1);

    for( int y = startIndexY - boarder; y < startIndexY + boarder; y++ )
    {
      for( int x = startIndexX -boarder; x < startIndexX + boarder; x++ )
      {
        float curDiff = 0.0;

        CImg<float> target4x4 = m_resultMatrix.get_crop( x, y, x + Size - 1, y + Size -1);
        CImg<float> diff4x4 = ( target4x4 - refer4x4 ).channels( 1, 2 );
        
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

	m_referenceAssocX = matchingX;
	m_referenceAssocY = matchingY;

  return bOk;
}


template<typename T>
bool  HSC::ReferenceAssociationNew(CImg<T> refStationary, CImg<T> targetImg, int& associationX, int& associationY)
{
  int refH = refStationary.height();
  int refW = refStationary.width();

  if( refH < 2 || refW < 2 )
  {
    return false;
  }

  bool bOk = false;
  float Centroid_Value = CENTROID_FLAG;

  // pitch size of reference stationary centroids
  float refStepX = (refStationary(refW-1, refH-1, 1) - refStationary(0, 0, 1))/(refW-1);
  float refStepY = (refStationary(refW-1, refH-1, 2) - refStationary(0, 0, 2))/(refH-1);

  int closedX = 0, closedY = 0;

  for(int i = 0; !bOk && i < targetImg.height(); i++)
  {
    for(int j = 0; !bOk && j < targetImg.width(); j++)
    {
      if( abs(refStationary(0,0,1) - targetImg(j, i, 1)) < 1.5*refStepX &&
          abs(refStationary(0,0,2) - targetImg(j, i, 2)) < 1.5*refStepY)
      {
        closedX = j;
        closedY = i;
        bOk = true;
      }
    }
  }

  if(bOk)
  {
    int boarder = 3;
    int startX = max(closedX - boarder, 0);
    int startY = max(closedY - boarder, 0);
    int endX = min(closedX + boarder, targetImg.width()-1);
    int endY = min(closedY + boarder, targetImg.height()-1);
    float minDiff = 1.0e99;

    for(int y = startY; y <= endY; y++)
    {
      for(int x = startX; x <= endX; x++)
      {
        float curDiff = 0.0;
        float validCount = refW * refH;

        CImg<float> target = targetImg.get_crop( x, y, x + refW - 1, y + refH -1);
        CImg<float> diffImg = ( target - refStationary ).channels( 1, 2 );
        
        cimg_forXY(target, i, j)
        {
          if( target( i, j, 0 ) < Centroid_Value )
          {
            diffImg( i, j, 0 ) = diffImg( i, j, 1 ) = 0.0;
            validCount -= 1.0;
          } 
        }

        if( validCount >= 1 )
        {
          // average difference square
          float curDiff = diffImg.pow( 2 ).sum() / validCount;

          if( curDiff < minDiff )
          {
            minDiff = curDiff;
            associationX = x;
            associationY = y;
          }
        }

      }
    }

  }// if(bOk)

  return bOk;
}


bool HSC::ReferenceAssociation()
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


void HSC::DrawCentroidDisplacement()
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

//	m_resultImage.display("result image", true);
}


void HSC::OutputResult( std::string outputName )
{
	FILE* pf = fopen(outputName.c_str(), "w");

	if( m_stationaryX < 0 || m_stationaryY < 0 ||
		m_discCenterX < 0 || m_discCenterY < 0 )
	{
		fclose(pf);
		return;
	}

	double Normal_Radius_mm = m_NORMAL_RADIUS_MM;
	double Pixel_Size_mm = m_PIXEL_SIZE_MM;
	double Focus_Length_mm = m_FOCUS_LENGTH_MM;

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

void HSC::GetCentroidResult( float* locX, float* locY, float* dx, float* dy, int& dataSize )
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
bool  HSC::SetAndCalcCentroid(
							T* data, int width, int height,		// image data
							float* pixelX, float* pixelY,			// reference data
							int* gridX, int* gridY, int* flag, int refSize,		// reference data
							float* locX, float* locY, float* dx, float* dy, int& dataSize	// output data size
						  )
{
	bool bOk = true;
	bOk &= SetStreamData( data, width, height );

	if(bOk)
	{
		SetReferenceData( pixelX, pixelY, gridX, gridY,  flag, refSize );

		bOk &= this->OnExcute();

		if(bOk)
		{
			GetCentroidResult(locX, locY, dx, dy, dataSize );
		}
	}

	DrawCentroidDisplacement();

  int  qSize = width * height;

  for (int i = 0; i <  qSize; i++) {
    data[i] = this->m_resultImage.data[i];
  }

  CImg<T> tempImage( data, width, height );
  tempImage.display("Result Image", true);

	return bOk;
}
  
  template<typename T>
  bool  HSC::SetAndCalcCentroidWithConfig(
							T* data, int width, int height,		// image data
							float* pixelX, float* pixelY,			// reference data
							int* gridX, int* gridY, int* flag, int refSize,		// reference data
							float* locX, float* locY, float* dx, float* dy, int& dataSize,	// output data size
							float Pixel_Size_MM = 4.65E-3, float Focus_Length_MM = 5.2, float Normal_Radius_MM	= 1.3,	// hardware parameters
							int brightSpotRadius = 75 // circle radius near center to be erased
						  )
  {
		this->m_PIXEL_SIZE_MM = Pixel_Size_MM;
		this->m_FOCUS_LENGTH_MM = Focus_length_MM;
		this->m_NORMAL_RADIUS_MM = Normal_Radius_MM;

		this->m_brightSpotRadius = brightSpotRadius;

		return  SetAndCalcCentroid<T>(
										data, width, height,
										pixelX, pixelY,
										gridX, gridY, flag, refSize,
										locX, locY, dx, dy, dataSize,
									  );
  }

bool  HSC::CalcKeratometerRingSpotsWithReference(float* data, int width, int height,
												 int refRadius, int refLocX, int refLocY,
												float* locX, float* locY, int& dataSize	// output data size
											  )
{
	if( width < 1 || height < 1 || data == NULL )
	{
		return false;
	}

	CImg<float> ringImg(data, width, height);
	//ringImg.blur(2);

	int R = 70;
	int CX, CY;

	this->CalcKeratometerRingCenter(data, width, height, refRadius, refLocX, refLocY, CX, CY );

	CImg<float> centerImg = ringImg.get_crop(CX-R, CY-R, CX+R, CY+R);

	int shiftX = CX - R,
		shiftY = CY - R;
	CX = R;
	CY = R;


	float dotProd = centerImg.dot(centerImg);

	// the radius of the ring
	float Radius = refRadius;
	int sideSize = 5;
	float PIE = cimg::PI;
	
	int maxX, maxY;

	// find one centroid from the 16 spots
	int vertX[5], vertY[5];
	float vertexMaxVal[4];
	float maxVal = 0;
	// start searching spots around the ring
	int x[3];
	int y[3];

	for(int k = 0; k<4; k++)
	{
		x[0] = CX +(int)(Radius*cos(0.5*k*PIE));
		y[0] = CY +(int)(Radius*sin(0.5*k*PIE));
		x[1] = x[0] - sideSize;
		y[1] = y[0] - sideSize;
		x[2] = x[0] + sideSize;
		y[2] = y[0] + sideSize;

		int tempX, tempY;
		float curVal = FindLocalMaxPixelPos1( centerImg, x[1], y[1], x[2], y[2], tempX, tempY );

		if( curVal > maxVal ) {
			maxVal = curVal;
			maxX = tempX;
			maxY = tempY;
		}
	}


	int ringSpotNum = 16;
	double startAngle = 0.0;
	if(abs(maxX-CX) > 1.0)
	{
		startAngle = atan(double(maxY-CY)/(maxX-CX));
	}

	double stepAngle = 2.0*cimg_library::cimg::PI/ringSpotNum;
	int boxHalfSize = (int)(0.5*Radius * sin(2.0*PIE/ringSpotNum));


	int count = 0;
	float ringSpotX[16];
	float ringSpotY[16];
	
	double angle = startAngle;

	x[0] = maxX;
	y[0] = maxY;

	for(int i=0; i<ringSpotNum; i++)
	{
		x[1] = x[0] - boxHalfSize;
		x[2] = x[0] + boxHalfSize;
		y[1] = y[0] - boxHalfSize;
		y[2] = y[0] + boxHalfSize;

		maxVal = FindLocalMaxPixelPos1(centerImg, x[1], y[1], x[2], y[2], maxX, maxY);

		bool bFound = false;

		if(abs(x[0]-maxX) < boxHalfSize && abs(y[0]-maxY) < boxHalfSize )
		{
			bFound = true;
		}
		else
		{
			x[1] = maxX - boxHalfSize;
			x[2] = maxX + boxHalfSize;
			y[1] = maxY - boxHalfSize;
			y[2] = maxY + boxHalfSize;
			maxVal = FindLocalMaxPixelPos1(centerImg, x[1], y[1], x[2], y[2], maxX, maxY);
			if(abs(x[0]-maxX) < boxHalfSize && abs(y[0]-maxY) < boxHalfSize )
			{
				bFound = true;	
			}
		}

		// move to the next spot
		if(bFound)	// move relative to the new found (maxX, maxY)
		{
			float subX, subY;
			this->CalcSubpixel(centerImg, maxX, maxY, subX, subY);
			ringSpotX[count] = subX;
			ringSpotY[count] = subY;

			x[0] = CX + (maxX - CX)*cos(stepAngle) - (maxY - CY)*sin(stepAngle);
			y[0] = CY + (maxY - CY)*cos(stepAngle) + (maxX - CX)*sin(stepAngle);
			count++;
		}
		else // move relative to the reference location (x[0], y[0])
		{
			int refX = x[0];
			int refY = y[0];
			x[0] = CX + (refX - CX)*cos(stepAngle) - (refY - CY)*sin(stepAngle);
			y[0] = CY + (refY - CY)*cos(stepAngle) + (refX - CX)*sin(stepAngle);
		}
	}

	dataSize = count;
	for( int i = 0; i < count; i++)
	{
		locX[i] = shiftX+ringSpotX[i];
		locY[i] = shiftY+ringSpotY[i];
	}

	for(int i=0; i<count; i++)
	{
		ringImg((int)(locX[i]+0.5), (int)(locY[i]+0.5)) = 255;
	}

	ringImg.display();	// erase
	return true;
}


bool  HSC::CalcKeratometerRingSpots(float* data, int width, int height,
								float* locX, float* locY, int& dataSize,	// output data size
								int refRadius, int refLocX, int refLocY		// reference varibles
							  )
{
	int centerX, centerY;
	
	CalcKeratometerRingMirrorCenter( data, width, height, 
									 refRadius, refLocX, refLocY, 
									 centerX, centerY );
	
	float radius;
	CalcKeratometerRingRadius(data, width, height, centerX, centerY, radius );

	float contrast = this->CalcContrast( data, width, height, radius, centerX, centerY );

	return CalcKeratometerRingSpotsWithReference(data, width, height, 
												 radius, centerX, centerY, 
												 locX, locY, dataSize );

}

void HSC::SetHardwareConfig(float Pixel_Size_MM, float Focus_Length_MM, float Normal_Radius_MM)
{
	this->m_PIXEL_SIZE_MM = Pixel_Size_MM;
	this->m_FOCUS_LENGTH_MM = Focus_Length_MM;
	this->m_NORMAL_RADIUS_MM = Normal_Radius_MM;
}

void  HSC::SetBrightSpotRadius( int Bright_Spot_Radius )
{
	this->m_brightSpotRadius = Bright_Spot_Radius;
}

void  HSC::SetReferenceData( float* pixelX, float* pixelY, 
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
void HSC::ConvertResultToMM()
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


	double Normal_Radius_mm = m_NORMAL_RADIUS_MM;
	double Pixel_Size_mm = m_PIXEL_SIZE_MM;
	double Focus_Length_mm = m_FOCUS_LENGTH_MM;

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
						/*
						fprintf(pf,"%f,%f,%f,%f\n", 
									distToCenterX / Normal_Radius, 
									distToCenterY / Normal_Radius,
									(m_referenceMatrix(x+DX, y+DY, 1) - m_resultMatrix(x, y, 1))/Focus_Length, 
									(m_referenceMatrix(x+DX, y+DY, 2) - m_resultMatrix(x, y, 2))/Focus_Length
									);
						*/
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

float HSC::CalcContrast( float* imageData, int imageWidth, int imageHeight, int refRadius, int refX, int refY )
{

	int innerR = (int)( refRadius - 10 );
	int BOUNDARY_STRIP = 10;

	if( imageData == NULL || 
		refX - refRadius - BOUNDARY_STRIP < 0 ||
		refY - refRadius - BOUNDARY_STRIP < 0 ||
		refX + refRadius + BOUNDARY_STRIP > imageWidth ||
		refY + refRadius + BOUNDARY_STRIP > imageHeight )
	{
		return 0.0;
	}

	CImg<float> focusImg( imageData, imageWidth, imageHeight );

	focusImg.crop(  refX - refRadius - BOUNDARY_STRIP,
					refY - refRadius - BOUNDARY_STRIP,
					refX + refRadius + BOUNDARY_STRIP,
					refY + refRadius + BOUNDARY_STRIP
				 );

	int one[3] = {1, 1, 1};
	CImg<float> labelImg0( focusImg );
	labelImg0.fill( 0 );
	CImg<float> labelImg1( focusImg );
	labelImg1.fill( 0 );
	labelImg0.draw_circle(refRadius + BOUNDARY_STRIP, refRadius + BOUNDARY_STRIP, refRadius + BOUNDARY_STRIP, one ); 
	labelImg1.draw_circle(refRadius + BOUNDARY_STRIP, refRadius + BOUNDARY_STRIP, innerR, one ); 

	CImg<float> labelImg = labelImg0 - labelImg1;

	labelImg.display();

	float mean = focusImg.mean();

	CImgList<float> grad = focusImg.get_gradient();
	grad[0].pow(2);
	grad[1].pow(2);

	grad[0] = grad[0]*labelImg0;
	grad[1] = grad[1]*labelImg0;

	float contrast = grad[0].sum() + grad[1].sum();

	// normalize
	contrast /= ( mean*mean * labelImg0.sum() );
	return contrast;
}

void HSC::CalcKeratometerRingCenter( float* data,int width, int height, int refRadius, int refLocX, int refLocY, int& centerX, int& centerY ) {

   CImg<float> img( data, width, height );
   img.blur(2);

   int x0 = refLocX;
   int y0 = refLocY;

   if( refLocX == 0 || refLocY==0 ) {
	   x0 = static_cast<int>(0.5*width);
	   y0 = static_cast<int>(0.5*height);
   }

   int Boundary_Margin = 10;
   int R = refRadius + Boundary_Margin;

   CImg<float> regionImg = img.get_crop(x0-R, y0-R, x0+R, y0+R);
   regionImg.blur(2);

   CImg<float> resultImg(regionImg.width(), regionImg.height() );
   resultImg.fill(0);
  

   int boxSize = 5;

   int xLow  = refRadius;
   int xHigh = xLow + 2*Boundary_Margin;
   int yLow  = refRadius;
   int yHigh = yLow + 2*Boundary_Margin;

   int Div = 8;
   double theta = cimg::PI/(double)Div;

   for( int i = yLow; i < yHigh; i++ ) {
     for( int j = xLow; j < xHigh; j++ ) {

		 for(int k=0; k<Div; k++)
		 {
			 int cx = j + refRadius*cos(k*theta);
			 int cy = i + refRadius*sin(k*theta);

			 // the reflection of (cx, cy) about (j, i)
			 int acx = j + refRadius*cos(k*theta+cimg::PI);
			 int acy = i + refRadius*sin(k*theta+cimg::PI);

			 /*resultImg(j, i) += __min( regionImg.get_crop(cx-boxSize, cy-boxSize, cx+boxSize, cy+boxSize).sum(),
									   regionImg.get_crop(acx-boxSize, acy-boxSize, acx+boxSize, acy+boxSize).sum());
			*/
			 resultImg(j, i) += regionImg.get_crop(cx-boxSize, cy-boxSize, cx+boxSize, cy+boxSize).sum()*
								regionImg.get_crop(acx-boxSize, acy-boxSize, acx+boxSize, acy+boxSize).sum();
		 }
     }
   }

   float max = this->FindLocalMaxPixelPos1( resultImg, 0, 0, 2*R, 2*R, centerX, centerY );

   centerX += x0-R;
   centerY += y0-R;
 }


void HSC::CalcKeratometerRingMirrorCenter( float* data,int width, int height, int refRadius, int refLocX, int refLocY , int& centerX, int& centerY )
{
	CImg<float> ringImg(data, width, height );

	int x0 = refLocX;
	int y0 = refLocY;

	if( refLocX == 0 || refLocY == 0)
	{
		x0 = static_cast<int>(0.5*width);
		y0 = static_cast<int>(0.5*height);
	}


	int halfSearchRange = static_cast<int>(1.2*refRadius);
	if( refRadius == 0 ) {
		halfSearchRange = static_cast<int>(0.08*__min(width, height));
	}

	int x1 = x0 - halfSearchRange, x2 = x0 + halfSearchRange;
	int y1 = y0 - halfSearchRange, y2 = y0 + halfSearchRange;
	int R = halfSearchRange; //refRadius + 10;

	CImg<float> leftImg;
	CImg<float> rightImg;
	CImg<float> topImg;
	CImg<float> bottomImg;

	int maxX, maxY;
	float dot = 0;
	float maxDot = 0;

	for(int x = x1; x < x2; x++) {
		leftImg = ringImg.get_crop(x-R,y0-R,x,y0+R);
		leftImg.mirror('x');

		rightImg = ringImg.get_crop(x, y0-R, x+R, y0+R);

		dot = leftImg.dot( rightImg );
		
		if( dot > maxDot )
		{
			maxDot = dot;
			maxX = x;
		}
	}

	maxDot = 0;
	for( int y = y1; y < y2; y++ ) {
		topImg = ringImg.get_crop(x0-R, y-R, x0+R, y).mirror('y');
		bottomImg = ringImg.get_crop(x0-R, y, x0+R, y+R);
		dot = topImg.dot( bottomImg );
		if( dot > maxDot ) {
			maxDot = dot;
			maxY = y;
		}
	}

	centerX = maxX;
	centerY = maxY;


	int white[3] = {255, 255, 255};
	ringImg.get_crop(maxX-R, maxY-R, maxX+R, maxY+R);		// erase
	ringImg.draw_point(maxX, maxY, white, 1).display();		// erase

}

void HSC::CalcKeratometerRingRadius( float* data,int width, int height, int centerX, int centerY, float& radius ) {

	CImg<float> ringImg(data, width, height);
	ringImg.blur(2);

	float radiusUpper = 50;
	float radiusLower = 15;

	int boxSize = 4;

	// don't go beyond the range
	radiusUpper = __min(centerX, radiusUpper + boxSize );
	radiusUpper = __min(width - centerX, radiusUpper + boxSize );
	radiusUpper = __min(centerY, radiusUpper + boxSize );
	radiusUpper = __min(height - centerY, radiusUpper + boxSize );

	if( centerX < radiusUpper + boxSize			||
		centerX + radiusUpper + boxSize > width ||
		centerY < radiusUpper + boxSize			||
		centerY + radiusUpper + boxSize > height ) {
		return;
	}

	int maxLocX[2] = { centerX, centerX },
		maxLocY[2] = { centerY, centerY };
	float maxDot = 0;
	int r = boxSize;

	for( int k = 0; k < 2; k++ ) {
		maxDot = 0;
		float dx = cos(0.5*cimg::PI*k);
		float dy = sin(0.5*cimg::PI*k);

		for( int i = radiusLower; i<radiusUpper; i++) {
			int x[2] = {(int)(centerX - dx*i), (int)(centerX + dx*i) };	
			int y[2] = {(int)(centerY - dy*i), (int)(centerY + dy*i) };
			CImg<float> orgImg = ringImg.get_crop(x[0]-r, y[0]-r, x[0]+r, y[0]+r);
			CImg<float> mirImg = ringImg.get_crop(x[1]-r, y[1]-r, x[1]+r, y[1]+r);

			if( k==0 )
			{
				mirImg.mirror('x');
			}
			else
			{
				mirImg.mirror('y');
			}
			float dot = orgImg.dot( mirImg );

			if(dot > maxDot )
			{
				maxDot = dot;
				maxLocX[k] = x[0];
				maxLocY[k] = y[0];
			}
		}		
	}

	float radius1 = sqrt((float)((centerX-maxLocX[0])*(centerX-maxLocX[0])) + 
						(float)((centerY-maxLocY[0])*(centerY-maxLocY[0])));

	float radius2 = sqrt((float)((centerX-maxLocX[1])*(centerX-maxLocX[1])) + 
						 (float)((centerY-maxLocY[1])*(centerY-maxLocY[1])));

	radius = 0.5*(radius1 + radius2);

}


void HSC::Test()
{



			


	 // simulate function SetAndCalcCentroid()

	//m_imgInput

	 int width =  m_referenceMatrix.width();
	 int height = m_referenceMatrix.height();

	 int size = width * height;
	 float CentroidVal = CENTROID_FLAG;
	 float StationaryVal = HS_STATIONARY_CENTROID;

	 float* refX = new float[size];
	 float* refY = new float[size];
	 int*   refGridX  = new int[size];
	 int*   refGridY  = new int[size];
	 int*   flag   = new int[size];

	 int count = 0;

	 for(int j = 0; j < height; j++)
	 {
		 for(int i = 0; i < width; i++)
		 {
			if(m_referenceMatrix(i, j, 0) >= CentroidVal )
			{
				refX[count] = m_referenceMatrix(i, j, 1);
				refY[count] = m_referenceMatrix(i, j, 2);
				refGridX[count] = i;
				refGridY[count] = j;
				flag[count] = 0;
				if(m_referenceMatrix(i, j, 0) == StationaryVal)
				{
					flag[count] = 0;
				}

				++count;
			 }
		 }
	 }
/*
	 int refSize = count;
	 float* locX = new float[2500];
	 float* locY = new float[2500];
	 float*   dx  = new float[2500];
	 float*   dy  = new float[2500];

	 this->SetAndCalcCentroid(  m_imgInput.data(), m_imgInput.width(), m_imgInput.height(), 
								refX, refY, refGridX, refGridY, flag, refSize, 
								locX, locY, dx, dy, dataSize );


	 DrawCentroidDisplacement();
	 	OutputResult( "CentroidResult.txt" );


	 delete[] refX;
	 delete[] refY;
	 delete[] refGridX;
	 delete[] refGridY;
	 delete[] flag;
	 delete[] locX;
	 delete[] locY;
	 delete[] dx;
	 delete[] dy;
*/
 }