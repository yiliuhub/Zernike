#include "StdAfx.h"
#include "HSC.h"

HSC::HSC(void)
:m_stepX(0)
,m_stepY(0)
,m_LSRadius(5)
,m_minGridSize(10)
,m_centroidThreshold(8.0f)
{
  m_imgInput.assign();
  m_image.assign();
  m_resultMatrix.assign();
  m_resultImage.assign();
}

HSC::~HSC(void)
{
  m_imgInput.assign();
  m_image.assign();
  m_resultMatrix.assign();
  m_resultImage.assign();
}

bool 
HSC::OpenImageFromFile(std::string filename)
{
  m_imgInput = CImg<float>(filename.c_str());
 
  // black-white image
  m_image = m_imgInput.get_channel(0);
  m_image.blur(4.0);

  // color image
  m_resultImage = CImg<float>(m_imgInput);
  m_resultImage.fill(0);

  bool bOk = FindGridPicth();

  long t0 = cimg::time();
//  SearchInDirection(m_centerX, m_centerY, m_startIndexX, m_startIndexY, 1, 0);
//  SearchInDirection(m_centerX, m_centerY, m_startIndexX, m_startIndexY,-1, 0);

  RadiantSearch();

  long t1 = cimg::time();
  long dt = t1 - t0;

  float minVal, maxVal;
  minVal = m_resultMatrix.min_max(maxVal);

  minVal = m_resultImage.min_max(maxVal);

  m_resultImage.display("Result Image", true);
  m_resultMatrix.display("Result Matrix", true);

  std::string locStr;
  locStr += "(";
  char* str = new char[10];
  float locX = 400.65f, locY = 643.77f;
  sprintf(str, "%.3f", locX);
  locStr = locStr + str + ", ";
  sprintf(str, "%.3g", locY);
  locStr = locStr + str + ")";

  unsigned char purple[] = {255, 0, 255};

  m_imgInput.draw_text(100, 500, locStr.c_str(), purple);
  m_imgInput.display("Input Image", true);

  delete[] str;

  return true;
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
  int y[2] = {(int)( 0.4*m_image.height() ), (int)( 0.6*m_image.height() ) };

  int count = 0;
  UInt localMaxX, localMaxY;

  float maxVal = FindLocalMaxPixelPos1(x[0], y[0], x[1], y[1], m_centerX, m_centerY);

  // calc the sub-pixel at (m_centerX, m_centerY)
  float sum3x3 = m_image(m_centerX-1, m_centerY-1) + m_image(m_centerX, m_centerY-1) + m_image(m_centerX+1, m_centerY-1) +
                 m_image(m_centerX-1, m_centerY  ) + m_image(m_centerX, m_centerY  ) + m_image(m_centerX+1, m_centerY  ) +
                 m_image(m_centerX-1, m_centerY+1) + m_image(m_centerX, m_centerY+1) + m_image(m_centerX+1, m_centerY+1);

  float avgX = (m_centerX*sum3x3  - (m_image(m_centerX-1, m_centerY-1) + m_image(m_centerX-1, m_centerY)+m_image(m_centerX-1, m_centerY+1))
                                  + (m_image(m_centerX+1, m_centerY-1) + m_image(m_centerX+1, m_centerY)+m_image(m_centerX+1, m_centerY+1))
               )/sum3x3;

  float avgY = (m_centerY*sum3x3  - (m_image(m_centerX-1, m_centerY-1) + m_image(m_centerX, m_centerY-1)+m_image(m_centerX+1, m_centerY-1)) 
                                  + (m_image(m_centerX-1, m_centerY+1) + m_image(m_centerX, m_centerY+1)+m_image(m_centerX+1, m_centerY+1))
               )/sum3x3;

  m_resultImage(m_centerX, m_centerY, 0, 0) = 3000;
  m_resultImage(m_centerX, m_centerY, 0, 1) = avgX;
  m_resultImage(m_centerX, m_centerY, 0, 2) = avgY;


  // stepX
  int deltaX=10, deltaY = 10;
  x[0] = m_centerX;           y[0] = m_centerY - deltaY;
  x[1] = x[0] + 4*deltaX;     y[1] = m_centerY + deltaY;

  bOk = false;
  while( !bOk && x[0] < this->m_image.width() - 4*(int)deltaX)
  {
    x[0] += deltaX;
    x[1] = x[0] + 4*deltaX;

    maxVal = FindLocalMaxPixelPos(x[0], y[0], x[1], y[1], localMaxX, localMaxY);

    if( (maxVal < m_centroidThreshold) || ((int)localMaxX - x[0] < deltaX) || (x[1] - (int)localMaxX < deltaX) )
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

    if( (maxVal < m_centroidThreshold) ||((int)localMaxY - y[0] < deltaY) || (y[1] - (int)localMaxY < deltaY) )
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

  // calc the index of centroid @ (m_centerX, m_centerY) in m_resultMatrix
  m_startIndexX = cimg::max<UInt, int>(m_centerX, m_image.width()-m_centerX)/m_stepX + 1;
  m_startIndexY = cimg::max<UInt, int>(m_centerY, m_image.height()-m_centerY)/m_stepY + 1;

  m_resultMatrix = CImg<float>(2*m_startIndexX +1, 2*m_startIndexY +1, 1, 3);
  m_resultMatrix.fill(0);

  float moduledValue = 500*((m_startIndexX+m_startIndexY)%2);
  m_resultMatrix(m_startIndexX, m_startIndexY, 0, 0) = 2000 + moduledValue;
  m_resultMatrix(m_startIndexX, m_startIndexY, 0, 1) = avgX;
  m_resultMatrix(m_startIndexX, m_startIndexY, 0, 2) = avgY;

  return bOk;
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
    int Radius = 5;
    int Minimum_Grid = 10;
    int dx = cimg::sign(dirX);
    int dy = cimg::sign(dirY);

    if( (dirX == 0 && dirY == 0) || (abs(dirX) > 1) || (abs(dirY) > 1) ||
      (m_stepX < m_minGridSize) || (m_stepY < m_minGridSize)  )
    {
        return;
    }

    int stepX = dirX*m_stepX;
    int stepY = dirY*m_stepY;

    int centerX = startX + stepX;
    int centerY = startY + stepY;

    UInt resultIndexX = startIndexX + dirX;
    UInt resultIndexY = startIndexY + dirY;

    UInt localMaxX, localMaxY;

    while(centerX>Radius && centerX < m_image.width()-Radius &&
          centerY>Radius && centerY < m_image.height()-Radius)
    {
        float maxVal = FindLocalMaxPixelPos1( centerX-m_LSRadius, centerY-m_LSRadius, 
                                              centerX+m_LSRadius, centerY+m_LSRadius,
                                              localMaxX, localMaxY );
        if( maxVal > m_centroidThreshold )
        {
          int x = localMaxX;
          int y = localMaxY;
          float sum3x3 = m_image(x-1, y-1) + m_image(x, y-1) + m_image(x+1, y-1) +
                         m_image(x-1, y  ) + m_image(x, y  ) + m_image(x+1, y  ) +
                         m_image(x-1, y+1) + m_image(x, y+1) + m_image(x+1, y+1);

          float avgX =  (x*sum3x3  - (m_image(x-1, y-1) + m_image(x-1, y)+m_image(x-1, y+1))
                                   + (m_image(x+1, y-1) + m_image(x+1, y)+m_image(x+1, y+1))
                        )/sum3x3;

           float avgY = (y*sum3x3  - (m_image(x-1, y-1) + m_image(x, y-1)+m_image(x+1, y-1)) 
                                   + (m_image(x-1, y+1) + m_image(x, y+1)+m_image(x+1, y+1))
                        )/sum3x3;

            m_resultImage(localMaxX, localMaxY, 0, 0) = 3000;
            m_resultImage(localMaxX, localMaxY, 0, 1) = avgX;
            m_resultImage(localMaxX, localMaxY, 0, 2) = avgY;

            if( (int)resultIndexX > 0 && resultIndexY > 0 &&
                (int)resultIndexX < m_resultMatrix.width() &&
                (int)resultIndexY < m_resultMatrix.height() )
            {
                float moduledValue = 500*((resultIndexX+resultIndexY)%2);
                m_resultMatrix(resultIndexX, resultIndexY, 0, 0) = 2000 + moduledValue;
                m_resultMatrix(resultIndexX, resultIndexY, 0, 1) = avgX;
                m_resultMatrix(resultIndexX, resultIndexY, 0, 2) = avgY;
            }

            centerX = localMaxX + stepX;
            centerY = localMaxY + stepY;
        }
        else
        {
            centerX += stepX;
            centerY += stepY;
        }
        // update centroid index in result matrix along searching direction
        resultIndexX += dirX;
        resultIndexY += dirY;
    }
}

void HSC::RadiantSearch()
{
  	SearchInDirection(m_centerX, m_centerY, m_startIndexX, m_startIndexY, 0,  1);
    SearchInDirection(m_centerX, m_centerY, m_startIndexX, m_startIndexY, 0, -1);

    long t1 = cimg::time();

	  for(int i = 0; i<m_resultMatrix.height(); i++)
	  {
		  if(m_resultMatrix(m_startIndexX, i)>0)
		  {
			  int cx = (int)(m_resultMatrix(m_startIndexX, i, 0, 1)+0.5);
			  int cy = (int)(m_resultMatrix(m_startIndexX, i, 0, 2)+0.5);

        this->SearchInDirection(cx, cy, m_startIndexX, i,  1, 0); // to the right
        this->SearchInDirection(cx, cy, m_startIndexX, i, -1, 0); // to the left
		  }
	  }
}