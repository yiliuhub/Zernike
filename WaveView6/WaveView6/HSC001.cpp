#include "StdAfx.h"
#include "HSC.h"

HSC::HSC(void)
:m_stepX(0)
,m_stepY(0)
,m_LSRadius(5)
,m_minGridSize(10)
,m_CENTROID_THRESHOLD(8.0f)
{
  m_imgInput.assign();
  m_image.assign();
  m_resultMatrix.assign();
  m_resultImage.assign();
}

HSC::~HSC(void)
{
  Clear();
}

void HSC::Clear()
{
  m_stepX = 0.0;
  m_stepY = 0.0;
  m_LSRadius = 5;
  m_minGridSize = 10;
  m_CENTROID_THRESHOLD = 8.0f;

  m_imgInput.assign();
  m_image.assign();
  m_resultMatrix.assign();
  m_resultImage.assign();
}

bool 
HSC::OpenImageFromFile(std::string filename)
{
	this->Clear();
  m_imgInput = CImg<float>(filename.c_str());
 
  // black-white image
  m_image = m_imgInput.get_channel(0);
  m_image.blur(4.0);

  m_CENTROID_THRESHOLD = m_image.mean();

  // color image
  m_resultImage = CImg<float>(m_imgInput);
  m_resultImage.fill(0);

  bool bOk = FindGridPicth();

  RadiantSearch();

  UpdateResultImage();
  //ResultDisplay();

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
  int y[2] = {(int)( 0.6*m_image.height() ), (int)( 0.8*m_image.height() ) };

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

  // calc the index of centroid @ (m_centerX, m_centerY) in m_resultMatrix
  m_startIndexX = cimg::max<UInt, int>(m_centerX, m_image.width()-m_centerX)/m_stepX + 1;
  m_startIndexY = cimg::max<UInt, int>(m_centerY, m_image.height()-m_centerY)/m_stepY + 1;

  m_resultMatrix = CImg<float>(2*m_startIndexX +1, 2*m_startIndexY +1, 1, 3);
  m_resultMatrix.fill(0);

  float moduledValue = 1000*((m_startIndexX+m_startIndexY)%2);
  m_resultMatrix(m_startIndexX, m_startIndexY, 0, 0) = 3000 + moduledValue;
  m_resultMatrix(m_startIndexX, m_startIndexY, 0, 1) = avgX;
  m_resultMatrix(m_startIndexX, m_startIndexY, 0, 2) = avgY;

  return bOk;
}

bool HSC::LocateMissingCentroid()
{
  bool bOk = false;
  int x[2] = {(int)( 0.5*m_image.width() - 4*m_stepX ),  (int)( 0.5*m_image.width() - 2*m_stepX) };
  int y[2] = {(int)( 0.5*m_image.height() - 4*m_stepY ), (int)( 0.5*m_image.height() - 2*m_stepY) };

  UInt localMaxX, localMaxY;
  float REGION_MAX = FindLocalMaxPixelPos1(x[0], y[0], x[1], y[1], localMaxX, localMaxY);

  float startX = localMaxX;
  float startY = localMaxY;
  float maxDist = 999999999.0;
  float prevX, prevY;
  UINT  closestToCenterX = localMaxX,
        closestToCenterY = localMaxY;

  for( y[0] = startY - m_LSRadius; !bOk && y[0] < 0.5*m_image.height() + 2*m_stepY; y[0] += m_stepY )
  {
    y[1] = y[0] + 2*m_LSRadius;

    for( x[0] = startX - m_LSRadius; x[0] < 0.5*m_image.width() + 2*m_stepX; x[0] += m_stepX )
    {
      x[1] = x[0] + 2*m_LSRadius;
      float local_max = FindLocalMaxPixelPos1(x[0], y[0], x[1], y[1], localMaxX, localMaxY);

      if( maxDist > abs( localMaxX - 0.5*m_image.width()) + abs( localMaxY - 0.5*m_image.height()) )
      {
        maxDist = abs( localMaxX - 0.5*m_image.width()) + abs( localMaxY - 0.5*m_image.height());
      }

      if( local_max < 0.5*REGION_MAX )
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
        if( abs(localMaxX-0.5*m_image.width()) + abs(localMaxY-0.5*m_image.height()) <
          abs(closestToCenterX-0.5*m_image.width()) + abs(closestToCenterY-0.5*m_image.height()) )
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
    local_max = FindLocalMaxPixelPos1(m_centerX-m_stepX-m_LSRadius, m_centerY-m_LSRadius, 
                                      m_centerX-m_stepX+m_LSRadius, m_centerY+m_LSRadius, locMaxX[0], locMaxY[0]);

    local_max = FindLocalMaxPixelPos1(m_centerX+m_stepX-m_LSRadius, m_centerY-m_LSRadius, 
                                      m_centerX+m_stepX+m_LSRadius, m_centerY+m_LSRadius, locMaxX[1], locMaxY[1]);

    local_max = FindLocalMaxPixelPos1(m_centerX-m_LSRadius, m_centerY-m_stepY-m_LSRadius, 
                                      m_centerX+m_LSRadius, m_centerY-m_stepY+m_LSRadius, locMaxX[2], locMaxY[2]);

    local_max = FindLocalMaxPixelPos1(m_centerX-m_LSRadius, m_centerY+m_stepY-m_LSRadius,
                                      m_centerX+m_LSRadius, m_centerY+m_stepY+m_LSRadius, locMaxX[3], locMaxY[3]);

    m_centerX = UINT( 0.25*( locMaxX[0] + locMaxX[1] + locMaxX[2] + locMaxX[3] ) + 0.5 );
    m_centerY = UINT( 0.25*( locMaxY[0] + locMaxY[1] + locMaxY[2] + locMaxY[3] ) + 0.5 );
  }
  else
  {
    m_centerX = closestToCenterX;
    m_centerY = closestToCenterY;
  }


  // calc the index of centroid @ (m_centerX, m_centerY) in m_resultMatrix
  m_startIndexX = cimg::max<UInt, int>(m_centerX, m_image.width()-m_centerX)/m_stepX + 1;
  m_startIndexY = cimg::max<UInt, int>(m_centerY, m_image.height()-m_centerY)/m_stepY + 1;

  m_resultMatrix = CImg<float>(2*m_startIndexX +1, 2*m_startIndexY +1, 1, 3);
  m_resultMatrix.fill(0);

  m_resultMatrix(m_startIndexX, m_startIndexY) = 1;
  m_resultMatrix(m_startIndexX, m_startIndexY, 1) = m_centerX;
  m_resultMatrix(m_startIndexX, m_startIndexY, 2) = m_centerY;

  for(int i = -2; i <= 2; i++)
  {
    m_resultImage(m_centerX+i, m_centerY, 0 ) = m_resultImage(m_centerX, m_centerY+i, 0 ) = CENTROID_FLAG;
    m_resultImage(m_centerX+i, m_centerY, 1 ) = m_resultImage(m_centerX, m_centerY+i, 1 ) = CENTROID_FLAG;
    m_resultImage(m_centerX+i, m_centerY, 2 ) = m_resultImage(m_centerX, m_centerY+i, 2 ) = 55;
  }

  if(!bOk)  // no artifact (missing centroid)
  {
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

    float moduledValue = 1000*( ( m_startIndexX+m_startIndexY )%2 );
    m_resultMatrix( m_startIndexX, m_startIndexY, 0, 0 ) = 3000 + moduledValue;
    m_resultMatrix( m_startIndexX, m_startIndexY, 0, 1 ) = avgX;
    m_resultMatrix( m_startIndexX, m_startIndexY, 0, 2 ) = avgY;

    for(int i = -2; i <= 2; i++)
    {
      m_resultImage(m_centerX+i, m_centerY, 0 ) = m_resultImage(m_centerX, m_centerY+i, 0 ) = CENTROID_FLAG;
      m_resultImage(m_centerX+i, m_centerY, 1 ) = m_resultImage(m_centerX, m_centerY+i, 1 ) = avgX;
      m_resultImage(m_centerX+i, m_centerY, 2 ) = m_resultImage(m_centerX, m_centerY+i, 2 ) = avgY;
    }
  }

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
        if( maxVal > m_CENTROID_THRESHOLD )
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
                float moduledValue = 1000*((resultIndexX+resultIndexY)%2);
                m_resultMatrix(resultIndexX, resultIndexY, 0, 0) = 3000 + moduledValue;
                m_resultMatrix(resultIndexX, resultIndexY, 0, 1) = avgX;
                m_resultMatrix(resultIndexX, resultIndexY, 0, 2) = avgY;
            }

            centerX = localMaxX + stepX;
            centerY = localMaxY + stepY;
        }
        else
        {
          if( (int)resultIndexX > 0 && resultIndexY > 0 &&
                (int)resultIndexX < m_resultMatrix.width() &&
                (int)resultIndexY < m_resultMatrix.height() )
          {
              m_resultMatrix(resultIndexX, resultIndexY, 0, 0) = 100;
              m_resultMatrix(resultIndexX, resultIndexY, 0, 1) = centerX;
              m_resultMatrix(resultIndexX, resultIndexY, 0, 2) = centerY;
          }

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
  	SearchInDirection(m_centerX, m_centerY, m_startIndexX, m_startIndexY, 0,  1); // up
    SearchInDirection(m_centerX, m_centerY, m_startIndexX, m_startIndexY, 0, -1); // down

	  for(int i = 0; i<m_resultMatrix.height(); i++)
	  {
		  if(m_resultMatrix(m_startIndexX, i) > 0)
		  {
			  int cx = (int)(m_resultMatrix(m_startIndexX, i, 0, 1)+0.5);
			  int cy = (int)(m_resultMatrix(m_startIndexX, i, 0, 2)+0.5);

        this->SearchInDirection(cx, cy, m_startIndexX, i,  1, 0); // to the right
        this->SearchInDirection(cx, cy, m_startIndexX, i, -1, 0); // to the left
		  }
	  }
}

bool HSC::LocateMissingCentroid()
{
  bool bOk = true;
  int boundX[2] = {(int)( 0.4*m_image.width() ),  (int)( 0.6*m_image.width() ) };
  int boundY[2] = {(int)( 0.4*m_image.height() ), (int)( 0.6*m_image.height() ) };

  UInt localMaxX, localMaxY;
  float REGION_MAX = FindLocalMaxPixelPos1(boundX[0], boundY[0], boundX[1], boundY[1], localMaxX, localMaxY);

  m_centerX = localMaxX;
  m_centerY = localMaxY;

  // stepX
  int dx=10, dy = 10;
  int x[2] = { localMaxX, localMaxX + 4*dx };
  int y[2] = { localMaxY-dy, localMaxY+dy };

  bOk = false;
  while( !bOk && x[0] < this->m_image.width() - 4*(int)dx)
  {
    x[0] += dx;
    x[1] = x[0] + 4*dx;

    float maxVal = FindLocalMaxPixelPos(x[0], y[0], x[1], y[1], localMaxX, localMaxY);
    if( (maxVal > 0.75*REGION_MAX) && ((int)localMaxX - x[0] < dx) && (x[1]-(int)localMaxX < dx) )
    {
      bOk = true;
    }
  }

  if(bOk)
    m_stepX = localMaxX - m_centerX;

  // reset bOk
  bOk = false;

  // stepY
  x[0] = m_centerX - dx;      y[0] = m_centerY;
  x[1] = m_centerX + dx;      y[1] = m_centerY + 4*dy;

  while(!bOk && y[0] < m_image.height() - 4*dy)
  {
    y[0] += dx;
    y[1] = y[0] + 4*dx;

    float maxVal = FindLocalMaxPixelPos(x[0], y[0], x[1], y[1], localMaxX, localMaxY);

    if( (maxVal > 0.75*REGION_MAX) &&((int)localMaxY - y[0] < dy) && (y[1] - (int)localMaxY < dy) )
    {
      bOk = true;
    }
  }

  if(bOk)
  {
    m_stepY = localMaxY - m_centerY;
  }

  // calc the index of centroid @ (m_centerX, m_centerY) in m_resultMatrix
  m_startIndexX = cimg::max<UInt, int>(m_centerX, m_image.width()-m_centerX)/m_stepX + 1;
  m_startIndexY = cimg::max<UInt, int>(m_centerY, m_image.height()-m_centerY)/m_stepY + 1;

  m_resultMatrix = CImg<float>(2*m_startIndexX +1, 2*m_startIndexY +1, 1, 3);
  m_resultMatrix.fill(0);

  return bOk;
}

void HSC::ResultDisplay()
{
	if(this->m_imgInput.data() == NULL)
		return;

  CImg<float> result2 = m_resultImage.get_channels(1, 2);

    result2.display("Result Image", true);
    m_resultMatrix.display("Result Matrix", true);

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
//	  if( m_resultImage(x, y, 0, 0) == 3000)
	  {
		  int fontSize = 24;(int)(13*__max(ratioX, ratioY));

		  locStr.clear();
		  sprintf(strX, "%.3f", (float)m_resultImage(x, y, 0, 1));
		  sprintf(strY, "%.3f", (float)m_resultImage(x, y, 0, 2));
		  locStr = locStr + "(" + strX + ", " + strY + ")";

		  CImg<float> img(m_imgInput);

		  if(x> 0.5*disp.width())
			x = x - locStr.size();
		  else
			x = x+5;
		  if(y<0.5*disp.height())
			y = y + 5;
		  else
			y = y-20;

		  img.draw_text(x,y,locStr.c_str(),foreground_color,background_color,1,fontSize);
	  }

      img.display(disp);
    }
    if( disp.is_resized()) disp.resize();
    disp.wait(25);
  }
  //delete[] str;

  return;
}

void HSC::OutputResult()
{
	FILE* pf = fopen("ResultOutput.txt", "w");

	for(int y = 0; y < m_resultMatrix.height(); y++)
	{
		for( int x = 0; x < m_resultMatrix.width(); x++)
		{
			if( m_resultMatrix(x, y, 0, 0) > 100 )
			{
				fprintf(pf,"%7.3f %7.3f\n", m_resultMatrix(x, y, 0, 1), m_resultMatrix(x, y, 0, 2));
			}
		}
	}

	fprintf(pf, "\n");
	fclose(pf);
}

void HSC::UpdateResultImage()
{
	for(int y = 0; y < m_resultMatrix.height(); y++)
	{
		for(int x = 0; x < m_resultMatrix.width(); x++)
		{
			if(m_resultMatrix(x, y, 0, 0) > 100)
			{
				int lx = (int)(m_resultMatrix(x, y, 0, 1) + 0.5);
				int ly = (int)(m_resultMatrix(x, y, 0, 2) + 0.5);
				for(int i = -5; i <= 5; i++)
					for(int j = -5; j <=5; j++)
					{
						m_resultImage(lx+j, ly+i, 0, 0) = CENTROID_REGION;
						m_resultImage(lx+j, ly+i, 0, 1) = m_resultMatrix(x, y, 0, 1);
						m_resultImage(lx+j, ly+i, 0, 2) = m_resultMatrix(x, y, 0, 2);
					}
			}
		}
	}
}