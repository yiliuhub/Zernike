

/*****************************************************************************
 FILE error.h
      define the TYPE of errors
******************************************************************************/

#ifndef ERROR_H
#define ERROR_H

#define FAILURE                 0
#define SUCCESS                 1

#define MAX_ERROR_NUM          14

#define GEOM_OK                 0
#define OUTOF_MEMORY            1
#define TOO_FEW_POINTS          2
#define THREE_POINTS_COLLINEAR  3
#define FOUR_POINTS_COLPLANER   4
#define TWO_LINE_PARALLEL       5
#define TWO_POINTS_SAME         6
#define ZERO_DIRECTION          7
#define ZERO_RADIUS             8
#define ZERO_LENGTH             9
#define ZERO_ANGLE             10
#define INVALID_ANGLE          11
#define POINTS_COLPLANER       12
#define POINTS_COLLINEAR       13

#endif
