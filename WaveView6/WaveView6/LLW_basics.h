/**************************************************************
*	File: LLW_basics.h
*	Author:	Yi Liu
*	Date:	October 2010
*
****************************************************************/
#ifndef LLW_Basics_H
#define LLW_Basics_H

#pragma once

template< class Type >
struct LLWPoint2D{
    Type x;
    Type y;
};

template< class Type >
struct LLWPoint3D{
    Type x;
    Type y;
    Type z;
};

#endif