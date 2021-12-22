// WaveView6.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// CWaveView6App:
// See WaveView6.cpp for the implementation of this class
//

class CWaveView6App : public CWinApp
{
public:
	CWaveView6App();

// Overrides
	public:
	virtual BOOL InitInstance();

// Implementation

	DECLARE_MESSAGE_MAP()
};

extern CWaveView6App theApp;