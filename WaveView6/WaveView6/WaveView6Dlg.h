// WaveView6Dlg.h : header file
//

#pragma once
#include "afxwin.h"
#include "HSC.h"

//#include "HSCentroid.h"
//using namespace LLW;


// CWaveView6Dlg dialog
class CWaveView6Dlg : public CDialog
{
// Construction
public:
	CWaveView6Dlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_WAVEVIEW6_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
private:
	// file names of reference and target images
	CStatic m_fileRef;
	CStatic m_fileTar;
	//LLW::Centroid::HSCentroid m_HSCentroid;
	HSC m_HSCentroid;
	
public:
	afx_msg void OnBnClickedBtnRef();
	afx_msg void OnBnClickedBtnTar();
	// static region where image will be displayed
	CStatic m_ImageRegion;
	afx_msg void OnStnClickedRefFile();
	afx_msg void OnTcnSelchangeTab1(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedButton1();
	afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedButton3();
	afx_msg void OnBnClickedButton4();
	afx_msg void OnBnClickedButton5();
};
