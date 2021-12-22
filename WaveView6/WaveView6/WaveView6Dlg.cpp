// WaveView6Dlg.cpp : implementation file
//

#include "stdafx.h"
#include "WaveView6.h"
#include "WaveView6Dlg.h"
#include <string>
using namespace std;
//#include "HSCentroid.h"
#include "Zernike.h"

// openCV
#include "cv.h"
#include "highgui.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CWaveView6Dlg dialog




CWaveView6Dlg::CWaveView6Dlg(CWnd* pParent /*=NULL*/)
	: CDialog(CWaveView6Dlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CWaveView6Dlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_REF_FILE, m_fileRef);
	DDX_Control(pDX, IDC_TAR_FILE, m_fileTar);
	DDX_Control(pDX, IDC_IMAGE_REGION, m_ImageRegion);
}

BEGIN_MESSAGE_MAP(CWaveView6Dlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BTN_REF, &CWaveView6Dlg::OnBnClickedBtnRef)
	ON_BN_CLICKED(IDC_BTN_TAR, &CWaveView6Dlg::OnBnClickedBtnTar)
	ON_STN_CLICKED(IDC_REF_FILE, &CWaveView6Dlg::OnStnClickedRefFile)
//	ON_NOTIFY(TCN_SELCHANGE, IDC_TAB1, &CWaveView6Dlg::OnTcnSelchangeTab1)
	ON_BN_CLICKED(IDOK, &CWaveView6Dlg::OnBnClickedOk)
	ON_BN_CLICKED(IDC_BUTTON1, &CWaveView6Dlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_BUTTON2, &CWaveView6Dlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_BUTTON3, &CWaveView6Dlg::OnBnClickedButton3)
	ON_BN_CLICKED(IDC_BUTTON4, &CWaveView6Dlg::OnBnClickedButton4)
	ON_BN_CLICKED(IDC_BUTTON5, &CWaveView6Dlg::OnBnClickedButton5)
END_MESSAGE_MAP()


// CWaveView6Dlg message handlers

BOOL CWaveView6Dlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CWaveView6Dlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CWaveView6Dlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CWaveView6Dlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void CWaveView6Dlg::OnBnClickedBtnRef()
{
	//CFileDialog fileDlg( TRUE, NULL, NULL, OFN_ALLOWMULTISELECT | OFN_HIDEREADONLY, L"Image Files (*.jpg, *.TIFF, *.PNG, *.BMP)|*.jpg; *.TIFF; *.PNG; *.BMP||", this);
	//fileDlg.m_ofn.lpstrTitle = L"<- click to open target image file";
	CFileDialog fileDlg( TRUE, NULL, NULL, OFN_ALLOWMULTISELECT | OFN_HIDEREADONLY, L"Image Files (*.BMP)|*.BMP||", this);
	fileDlg.m_ofn.lpstrTitle = L"";
	
	CString szlstfile;

	if ( fileDlg.DoModal() == IDOK)
	{
		szlstfile = fileDlg.GetPathName(); // This is your selected file name with path
		this->m_fileRef.SetWindowTextW( szlstfile );
	}

	LPPICTURE m_pPicture;

	USES_CONVERSION;
	HRESULT hr = ::OleLoadPicturePath(const_cast<LPOLESTR>(T2COLE(szlstfile)),
									  NULL,
									  0,
									  0,
									  IID_IPicture,
									  reinterpret_cast<LPVOID *>(&m_pPicture));

	TRACE(_T("CImgViewerDoc::OnOpenDocument - OleLoadPicturePath(\"%s\"): ")
		  _T("hr = 0x%08X, m_pPicture = 0x%08X\n"), szlstfile, hr, m_pPicture);

	SIZE m_sizeInPix;
	SIZE m_sizeInHiMetric;

	if (SUCCEEDED(hr) && m_pPicture != NULL)
	{
		// get width and height of picture
		m_pPicture->get_Width(&m_sizeInHiMetric.cx);
		m_pPicture->get_Height(&m_sizeInHiMetric.cy);

		const int HIMETRIC_PER_INCH = 2540;

		HDC hDCScreen = ::GetDC(NULL);
		ASSERT(hDCScreen != NULL);
		// Pixels per logical inch along width
		const int nPixelsPerInchX = ::GetDeviceCaps(hDCScreen, LOGPIXELSX);
		// Pixels per logical inch along height
		const int nPixelsPerInchY = ::GetDeviceCaps(hDCScreen, LOGPIXELSY);
		::ReleaseDC(NULL, hDCScreen);

		//// convert himetric to pixels
		//	m_sizeInPix.cx = (nPixelsPerInchX * m_sizeInHiMetric.cx +
		//					  HIMETRIC_PER_INCH / 2) / HIMETRIC_PER_INCH;
		//	m_sizeInPix.cy = (nPixelsPerInchY * m_sizeInHiMetric.cy +
		//					  HIMETRIC_PER_INCH / 2) / HIMETRIC_PER_INCH;
      
		long hmWidth  = 0;
        long hmHeight = 0;

        m_pPicture->get_Width (&hmWidth);
        m_pPicture->get_Height(&hmHeight);

		RECT rcDlg;
		this->GetClientRect(&rcDlg);

		RECT rcImgReg;
		this->m_ImageRegion.GetClientRect(&rcImgReg);

		long Margin = 10;
		CPoint ptDim(rcImgReg.right - rcImgReg.left - 2*Margin, rcImgReg.bottom - rcImgReg.top-3*Margin);
		CPoint ptUL(rcImgReg.left + Margin, rcImgReg.top  + 2*Margin);


		float ratio = __min((float)ptDim.x/hmWidth, (float)ptDim.y/hmHeight);
		ptUL.x +=  (long)(0.5*(ptDim.x - ratio*hmWidth));
		ptUL.y +=  (long)(0.5*(ptDim.y - ratio*hmHeight));

		ptDim.x = (long)(ratio*hmWidth);
		ptDim.y = (long)(ratio*hmHeight);


		CPoint start, dim;
		start.x = 0;
		start.y = hmHeight;
		dim.x = hmWidth;
		dim.y = -hmHeight;
		

		HRESULT hrP = m_pPicture->Render(this->m_ImageRegion.GetWindowDC()->m_hDC,
											ptUL.x,
											ptUL.y,
											ptDim.x,
											ptDim.y,
											start.x,
											start.y,		// the starting line (from bottom) to display in y
											dim.x,
											dim.y,			//-dim.y is the pixel size in y to display
											&rcImgReg);

		//HRESULT hrP = m_pPicture->Render(this->m_ImageRegion.GetWindowDC()->m_hDC,
		//									ptUL.x,
		//									ptUL.y,
		//									ptDim.x,
		//									ptDim.y,
		//									0,
		//									hmHeight,
		//									hmWidth,
		//									-hmHeight,
		//									&rcImgReg);
		this->m_ImageRegion.SetWindowTextW(/*L"Reference Image: " + */szlstfile);
		HBITMAP btmp_img = this->m_ImageRegion.GetBitmap();
	}
}

void CWaveView6Dlg::OnBnClickedBtnTar()
{
	CString szlstfile;
	this->m_ImageRegion.GetWindowTextW(szlstfile);

	CT2CA pszConvertedAnsiString (szlstfile);
	std::string fileName(pszConvertedAnsiString);

	CFileDialog fileDlg( TRUE, NULL, NULL, OFN_ALLOWMULTISELECT | OFN_HIDEREADONLY, L"Image Files (*.jpg, *.TIFF, *.PNG, *.BMP)|*.jpg; *.TIFF; *.PNG; *.BMP||", this);
	fileDlg.m_ofn.lpstrTitle = L"<- click to open target image file";

	if ( fileDlg.DoModal() == IDOK)
	{
		CString szlstfile = fileDlg.GetPathName(); // This is your selected file name with path
		this->m_fileTar.SetWindowTextW( szlstfile );

		this->m_ImageRegion.SetWindowTextW(L"Target Image");

		CT2CA pszConvertedAnsiString (szlstfile);
		const std::string fileName(pszConvertedAnsiString);
		this->m_fileTar.SetWindowTextW( szlstfile );

		m_HSCentroid.OpenImageFromFile(fileName);
	}

	//CDC* pDC = m_ImageRegion.GetDC();
	//pDC->Ellipse(50, 100, 400, 400);
}

void CWaveView6Dlg::OnStnClickedRefFile()
{
	// TODO: Add your control notification handler code here
}

void CWaveView6Dlg::OnTcnSelchangeTab1(NMHDR *pNMHDR, LRESULT *pResult)
{
	// TODO: Add your control notification handler code here
	*pResult = 0;
}

void CWaveView6Dlg::OnBnClickedOk()
{
	// TODO: Add your control notification handler code here
	OnOK();
}

void CWaveView6Dlg::OnBnClickedButton1()
{
	// TODO: Add your control notification handler code here
	CFileDialog fileDlg( TRUE, NULL, NULL, OFN_ALLOWMULTISELECT | OFN_HIDEREADONLY, L"Text file (*.txt)|*.txt||", this);

	if ( fileDlg.DoModal() == IDOK)
	{
		bool bOk = true;
	}
}

void CWaveView6Dlg::OnBnClickedButton2()
{
	// TODO: Add your control notification handler code here
	m_HSCentroid.ResultDisplay();
}

void CWaveView6Dlg::OnBnClickedButton3()
{
	Zernike zernike;
	zernike.Excute();
}

void CWaveView6Dlg::OnBnClickedButton4()
{
	m_HSCentroid.SaveAsReference();
}

void CWaveView6Dlg::OnBnClickedButton5()
{
	// TODO: Add your control notification handler code here
	m_HSCentroid.Test();
}
