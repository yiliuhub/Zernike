

#ifndef DllExport_h
#define DllExport_h 1

	#ifdef __DLL_EXPORT__
		#define DllExport __declspec(dllexport)
	#elif __DLL_IMPORT__
		#define DllExport  __declspec(dllimport)
	#else
		#define DllExport  
	#endif

#endif
