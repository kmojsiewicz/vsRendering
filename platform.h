#ifndef PLATFORM_H
#define PLATFORM_H

#include <stdlib.h>
#include <fcntl.h>
#include <io.h>
#include <iostream>
#include <io.h>
#include <algorithm>
#include <cstdint>
#include <string>
#include <atomic>
#include <thread>

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

#if defined(_WIN32) && !defined(__MINGW32__)
	#pragma comment(lib, "user32.lib")										// Visual Studio Only
	#pragma comment(lib, "gdi32.lib")										// For other Windows Compilers please add
	#pragma comment(lib, "opengl32.lib")									// these libs to your linker input
#endif

#include <dwmapi.h>
#include <GL/gl.h>
#if !defined(__MINGW32__)
	#pragma comment(lib, "Dwmapi.lib")
#endif

#if defined max
	#undef max
#endif
#if defined min
	#undef min
#endif

class CPlatform;
static std::unique_ptr<CPlatform> engine = nullptr;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#if defined(UNICODE) || defined(_UNICODE)
	#define strT(s) L##s
#else
	#define strT(s) s
#endif

std::wstring ConvertS2W(std::string s)
{
#ifdef __MINGW32__
	wchar_t* buffer = new wchar_t[s.length() + 1];
	mbstowcs(buffer, s.c_str(), s.length());
	buffer[s.length()] = L'\0';
#else
	int count = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, NULL, 0);
	wchar_t* buffer = new wchar_t[count];
	MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, buffer, count);
#endif
	std::wstring w(buffer);
	delete[] buffer;
	return w;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Pixel - Represents a 32-Bit RGBA colour

struct Pixel
{
	union {
		uint32_t n = (0xFF << 24);
		struct { uint8_t r; uint8_t g; uint8_t b; uint8_t a; };
	};

	enum Mode { NORMAL, MASK, ALPHA, CUSTOM };

	Pixel() {
		r = 0; g = 0; b = 0; a = 0xFF;
	}
	Pixel(uint8_t red, uint8_t green, uint8_t blue, uint8_t alpha = 0xFF) {
		n = red | (green << 8) | (blue << 16) | (alpha << 24);
	}
	Pixel(uint32_t p) {
		n = p;
	}

	Pixel  operator * (const float i) const;
};

static const Pixel
	GREY(192, 192, 192), DARK_GREY(128, 128, 128), VERY_DARK_GREY(64, 64, 64),
	RED(255, 0, 0), DARK_RED(128, 0, 0), VERY_DARK_RED(64, 0, 0),
	YELLOW(255, 255, 0), DARK_YELLOW(128, 128, 0), VERY_DARK_YELLOW(64, 64, 0),
	GREEN(0, 255, 0), DARK_GREEN(0, 128, 0), VERY_DARK_GREEN(0, 64, 0),
	CYAN(0, 255, 255), DARK_CYAN(0, 128, 128), VERY_DARK_CYAN(0, 64, 64),
	BLUE(0, 0, 255), DARK_BLUE(0, 0, 128), VERY_DARK_BLUE(0, 0, 64),
	MAGENTA(255, 0, 255), DARK_MAGENTA(128, 0, 128), VERY_DARK_MAGENTA(64, 0, 64),
	WHITE(255, 255, 255), BLACK(0, 0, 0), BLANK(0, 0, 0, 0);

Pixel  Pixel::operator * (const float i) const
{
	float fR = std::min(255.0f, std::max(0.0f, float(r) * i));
	float fG = std::min(255.0f, std::max(0.0f, float(g) * i));
	float fB = std::min(255.0f, std::max(0.0f, float(b) * i));
	return Pixel(uint8_t(fR), uint8_t(fG), uint8_t(fB), a);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class CLayer
{
public:
	float fOffsetX = 0.0;
	float fOffsetY = 0.0;
	float fScaleX = 1.0;
	float fScaleY = 1.0;
	bool bShow = false;
	bool bUpdate = false;
	Pixel tint = WHITE;

	CLayer();
	~CLayer();
	void	Create(uint32_t width, uint32_t height, bool filter = false, bool clamp = true);
	Pixel*	GetData();
	Pixel	GetPixel(int32_t x, int32_t y) const;
	bool	SetPixel(int32_t x, int32_t y, Pixel p);
	void	Clear(Pixel p);
	void	Update();

public:
	int32_t id = -1;
	int32_t width = 0;
	int32_t height = 0;
	std::vector<Pixel> pColData;
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class CPlatform
{
public:
	CPlatform();
	virtual ~CPlatform();

	int Construct(int32_t screen_w, int32_t screen_h, bool full_screen = false, bool vsync = false, bool console = false);
	int Start();

private:
	HWND		hWnd = nullptr;

	int			iWindowSizeX;
	int			iWindowSizeY;
	int			iViewPosX = 0;
	int			iViewPosY = 0;
	int			iViewSizeX = 0;
	int			iViewSizeY = 0;
	int			iScreenSizeX;
	int			iScreenSizeY;

	bool		bFullScreen;
	bool		bEnableVSYNC;
	bool		bEnableConsole;

	uint32_t	nLastFPS = 0;
	int			nFrameCount = 0;
	float		fFrameTimer = 1.0f;
	float		fLastElapsed = 0.0f;
	std::chrono::time_point<std::chrono::system_clock> m_tp1, m_tp2;

	HDC			glDeviceContext = 0;
	HGLRC		glRenderContext = 0;
	bool		b_glSync = false;

	CLayer		layer;
private:
	static std::atomic<bool> bAtomActive;									// If anything sets this flag to false, the engine "should" shut down gracefully
	void		EngineThread();
	int			ThreadStartUp()					{ return 1; }				// The main engine thread
	int			ThreadCleanUp()					{ RenderDestroyDevice(); PostMessage(hWnd, WM_DESTROY, 0, 0); return 1; }
	void		RedirectIOToConsole();
	int			CreateWindowPane(int iWindowPosX, int iWindowPosY);
	void		UpdateWindowSize(int32_t x, int32_t y) { iWindowSizeX = x; iWindowSizeY = y; UpdateViewport(); }
	void		UpdateViewport();
	int			CreateGraphics();
	int			ApplicationStartUp()			{ return 1; }
	int			ApplicationCleanUp()			{ return 1; }
	int			StartSystemEventLoop();
	int			RenderCreateDevice(std::vector<void*> params);
	int			RenderDestroyDevice();
	void		RenderUpdateViewport();
	void		RenderClearBuffer(Pixel p, bool bDepth);
	void		RenderDisplayFrame();
	void		RenderPrepareDrawing();
public:
	uint32_t	RenderCreateTexture(const uint32_t width, const uint32_t height, const bool filtered, const bool clamp);
	void		RenderApplyTexture(uint32_t id);
	uint32_t	RenderDeleteTexture(const uint32_t id);
	void		RenderUpdateTexture(uint32_t id, int width, int height, const Pixel* pixels);
	void		RenderReadTexture(uint32_t id, int width, int height, Pixel* pixels);
	void		RenderDrawLayerQuad(const float offsetX, const float offsetY, const float scaleX, const float scaleY, const Pixel tint);

public: // User Override Interfaces
	// User must override these functions as required. I have not made
	// them abstract because I do need a default behaviour to occur if
	// they are not overwritten

	virtual bool OnUserCreate()					{ return false; }			// Called once on application startup, use to load your resources
	virtual bool OnUserUpdate(float fElapsedTime)	{ (void)(fElapsedTime); return false;  }	// Called every frame, and provides you with a time per frame value
	virtual bool OnUserDestroy()				{ return true; }			// Called once on application termination, so you can be one clean coder
	void		PrepareEngine();
	void		CoreUpdate();
	void		Terminate() { bAtomActive = false; }

public:
	int32_t ScreenWidth() const					{ return iScreenSizeX; }	// Returns the width of the screen in "pixels"
	int32_t ScreenHeight() const				{ return iScreenSizeY; }	// Returns the height of the screen in "pixels"
	void Clear(Pixel p)							{ layer.Clear(p); }			// Clears entire draw target to Pixel
	CLayer* GetDrawTarget()						{ return &layer; }

public:																		// Branding
	std::string sAppName;

};

std::atomic<bool> CPlatform::bAtomActive{ false };							// singleton instances
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// CLayer IMPLEMENTATION

CLayer::CLayer()
{
	width = 0;
	height = 0;
}

CLayer::~CLayer()
{
	pColData.clear();
}

void CLayer::Create(uint32_t w, uint32_t h, bool filter, bool clamp)
{
	width = w;
	height = h;
	pColData.resize(width * height, (0xFF << 24));
	if (pColData.data() != nullptr) {
		id = engine->RenderCreateTexture(width, height, filter, clamp);
	}
}

Pixel* CLayer::GetData()
{
	return pColData.data();
}

Pixel CLayer::GetPixel(int32_t x, int32_t y) const
{
	if (x >= 0 && x < width && y >= 0 && y < height)
		return pColData[y * width + x];
	else
		return Pixel(0, 0, 0, 0);
}

bool CLayer::SetPixel(int32_t x, int32_t y, Pixel p)
{
	if (x >= 0 && x < width && y >= 0 && y < height) {
		pColData[y * width + x] = p;
		return true;
	}
	else
		return false;
}

void CLayer::Clear(Pixel p)
{
	std:fill(pColData.begin(), pColData.end(), p);
}

void CLayer::Update()
{
	if (pColData.data() == nullptr) return;
	engine->RenderApplyTexture(id);
	engine->RenderUpdateTexture(id, width, height, pColData.data());
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// CPlatform IMPLEMENTATION

CPlatform::CPlatform()
{
	sAppName = "Undefined";
	engine = std::unique_ptr<CPlatform>(this);
}

CPlatform::~CPlatform()
{
	engine.release();
}

int CPlatform::StartSystemEventLoop()
{
	MSG msg;
	while (GetMessage(&msg, NULL, 0, 0) > 0) {
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return 1;
}

// Windows Event Handler - this is statically connected to the windows event system
static LRESULT CALLBACK WindowEvent(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	switch (uMsg) {
	case WM_CLOSE:		engine->Terminate(); return 0;
	case WM_DESTROY:	PostQuitMessage(0);	return 0;
	}
	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}

int CPlatform::Construct(int32_t screen_w, int32_t screen_h, bool full_screen, bool vsync, bool console)
{
	iScreenSizeX = screen_w;
	iScreenSizeY = screen_h;
	iWindowSizeX = screen_w * 1;
	iWindowSizeY = screen_h * 1;
	bFullScreen = full_screen;
	bEnableVSYNC = vsync;
	bEnableConsole = console;

	if (iScreenSizeX <= 0 || iScreenSizeY <= 0)
		return 0;
	return 1;
}

void CPlatform::UpdateViewport()
{
	int32_t ww = iScreenSizeX * 1;
	int32_t wh = iScreenSizeY * 1;
	float wasp = (float)ww / (float)wh;

	iViewSizeX = (int32_t)iWindowSizeX;
	iViewSizeY = (int32_t)((float)iViewSizeX / wasp);

	if (iViewSizeY > iWindowSizeY) {
		iViewSizeY = iWindowSizeY;
		iViewSizeX = (int32_t)((float)iViewSizeY * wasp);
	}

	iViewPosX = (iWindowSizeX - iViewSizeX) >> 1;
	iViewPosY = (iWindowSizeY - iViewSizeY) >> 1;
}

int CPlatform::Start()
{
	if (ApplicationStartUp() != 1) return 0;
	if (CreateWindowPane(30, 30) != 1) return 0;							// Construct the window
	UpdateWindowSize(iWindowSizeX, iWindowSizeY);

	bAtomActive = true;
	std::thread t = std::thread(&CPlatform::EngineThread, this);				// Start the thread
	StartSystemEventLoop();													// Some implementations may form an event loop here
	t.join();																// Wait for thread to be exited

	if (ApplicationCleanUp() != 1) return 0;

	return 1;
}

int CPlatform::CreateWindowPane(int iWindowPosX, int iWindowPosY)
{
	WNDCLASS wc;
	wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	wc.hInstance = GetModuleHandle(nullptr);
	wc.lpfnWndProc = WindowEvent;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.lpszMenuName = nullptr;
	wc.hbrBackground = nullptr;
	wc.lpszClassName = strT("ENGINE");
	RegisterClass(&wc);

	DWORD dwExStyle = WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;					// Define window furniture
	DWORD dwStyle = WS_CAPTION | WS_SYSMENU | WS_VISIBLE | WS_THICKFRAME;

	int iWindowTopLeftX = iWindowPosX;
	int iWindowTopLeftY = iWindowPosY;

	if (bFullScreen) {														// Handle Fullscreen
		dwExStyle = 0;
		dwStyle = WS_VISIBLE | WS_POPUP;
		HMONITOR hmon = MonitorFromWindow(hWnd, MONITOR_DEFAULTTONEAREST);
		MONITORINFO mi = { sizeof(mi) };
		if (!GetMonitorInfo(hmon, &mi)) return 0;
		iWindowSizeX = mi.rcMonitor.right;
		iWindowSizeY = mi.rcMonitor.bottom;
		iWindowTopLeftX = 0;
		iWindowTopLeftY = 0;
	}

	RECT rWndRect = { 0, 0, iWindowSizeX, iWindowSizeY };					// Keep client size as requested
	AdjustWindowRectEx(&rWndRect, dwStyle, FALSE, dwExStyle);
	int width = rWndRect.right - rWndRect.left;
	int height = rWndRect.bottom - rWndRect.top;

	if (bEnableConsole) {
		RedirectIOToConsole();
	}

	hWnd = CreateWindowEx(dwExStyle, strT("ENGINE"), strT(""), dwStyle,
		iWindowTopLeftX, iWindowTopLeftY, width, height, NULL, NULL, GetModuleHandle(nullptr), NULL);

	return 1;
}

void CPlatform::RedirectIOToConsole() 
{
	AllocConsole();															// Create a console for this application

	HANDLE ConsoleOutput = GetStdHandle(STD_OUTPUT_HANDLE);					// Get STDOUT handle
	int SystemOutput = _open_osfhandle(intptr_t(ConsoleOutput), _O_TEXT);
	FILE* COutputHandle = _fdopen(SystemOutput, "w");

	HANDLE ConsoleError = GetStdHandle(STD_ERROR_HANDLE);					// Get STDERR handle
	int SystemError = _open_osfhandle(intptr_t(ConsoleError), _O_TEXT);
	FILE* CErrorHandle = _fdopen(SystemError, "w");

	HANDLE ConsoleInput = GetStdHandle(STD_INPUT_HANDLE);					// Get STDIN handle
	int SystemInput = _open_osfhandle(intptr_t(ConsoleInput), _O_TEXT);
	FILE* CInputHandle = _fdopen(SystemInput, "r");

	std::ios::sync_with_stdio(true);										// make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog point to console as well

	freopen_s(&CInputHandle, "CONIN$", "r", stdin);							// Redirect the CRT standard input, output, and error handles to the console
	freopen_s(&COutputHandle, "CONOUT$", "w", stdout);
	freopen_s(&CErrorHandle, "CONOUT$", "w", stderr);

	//Clear the error state for each of the C++ standard stream objects. We need to do this, as
	//attempts to access the standard streams before they refer to a valid target will cause the
	//iostream objects to enter an error state. In versions of Visual Studio after 2005, this seems
	//to always occur during startup regardless of whether anything has been read from or written to
	//the console or not.
	std::wcout.clear();
	std::cout.clear();
	std::wcerr.clear();
	std::cerr.clear();
	std::wcin.clear();
	std::cin.clear();
}

int CPlatform::CreateGraphics()
{
	if (RenderCreateDevice({ hWnd }) == 1) {
		RenderUpdateViewport();
		return 1;
	}
	else
		return 0;
}

void CPlatform::PrepareEngine()
{
	if (CreateGraphics() == 0) return;										// Start OpenGL, the context is owned by the game thread

	layer.Create(iScreenSizeX, iScreenSizeY);
	layer.bUpdate = true;
	layer.bShow = true;

	m_tp1 = std::chrono::system_clock::now();
	m_tp2 = std::chrono::system_clock::now();
}

void CPlatform::CoreUpdate()
{
	m_tp2 = std::chrono::system_clock::now();								// Handle Timing
	std::chrono::duration<float> elapsedTime = m_tp2 - m_tp1;
	m_tp1 = m_tp2;

	float fElapsedTime = elapsedTime.count();								// Our time per frame coefficient
	fLastElapsed = fElapsedTime;

	try {
		if (!OnUserUpdate(fElapsedTime)) bAtomActive = false;
	}
	catch (const std::exception& e) {
		std::cout << "CoreUpdate Exception : " << e.what() << std::endl;
	}

	RenderUpdateViewport();													// Display Frame
	RenderClearBuffer(BLACK, true);
	RenderPrepareDrawing();

	if (layer.bShow) {
		RenderApplyTexture(layer.id);
		if (layer.bUpdate) {
			layer.Update();
		}
		RenderDrawLayerQuad(layer.fOffsetX, layer.fOffsetY, layer.fScaleX, layer.fScaleY, layer.tint);
	}

	RenderDisplayFrame();													// Present Graphics to screen

	fFrameTimer += fElapsedTime;											// Update Title Bar
	nFrameCount++;
	if (fFrameTimer >= 1.0f) {
		nLastFPS = nFrameCount;
		fFrameTimer -= 1.0f;
		std::string sTitle = "Engine - " + sAppName + " - FPS: " + std::to_string(nFrameCount);
		#ifdef UNICODE
			SetWindowText(hWnd, ConvertS2W(sTitle).c_str());
		#else
			SetWindowText(hWnd, sTitle.c_str());
		#endif
		nFrameCount = 0;
	}

	Sleep(1);
}

void CPlatform::EngineThread()
{
	if (ThreadStartUp() == 0)	return;										// Allow platform to do stuff here if needed, since its now in the context of this thread
	try {
		PrepareEngine();													// Do engine context specific initialisation
	}
	catch (const std::exception& e) {
		std::cout << "EngineThread Exception : " << e.what() << std::endl;
	}

	if (!OnUserCreate()) bAtomActive = false;								// Create user resources as part of this thread

	while (bAtomActive) {
		while (bAtomActive) { CoreUpdate(); }								// Run as fast as possible
		if (!OnUserDestroy()) {												// Allow the user to free resources if they have overrided the destroy function
			bAtomActive = true;												// User denied destroy for some reason, so continue running
		}
	}
	ThreadCleanUp();
}

typedef BOOL(WINAPI wglSwapInterval_t) (int interval);
static wglSwapInterval_t* wglSwapInterval = nullptr;

int CPlatform::RenderCreateDevice(std::vector<void*> params)
{
	glDeviceContext = GetDC((HWND)(params[0]));								// Create Device Context
	PIXELFORMATDESCRIPTOR pfd =
	{
		sizeof(PIXELFORMATDESCRIPTOR), 1,
		PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,
		PFD_TYPE_RGBA, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		PFD_MAIN_PLANE, 0, 0, 0, 0
	};

	int pf = 0;
	if (!(pf = ChoosePixelFormat(glDeviceContext, &pfd))) return 0;
	SetPixelFormat(glDeviceContext, pf, &pfd);

	if (!(glRenderContext = wglCreateContext(glDeviceContext))) return 0;
	wglMakeCurrent(glDeviceContext, glRenderContext);

	wglSwapInterval = (wglSwapInterval_t*)wglGetProcAddress("wglSwapIntervalEXT"); 	// Remove Frame cap
	if (wglSwapInterval && !bEnableVSYNC) wglSwapInterval(0);

	glEnable(GL_TEXTURE_2D);												// Turn on texturing
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	return 1;
}

int CPlatform::RenderDestroyDevice()
{
	wglDeleteContext(glRenderContext);
	return 1;
}

void CPlatform::RenderUpdateViewport()
{
	glViewport(iViewPosX, iViewPosY, iViewSizeX, iViewSizeY);
}

void CPlatform::RenderClearBuffer(Pixel p, bool bDepth)
{
	glClearColor(float(p.r) / 255.0f, float(p.g) / 255.0f, float(p.b) / 255.0f, float(p.a) / 255.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	if (bDepth) glClear(GL_DEPTH_BUFFER_BIT);
}

void CPlatform::RenderDisplayFrame()
{
	SwapBuffers(glDeviceContext);
	if (b_glSync) DwmFlush();												// Woooohooooooo!!!! SMOOOOOOOTH!
}

void CPlatform::RenderPrepareDrawing()
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

uint32_t CPlatform::RenderCreateTexture(const uint32_t width, const uint32_t height, const bool filtered, const bool clamp)
{
	(void)(width);
	(void)(height);
	uint32_t id = 0;
	glGenTextures(1, &id);
	glBindTexture(GL_TEXTURE_2D, id);
	
	if (filtered) {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	}
	else {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	}

	if (clamp) {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	}
	else {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	}

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	return id;
}

void CPlatform::RenderApplyTexture(uint32_t id)
{
	glBindTexture(GL_TEXTURE_2D, id);
}

uint32_t CPlatform::RenderDeleteTexture(const uint32_t id)
{
	glDeleteTextures(1, &id);
	return id;
}

void CPlatform::RenderUpdateTexture(uint32_t id, int width, int height, const Pixel* pixels) 
{
	(void)(id);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
}

void CPlatform::RenderReadTexture(uint32_t id, int width, int height, Pixel* pixels)
{
	glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
}

void CPlatform::RenderDrawLayerQuad(const float offsetX, const float offsetY, const float scaleX, const float scaleY, const Pixel tint)
{
	glBegin(GL_QUADS);
	glColor4ub(tint.r, tint.g, tint.b, tint.a);
	glTexCoord2f(0.0f * scaleX + offsetX, 1.0f * scaleY + offsetY);
	glVertex3f(-1.0f, -1.0f, 0.0f);
	glTexCoord2f(0.0f * scaleX + offsetX, 0.0f * scaleY + offsetY);
	glVertex3f(-1.0f, 1.0f , 0.0f);
	glTexCoord2f(1.0f * scaleX + offsetX, 0.0f * scaleY + offsetY);
	glVertex3f(1.0f , 1.0f , 0.0f);
	glTexCoord2f(1.0f * scaleX + offsetX, 1.0f * scaleY + offsetY);
	glVertex3f(1.0f, -1.0f , 0.0f);
	glEnd();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#endif // PLATFORM_H
