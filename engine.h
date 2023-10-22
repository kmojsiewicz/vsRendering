#ifndef ENGINE_H
#define ENGINE_H

#if __cplusplus >= 201703L                                                      // code for C++17 and later
    #include <filesystem>
#else
    #ifdef _WIN32
        #include <direct.h>
        #define getcwd _getcwd                                                  // stupid MSFT "deprecation" warning
    #else
        #include <unistd.h>
    #endif
#endif
#include <vector>
#include <fstream>
#include <iostream>
#include "platform.h"

#define WND_WIDTH   800
#define WND_HEIGHT  600

static float z_buffer[WND_WIDTH][WND_HEIGHT];                               // z-buffer uses heap data

struct TTriangle;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class CEngine : public CPlatform 
{
public:
    static void ClrZBuffer(void);
    inline void PutPixel(int x, int y, float z, Pixel p, float light);
    bool Draw(int32_t x, int32_t y, Pixel p);
    void DrawLine(int32_t x1, int32_t y1, int32_t x2, int32_t y2, Pixel p = WHITE, uint32_t pattern = 0xFFFFFFFF);
    void DrawTriangle(TTriangle t);
    void FillTriangle(TTriangle t, Pixel p);
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

float q_rsqrt(float number)
{
    long i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = number * 0.5F;
    y = number;
    i = *(long*)&y;                                                         // evil floating point bit level hacking
    i = 0x5f3759df - (i >> 1);                                              // what the fuck?
    y = *(float*)&i;
    y = y * (threehalfs - (x2 * y * y));                                    // 1st iteration
    //y = y * ( threehalfs - ( x2 * y * y ) );                              // 2nd iteration, this can be removed

    return y;
}

static inline unsigned int BitCountByMask(unsigned int Mask)
{
    unsigned int BitCount = 0;
    while (Mask) {
        Mask &= Mask - 1;
        BitCount++;
    }
    return BitCount;
}

static inline unsigned int BitCountToMask(unsigned int BitCount) {
    return (BitCount == 32) ? 0xFFFFFFFF : (1 << BitCount) - 1;
}

static inline unsigned int BitPositionByMask(unsigned int Mask) {
    return BitCountByMask((Mask & (~Mask + 1)) - 1);
}

static inline unsigned int ComponentByMask(unsigned int Color, unsigned int Mask) {
    unsigned int Component = Color & Mask;
    return Component >> BitPositionByMask(Mask);
}

static unsigned int Convert(unsigned int Color, unsigned int FromBitCount, unsigned int ToBitCount) {
    if (ToBitCount < FromBitCount) {
        Color >>= (FromBitCount - ToBitCount);
    }
    else {
        Color <<= (ToBitCount - FromBitCount);
        if (Color > 0) {
            Color |= BitCountToMask(ToBitCount - FromBitCount);
        }
    }
    return Color;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class CTexture
{
private:
    static int id;
    int width;
    int height;
    std::vector<Pixel> pColData;

public:
    CTexture() {
        width = 0;
        height = 0;
    }
    ~CTexture() {
        pColData.clear();
    }

    void Create(uint32_t w, uint32_t h) 
    {
        width = w;
        height = h;
        pColData.resize(width * height, (0xFF << 24));
        if (pColData.data() != nullptr) {
            id++;
        }
    }

    void Clear(Pixel p) 
    {
        fill(pColData.begin(), pColData.end(), p);
    }

    Pixel GetPixel(int32_t x, int32_t y) const 
    {
        if (x >= 0 && x < width && y >= 0 && y < height)
            return pColData[y * width + x];
        else
            return Pixel(0, 0, 0, 0);
    }

    bool SetPixel(int32_t x, int32_t y, Pixel p) 
    {
        if (x >= 0 && x < width && y >= 0 && y < height) {
            pColData[y * width + x] = p;
            return true;
        }
        else
            return false;
    }

    void Draw(int x, int y) 
    {
        int i;
        int j;
        int offset;
        CLayer* pDrawTarget = engine->GetDrawTarget();

        offset = 0;
        for (j = 0; j < height; j++) {
            for (i = 0; i < width; i++) {
                pDrawTarget->SetPixel(i + x, j + y, pColData[i + offset]);
            }
            offset += width;
        }
    }

    int GetWidth()   { return width;  }
    int GetHeight()  { return height; }
    Pixel* GetData() { return pColData.data(); }
    int GetId()      { return id;  }

    bool LoadFromBitmap(std::string sfilename) 
    {
        #pragma pack(push,1)

        typedef struct _BITMAP_FILEHEADER {
            uint16_t Signature;
            uint32_t Size;
            uint32_t Reserved;
            uint32_t BitsOffset;
        } BITMAP_FILEHEADER;

        #define BITMAP_FILEHEADER_SIZE 14

        typedef struct _BITMAP_HEADER {
            uint32_t HeaderSize;
            int32_t Width;
            int32_t Height;
            uint16_t Planes;
            uint16_t BitCount;
            uint32_t Compression;
            uint32_t SizeImage;
            int32_t PelsPerMeterX;
            int32_t PelsPerMeterY;
            uint32_t ClrUsed;
            uint32_t ClrImportant;
            uint32_t RedMask;
            uint32_t GreenMask;
            uint32_t BlueMask;
            uint32_t AlphaMask;
            uint32_t CsType;
            uint32_t Endpoints[9];                                          // see http://msdn2.microsoft.com/en-us/library/ms536569.aspx
            uint32_t GammaRed;
            uint32_t GammaGreen;
            uint32_t GammaBlue;
        } BITMAP_HEADER;

        typedef struct _BGRA {
            uint8_t Blue;
            uint8_t Green;
            uint8_t Red;
            uint8_t Alpha;
        } BGRA;

        #pragma pack(pop)

        int x = 0;
        int y = 0;
        int offset;
        int paddedWidth;
        int index = 0;
        unsigned char c;
        std::vector<BGRA> ColorTable;
        unsigned int ColorTableSize = 0;
        Pixel pixel;
        BITMAP_FILEHEADER m_BitmapFileHeader;
        BITMAP_HEADER m_BitmapHeader;

        std::ifstream ifs(sfilename.c_str(), std::ifstream::in | std::ios::binary);

        if (!ifs.good()) {
            std::string sfileName2;

            #if __cplusplus >= 201703L                                      // code for C++17 and later
                sfileName2 = std::filesystem::current_path();
            #else                                                           // code for earlier versions of c++
                char buffer[FILENAME_MAX];
                if (getcwd(buffer, FILENAME_MAX) != nullptr) {
                    sfileName2 = buffer;
                }
            #endif

            sfileName2 += "\\BMP\\" + std::string(sfilename);
            ifs.open(sfileName2.c_str(), std::ifstream::in | std::ios::binary);

            if (!ifs.good()) {
                OutputDebugStringA(std::string("Error opening file " + sfilename).c_str());
                return false;
            }
        }
        ifs.read((char*)&m_BitmapFileHeader, BITMAP_FILEHEADER_SIZE);

        // We check whether the appropriate bitmap file can be displayed
        if (m_BitmapFileHeader.Signature != 19778 || m_BitmapFileHeader.Reserved != 0) {    // Not a valid bitmap file - don't display
            ifs.close();
            OutputDebugStringA("Not a valid bitmap file - don't display");
            return false;
        }

        ifs.read((char*)&m_BitmapHeader, sizeof(BITMAP_HEADER));

        if (m_BitmapHeader.Compression != 0 && m_BitmapHeader.Compression != 1
            && m_BitmapHeader.Compression != 3) {                                           // Compressed file - don't display
            ifs.close();
            OutputDebugStringA("Compressed file - don't display");
            return false;
        }

        if (m_BitmapHeader.BitCount == 1) {
            ColorTableSize = 2;
        }
        else if (m_BitmapHeader.BitCount == 4) {
            ColorTableSize = 16;
        }
        else if (m_BitmapHeader.BitCount == 8) {
            ColorTableSize = 256;
        }

        if (m_BitmapHeader.Compression == 1) {
            ColorTableSize = m_BitmapHeader.ClrUsed;
        }
        else if (m_BitmapHeader.Compression == 3) {
            ColorTableSize = m_BitmapHeader.ClrUsed;
        }

        if (ColorTableSize > 0) {
            try {
                ColorTable.resize(ColorTableSize);                              // std::bad_alloc exception should be thrown if memory is not available
            }
            catch (const std::exception& e) {
                (void)e;
                ifs.close();
                return false;
            }

            // Load Color Table
            ifs.seekg(BITMAP_FILEHEADER_SIZE + m_BitmapHeader.HeaderSize, std::ios::beg);
            ifs.read((char*)ColorTable.data(), sizeof(BGRA)* ColorTableSize);          // we read the color palette
        }

        ifs.seekg(m_BitmapFileHeader.BitsOffset, std::ios::beg);

        offset = m_BitmapHeader.Width * (m_BitmapHeader.Height - 1);            // Set appropriate position to display graphics data
        paddedWidth = m_BitmapHeader.Width & 0xFFFC;                            // Pad line length to 4bytes
        paddedWidth = m_BitmapHeader.Width & 0xFFFC;                            // Pad line length to 4bytes
        if (m_BitmapHeader.Width != paddedWidth)
            paddedWidth += 4;

        Create(m_BitmapHeader.Width, m_BitmapHeader.Height);
        y = 0;                                                                  // We've read no lines so far

        if (m_BitmapHeader.Compression == 0) {
            if (m_BitmapHeader.BitCount == 8) {

                while (y < m_BitmapHeader.Height) {                             // Decode and display graphics
                    for (x = 0; x < paddedWidth; x++) {                         // Read next line
                        c = ifs.get();
                        pixel.r = ColorTable[c].Red;
                        pixel.g = ColorTable[c].Green;
                        pixel.b = ColorTable[c].Blue;
                        pixel.a = 255;
                        pColData[offset + x] = pixel;
                    }
                    offset -= m_BitmapHeader.Width;                             // Move up one line on the screen
                    y++;                                                        // increase amount of lines read
                }
            }
            else {                                                              // Other than 8-bit colour
                ifs.seekg(m_BitmapFileHeader.BitsOffset, std::ios::beg);

                if (m_BitmapHeader.BitCount == 16) {
                    unsigned short col16bit;
                    while (y < m_BitmapHeader.Height) {                         // Decode and display graphics
                        for (x = 0; x < paddedWidth; x++) {                     // Read next line
                            col16bit = ifs.get();
                            col16bit |= (unsigned short)ifs.get() << 8;
                            pixel.b = (col16bit & 0x1f) << 3;
                            pixel.g = ((col16bit >> 5) & 0x1f) << 3;
                            pixel.r = ((col16bit >> 10) & 0x1f) << 3;
                            pixel.a = 255;
                            pColData[offset + x] = pixel;
                        }
                        offset -= m_BitmapHeader.Width;                         // Move up one line on the screen
                        y++;                                                    // increase amount of lines read
                    }
                }
                if (m_BitmapHeader.BitCount == 24) {
                    while (y < m_BitmapHeader.Height) {                         // Decode and display graphics
                        for (x = 0; x < paddedWidth; x++) {                     // Read next line
                            pixel.b = ifs.get();
                            pixel.g = ifs.get();
                            pixel.r = ifs.get();
                            pixel.a = 255;
                            pColData[offset + x] = pixel;
                        }
                        offset -= m_BitmapHeader.Width;                         // Move up one line on the screen
                        y++;                                                    // increase amount of lines read
                    }
                }
                else if (m_BitmapHeader.BitCount == 32) {
                    while (y < m_BitmapHeader.Height) {                         // Decode and display graphics
                        for (x = 0; x < paddedWidth; x++) {                     // Read next line
                            pixel.b = ifs.get();
                            pixel.g = ifs.get();
                            pixel.r = ifs.get();
                            pixel.a = ifs.get();
                            pColData[offset + x] = pixel;
                        }
                        offset -= m_BitmapHeader.Width;                         // Move up one line on the screen
                        y++;                                                    // increase amount of lines read
                    }
                }
            }
        }
        else if (m_BitmapHeader.Compression == 1) {                             // RLE 8
            uint8_t Count = 0;
            c = 0;

            while (ifs.eof() == false) {
                ifs.read((char*)&Count, sizeof(uint8_t));
                ifs.read((char*)&c, sizeof(uint8_t));

                if (Count > 0) {
                    offset = x + y * width;
                    for (int k = 0; k < Count; k++) {
                        pixel.r = ColorTable[c].Red;
                        pixel.g = ColorTable[c].Green;
                        pixel.b = ColorTable[c].Blue;
                        pixel.a = ColorTable[c].Alpha;
                        pColData[offset + k] = pixel;
                    }
                    x += Count;
                }
                else if (Count == 0) {
                    int Flag = c;
                    if (Flag == 0) {
                        x = 0;
                        y++;
                    }
                    else if (Flag == 1) {
                        break;
                    }
                    else if (Flag == 2) {
                        char rx = 0;
                        char ry = 0;
                        ifs.read((char*)&rx, sizeof(char));
                        ifs.read((char*)&ry, sizeof(char));
                        x += rx;
                        y += ry;
                    }
                    else {
                        Count = Flag;
                        offset = x + y * width;
                        for (int k = 0; k < Count; k++) {
                            ifs.read((char*)&c, sizeof(uint8_t));
                            pixel.r = ColorTable[c].Red;
                            pixel.g = ColorTable[c].Green;
                            pixel.b = ColorTable[c].Blue;
                            pixel.a = ColorTable[c].Alpha;
                            pColData[offset + k] = pixel;
                        }
                        x += Count;
                        // Attention: Current Microsoft STL implementation seems to be buggy, tellg() always returns 0.
                        if (ifs.tellg() & 1) {
                            ifs.seekg(1, std::ios::cur);
                        }
                    }
                }
            }
        }
        else if (m_BitmapHeader.Compression == 3) {
            /* We assumes that mask of each color component can be in any order */

            uint32_t BitCountRed = BitCountByMask(m_BitmapHeader.RedMask);
            uint32_t BitCountGreen = BitCountByMask(m_BitmapHeader.GreenMask);
            uint32_t BitCountBlue = BitCountByMask(m_BitmapHeader.BlueMask);
            uint32_t BitCountAlpha = BitCountByMask(m_BitmapHeader.AlphaMask);

            for (y = 0; y < height; y++) {

                for (x = 0; x < width; x++) {

                    uint32_t Color = 0;

                    if (m_BitmapHeader.BitCount == 16) {
                        Color = ifs.get();
                        Color |= (unsigned short)ifs.get() << 8;
                    }
                    else if (m_BitmapHeader.BitCount == 32) {
                        Color = ifs.get();
                        Color |= (unsigned short)ifs.get() << 8;
                        Color |= (unsigned short)ifs.get() << 16;
                        Color |= (unsigned short)ifs.get() << 24;
                    }
                    else {
                        // Other formats are not valid
                    }
                    pixel.r = Convert(ComponentByMask(Color, m_BitmapHeader.RedMask), BitCountRed, 8);
                    pixel.g = Convert(ComponentByMask(Color, m_BitmapHeader.GreenMask), BitCountGreen, 8);
                    pixel.b = Convert(ComponentByMask(Color, m_BitmapHeader.BlueMask), BitCountBlue, 8);
                    pixel.a = Convert(ComponentByMask(Color, m_BitmapHeader.AlphaMask), BitCountAlpha, 8);

                    pColData[offset + x] = pixel;
                }
                offset -= m_BitmapHeader.Width;
            }
        }

        ifs.close();
        return true;
    }
};

int CTexture::id = 0;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template <class T>
void swap_data(T& x, T& y)
{
    T temp;
    temp = x;
    x = y;
    y = temp;
}

struct TVec3d {
    float x;
    float y;
    float z;
};

struct TVertex {
    float x;
    float y;
    float z;
    float u;
    float v;
};

struct TTriangle {
    TVertex     V1;
    TVertex     V2;
    TVertex     V3;
    CTexture*   texture;
    float       light;

    void SetTexture(CTexture* t)
    {
        texture = t;
    }

    void Translate(float x, float y, float z) noexcept
    {
        V1.x += x; V2.x += x; V3.x += x;
        V1.y += y; V2.y += y; V3.y += y;
        V1.z += z; V2.z += z; V3.z += z;
    }

    void Scale(float x, float y, float z) noexcept
    {
        V1.x *= x; V2.x *= x; V3.x *= x;
        V1.y *= y; V2.y *= y; V3.y *= y;
        V1.z *= z; V2.z *= z; V3.z *= z;
    }

    TVec3d NormalVector() noexcept
    {
        TVec3d normal, line1, line2;

        line1.x = V2.x - V1.x;
        line1.y = V2.y - V1.y;
        line1.z = V2.z - V1.z;

        line2.x = V3.x - V1.x;
        line2.y = V3.y - V1.y;
        line2.z = V3.z - V1.z;

        normal.x = line1.y * line2.z - line1.z * line2.y;
        normal.y = line1.z * line2.x - line1.x * line2.z;
        normal.z = line1.x * line2.y - line1.y * line2.x;

        float inv_sqrt_nl = q_rsqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        normal.x *= inv_sqrt_nl;
        normal.y *= inv_sqrt_nl;
        normal.z *= inv_sqrt_nl;

        return normal;
    }
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

struct TMesh
{
    std::vector<TTriangle> triangles;

    void MakeQube(CTexture* texture, float light) 
    {
        triangles.clear();
        // FRONT
        triangles.push_back({ { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 0.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 0.0f }, texture, light });
        triangles.push_back({ { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 0.0f, 1.0f, 1.0f }, texture, light });
        // RIGHT
        triangles.push_back({ { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 0.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f }, texture, light });
        triangles.push_back({ { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 1.0f, 1.0f, 1.0f }, texture, light });
        // BACK
        triangles.push_back({ { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 1.0f, 1.0f, 1.0f, 0.0f }, texture, light });
        triangles.push_back({ { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 0.0f, 0.0f, 1.0f, 1.0f, 1.0f }, texture, light });
        // LEFT
        triangles.push_back({ { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 1.0f, 0.0f, 1.0f, 0.0f }, texture, light });
        triangles.push_back({ { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 0.0f, 1.0f, 0.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f }, texture, light });
        // TOP
        triangles.push_back({ { 0.0f, 1.0f, 0.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f }, texture, light });
        triangles.push_back({ { 0.0f, 1.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 1.0f }, texture, light });
        // BOTTOM
        triangles.push_back({ { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 0.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 0.0f }, texture, light });
        triangles.push_back({ { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 0.0f, 1.0f, 1.0f }, texture, light });
    }

    bool LoadFromObjectFile(std::string sfilename, CTexture *texture, float light) {

        std::vector<TVertex> verts;
        std::string line, key;
        std::ifstream ifs(sfilename.c_str(), std::ifstream::in);

        if (!ifs.good()) {
            std::string sfileName2;

            #if __cplusplus >= 201703L                                      // code for C++17 and later
                sfileName2 = std::filesystem::current_path();
            #else                                                           // code for earlier versions of c++
            char buffer[FILENAME_MAX];
            if (getcwd(buffer, FILENAME_MAX) != nullptr) {
                sfileName2 = buffer;
            }
            #endif

            sfileName2 += "\\OBJ\\" + std::string(sfilename);
            ifs.open(sfileName2.c_str(), std::ifstream::in);

            if (!ifs.good()) {
                OutputDebugStringA(std::string("Error opening file " + sfilename).c_str());
                return false;
            }
        }

        while (ifs.good() && !ifs.eof() && std::getline(ifs, line)) {
            key = "";
            std::strstream s;
            s << line;
            s >> key >> std::ws;

            if (key == "v") {                                               // vertex
                TVertex v;
                s >> v.x >> v.y >> v.z;
                v.u = 0.0;
                v.v = 0.0;
                verts.push_back(v);
            }
            else if (key == "vp") {                                         // parameter
            }
            else if (key == "vt") {                                         // texture coordinate
            }
            else if (key == "f") {                                          // face
                int v[3], t[3], n[3];
                while (!s.eof()) {
                    s >> v[0] >> std::ws;
                    if (s.peek() == '/') {
                        s.get();
                        if (s.peek() == '/') {
                            s.get();
                            s >> n[0] >> std::ws;
                        }
                        else {
                            s >> t[0] >> std::ws;
                            if (s.peek() == '/') {
                                s.get();
                                s >> n[0] >> std::ws;
                            }
                        }
                    }
                    while (!s.eof()) {
                        s >> v[1] >> std::ws;
                        if (s.peek() == '/') {
                            s.get();
                            if (s.peek() == '/') {
                                s.get();
                                s >> n[1] >> std::ws;
                            }
                            else {
                                s >> t[1] >> std::ws;
                                if (s.peek() == '/') {
                                    s.get();
                                    s >> n[1] >> std::ws;
                                }
                            }
                        }
                        while (!s.eof()) {
                            s >> v[2] >> std::ws;
                            if (s.peek() == '/') {
                                s.get();
                                if (s.peek() == '/') {
                                    s.get();
                                    s >> n[2] >> std::ws;
                                }
                                else {
                                    s >> t[2] >> std::ws;
                                    if (s.peek() == '/') {
                                        s.get();
                                        s >> n[2] >> std::ws;
                                    }
                                }
                            }
                            triangles.push_back({ verts[v[0] - 1], verts[v[1] - 1], verts[v[2] - 1], texture , light });
                        }
                    }

                }
            }
        }
        ifs.close();

        return true;
    }
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void CEngine::ClrZBuffer(void)
{
    int x, y;
    for (x = 0; x < WND_WIDTH; x++) {
        for (y = 0; y < WND_HEIGHT; y++) {
            z_buffer[x][y] = (float)INT_MAX;
        }
    }
}

inline void CEngine::PutPixel(int x, int y, float z, Pixel p, float light)
{
    if (x < 0) x = 0;
    if (x > (WND_WIDTH - 1)) x = WND_WIDTH - 1;
    if (y < 0) y = 0;
    if (y > (WND_HEIGHT - 1)) y = WND_HEIGHT - 1;

    //Pixel p(col.R * light, col.G * light, col.B * light);
    //Pixel p((uint8_t)(255.0 * light), (uint8_t)(255.0 * light), (uint8_t)(255.0 * light));
    CLayer* pDrawTarget = GetDrawTarget();

    //if (z_buffer[x][y] > z) {
    //    z_buffer[x][y] = z;
    pDrawTarget->SetPixel(x, y, p);
    //}
}

bool CEngine::Draw(int32_t x, int32_t y, Pixel p)
{
    CLayer* pDrawTarget = GetDrawTarget();
    pDrawTarget->SetPixel(x, y, p);
    return true;
}

void CEngine::DrawLine(int32_t x1, int32_t y1, int32_t x2, int32_t y2, Pixel p, uint32_t pattern)
{
    auto rol = [&](void) { pattern = (pattern << 1) | (pattern >> 31); return pattern & 1; };

    Pixel* pixels = GetDrawTarget()->GetData();

    if (y1 > y2) {                                                          // sorting point by their Y values
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    int screenx = ScreenWidth();
    int screeny = ScreenHeight();

    int dx21 = (x2 > x1) ? (x2 - x1) : (x1 - x2);
    int dy21 = -(y2 - y1);                                                  // dy21 will always be negative 
    int ix21;
    int ex;
    int e21 = dx21 + dy21;
    int ne21;

    if (x1 <= x2) {
        ix21 = 1;
        ex = screenx - 1;
    }
    else {
        ix21 = -1;
        ex = 0;
    }

    while (1) {
        if ((x1 >= 0) && (y1 >= 0) && (x1 < screenx) && (y1 < screeny) && rol()) {
            pixels[x1 + y1 * screenx] = p;
        }
        if ((x1 == x2) && (y1 == y2)) break;
        ne21 = 2 * e21;
        if (ne21 >= dy21) {
            e21 += dy21;
            if ((x1 == x2) || (x1 == ex)) break;
            x1 += ix21;
        }
        if (ne21 <= dx21) {
            e21 += dx21;
            if (y1 == y2) break;
            y1++;
            if (y1 == screeny) break;
        }
    }
}

void CEngine::DrawTriangle(TTriangle t)
{
    DrawLine((int32_t)t.V1.x, (int32_t)t.V1.y, (int32_t)t.V2.x, (int32_t)t.V2.y);
    DrawLine((int32_t)t.V2.x, (int32_t)t.V2.y, (int32_t)t.V3.x, (int32_t)t.V3.y);
    DrawLine((int32_t)t.V1.x, (int32_t)t.V1.y, (int32_t)t.V3.x, (int32_t)t.V3.y);
}

void CEngine::FillTriangle(TTriangle t, Pixel p)
{
    Pixel* pixels = GetDrawTarget()->GetData();
    int screenx = ScreenWidth();
    int screeny = ScreenHeight();

    auto drawline = [&](int sx, int ex, int ny) {
        if (ny < 0) return;
        if (ny >= screeny) return;
        if (sx > ex) std::swap(sx, ex);
        if (sx < 0) sx = 0;
        if (ex >= screenx) ex = screenx - 1;
        int offset = ny * screenx + sx;
        for (int i = sx; i <= ex; i++) {
            pixels[offset++] = p;
        }
    };

    int x1, y1, x2, y2, x3, y3;
    x1 = (int)t.V1.x;
    y1 = (int)t.V1.y;
    x2 = (int)t.V2.x;
    y2 = (int)t.V2.y;
    x3 = (int)t.V3.x;
    y3 = (int)t.V3.y;

    if (y1 > y2) {                                                          // sort the vertices (V1,V2,V3) by their Y values
        std::swap(x1, x2);
        std::swap(y1, y2);
    }
    if (y1 > y3) {
        std::swap(x1, x3);
        std::swap(y1, y3);
    }
    if (y2 > y3) {
        std::swap(x2, x3);
        std::swap(y2, y3);
    }

    int dx21 = (x2 > x1) ? (x2 - x1) : (x1 - x2);
    int dx31 = (x3 > x1) ? (x3 - x1) : (x1 - x3);
    int dx32 = (x3 > x2) ? (x3 - x2) : (x2 - x3);
    int dy21 = -(y2 - y1);                                                  // dy21 will always be negative 
    int dy31 = -(y3 - y1);                                                  // dy31 will always be negative 
    int dy32 = -(y3 - y2);                                                  // dy31 will always be negative 
    int ix21 = (x1 <= x2) ? 1 : -1;
    int ix31 = (x1 <= x3) ? 1 : -1;
    int ix32 = (x2 <= x3) ? 1 : -1;
    int e21 = dx21 + dy21;
    int e31 = dx31 + dy31;
    int e32 = dx32 + dy32;
    int ne21, ne31, ne32;
    int xa = x1;
    int xb = x1;
    int y = y1;

    while (1) {
        if ((xa == x2) && (y == y2)) break;
        ne21 = 2 * e21;
        if (ne21 >= dy21) {
            e21 += dy21;
            if (xa == x2) break;
            xa += ix21;
        }
        if (ne21 <= dx21) {
            e21 += dx21;
            
            while (1) {
                ne31 = 2 * e31;
                if (ne31 >= dy31) {
                    e31 += dy31;
                    if (xb == x3) break;
                    xb += ix31;
                }
                if (ne31 <= dx31) {
                    e31 += dx31;
                    break;
                }
            }

            drawline(xa, xb, y);                                            // Draw line from min to max points found on the y

            if (y == y2) break;
            y++;
            if (y == screeny) return;
        }
    }

    while (1) {
        if ((xb >= 0) && (y >= 0) && (xb < screenx) && (y < screeny)) {
            pixels[xb + y * screenx] = p;
        }
        if ((xb == x3) && (y == y3)) break;
        ne31 = 2 * e31;
        if (ne31 >= dy31) {
            e31 += dy31;
            if (xb == x3) break;
            xb += ix31;
        }
        if (ne31 <= dx31) {
            e31 += dx31;

            while (1) {
                ne32 = 2 * e32;
                if (ne32 >= dy32) {
                    e32 += dy32;
                    if (xa == x3) break;
                    xa += ix32;
                }
                if (ne32 <= dx32) {
                    e32 += dx32;
                    break;
                }
            }

            drawline(xa, xb, y);                                            // Draw line from min to max points found on the y

            if (y == y3) break;
            y++;
            if (y == screeny) break;
        }
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#endif // ENGINE_H
