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
    int x, y, xi, xe, ye;
    int screenx, screeny;
    int D21;
    int dx21, dx21tmp;
    int dy21, dy21tmp;

    screenx = ScreenWidth();
    screeny = ScreenHeight();

    if (y1 > y2) {                                                      // sorting point by their Y values
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    dx21 = (x2 > x1) ? (x2 - x1) : (x1 - x2);
    dy21 = (y2 - y1);                                                   // dy21 will always be positive 

    if (dx21 == 0) {                                                    // Line is vertical
        if (x1 < 0) return;
        if (x1 >= screenx) return;
        if (y1 < 0) y1 = 0;
        if (y2 >= screeny) y2 = screeny - 1;
        int offset = y1 * screenx + x1;
        for (y = y1; y <= y2; y++) {
            if (rol()) {
                pixels[offset] = p;
            }
            offset += screenx;
        }
        return;
    }

    if (dy21 == 0) {                                                    // Line is horizontal
        if (y1 < 0) return;
        if (y1 >= screeny) return;
        if (x2 < x1) std::swap(x1, x2);
        if (x1 < 0) x1 = 0;
        if (x2 >= screenx) x2 = screenx - 1;
        int offset = y1 * screenx + x1;
        for (x = x1; x <= x2; x++) {
            if (rol()) {
                pixels[offset] = p;
            }
            offset++;
        }
        return;
    }

    x = x1;
    y = y1;
    xe = x2;
    ye = y2;

    if (y2 < 0) return;                                                 // line off screen
    if (y1 >= screeny) return;                                          // line off screen
    if (x1 <= x2) {
        if (x2 < 0) return;                                             // line off screen
        if (x1 >= screenx) return;                                      // line off screen
        xi = 1;
        if (xe >= screenx) xe = screenx - 1;
    }
    else {
        if (x1 < 0) return;                                             // line off screen
        if (x2 >= screenx) return;                                      // line off screen
        xi = -1;
        if (xe < 0) xe = 0;
    }
    if (ye >= screeny) ye = screeny - 1;

    dx21tmp = dx21;
    dy21tmp = dy21;
    dx21 += dx21;                                                       // dx21 = 2 * dx21
    dy21 += dy21;                                                       // dy21 = 2 * dy21

    if (dy21 < dx21) {
        D21 = dy21 - dx21tmp;                                           // D21 = 2 * dy21 - dx21;
        for (; dx21tmp; dx21tmp--) {
            if ((x >= 0) && (y >= 0) && (x < screenx) && (y < screeny) && rol()) {
                pixels[x + y * screenx] = p;
            }
            if (x == xe) return;
            x += xi;
            if (D21 > 0) {
                if (y == ye) return;
                y++;
                D21 -= dx21;
            }
            D21 += dy21;
        }
    }
    else {
        D21 = dx21 - dy21tmp;                                           // D21 = 2 * dx21 - dy21;
        for (; dy21tmp; dy21tmp--) {
            if ((x >= 0) && (y >= 0) && (x < screenx) && (y < screeny) && rol()) {
                pixels[x + y * screenx] = p;
            }
            if (y == ye) return;
            y++;
            if (D21 > 0) {
                if (x == xe) return;
                x += xi;
                D21 -= dy21;
            }
            D21 += dx21;
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

    auto drawline = [&](int sx, int ex, int ny) {
        if (ny < 0) return;
        if (ny >= ScreenHeight()) return;
        if (sx < 0) sx = 0;
        if (ex >= ScreenWidth()) ex = ScreenWidth() - 1;
        int offset = ny * ScreenWidth() + sx;
        for (int i = sx; i <= ex; i++) {
            pixels[offset++] = p;
        }
    };

    int t1x, t2x, y, minx, maxx, t1xp, t2xp;
    bool changed1 = false;
    bool changed2 = false;
    int signx1, signx2, dx1, dy1, dx2, dy2;
    int e1, e2;

    if (t.V1.y > t.V2.y) { std::swap(t.V1, t.V2); }                         // Sort vertices
    if (t.V1.y > t.V3.y) { std::swap(t.V1, t.V3); }
    if (t.V2.y > t.V3.y) { std::swap(t.V2, t.V3); }

    int x1, y1, x2, y2, x3, y3;
    x1 = (int)t.V1.x; y1 = (int)t.V1.y;
    x2 = (int)t.V2.x; y2 = (int)t.V2.y;
    x3 = (int)t.V3.x; y3 = (int)t.V3.y;

    t1x = t2x = x1;                                                         // Starting points 
    y = y1;       
    dx1 = (int)(x2 - x1);
    if (dx1 < 0) { dx1 = -dx1; signx1 = -1; }
    else signx1 = 1;
    dy1 = (int)(y2 - y1);

    dx2 = (int)(x3 - x1);
    if (dx2 < 0) { dx2 = -dx2; signx2 = -1; }
    else signx2 = 1;
    dy2 = (int)(y3 - y1);

    if (dy1 > dx1) { std::swap(dx1, dy1); changed1 = true; }
    if (dy2 > dx2) { std::swap(dy2, dx2); changed2 = true; }

    e2 = (int)(dx2 >> 1);

    if (y1 == y2) goto next;                                                // Flat top, just process the second half
    e1 = (int)(dx1 >> 1);

    for (int i = 0; i < dx1;) {
        t1xp = 0; t2xp = 0;
        if (t1x < t2x) { minx = t1x; maxx = t2x; }
        else { minx = t2x; maxx = t1x; }

        while (i < dx1) {                                                   // process first line until y value is about to change
            i++;
            e1 += dy1;
            while (e1 >= dx1) {
                e1 -= dx1;
                if (changed1) t1xp = signx1;                                // t1x += signx1;
                else          goto next1;
            }
            if (changed1) break;
            else t1x += signx1;
        }

    next1:                                                                  // Move line
        while (1) {                                                         // process second line until y value is about to change
            e2 += dy2;
            while (e2 >= dx2) {
                e2 -= dx2;
                if (changed2) t2xp = signx2;                                // t2x += signx2;
                else          goto next2;
            }
            if (changed2)     break;
            else              t2x += signx2;
        }
    next2:
        if (minx > t1x) minx = t1x;
        if (minx > t2x) minx = t2x;
        if (maxx < t1x) maxx = t1x;
        if (maxx < t2x) maxx = t2x;
        drawline(minx, maxx, y);                                            // Draw line from min to max points found on the y

        if (!changed1) t1x += signx1;                                       // Now increase y
        t1x += t1xp;
        if (!changed2) t2x += signx2;
        t2x += t2xp;
        y += 1;
        if (y == y2) break;
    }
next:
    dx1 = (int)(x3 - x2); if (dx1 < 0) { dx1 = -dx1; signx1 = -1; }         // Second half
    else signx1 = 1;
    dy1 = (int)(y3 - y2);
    t1x = x2;

    if (dy1 > dx1) {                                                        // swap values
        std::swap(dy1, dx1);
        changed1 = true;
    }
    else changed1 = false;

    e1 = (int)(dx1 >> 1);
    for (int i = 0; i <= dx1; i++) {
        t1xp = 0; t2xp = 0;
        if (t1x < t2x) { minx = t1x; maxx = t2x; }
        else { minx = t2x; maxx = t1x; }
        while (i < dx1) {                                                   // process first line until y value is about to change
            e1 += dy1;
            while (e1 >= dx1) {
                e1 -= dx1;
                if (changed1) { t1xp = signx1; break; }                     // t1x += signx1;
                else          goto next3;
            }
            if (changed1) break;
            else   	   	  t1x += signx1;
            if (i < dx1) i++;
        }
    next3:
        while (t2x != x3) {                                                 // process second line until y value is about to change
            e2 += dy2;
            while (e2 >= dx2) {
                e2 -= dx2;
                if (changed2) t2xp = signx2;
                else          goto next4;
            }
            if (changed2)     break;
            else              t2x += signx2;
        }
    next4:
        if (minx > t1x) minx = t1x;
        if (minx > t2x) minx = t2x;
        if (maxx < t1x) maxx = t1x;
        if (maxx < t2x) maxx = t2x;
        drawline(minx, maxx, y);
        if (!changed1) t1x += signx1;
        t1x += t1xp;
        if (!changed2) t2x += signx2;
        t2x += t2xp;
        y += 1;
        if (y > y3) return;
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#endif // ENGINE_H
