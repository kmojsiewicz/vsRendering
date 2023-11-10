#ifndef ENGINE_H
#define ENGINE_H

#if __cplusplus >= 201703L                                                  // code for C++17 and later
    #include <filesystem>
#else
    #ifdef _WIN32
        #include <direct.h>
        #define getcwd _getcwd                                              // stupid MSFT "deprecation" warning
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

static float z_buffer[WND_WIDTH * WND_HEIGHT];                                // z-buffer uses heap data

struct TTriangle;
class CTexture;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class CEngine : public CPlatform 
{
public:
    static void ClrZBuffer(void);
    inline void PutPixel(int x, int y, float z, Pixel pixel, float light);
    bool Draw(int32_t x, int32_t y, Pixel pixel);
    void DrawLine(int32_t x1, int32_t y1, int32_t x2, int32_t y2, Pixel pixel = WHITE, uint32_t pattern = 0xFFFFFFFF);
    void DrawTriangle(TTriangle* t, Pixel pixel);
    void TexturedTriangle(TTriangle* t);
    void ColoredTriangle(TTriangle* t);
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

    Pixel GetPixel(int32_t offset) const
    {
        if (offset >= 0 && offset < pColData.size())
            return pColData[offset];
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

struct TVec2d {
    float u = 0.0f;
    float v = 0.0f;
    float w = 1.0f;
};

struct TVec3d {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    float w = 1.0f;
};

TVec3d Vector_Add(TVec3d& v1, TVec3d& v2)
{
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

TVec3d Vector_Sub(TVec3d& v1, TVec3d& v2)
{
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

TVec3d Vector_Mul(TVec3d& v, float k)
{
    return { v.x * k, v.y * k, v.z * k };
}

TVec3d Vector_Div(TVec3d& v, float k)
{
    return { v.x / k, v.y / k, v.z / k };
}

float Vector_DotProduct(TVec3d& v1, TVec3d& v2)
{
    return  v1.x * v2.x + v1.y * v2.y + v1.z * v2.z ;
}

float Vector_Length(TVec3d& v)
{
    return sqrtf(Vector_DotProduct(v, v));
}

TVec3d Vector_Normalise(TVec3d& v)
{
    float inv_sqrt_l = q_rsqrt(Vector_DotProduct(v, v));
    return { v.x * inv_sqrt_l, v.y * inv_sqrt_l, v.z * inv_sqrt_l };
}

TVec3d Vector_CrossProduct(TVec3d& v1, TVec3d& v2)
{
    TVec3d v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}

TVec3d Vector_IntersectPlane(TVec3d& plane_p, TVec3d& plane_n, TVec3d& lineStart, TVec3d& lineEnd, float& t)
{
    plane_n = Vector_Normalise(plane_n);
    float plane_d = -Vector_DotProduct(plane_n, plane_p);
    float ad = Vector_DotProduct(lineStart, plane_n);
    float bd = Vector_DotProduct(lineEnd, plane_n);
    t = (-plane_d - ad) / (bd - ad);
    TVec3d lineStartToEnd = Vector_Sub(lineEnd, lineStart);
    TVec3d lineToIntersect = Vector_Mul(lineStartToEnd, t);
    return Vector_Add(lineStart, lineToIntersect);
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

struct TMat4x4 {
    float m[4][4] = { 0 };
};

TVec3d Matrix_MultiplyVector(TMat4x4& m, TVec3d& i)
{
    TVec3d v;
    v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
    v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
    v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
    v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];

    return v;
}

TMat4x4 Matrix_MakeIdentity()
{
    TMat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    return matrix;
}

TMat4x4 Matrix_MakeRotationX(float fAngleRad)
{
    TMat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = cosf(fAngleRad);
    matrix.m[1][2] = sinf(fAngleRad);
    matrix.m[2][1] = -sinf(fAngleRad);
    matrix.m[2][2] = cosf(fAngleRad);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

TMat4x4 Matrix_MakeRotationY(float fAngleRad)
{
    TMat4x4 matrix;
    matrix.m[0][0] = cosf(fAngleRad);
    matrix.m[0][2] = sinf(fAngleRad);
    matrix.m[2][0] = -sinf(fAngleRad);
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = cosf(fAngleRad);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

TMat4x4 Matrix_MakeRotationZ(float fAngleRad)
{
    TMat4x4 matrix;
    matrix.m[0][0] = cosf(fAngleRad);
    matrix.m[0][1] = sinf(fAngleRad);
    matrix.m[1][0] = -sinf(fAngleRad);
    matrix.m[1][1] = cosf(fAngleRad);
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    return matrix;
}

TMat4x4 Matrix_MakeTranslation(float x, float y, float z)
{
    TMat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    matrix.m[3][0] = x;
    matrix.m[3][1] = y;
    matrix.m[3][2] = z;
    return matrix;
}

TMat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar, float fScale = 1.0f)
{
    float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
    TMat4x4 matrix;
    matrix.m[0][0] = fScale * fAspectRatio * fFovRad;
    matrix.m[1][1] = fScale * fFovRad;
    matrix.m[2][2] = fScale * fFar / (fFar - fNear);
    matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
    matrix.m[2][3] = 1.0f;
    matrix.m[3][3] = 0.0f;
    return matrix;
}

TMat4x4 Matrix_MultiplyMatrix(TMat4x4& m1, TMat4x4& m2)
{
    TMat4x4 matrix;
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
    return matrix;
}

TMat4x4 Matrix_PointAt(TVec3d& pos, TVec3d& target, TVec3d& up)
{
    TVec3d newForward = Vector_Sub(target, pos);                            // Calculate new forward direction
    newForward = Vector_Normalise(newForward);

    TVec3d a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));   // Calculate new Up direction
    TVec3d newUp = Vector_Sub(up, a);
    newUp = Vector_Normalise(newUp);

    TVec3d newRight = Vector_CrossProduct(newUp, newForward);               // New Right direction is easy, its just cross product

    TMat4x4 matrix;                                                         // Construct Dimensioning and Translation Matrix	
    matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
    return matrix;
}

TMat4x4 Matrix_QuickInverse(TMat4x4& m)                                     // Only for Rotation/Translation Matrices
{
    TMat4x4 matrix;
    matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
    matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
    matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

struct TVertex {
    TVec3d p;
    TVec2d t;
    Pixel pixel;
};

struct TTriangle {
    TVertex     V1;
    TVertex     V2;
    TVertex     V3;
    CTexture*   texture = nullptr;
    float       light;

    void SetTexture(CTexture* t)
    {
        texture = t;
    }

    void Translate(TVec3d v)
    {
        V1.p = Vector_Add(V1.p, v);
        V2.p = Vector_Add(V2.p, v);
        V3.p = Vector_Add(V3.p, v);
    }

    void Scale(float x, float y, float z)
    {
        V1.p.x *= x; V2.p.x *= x; V3.p.x *= x;
        V1.p.y *= y; V2.p.y *= y; V3.p.y *= y;
        V1.p.z *= z; V2.p.z *= z; V3.p.z *= z;
    }

    TVec3d NormalVector()
    {
        TVec3d line1 = Vector_Sub(V2.p, V1.p);
        TVec3d line2 = Vector_Sub(V3.p, V1.p);
        TVec3d normal = Vector_CrossProduct(line1, line2);
        return Vector_Normalise(normal);
    }

    int Triangle_ClipAgainstPlane(TVec3d plane_p, TVec3d plane_n, TTriangle& out_tri1, TTriangle& out_tri2)
    {
        plane_n = Vector_Normalise(plane_n);                                // Make sure plane normal is indeed normal

        auto dist = [&](TVec3d& p) {                                        // Return signed shortest distance from point to plane, plane normal must be normalised
            TVec3d n = Vector_Normalise(p);
            return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
            };

        TVec3d* inside_points[3];  int nInsidePointCount = 0;               // Create two temporary storage arrays to classify points either side of plane
        TVec3d* outside_points[3]; int nOutsidePointCount = 0;              // If distance sign is positive, point lies on "inside" of plane
        TVec2d* inside_tex[3];  int nInsideTexCount = 0;
        TVec2d* outside_tex[3]; int nOutsideTexCount = 0;

        float d0 = dist(V1.p);                                       // Get signed distance of each point in triangle to plane
        float d1 = dist(V2.p);
        float d2 = dist(V3.p);
        if (d0 >= 0) {
            inside_points[nInsidePointCount++] = &V1.p;
            inside_tex[nInsideTexCount++] = &V1.t;
        }
        else {
            outside_points[nOutsidePointCount++] = &V1.p;
            outside_tex[nOutsideTexCount++] = &V1.t;
        }
        if (d1 >= 0) {
            inside_points[nInsidePointCount++] = &V2.p;
            inside_tex[nInsideTexCount++] = &V2.t;
        }
        else {
            outside_points[nOutsidePointCount++] = &V2.p;
            outside_tex[nOutsideTexCount++] = &V2.t;
        }
        if (d2 >= 0) {
            inside_points[nInsidePointCount++] = &V3.p;
            inside_tex[nInsideTexCount++] = &V3.t;
        }
        else {
            outside_points[nOutsidePointCount++] = &V3.p;
            outside_tex[nOutsideTexCount++] = &V3.t;
        }

        // Now classify triangle points, and break the input triangle into 
        // smaller output triangles if required. There are four possible outcomes...
        if (nInsidePointCount == 0) {                                       // All points lie on the outside of plane, so clip whole triangle
            return 0;                                                       // No returned triangles are valid
        }

        if (nInsidePointCount == 3) {                                       // All points lie on the inside of plane, so do nothing
            out_tri1 = *this;                                               // and allow the triangle to simply pass through
            return 1;                                                       // Just the one returned original triangle is valid
        }

        if (nInsidePointCount == 1 && nOutsidePointCount == 2) {            // Triangle should be clipped. As two points lie outside the plane, the triangle simply becomes a smaller triangle
            out_tri1.light = light;                                         // Copy appearance info to new triangle
            out_tri1.texture = texture;
            out_tri1.V1.p = *inside_points[0];                              // The inside point is valid, so keep that...
            out_tri1.V1.t = *inside_tex[0];
            // but the two new points are at the locations where the 
            // original sides of the triangle (lines) intersect with the plane
            float t;
            out_tri1.V2.p = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
            out_tri1.V2.t.u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
            out_tri1.V2.t.v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
            out_tri1.V2.t.w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;
            out_tri1.V3.p = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1], t);
            out_tri1.V3.t.u = t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
            out_tri1.V3.t.v = t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;
            out_tri1.V3.t.w = t * (outside_tex[1]->w - inside_tex[0]->w) + inside_tex[0]->w;

            out_tri1.V1.pixel = V1.pixel;                                   // debugging out_tri1.V1.pixel = RED;
            out_tri1.V2.pixel = V2.pixel;                                   // debugging out_tri1.V2.pixel = RED;
            out_tri1.V3.pixel = V3.pixel;                                   // debugging out_tri1.V3.pixel = RED;
            return 1;                                                       // Return the newly formed single triangle
        }

        if (nInsidePointCount == 2 && nOutsidePointCount == 1) {            // Triangle should be clipped. As two points lie inside the plane, the clipped triangle becomes a "quad". Fortunately, we can
            out_tri1.light = light;                                         // represent a quad with two new triangles
            out_tri2.light = light;                                         // Copy appearance info to new triangles
            out_tri1.texture = texture;
            out_tri2.texture = texture;
            // The first triangle consists of the two inside points and a new
            // point determined by the location where one side of the triangle
            // intersects with the plane
            float t;
            out_tri1.V1.p = *inside_points[0];
            out_tri1.V2.p = *inside_points[1];
            out_tri1.V1.t = *inside_tex[0];
            out_tri1.V2.t = *inside_tex[1];
            out_tri1.V3.p = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
            out_tri1.V3.t.u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
            out_tri1.V3.t.v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
            out_tri1.V3.t.w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;
            // The second triangle is composed of one of he inside points, a
            // new point determined by the intersection of the other side of the 
            // triangle and the plane, and the newly created point above
            out_tri2.V1.p = *inside_points[1];
            out_tri2.V1.t = *inside_tex[1];
            out_tri2.V2.p = out_tri1.V3.p;
            out_tri2.V2.t = out_tri1.V3.t;
            out_tri2.V3.p = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0], t);
            out_tri2.V3.t.u = t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
            out_tri2.V3.t.v = t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;
            out_tri2.V3.t.w = t * (outside_tex[0]->w - inside_tex[1]->w) + inside_tex[1]->w;

            out_tri1.V1.pixel = V1.pixel;                                   // debugging out_tri1.V1.pixel = GREEN;
            out_tri1.V2.pixel = V2.pixel;                                   // debugging out_tri1.V2.pixel = GREEN;
            out_tri1.V3.pixel = V3.pixel;                                   // debugging out_tri1.V3.pixel = GREEN;
            out_tri2.V1.pixel = V1.pixel;                                   // debugging out_tri2.V1.pixel = BLUE;
            out_tri2.V2.pixel = V2.pixel;                                   // debugging out_tri2.V2.pixel = BLUE;
            out_tri2.V3.pixel = V3.pixel;                                   // debugging out_tri2.V3.pixel = BLUE;
            return 2;                                                       // Return two newly formed triangles which form a quad
        }

        return 0;
    }
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

struct TMesh
{
    std::vector<TTriangle> triangles;

    void MakeQube(CTexture* texture, Pixel defPixel, float light) 
    {
        triangles.clear();
        // FRONT
        triangles.push_back({ { { 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 0.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f }, defPixel },    { { 1.0f, 1.0f, 0.0f }, { 1.0f, 0.0f, 1.0f }, defPixel }, texture, light });
        triangles.push_back({ { { 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 1.0f, 1.0f, 0.0f }, { 1.0f, 0.0f, 1.0f }, defPixel },    { { 1.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f }, defPixel }, texture, light });
        // RIGHT
        triangles.push_back({ { { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 1.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f }, defPixel },    { { 1.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f }, defPixel }, texture, light });
        triangles.push_back({ { { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 1.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f }, defPixel },    { { 1.0f, 0.0f, 1.0f }, { 1.0f, 1.0f, 1.0f }, defPixel }, texture, light });
        // BACK
        triangles.push_back({ { { 1.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 1.0f, 1.0f, 1.0f }, { 0.0f, 0.0f, 1.0f }, defPixel },    { { 0.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f }, defPixel }, texture, light });
        triangles.push_back({ { { 1.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 0.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f }, defPixel },    { { 0.0f, 0.0f, 1.0f }, { 1.0f, 1.0f, 1.0f }, defPixel }, texture, light });
        // LEFT
        triangles.push_back({ { { 0.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 0.0f, 1.0f, 1.0f }, { 0.0f, 0.0f, 1.0f }, defPixel },    { { 0.0f, 1.0f, 0.0f }, { 1.0f, 0.0f, 1.0f }, defPixel }, texture, light });
        triangles.push_back({ { { 0.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 0.0f, 1.0f, 0.0f }, { 1.0f, 0.0f, 1.0f }, defPixel },    { { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f }, defPixel }, texture, light });
        // TOP
        triangles.push_back({ { { 0.0f, 1.0f, 0.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 0.0f, 1.0f, 1.0f }, { 0.0f, 0.0f, 1.0f } , defPixel},    { { 1.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f }, defPixel }, texture, light });
        triangles.push_back({ { { 0.0f, 1.0f, 0.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 1.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f }, defPixel },    { { 1.0f, 1.0f, 0.0f }, { 1.0f, 1.0f, 1.0f }, defPixel }, texture, light });
        // BOTTOM
        triangles.push_back({ { { 1.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 0.0f, 0.0f, 1.0f }, { 0.0f, 0.0f, 1.0f }, defPixel },    { { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 1.0f }, defPixel }, texture, light });
        triangles.push_back({ { { 1.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 1.0f }, defPixel },    { { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 1.0f }, defPixel },    { { 1.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f }, defPixel }, texture, light });
    }

    bool LoadFromObjectFile(std::string sfilename, CTexture *texture, Pixel defPixel, float light) {

        bool isQuad;
        std::vector<TVec3d> verts;
        std::vector<TVec2d> texs;
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
                TVec3d v;
                s >> v.x >> v.y >> v.z;
                verts.push_back(v);
            }
            else if (key == "vp") {                                         // parameter
            }
            else if (key == "vt") {                                         // texture coordinate
                TVec2d t;
                s >> t.u >> t.v;
                texs.push_back(t);
            }
            else if (key == "f") {                                          // face
                isQuad = false;
                int v[4], t[4] = { 0 }, n[4] = { 0 };
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
                            while (!s.eof()) {
                                s >> v[3] >> std::ws;
                                isQuad = true;
                                if (s.peek() == '/') {
                                    s.get();
                                    if (s.peek() == '/') {
                                        s.get();
                                        s >> n[3] >> std::ws;
                                    }
                                    else {
                                        s >> t[3] >> std::ws;
                                        if (s.peek() == '/') {
                                            s.get();
                                            s >> n[3] >> std::ws;
                                        }
                                    }
                                }
                            }
                            int v0_index, v1_index, v2_index, v3_index;
                            int t0_index = 0, t1_index = 0, t2_index = 0, t3_index = 0;

                            v0_index = v[0] - 1;
                            v1_index = v[1] - 1;
                            v2_index = v[2] - 1;
                            if (v[0] < 0) v0_index = (int)verts.size() + v[0];
                            if (v[1] < 0) v1_index = (int)verts.size() + v[1];
                            if (v[2] < 0) v2_index = (int)verts.size() + v[2];

                            if ((v0_index > 0) && (v1_index > 0) && (v2_index > 0)) {
                                TVertex v1, v2, v3, v4;
                                v1.p = verts[v0_index]; v1.pixel = defPixel;
                                v2.p = verts[v1_index]; v2.pixel = defPixel;
                                v3.p = verts[v2_index]; v3.pixel = defPixel;
                                if ((t[0] != 0) && (t[1] != 0) && (t[2] != 0)) {
                                    t0_index = t[0] - 1;
                                    t1_index = t[1] - 1;
                                    t2_index = t[2] - 1;
                                    if (t[0] < 0) t0_index = (int)texs.size() + t[0];
                                    if (t[1] < 0) t1_index = (int)texs.size() + t[1];
                                    if (t[2] < 0) t2_index = (int)texs.size() + t[2];
                                    v1.t = texs[t0_index];
                                    v2.t = texs[t1_index];
                                    v3.t = texs[t2_index];
                                }

                                triangles.push_back({ v1, v2, v3, texture , light });

                                if (isQuad) {
                                    v3_index = v[3] - 1;
                                    if (v[3] < 0) v3_index = (int)verts.size() + v[3];
                                    if ((v3_index > 0)) {
                                        v4.p = verts[v3_index]; v4.pixel = defPixel;
                                        if ((t[0] != 0) && (t[1] != 0) && (t[2] != 0) && (t[3] != 0)) {
                                            t3_index = t[3] - 1;
                                            if (t[3] < 0) t3_index = (int)texs.size() + t[3];
                                            v4.t = texs[t3_index];
                                        }
                                        triangles.push_back({ v1, v3, v4, texture, light });
                                    }
                                }
                            }
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
    int i;
    for (i = 0; i < WND_WIDTH * WND_HEIGHT; i++) {
        z_buffer[i] = 0.0f;
    }
}

inline void CEngine::PutPixel(int x, int y, float z, Pixel pixel, float light)
{
    if (x < 0) x = 0;
    if (x > (WND_WIDTH - 1)) x = WND_WIDTH - 1;
    if (y < 0) y = 0;
    if (y > (WND_HEIGHT - 1)) y = WND_HEIGHT - 1;

    CLayer* pDrawTarget = GetDrawTarget();
    pDrawTarget->SetPixel(x, y, pixel * light);
}

bool CEngine::Draw(int32_t x, int32_t y, Pixel pixel)
{
    CLayer* pDrawTarget = GetDrawTarget();
    pDrawTarget->SetPixel(x, y, pixel);
    return true;
}

void CEngine::DrawLine(int32_t x1, int32_t y1, int32_t x2, int32_t y2, Pixel pixel, uint32_t pattern)
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
            pixels[x1 + y1 * screenx] = pixel;
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

void CEngine::DrawTriangle(TTriangle *t, Pixel pixel)
{
    DrawLine((int32_t)t->V1.p.x, (int32_t)t->V1.p.y, (int32_t)t->V2.p.x, (int32_t)t->V2.p.y, pixel);
    DrawLine((int32_t)t->V2.p.x, (int32_t)t->V2.p.y, (int32_t)t->V3.p.x, (int32_t)t->V3.p.y, pixel);
    DrawLine((int32_t)t->V1.p.x, (int32_t)t->V1.p.y, (int32_t)t->V3.p.x, (int32_t)t->V3.p.y, pixel);
}

void CEngine::TexturedTriangle(TTriangle* t)
{
    if (t->texture == nullptr) return;

    Pixel* pixels = GetDrawTarget()->GetData();
    CTexture* tex = t->texture;

    struct intVertex {
        int x, y;
        int r, g, b;
        float u,v,w;
    } v1, v2, v3;

    v1.x = (int)(t->V1.p.x); v2.x = (int)(t->V2.p.x); v3.x = (int)(t->V3.p.x);
    v1.y = (int)(t->V1.p.y); v2.y = (int)(t->V2.p.y); v3.y = (int)(t->V3.p.y);
    v1.w = (t->V1.t.w); v2.w = (t->V2.t.w); v3.w = (t->V3.t.w);
    v1.u = ((t->V1.t.u * (t->texture->GetWidth() - 1)));
    v1.v = ((t->V1.t.v * (t->texture->GetHeight() - 1)));
    v2.u = ((t->V2.t.u * (t->texture->GetWidth() - 1)));
    v2.v = ((t->V2.t.v * (t->texture->GetHeight() - 1)));
    v3.u = ((t->V3.t.u * (t->texture->GetWidth() - 1)));
    v3.v = ((t->V3.t.v * (t->texture->GetHeight() - 1)));

    auto lineTextured = [&](int ax, int bx, int y, float tex_su, float tex_eu, float tex_sv, float tex_ev, float tex_sw, float tex_ew) 
    {
        if (ax > bx) {
            std::swap(ax, bx);
            std::swap(tex_su, tex_eu);
            std::swap(tex_sv, tex_ev);
            std::swap(tex_sw, tex_ew);
        }
        
        int dtx_div = (bx - ax);
        float tex_u = tex_su;
        float tex_v = tex_sv;
        float tex_w = tex_sw;
        float dtx_u = (tex_eu - tex_su) / (float)dtx_div;
        float dtx_v = (tex_ev - tex_sv) / (float)dtx_div;
        float dtx_w = (tex_ew - tex_sw) / (float)dtx_div;

        int offset = y * ScreenWidth() + ax;
        for (int j = ax; j <= bx; j++) {
            if (tex_w > z_buffer[offset]) {
                z_buffer[offset] = tex_w;
                pixels[offset] = tex->GetPixel((int)(tex_u / tex_w) + (int)(tex_v / tex_w) * t->texture->GetWidth()) * t->light;
            }
            offset++;
            tex_u += dtx_u;
            tex_v += dtx_v;
            tex_w += dtx_w;
        }
    };

    if (v1.y > v2.y) std::swap(v1, v2);                                     // sort the vertices (v1,v2,v3) by their Y values
    if (v1.y > v3.y) std::swap(v1, v3);
    if (v2.y > v3.y) std::swap(v2, v3);

    float tex_su = v1.u;
    float tex_sv = v1.v;
    float tex_sw = v1.w;
    float tex_eu = v1.u;
    float tex_ev = v1.v;
    float tex_ew = v1.w;

    float du1_step = 0, dv1_step = 0,
        du2_step = 0, dv2_step = 0,
        dw1_step = 0, dw2_step = 0;
    int xA, xB;
    int dxA_mul = 0, dxB_mul = 0;
    int dxA_div = 1, dxB_div = 1;

    int dy1 = abs(v2.y - v1.y);
    int dy2 = abs(v3.y - v1.y);

    if (dy1) {
        dxA_mul = (v2.x - v1.x);
        dxA_div = dy1;
        du1_step = (v2.u - v1.u) / (float)dy1;
        dv1_step = (v2.v - v1.v) / (float)dy1;
        dw1_step = (v2.w - v1.w) / (float)dy1;
    }
    if (dy2) {
        dxB_mul = (v3.x - v1.x);
        dxB_div = dy2;
        du2_step = (v3.u - v1.u) / (float)dy2;
        dv2_step = (v3.v - v1.v) / (float)dy2;
        dw2_step = (v3.w - v1.w) / (float)dy2;
    }

    xB = v1.x * dxB_div;
    if (dy1) {
        xA = v1.x * dxA_div;
        for (int y = v1.y; y <= v2.y; y++) {
            lineTextured(xA / dxA_div, xB / dxB_div, y, tex_su, tex_eu, tex_sv, tex_ev, tex_sw, tex_ew);
            xA += dxA_mul;
            xB += dxB_mul;
            tex_su += du1_step; tex_sv += dv1_step; tex_sw += dw1_step;
            tex_eu += du2_step; tex_ev += dv2_step; tex_ew += dw2_step;
        }
        xB -= dxB_mul;
        tex_eu -= du2_step; tex_ev -= dv2_step; tex_ew -= dw2_step;
    }

    du1_step = 0, 
    dv1_step = 0;
    dy1 = abs(v3.y - v2.y);
    if (dy1) {
        dxA_mul = (v3.x - v2.x);
        dxA_div = dy1;
        du1_step = (v3.u - v2.u) / (float)(dy1);
        dv1_step = (v3.v - v2.v) / (float)(dy1);
        dw1_step = (v3.w - v2.w) / (float)(dy1);
    }
    else { 
        dxA_mul = 0;
        dxA_div = 1;
    }

    if (dy1) {
        xA = v2.x * dxA_div;
        tex_su = v2.u;
        tex_sv = v2.v;
        tex_sw = v2.w;
        for (int y = v2.y; y <= v3.y; y++) {
            lineTextured(xA / dxA_div, xB / dxB_div, y, tex_su, tex_eu, tex_sv, tex_ev, tex_sw, tex_ew);
            xA += dxA_mul;
            xB += dxB_mul;
            tex_su += du1_step; tex_sv += dv1_step; tex_sw += dw1_step;
            tex_eu += du2_step; tex_ev += dv2_step; tex_ew += dw2_step;
        }
    }
}

void CEngine::ColoredTriangle(TTriangle* t)
{
    struct intVertex {
        int x, y;
        int r, g, b;
    } v1, v2, v3;

    v1.r = (int)(t->V1.pixel.r * t->light); v1.g = (int)(t->V1.pixel.g * t->light); v1.b = (int)(t->V1.pixel.b * t->light);
    v2.r = (int)(t->V2.pixel.r * t->light); v2.g = (int)(t->V2.pixel.g * t->light); v2.b = (int)(t->V2.pixel.b * t->light);
    v3.r = (int)(t->V3.pixel.r * t->light); v3.g = (int)(t->V3.pixel.g * t->light); v3.b = (int)(t->V3.pixel.b * t->light);
    v1.x = (int)(t->V1.p.x);
    v2.x = (int)(t->V2.p.x);
    v3.x = (int)(t->V3.p.x);
    v1.y = (int)(t->V1.p.y);
    v2.y = (int)(t->V2.p.y);
    v3.y = (int)(t->V3.p.y);

    if (v1.y > v2.y) std::swap(v1, v2);                                     // sort the vertices (v1,v2,v3) by their Y values
    if (v1.y > v3.y) std::swap(v1, v3);
    if (v2.y > v3.y) std::swap(v2, v3);

    Pixel* pixels = GetDrawTarget()->GetData();
    int screenx = ScreenWidth();
    int screeny = ScreenHeight();

    int ix21 = (v1.x <= v2.x) ? 1 : -1;
    int ix31 = (v1.x <= v3.x) ? 1 : -1;
    int ix32 = (v2.x <= v3.x) ? 1 : -1;
    int dx21 = (v1.x <= v2.x) ? (v2.x - v1.x) : (v1.x - v2.x);
    int dx31 = (v1.x <= v3.x) ? (v3.x - v1.x) : (v1.x - v3.x);
    int dx32 = (v2.x <= v3.x) ? (v3.x - v2.x) : (v2.x - v3.x);
    int dy21 = (v2.y - v1.y);                                              // dy21 will always be positive 
    int dy31 = (v3.y - v1.y);                                              // dy31 will always be positive 
    int dy32 = (v3.y - v2.y);                                              // dy32 will always be positive 

    int dRdX_fract, dGdX_fract, dBdX_fract;
    dRdX_fract = ((v3.r - v1.r) * dy21 + (v1.r - v2.r) * dy31);
    dGdX_fract = ((v3.g - v1.g) * dy21 + (v1.g - v2.g) * dy31);
    dBdX_fract = ((v3.b - v1.b) * dy21 + (v1.b - v2.b) * dy31);
    int dX_denom   = ((v3.x - v1.x) * dy21 + (v1.x - v2.x) * dy31);

    auto drawline = [&](int sx, int ex, int ny, int _r, int _g, int _b, int _dy_denom) {

        auto clipcolor = [](int c) {
            if (c < 0) return 0;
            if (c > 255) return 255;
            return c;
        };

        if (ny < 0) return;
        if (ny >= screeny) return;
        if (sx < 0) sx = 0;
        if (ex >= screenx) ex = screenx - 1;
        if (dX_denom == 0) dX_denom = 1;

        int red   = (_r / _dy_denom) * dX_denom;
        int green = (_g / _dy_denom) * dX_denom;
        int blue  = (_b / _dy_denom) * dX_denom;
        int offset = ny * screenx + sx;
        for (int i = sx; i <= ex; i++) {
            pixels[offset] = Pixel(clipcolor(red / dX_denom), clipcolor(green / dX_denom), clipcolor(blue / dX_denom), 255);
            offset++;
            red += dRdX_fract;
            green += dGdX_fract;
            blue += dBdX_fract;
        }
    };

    bool changed1, changed2, v2OnRightSide;
    int eA, eB, eAp, eBp;
    int minx, maxx;
    int xA = v1.x;
    int xB = v1.x;
    int y = v1.y;

    int r = 0, g = 0, b = 0;
    int dR_fract = 0, dG_fract = 0, dB_fract = 0;
    int dY_denom = dy21;

    dR_fract = (v2.r - v1.r);                                               // dR_fract = v2.r - v1.r + 4 * dRdX_fract / dX_denom;
    dG_fract = (v2.g - v1.g);                                               // dG_fract = v2.g - v1.g + 4 * dGdX_fract / dX_denom;
    dB_fract = (v2.b - v1.b);                                               // dB_fract = v2.b - v1.b + 4 * dBdX_fract / dX_denom;

    if (((v2.x - v1.x) * dy31) > ((v3.x - v1.x) * dy21)) {                  // V2 on the right side
        dR_fract = (v3.r - v1.r);
        dG_fract = (v3.g - v1.g);
        dB_fract = (v3.b - v1.b);
        dY_denom = dy31;
        v2OnRightSide = true; 
    }
    else {
        v2OnRightSide = false;
    }
    if (dY_denom == 0) dY_denom++;

    r = v1.r * dY_denom;                                                    // (v1.r - 4 * dRdX_fract / dX_denom) * dY_denom;
    g = v1.g * dY_denom;                                                    // (v1.g - 4 * dGdX_fract / dX_denom) * dY_denom;
    b = v1.b * dY_denom;                                                    // (v1.b - 4 * dBdX_fract / dX_denom) * dY_denom;

    if (dy21 > dx21) { std::swap(dx21, dy21); changed1 = true; }
    else { changed1 = false; }
    if (dy31 > dx31) { std::swap(dy31, dx31); changed2 = true; }
    else { changed2 = false; }

    eB = (int)(dx31 >> 1);
    if (v1.y != v2.y) {                                                     // Flat top, just process the second half
        eA = (int)(dx21 >> 1);

        for (int i = 0; i < dx21;) {
            eAp = 0;
            eBp = 0;
            if (xA < xB) { minx = xA; maxx = xB; }
            else { minx = xB; maxx = xA; }

            eA += dy21;
            if (changed1) {
                if (i < dx21) {
                    i++;
                    while (eA >= dx21) {
                        eA -= dx21;
                        eAp = ix21;
                    }
                }
            }
            else {
                while ((i < dx21) && (eA < dx21)) {
                    xA += ix21;
                    eA += dy21;
                    i++;
                }
                eA -= dx21;
                eAp += ix21;
                if (minx > xA) minx = xA;
                if (maxx < xA) maxx = xA;
            }

            eB += dy31;                                                     // process second line until y value is about to change
            if (changed2) {
                while (eB >= dx31) {
                    eB -= dx31;
                    eBp = ix31;
                }
            }
            else {
                while (eB < dx31) {
                    xB += ix31;
                    eB += dy31;
                }
                eB -= dx31;
                eBp += ix31;
                if (minx > xB) minx = xB;
                if (maxx < xB) maxx = xB;
            }

            drawline(minx, maxx, y, r, g, b, dY_denom);
            xA += eAp;
            xB += eBp;
            y++;                                                            // Now increase y
            r += dR_fract;
            g += dG_fract;
            b += dB_fract;
            if (y == v2.y) break;
            if (y == screeny) return;
        }
    }

    if (!v2OnRightSide) {
        dY_denom = dy32;
        if (dY_denom == 0) dY_denom++;
        dR_fract = (v3.r - v2.r);
        dG_fract = (v3.g - v2.g);
        dB_fract = (v3.b - v2.b);
        r = v2.r * dY_denom;
        g = v2.g * dY_denom;
        b = v2.b * dY_denom;
    }

    if (dy32 > dx32) { std::swap(dx32, dy32); changed1 = true; }
    else { changed1 = false; }

    eA = (int)(dx32 >> 1);
    xA = v2.x;

    for (int i = 0; i <= dx32; i++) {
        eAp = 0;
        eBp = 0;
        if (xA < xB) { minx = xA; maxx = xB; }
        else { minx = xB; maxx = xA; }

        eA += dy32;
        if (changed1) {
            if (i < dx32) {
                if (eA >= dx32) {
                    eA -= dx32;
                    eAp = ix32;
                }
            }
        }
        else {
            while ((i < dx32) && (eA < dx32)) {
                eA += dy32;
                xA += ix32;
                i++;
            }
            eA -= dx32;
            eAp += ix32;
            if (minx > xA) minx = xA;
            if (maxx < xA) maxx = xA;
        }

        if (changed2) {
            if (xB != v3.x) {
                eB += dy31;
                while (eB >= dx31) {
                    eB -= dx31;
                    eBp = ix31;
                }
            }
        }
        else {
            while (xB != v3.x) {
                eB += dy31;
                if (eB < dx31) xB += ix31;
                else {
                    eB -= dx31;
                    break;
                }
            }
            eBp += ix31;
            if (minx > xB) minx = xB;
            if (maxx < xB) maxx = xB;
        }

        drawline(minx, maxx, y, r, g, b, dY_denom);
        xA += eAp;
        xB += eBp;
        y++;                                                                // Now increase y
        r += dR_fract;
        g += dG_fract;
        b += dB_fract;
        if (y > v3.y) break;
        if (y == screeny) return;
    }

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#endif // ENGINE_H
