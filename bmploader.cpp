#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <cstdint>
#include <stdio.h>
#include <string.h>
#include <string>
#include "bmploader.h"

size_t strlcpy(char *dst, const char *src, size_t n) {
    char *p = dst;

    if (n != 0) {
        for (; --n != 0; p++, src++) {
            if ((*p = *src) == '\0')
                return p - dst;
        }
        *p = '\0';
    }
    return (p - dst) + strlen(src);
}

size_t strlcat(char *dst, const char *src, size_t n) {
    char *p = dst;

    while (n != 0 && *p != '\0') {
        p++;
        n--;
    }
    if (n != 0) {
        for (; --n != 0; p++, src++) {
            if ((*p = *src) == '\0')
                return p - dst;
        }
        *p = '\0';
    }
    return (p - dst) + strlen(src);
}

BMPLoader::BMPLoader()
{
}

TRGBColor *BMPLoader::loadTexture(const std::string fileName, int &width, int &height)
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
        uint32_t Endpoints[9];                                                  // see http://msdn2.microsoft.com/en-us/library/ms536569.aspx
        uint32_t GammaRed;
        uint32_t GammaGreen;
        uint32_t GammaBlue;
    } BITMAP_HEADER;

    #pragma pack(pop)

    int  i;
    char palette[256][4];                                                       // needed for the color palette
    int  offset;
    int  lines;
    int  paddedWidth;
    int  color;

    FILE *BMPFile;                                                              // file handle
    BITMAP_FILEHEADER m_BitmapFileHeader;
    BITMAP_HEADER m_BitmapHeader;

    fopen_s(&BMPFile, fileName.c_str(), "rb");
    if (BMPFile == NULL) {
        std::string fileName2;
        fileName2 = "..\\Texturing\\" + std::string(fileName);
        fopen_s(&BMPFile, fileName2.c_str(), "rb");
        if (BMPFile == NULL) {
            OutputDebugStringA(std::string("Error opening file " + fileName).c_str());
            return NULL;
        }
        return NULL;
    }
    fread(&m_BitmapFileHeader, BITMAP_FILEHEADER_SIZE, 1, BMPFile);             // we read the file header

    // We check whether the appropriate bitmap file can be displayed
    if (m_BitmapFileHeader.Signature != 19778 || m_BitmapFileHeader.Reserved != 0) {    // Not a valid bitmap file - don't display
        fclose(BMPFile);
        OutputDebugStringA("Not a valid bitmap file - don't display");
        return NULL;
    }

    fread(&m_BitmapHeader, sizeof(BITMAP_HEADER), 1, BMPFile);

    if (m_BitmapHeader.Compression != 0) {                                      // Compressed file - don't display
        fclose(BMPFile);
        OutputDebugStringA("Compressed file - don't display");
        return NULL;
    }

    fseek(BMPFile, BITMAP_FILEHEADER_SIZE + m_BitmapHeader.HeaderSize, SEEK_SET);

    offset = m_BitmapHeader.Width * (m_BitmapHeader.Height-1);                  // Set appropriate position to display graphics data
    paddedWidth = m_BitmapHeader.Width & 0xFFFC;                                // Pad line length to 4bytes
    if (m_BitmapHeader.Width != paddedWidth)
        paddedWidth += 4;

    TRGBColor * texture = new TRGBColor[m_BitmapHeader.Width * m_BitmapHeader.Height];
    width = m_BitmapHeader.Width;
    height = m_BitmapHeader.Height;
    lines = 0;                                                                  // We've read no lines so far

    if (m_BitmapHeader.BitCount == 8) {
        fread (&palette, 1024, 1, BMPFile);                                     // we read the color palette

        while (lines < m_BitmapHeader.Height)	{                               // Decode and display graphics
            for (i=0; i < paddedWidth; i++) {                                   // Read next line
                color =  getc(BMPFile);
                texture[offset + i].R = palette[color & 0xFF][2];
                texture[offset + i].G = palette[color & 0xFF][1];
                texture[offset + i].B = palette[color & 0xFF][0];
                texture[offset + i].alpha = palette[color & 0xFF][3];
            }
            offset -= m_BitmapHeader.Width;                                     // Move up one line on the screen
            lines++;                                                            // increase amount of lines read
        }

        fclose(BMPFile);
        return (TRGBColor*) texture;
    }
    else {                                                                      // Other than 8-bit colour - don't display
        fseek(BMPFile, m_BitmapFileHeader.BitsOffset, SEEK_SET);

        if (m_BitmapHeader.BitCount == 16) {
            unsigned short col16bit;
            while (lines < m_BitmapHeader.Height)	{                           // Decode and display graphics
                for (i=0; i < paddedWidth; i++) {                               // Read next line
                    col16bit  = getc(BMPFile);
                    col16bit |= (unsigned short)getc(BMPFile) << 8;
                    texture[offset + i].B = (col16bit & 0x1f) << 3;
                    texture[offset + i].G = ((col16bit >> 5) & 0x1f) << 3;
                    texture[offset + i].R = ((col16bit >> 10) & 0x1f) << 3;
                    texture[offset + i].alpha = 255;
                }
                offset -= m_BitmapHeader.Width;                                 // Move up one line on the screen
                lines++;                                                        // increase amount of lines read
            }
        }
        if (m_BitmapHeader.BitCount == 24) {
            while (lines < m_BitmapHeader.Height)	{                           // Decode and display graphics
                for (i=0; i < paddedWidth; i++) {                               // Read next line
                    texture[offset + i].B = getc(BMPFile);
                    texture[offset + i].G = getc(BMPFile);
                    texture[offset + i].R = getc(BMPFile);
                    texture[offset + i].alpha = 255;
                }
                offset -= m_BitmapHeader.Width;                                 // Move up one line on the screen
                lines++;                                                        // increase amount of lines read
            }
        }
        else if (m_BitmapHeader.BitCount == 32) {
            while (lines < m_BitmapHeader.Height)	{                           // Decode and display graphics
                for (i=0; i < paddedWidth; i++) {                               // Read next line
                    texture[offset + i].B = getc(BMPFile);
                    texture[offset + i].G = getc(BMPFile);
                    texture[offset + i].R = getc(BMPFile);
                    texture[offset + i].alpha = getc(BMPFile);
                }
                offset -= m_BitmapHeader.Width;                                 // Move up one line on the screen
                lines++;                                                        // increase amount of lines read
            }
        }

        fclose(BMPFile);
        return (TRGBColor*) texture;
    }

    fclose(BMPFile);
    return NULL;
};
