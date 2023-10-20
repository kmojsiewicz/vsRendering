#ifndef BMPLOADER_H
#define BMPLOADER_H

#include <stdlib.h>

union TRGBColor {
    struct {
        unsigned char B;
        unsigned char G;
        unsigned char R;
        unsigned char alpha;
    };
    unsigned int RGBa;

    TRGBColor() : R(0), G(0), B(0), alpha(255) {}
    TRGBColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a) : R(r), G(g), B(b), alpha(a) {}
};

class BMPLoader
{
public:
    BMPLoader();
    static TRGBColor* loadTexture(const std::string fileName, int& width, int& height);
};

#endif // BMPLOADER_H
