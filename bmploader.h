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
};

class BMPLoader
{
public:
    BMPLoader();
    static TRGBColor *loadTexture(const std::string fileName, int &width, int &height);
};

#endif // BMPLOADER_H
