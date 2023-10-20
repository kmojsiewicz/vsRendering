#ifndef TEXTURE_H
#define TEXTURE_H

#include <string>
#include "bmploader.h"

class Texture
{
public:
    TRGBColor *data;
    int width;
    int height;

    Texture() { width = 0; height = 0; data = NULL; }
    ~Texture() { if (data) delete data; }

    void draw();
    int loadFromBitmap(const std::string fileName);
    TRGBColor getColor(int x, int y);
};

#endif // TEXTURE_H
