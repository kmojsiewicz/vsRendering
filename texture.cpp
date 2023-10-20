#include "texture.h"

void Texture::draw()
{
    int x;
    int y;
    int offset;
    TRGBColor Color;

    if (data == NULL) {
        return;
    }

    offset = 0;
    for (y=0; y<height; y++) {
        for (x=0; x<width; x++) {
            //painter->setPen(QColor::fromRgb((QRgb)data[x + offset]));
            //painter->drawPoint(x, y);
        }
        offset += width;
    }
}

int Texture::loadFromBitmap(const std::string fileName)
{
    if ((data = BMPLoader::loadTexture(fileName, width, height)) != NULL) {
        return 1;
    }

    return 0;
}

TRGBColor Texture::getColor(int x, int y)
{
    if (x < 0) x = 0;
    if (x > (width - 1)) x = width -1;
    if (y < 0) y = 0;
    if (y > (height - 1)) y = height -1;
    return data[x + y*width];
}
