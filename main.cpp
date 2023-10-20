#define _USE_MATH_DEFINES
#include <vector>
#include <math.h>
#include <string>

#include "platform.h"
#include "bmploader.h"
#include "texture.h"

#define WND_WIDTH   800
#define WND_HEIGHT  600

using namespace std;

static float z_buffer[WND_WIDTH][WND_HEIGHT];                  // z-buffer uses heap data

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
    TVertex V1;
    TVertex V2;
    TVertex V3;
    Texture* texture;
    float   light;
};

struct TMesh
{
    vector<TTriangle> triangles;
};

struct TMat4x4 {
    float m[4][4] = { 0 };
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class Engine3D : public CPlatform
{
public:
    Engine3D() {
        sAppName = "Texturing";
    }

private:
    Texture leftTexture, topTexture, rightTexture, bottomTexture, frontTexture, backTexture;

    TMesh meshCube;
    TMat4x4 matProj;

    TVec3d vCamera;

    float fTheta;
    float fFar = 1000.0f;

    void clrZBuffer(void)
    {
        int x, y;
        for (x = 0; x < WND_WIDTH; x++) {
            for (y = 0; y < WND_HEIGHT; y++) {
                z_buffer[x][y] = (float)INT_MAX;
            }
        }
    }

    void MultiplyMatrixVector(TVertex& i, TVertex& o, TMat4x4& m)
    {
        o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
        o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
        o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
        float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

        if (w != 0.0f)
        {
            o.x /= w; o.y /= w; o.z /= w;
        }
    }

    inline void putPixel(int x, int y, float z, TRGBColor col, float light)
    {
        if (x < 0) x = 0;
        if (x > (WND_WIDTH - 1)) x = WND_WIDTH - 1;
        if (y < 0) y = 0;
        if (y > (WND_HEIGHT - 1)) y = WND_HEIGHT - 1;

        Pixel p(col.R * light, col.G * light, col.B * light);
        CLayer* pDrawTarget = GetDrawTarget();

        //if (z_buffer[x][y] > z) {
        //    z_buffer[x][y] = z;
        pDrawTarget->SetPixel(x, y, p);
        //}
    }

#define SUB_PIX(a) (ceil(a)-a)
    void drawTriangle(TTriangle t)
    {
        if ((t.texture->width == 0) || (t.texture->height == 0)) return;

        if (t.V1.y > t.V2.y) {                                              // sort the vertices (V1,V2,V3) by their Y values
            swap_data(t.V1, t.V2);
        }
        if (t.V1.y > t.V3.y) {
            swap_data(t.V1, t.V3);
        }
        if (t.V2.y > t.V3.y) {
            swap_data(t.V2, t.V3);
        }

        if ((int)t.V1.y == (int)t.V3.y) return;                             // check if we have more than a zero height triangle

        // We have to decide whether V2 is on the left side or the right one. We could do that by findng the V4, and
        // check the disatnce from V4 to V2 (V4.y = V2.y). V4 is one the edge (V1V3)
        // Using formula : I(t) = A + t(B-A) we find y = V1.y + t(V3.y - V1.y) and y = V2.y
        // V2.y - V1.y = t(V3.y - V1.y) -->  t = (V2.y - V1.y)/(V3.y - V1.y)
        // V4.x = V1.x + t(V3.x - V1.x) and V4.y = V2.y
        // float distance = V4.x - V2.x
        // if (distance > 0) then the middle vertex is on the left side (V1V3 is the longest edge on the right side)

        float dY21 = (float)1.0 / ceil(t.V2.y - t.V1.y);
        float dY31 = (float)1.0 / ceil(t.V3.y - t.V1.y);
        float dY32 = (float)1.0 / ceil(t.V3.y - t.V2.y);

        float dXdY21 = (float)(t.V2.x - t.V1.x) * dY21;                      // dXdY means deltaX/deltaY
        float dXdY31 = (float)(t.V3.x - t.V1.x) * dY31;
        float dXdY32 = (float)(t.V3.x - t.V2.x) * dY32;
        float dXdY31tmp = dXdY31;

        float dX = (float)1.0 / ((t.V3.x - t.V1.x) * ceil(t.V2.y - t.V1.y) + (t.V1.x - t.V2.x) * ceil(t.V3.y - t.V1.y));

        // we calculate delta values ​​to find the z value
        float dZdY21 = (float)(t.V2.z - t.V1.z) * dY21;
        float dZdY31 = (float)(t.V3.z - t.V1.z) * dY31;
        float dZdY32 = (float)(t.V3.z - t.V2.z) * dY32;
        float dZdX = (float)((t.V3.z - t.V1.z) * ceil(t.V2.y - t.V1.y) + (t.V1.z - t.V2.z) * ceil(t.V3.y - t.V1.y)) * dX;

        // we calculate delta values ​​to find the u-value of the texture
        float dUdY21 = (float)(t.V2.u - t.V1.u) * dY21;
        float dUdY31 = (float)(t.V3.u - t.V1.u) * dY31;
        float dUdY32 = (float)(t.V3.u - t.V2.u) * dY32;
        float dUdX = (float)((t.V3.u - t.V1.u) * ceil(t.V2.y - t.V1.y) + (t.V1.u - t.V2.u) * ceil(t.V3.y - t.V1.y)) * dX;
        if ((t.texture->width > 0) && (t.texture->height > 0)) {
            dUdY21 *= (t.texture->width - 1);
            dUdY31 *= (t.texture->width - 1);
            dUdY32 *= (t.texture->width - 1);
            dUdX *= (t.texture->width - 1);
        }

        // we calculate delta values ​​to find the v-value of the texture
        float dVdY21 = (float)(t.V2.v - t.V1.v) * dY21;
        float dVdY31 = (float)(t.V3.v - t.V1.v) * dY31;
        float dVdY32 = (float)(t.V3.v - t.V2.v) * dY32;
        float dVdX = (float)((t.V3.v - t.V1.v) * ceil(t.V2.y - t.V1.y) + (t.V1.v - t.V2.v) * ceil(t.V3.y - t.V1.y)) * dX;
        if ((t.texture->width > 0) && (t.texture->height > 0)) {
            dVdY21 *= (t.texture->height - 1);
            dVdY31 *= (t.texture->height - 1);
            dVdY32 *= (t.texture->height - 1);
            dVdX *= (t.texture->height - 1);
        }

        if (dXdY21 > dXdY31) {
            swap_data(dXdY21, dXdY31);
            dZdY21 = dZdY31;
            dUdY21 = dUdY31;
            dVdY21 = dVdY31;
        }

        int prestep = (int)SUB_PIX(t.V1.y);
        int x;
        int x_end;
        float x_left = t.V1.x + prestep * dXdY21;
        float x_right = t.V1.x + prestep * dXdY31;
        int y = (int)ceil(t.V1.y);
        float z, u, v;
        float zp = (t.V1.z + prestep * dZdY21);
        float up = (t.V1.u + prestep * dUdY21);
        float vp = (t.V1.v + prestep * dVdY21);
        if ((t.texture->width > 0) && (t.texture->height > 0)) {
            up *= (t.texture->width - 1);
            vp *= (t.texture->height - 1);
        }

        while (y < t.V2.y) {
            x_end = (int)ceil(x_right);
            z = ceil(zp);
            u = ceil(up);
            v = ceil(vp);
            if (x_left < x_right) {
                for (x = (int)ceil(x_left); x < x_end; x++) {
                    z += dZdX;
                    u += dUdX;
                    v += dVdX;
                    putPixel(x, y, z, t.texture->getColor((int)u, (int)v), t.light);
                }
            }
            else {
                x_end = (int)x_left;
                for (x = (int)ceil(x_right); x < x_end; x++) {
                    z += dZdX;
                    u += dUdX;
                    v += dVdX;
                    putPixel(x, y, z, t.texture->getColor((int)u, (int)v), t.light);
                }

            }
            x_left += dXdY21;
            x_right += dXdY31;
            zp += dZdY21;
            up += dUdY21;
            vp += dVdY21;
            y++;
        }

        dXdY31 = dXdY31tmp;
        if (dXdY32 < dXdY31) {
            swap_data(dXdY31, dXdY32);
            dZdY32 = dZdY31;
            dUdY32 = dUdY31;
            dVdY32 = dVdY31;
        }

        prestep = (int)SUB_PIX(t.V2.y);
        if (t.V2.x > t.V1.x) {
            x_right = t.V2.x + prestep * dXdY31;
        }
        else {
            x_left = t.V2.x + SUB_PIX(t.V2.y) * dXdY32;
            zp = (t.V2.z + prestep * dZdY32);
            up = (t.V2.u + prestep * dUdY32);
            vp = (t.V2.v + prestep * dVdY32);
            if ((t.texture->width > 0) && (t.texture->height > 0)) {
                up *= (t.texture->width - 1);
                vp *= (t.texture->height - 1);
            }
        }

        while (y < t.V3.y) {
            x_end = (int)ceil(x_right);
            z = ceil(zp);
            u = ceil(up);
            v = ceil(vp);

            for (x = (int)ceil(x_left); x < x_end; x++) {
                z += dZdX;
                u += dUdX;
                v += dVdX;
                putPixel(x, y, z, t.texture->getColor((int)u, (int)v), t.light);
            }
            x_left += dXdY32;
            x_right += dXdY31;
            zp += dZdY32;
            up += dUdY32;
            vp += dVdY32;
            y++;

        }

        return;
    }

public:
    bool OnUserCreate() override
    {
        leftTexture.loadFromBitmap("negx.bmp");
        topTexture.loadFromBitmap("posy.bmp");
        rightTexture.loadFromBitmap("posx.bmp");
        bottomTexture.loadFromBitmap("negy.bmp");
        frontTexture.loadFromBitmap("negz.bmp");
        backTexture.loadFromBitmap("posz.bmp");

        meshCube.triangles = {

            // FRONT
            { { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 0.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 0.0f },  &frontTexture },
            { { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 0.0f, 1.0f, 1.0f },  &frontTexture },

            // RIGHT
            { { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 0.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },  &rightTexture },
            { { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 1.0f, 1.0f, 1.0f },  &rightTexture },

            // BACK
            { { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 1.0f, 1.0f, 1.0f, 0.0f },  &backTexture },
            { { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 0.0f, 0.0f, 1.0f, 1.0f, 1.0f },  &backTexture },

            // LEFT
            { { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 1.0f, 0.0f, 1.0f, 0.0f },  &leftTexture },
            { { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 0.0f, 1.0f, 0.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f },  &leftTexture },

            // TOP
            { { 0.0f, 1.0f, 0.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },  &topTexture },
            { { 0.0f, 1.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 1.0f },  &topTexture },

            // BOTTOM
            { { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 0.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 0.0f },  &bottomTexture },
            { { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 0.0f, 1.0f, 1.0f },  &bottomTexture },
        };

        // Projection matrix
        float fScale = 1.0f;
        float fNear = 0.1f;
        float fFov = 90.0f;
        float fAspectRatio = (float)WND_HEIGHT / (float)WND_WIDTH;
        float fFovRad = (int)1.0f / tanf(fFov * 0.5f / 180.0f * (int)M_PI);

        matProj.m[0][0] = fScale * fAspectRatio * fFovRad;
        matProj.m[1][1] = fScale * fFovRad;
        matProj.m[2][2] = fScale * fFar / (fFar - fNear);
        matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
        matProj.m[2][3] = 1.0f;
        matProj.m[3][3] = 0.0f;


        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override
    {

        Clear(BLACK);
        clrZBuffer();                                                       // Clear Screen and Z buffer

        TMat4x4 matRotX, matRotZ;                                           // Set up rotation matrices
        fTheta += 1.0f * fElapsedTime;

        matRotZ.m[0][0] = cosf(fTheta);                                     // Rotation Z
        matRotZ.m[0][1] = sinf(fTheta);
        matRotZ.m[1][0] = -sinf(fTheta);
        matRotZ.m[1][1] = cosf(fTheta);
        matRotZ.m[2][2] = 1;
        matRotZ.m[3][3] = 1;

        matRotX.m[0][0] = 1;                                                // Rotation X
        matRotX.m[1][1] = cosf(fTheta * 0.5f);
        matRotX.m[1][2] = sinf(fTheta * 0.5f);
        matRotX.m[2][1] = -sinf(fTheta * 0.5f);
        matRotX.m[2][2] = cosf(fTheta * 0.5f);
        matRotX.m[3][3] = 1;

        for (auto tri : meshCube.triangles)                                 // Draw Triangles
        {
            TTriangle triProjected, triTranslated, triRotatedZ, triRotatedZX;

            triProjected = tri;                                            // get u and v data into projected triangles

            // Rotate in Z-Axis
            MultiplyMatrixVector(tri.V1, triRotatedZ.V1, matRotZ);
            MultiplyMatrixVector(tri.V2, triRotatedZ.V2, matRotZ);
            MultiplyMatrixVector(tri.V3, triRotatedZ.V3, matRotZ);

            // Rotate in X-Axis
            MultiplyMatrixVector(triRotatedZ.V1, triRotatedZX.V1, matRotX);
            MultiplyMatrixVector(triRotatedZ.V2, triRotatedZX.V2, matRotX);
            MultiplyMatrixVector(triRotatedZ.V3, triRotatedZX.V3, matRotX);

            // Offset into the screen
            triTranslated = triRotatedZX;
            triTranslated.V1.z = triRotatedZX.V1.z + 3.0f;
            triTranslated.V2.z = triRotatedZX.V2.z + 3.0f;
            triTranslated.V3.z = triRotatedZX.V3.z + 3.0f;

            TVec3d normal, line1, line2;

            line1.x = triTranslated.V2.x - triTranslated.V1.x;
            line1.y = triTranslated.V2.y - triTranslated.V1.y;
            line1.z = triTranslated.V2.z - triTranslated.V1.z;

            line2.x = triTranslated.V3.x - triTranslated.V1.x;
            line2.y = triTranslated.V3.y - triTranslated.V1.y;
            line2.z = triTranslated.V3.z - triTranslated.V1.z;

            normal.x = line1.y * line2.z - line1.z * line2.y;
            normal.y = line1.z * line2.x - line1.x * line2.z;
            normal.z = line1.x * line2.y - line1.y * line2.x;

            float l = sqrtf(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
            normal.x /= l;
            normal.y /= l;
            normal.z /= l;

            //if (normal.z < 0) {
            // 
            if (normal.x * (triTranslated.V1.x - vCamera.x) +
                normal.y * (triTranslated.V1.y - vCamera.y) +
                normal.z * (triTranslated.V1.z - vCamera.z) < 0.0)
            {
                // Illumination
                TVec3d light_direction = { 0.0f, 0.0f, -1.0f };

                float ll = sqrtf(light_direction.x * light_direction.x + light_direction.y * light_direction.y + light_direction.z * light_direction.z);
                light_direction.x /= ll;
                light_direction.y /= ll;
                light_direction.z /= ll;

                float dp = normal.x * light_direction.x +
                    normal.y * light_direction.y +
                    normal.z * light_direction.z;

                triProjected.light = dp;

                // Project triangles from 3D --> 2D
                MultiplyMatrixVector(triTranslated.V1, triProjected.V1, matProj);
                MultiplyMatrixVector(triTranslated.V2, triProjected.V2, matProj);
                MultiplyMatrixVector(triTranslated.V3, triProjected.V3, matProj);

                // Scale into view
                triProjected.V1.x += 1.0f; triProjected.V1.y += 1.0f;
                triProjected.V2.x += 1.0f; triProjected.V2.y += 1.0f;
                triProjected.V3.x += 1.0f; triProjected.V3.y += 1.0f;
                triProjected.V1.x *= 0.5f * (float)ScreenWidth();
                triProjected.V1.y *= 0.5f * (float)ScreenHeight();
                triProjected.V1.z *= 0.5f * fFar;
                triProjected.V2.x *= 0.5f * (float)ScreenWidth();
                triProjected.V2.y *= 0.5f * (float)ScreenHeight();
                triProjected.V2.z *= 0.5f * fFar;
                triProjected.V3.x *= 0.5f * (float)ScreenWidth();
                triProjected.V3.y *= 0.5f * (float)ScreenHeight();
                triProjected.V3.z *= 0.5f * fFar;

                drawTriangle(triProjected);                                     // Rasterize triangle
            }

        }

        return true;
    }
};

int WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR CmdLine, int nShowCmd)
{
    Engine3D demo;

    if (demo.Construct(WND_WIDTH, WND_HEIGHT, false, false, false)) {
        //if (demo.Construct(WND_WIDTH, WND_HEIGHT, 1, 1, false, false)) {
        demo.Start();
    }

    return 0;
}

