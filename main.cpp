#define _USE_MATH_DEFINES
#include <vector>
#include <math.h>
#include <string>
#include <fstream>
#include <strstream>
#include <algorithm>

#include "platform.h"
#include "engine.h"

using namespace std;



struct TMat4x4 {
    float m[4][4] = { 0 };
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class Engine3D : public CEngine
{
public:
    Engine3D() {
        sAppName = "Rendering";
    }

private:
    bool printDurTimeOnce = true;
    std::chrono::high_resolution_clock::time_point tp_start, tp_end;
    std::chrono::duration<float> elapsed_sec;
    #define get_timepoint()           { if (printDurTimeOnce) { tp_start = std::chrono::high_resolution_clock::now(); } }
    #define print_diff_timepoint(...) { if (printDurTimeOnce) { tp_end = std::chrono::high_resolution_clock::now(); \
                                      elapsed_sec = tp_end - tp_start; \
                                      cout << __VA_ARGS__ << std::chrono::duration_cast<std::chrono::microseconds>(elapsed_sec).count() << " us" << endl; } }


    CTexture leftTexture, topTexture, rightTexture, bottomTexture, frontTexture, backTexture;
    TMesh mesh;
    TMat4x4 matProj;
    TVec3d vCamera;
    float fTheta;
    float fFar = 1000.0f;

    void MultiplyMatrixVector(TVertex& i, TVertex& o, TMat4x4& m)
    {
        o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
        o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
        o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
        float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

        if (w != 0.0f) {
            o.x /= w; o.y /= w; o.z /= w;
        }
    }

    void MultiplyMatrixByMatrix(TMat4x4& o, TMat4x4& m1, TMat4x4& m2)
    {
        int i, j, k;
        for (i = 0; i < 4; ++i) {
            for (j = 0; j < 4; ++j) {
                for (k = 0; k < 4; ++k) {
                    o.m[i][j] += m1.m[i][k] * m2.m[k][j];
                }
            }
        }
    }

#define SUB_PIX(a) (ceil(a)-a)
    void fillTriangleT(TTriangle t)
    {
        int u_width = t.texture->GetWidth();
        int v_height = t.texture->GetHeight();

        if ((u_width == 0) || (v_height == 0)) return;

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
        if ((u_width > 0) && (v_height > 0)) {
            dUdY21 *= (u_width - 1);
            dUdY31 *= (u_width - 1);
            dUdY32 *= (u_width - 1);
            dUdX *= (u_width - 1);
        }

        // we calculate delta values ​​to find the v-value of the texture
        float dVdY21 = (float)(t.V2.v - t.V1.v) * dY21;
        float dVdY31 = (float)(t.V3.v - t.V1.v) * dY31;
        float dVdY32 = (float)(t.V3.v - t.V2.v) * dY32;
        float dVdX = (float)((t.V3.v - t.V1.v) * ceil(t.V2.y - t.V1.y) + (t.V1.v - t.V2.v) * ceil(t.V3.y - t.V1.y)) * dX;
        if ((u_width > 0) && (v_height > 0)) {
            dVdY21 *= (v_height - 1);
            dVdY31 *= (v_height - 1);
            dVdY32 *= (v_height - 1);
            dVdX *= (v_height - 1);
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
        if ((u_width > 0) && (v_height > 0)) {
            up *= (u_width - 1);
            vp *= (v_height - 1);
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
                    PutPixel(x, y, z, t.texture->GetPixel((int)u, (int)v), t.light);
                }
            }
            else {
                x_end = (int)x_left;
                for (x = (int)ceil(x_right); x < x_end; x++) {
                    z += dZdX;
                    u += dUdX;
                    v += dVdX;
                    PutPixel(x, y, z, t.texture->GetPixel((int)u, (int)v), t.light);
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
            if ((u_width > 0) && (v_height > 0)) {
                up *= (u_width - 1);
                vp *= (v_height - 1);
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
                PutPixel(x, y, z, t.texture->GetPixel((int)u, (int)v), t.light);
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
        frontTexture.LoadFromBitmap("negz.bmp");
        mesh.LoadFromObjectFile("test.obj", &frontTexture, 1.0);
        //mesh.MakeQube(nullptr, 1.0);
        //frontTexture.LoadFromBitmap("negz.bmp");   mesh.triangles[0].SetTexture(&frontTexture);   mesh.triangles[1].SetTexture(&frontTexture);      // front
        //rightTexture.LoadFromBitmap("posx.bmp");   mesh.triangles[2].SetTexture(&rightTexture);   mesh.triangles[3].SetTexture(&rightTexture);      // right
        //backTexture.LoadFromBitmap("posz.bmp");    mesh.triangles[4].SetTexture(&backTexture);    mesh.triangles[5].SetTexture(&backTexture);       // back
        //leftTexture.LoadFromBitmap("negx.bmp");    mesh.triangles[6].SetTexture(&leftTexture);    mesh.triangles[7].SetTexture(&leftTexture);       // left
        //topTexture.LoadFromBitmap("posy.bmp");     mesh.triangles[8].SetTexture(&topTexture);     mesh.triangles[9].SetTexture(&topTexture);        // top
        //bottomTexture.LoadFromBitmap("negy.bmp");  mesh.triangles[10].SetTexture(&bottomTexture); mesh.triangles[11].SetTexture(&bottomTexture);    // bottom

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
        vector<TTriangle> vTrianglesToRaster;
        TMat4x4 matRotX, matRotZ, matRotZX;                                 // Set up rotation matrices

        Clear(BLACK);
        ClrZBuffer();                                                       // Clear Screen and Z buffer
        fTheta += 1.0f * fElapsedTime;
        //fTheta += 0.1f * fElapsedTime;

        get_timepoint();
            matRotZ.m[0][0] = cosf(fTheta);                                 // Rotation Z
            matRotZ.m[0][1] = sinf(fTheta);
            matRotZ.m[1][0] = -sinf(fTheta);
            matRotZ.m[1][1] = cosf(fTheta);
            matRotZ.m[2][2] = 1;
            matRotZ.m[3][3] = 1;
            matRotX.m[0][0] = 1;                                            // Rotation X
            matRotX.m[1][1] = cosf(fTheta * 0.5f);
            matRotX.m[1][2] = sinf(fTheta * 0.5f);
            matRotX.m[2][1] = -sinf(fTheta * 0.5f);
            matRotX.m[2][2] = cosf(fTheta * 0.5f);
            matRotX.m[3][3] = 1;
            MultiplyMatrixByMatrix(matRotZX, matRotZ, matRotX);             // Rotation in Z-Axis and X-Axis
        print_diff_timepoint("Rotation time     : ");

        get_timepoint();
        for (auto tri : mesh.triangles) 
        {
            TVec3d normal;
            TTriangle triProjected, triTranslated, triRotatedZX;

            triProjected = tri;                                             // get u and v data into projected triangles

            MultiplyMatrixVector(tri.V1, triRotatedZX.V1, matRotZX);
            MultiplyMatrixVector(tri.V2, triRotatedZX.V2, matRotZX);
            MultiplyMatrixVector(tri.V3, triRotatedZX.V3, matRotZX);

            triTranslated = triRotatedZX;
            triTranslated.V1.x += 3.0f;                                     // Offset into the screen
            triTranslated.V2.x += 3.0f;
            triTranslated.V3.x += 3.0f;
            triTranslated.V1.z += 8.0f;                                     // Offset into the screen
            triTranslated.V2.z += 8.0f;
            triTranslated.V3.z += 8.0f;

            normal = triTranslated.NormalVector();
            if (normal.x * (triTranslated.V1.x - vCamera.x) +
                normal.y * (triTranslated.V1.y - vCamera.y) +
                normal.z * (triTranslated.V1.z - vCamera.z) < 0.0)
            {
                TVec3d light_direction = { 0.0f, 0.0f, -1.0f };             // Illumination

                float inv_sqrt_ll = q_rsqrt(light_direction.x * light_direction.x + light_direction.y * light_direction.y + light_direction.z * light_direction.z);
                light_direction.x *= inv_sqrt_ll;
                light_direction.y *= inv_sqrt_ll;
                light_direction.z *= inv_sqrt_ll;
                float dp = normal.x * light_direction.x + normal.y * light_direction.y + normal.z * light_direction.z;

                // Project triangles from 3D --> 2D
                MultiplyMatrixVector(triTranslated.V1, triProjected.V1, matProj);
                MultiplyMatrixVector(triTranslated.V2, triProjected.V2, matProj);
                MultiplyMatrixVector(triTranslated.V3, triProjected.V3, matProj);

                triProjected.light = dp;
                triProjected.Translate(1.0f, 1.0f, 0.0f);
                triProjected.Scale(0.5f * (float)ScreenWidth(), 0.5f * (float)ScreenHeight(), 0.5f * fFar);
                vTrianglesToRaster.push_back(triProjected);
            }
        }
        print_diff_timepoint("Matrix projection : ");

        // sort triangles from back to front
        get_timepoint();
            sort(vTrianglesToRaster.begin(), vTrianglesToRaster.end(), [](TTriangle& t1, TTriangle t2) {
                float z1 = (t1.V1.z + t1.V2.z + t1.V3.z) / 3.0f;
                float z2 = (t2.V1.z + t2.V2.z + t2.V3.z) / 3.0f;
                return z1 > z2;
                });
        print_diff_timepoint("Sorting time      : ");

        //get_timepoint();
        //for (auto triProjected : vTrianglesToRaster)                        // Draw Triangles sorted
        //{
        //    fillTriangle(triProjected, WHITE * triProjected.light);         // Rasterize triangle
        //    //fillTriangleT(triProjected);
        //}
        //print_diff_timepoint("Rendering time 1  : ");

        get_timepoint();
        for (const auto &triProjected : vTrianglesToRaster)                 // Draw Triangles sorted
        {
            //fillTriangle(triProjected, WHITE * triProjected.light);       // Rasterize triangle
            //fillTriangleT(triProjected);
            DrawTriangle(triProjected);
        }
        print_diff_timepoint("Rendering time 2  : ");

        printDurTimeOnce = false;
        return true;
    }
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR CmdLine, _In_opt_ int nShowCmd)
{
    Engine3D demo;

    if (demo.Construct(WND_WIDTH, WND_HEIGHT, false, false, true)) {
        //if (demo.Construct(WND_WIDTH, WND_HEIGHT, 1, 1, false, false)) {
        demo.Start();
    }

    return 0;
}

