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

public:
    bool OnUserCreate() override
    {
        //mesh.LoadFromObjectFile("test.obj", nullptr, WHITE, 1.0);
        mesh.MakeQube(nullptr, WHITE, 1.0);
        frontTexture.LoadFromBitmap("negz.bmp");   mesh.triangles[0].SetTexture(&frontTexture);   mesh.triangles[1].SetTexture(&frontTexture);      // front
        rightTexture.LoadFromBitmap("posx.bmp");   mesh.triangles[2].SetTexture(&rightTexture);   mesh.triangles[3].SetTexture(&rightTexture);      // right
        backTexture.LoadFromBitmap("posz.bmp");    mesh.triangles[4].SetTexture(&backTexture);    mesh.triangles[5].SetTexture(&backTexture);       // back
        leftTexture.LoadFromBitmap("negx.bmp");    mesh.triangles[6].SetTexture(&leftTexture);    mesh.triangles[7].SetTexture(&leftTexture);       // left
        topTexture.LoadFromBitmap("posy.bmp");     mesh.triangles[8].SetTexture(&topTexture);     mesh.triangles[9].SetTexture(&topTexture);        // top
        bottomTexture.LoadFromBitmap("negy.bmp");  mesh.triangles[10].SetTexture(&bottomTexture); mesh.triangles[11].SetTexture(&bottomTexture);    // bottom

        float fScale = 1.0f;
        float fNear = 0.1f;
        float fFov = 90.0f;
        float fAspectRatio = (float)WND_HEIGHT / (float)WND_WIDTH;
        float fFovRad = (int)1.0f / tanf(fFov * 0.5f / 180.0f * (int)M_PI);

        matProj.m[0][0] = fScale * fAspectRatio * fFovRad;                  // Projection matrix
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

        Clear(DARK_GREY);
        ClrZBuffer();                                                       // Clear Screen and Z buffer
        fTheta += 1.0f * fElapsedTime;

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
            triTranslated.V1.x += 0.0f;                                     // Offset into the screen
            triTranslated.V2.x += 0.0f;
            triTranslated.V3.x += 0.0f;
            triTranslated.V1.z += 2.0f;                                     // Offset into the screen
            triTranslated.V2.z += 2.0f;
            triTranslated.V3.z += 2.0f;

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

        get_timepoint();                                                    // sort triangles from back to front
            sort(vTrianglesToRaster.begin(), vTrianglesToRaster.end(), [](TTriangle& t1, TTriangle t2) {
                float z1 = (t1.V1.z + t1.V2.z + t1.V3.z) / 3.0f;
                float z2 = (t2.V1.z + t2.V2.z + t2.V3.z) / 3.0f;
                return z1 > z2;
                });
        print_diff_timepoint("Sorting time      : ");

        get_timepoint();
        for (const auto &triProjected : vTrianglesToRaster) {               // Draw Triangles sorted
            FillTriangle(triProjected);                                     // Rasterize triangle
        }
        print_diff_timepoint("Rendering time  : ");

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

