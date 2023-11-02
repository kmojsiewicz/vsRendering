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
    TVec3d vCamera;                                                         // Location of camera in world space
    TVec3d vLookDir;	                                                    // Direction vector along the direction camera points
    float fYaw;		                                                        // FPS Camera rotation in XZ plane
    float fTheta;                                                           // Spins World transform

public:
    bool OnUserCreate() override
    {
        matProj = Matrix_MakeProjection(90.0f, (float)WND_HEIGHT / (float)WND_WIDTH, 0.1f, 1000.0f, 1.0f);      // field of view = 90 degree, fNear = 0.1f, fFar = 1000.0f, fScale = 1.0f

        mesh.LoadFromObjectFile("axis.obj", nullptr, WHITE, 1.0);
        //mesh.MakeQube(nullptr, WHITE, 1.0);
        //frontTexture.LoadFromBitmap("negz.bmp");   mesh.triangles[0].SetTexture(&frontTexture);   mesh.triangles[1].SetTexture(&frontTexture);      // front
        //rightTexture.LoadFromBitmap("posx.bmp");   mesh.triangles[2].SetTexture(&rightTexture);   mesh.triangles[3].SetTexture(&rightTexture);      // right
        //backTexture.LoadFromBitmap("posz.bmp");    mesh.triangles[4].SetTexture(&backTexture);    mesh.triangles[5].SetTexture(&backTexture);       // back
        //leftTexture.LoadFromBitmap("negx.bmp");    mesh.triangles[6].SetTexture(&leftTexture);    mesh.triangles[7].SetTexture(&leftTexture);       // left
        //topTexture.LoadFromBitmap("posy.bmp");     mesh.triangles[8].SetTexture(&topTexture);     mesh.triangles[9].SetTexture(&topTexture);        // top
        //bottomTexture.LoadFromBitmap("negy.bmp");  mesh.triangles[10].SetTexture(&bottomTexture); mesh.triangles[11].SetTexture(&bottomTexture);    // bottom
        
        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override
    {
        if (GetKey(VK_UP).bHeld)    vCamera.y += 8.0f * fElapsedTime;
        if (GetKey(VK_DOWN).bHeld)  vCamera.y -= 8.0f * fElapsedTime;
        if (GetKey(VK_LEFT).bHeld)  vCamera.x -= 8.0f * fElapsedTime;
        if (GetKey(VK_RIGHT).bHeld) vCamera.x += 8.0f * fElapsedTime;
        
        TVec3d vForward = Vector_Mul(vLookDir, 8.0f * fElapsedTime);
        if (GetKey(L'W').bHeld) vCamera = Vector_Add(vCamera, vForward);
        if (GetKey(L'S').bHeld) vCamera = Vector_Sub(vCamera, vForward);
        if (GetKey(L'A').bHeld) fYaw -= 2.0f * fElapsedTime;
        if (GetKey(L'D').bHeld) fYaw += 2.0f * fElapsedTime;
        
        get_timepoint();
            //fTheta += 1.0f * fElapsedTime;
            TMat4x4 matRotZ = Matrix_MakeRotationZ(fTheta);                 // Rotation Z
            TMat4x4 matRotX = Matrix_MakeRotationX(fTheta * 0.5f);          // Rotation X
            TMat4x4 matTrans = Matrix_MakeTranslation(0.0f, 0.0f, 18.0f);
            TMat4x4 matWorld = Matrix_MakeIdentity();                       // Form World Matrix
            matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);             // Transform by rotation (Rotation in Z-Axis and X-Axis)
            matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);           // Transform by translation

            TVec3d vUp = { 0, 1, 0 };
            TVec3d vTarget = { 0, 0, 1 };
            TMat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
            vLookDir = Matrix_MultiplyVector(matCameraRot, vTarget);
            vTarget = Vector_Add(vCamera, vLookDir);
            TMat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);      // Create "Point At" Matrix for camera
            TMat4x4 matView = Matrix_QuickInverse(matCamera);               // Make view matrix from camera
        print_diff_timepoint("Rotation and translation time     : ");


        get_timepoint();
        vector<TTriangle> vTrianglesToRaster;
        for (auto tri : mesh.triangles) 
        {
            TTriangle triTransformed, triVieved, triProjected = tri;        // get u and v data into projected triangles
            
            // World Matrix Transform
            triTransformed.V1.p = Matrix_MultiplyVector(matWorld, tri.V1.p);    
            triTransformed.V2.p = Matrix_MultiplyVector(matWorld, tri.V2.p);
            triTransformed.V3.p = Matrix_MultiplyVector(matWorld, tri.V3.p);

            TVec3d normal = triTransformed.NormalVector();
            TVec3d vCameraRay = Vector_Sub(triTransformed.V1.p, vCamera);
            if (Vector_DotProduct(normal, vCameraRay) < 0.0)
            {
                TVec3d light_direction = { 0.0f, 1.0f, -1.0f };             // Illumination
                light_direction = Vector_Normalise(light_direction);
                triProjected.light = max(0.1f, Vector_DotProduct(light_direction, normal));

                // Convert World Space --> View Space
                triVieved.V1.p = Matrix_MultiplyVector(matView, triTransformed.V1.p);
                triVieved.V2.p = Matrix_MultiplyVector(matView, triTransformed.V2.p);
                triVieved.V3.p = Matrix_MultiplyVector(matView, triTransformed.V3.p);

                // Project triangles from 3D --> 2D with normalising into cartesian space
                triProjected.V1.p = Matrix_MultiplyVector(matProj, triVieved.V1.p);
                triProjected.V2.p = Matrix_MultiplyVector(matProj, triVieved.V2.p);
                triProjected.V3.p = Matrix_MultiplyVector(matProj, triVieved.V3.p);
                triProjected.V1.p = Vector_Div(triProjected.V1.p, triProjected.V1.p.w);
                triProjected.V2.p = Vector_Div(triProjected.V2.p, triProjected.V2.p.w);
                triProjected.V3.p = Vector_Div(triProjected.V3.p, triProjected.V3.p.w);

                // X/Y are inverted so put them back
                triProjected.V1.p.x *= -1.0f;
                triProjected.V2.p.x *= -1.0f;
                triProjected.V3.p.x *= -1.0f;
                triProjected.V1.p.y *= -1.0f;
                triProjected.V2.p.y *= -1.0f;
                triProjected.V3.p.y *= -1.0f;

                triProjected.Translate({ 1.0f, 1.0f, 0.0f });               // Offset verts into visible normalised space
                triProjected.Scale(0.5f * (float)ScreenWidth(), 0.5f * (float)ScreenHeight(), 0.5f * 100.0);
                
                vTrianglesToRaster.push_back(triProjected);
            }
        }
        print_diff_timepoint("Matrix projection : ");


        get_timepoint();                                                    // sort triangles from back to front
            sort(vTrianglesToRaster.begin(), vTrianglesToRaster.end(), [](TTriangle& t1, TTriangle t2) {
                float z1 = (t1.V1.p.z + t1.V2.p.z + t1.V3.p.z) / 3.0f;
                float z2 = (t2.V1.p.z + t2.V2.p.z + t2.V3.p.z) / 3.0f;
                return z1 > z2;
             });
        print_diff_timepoint("Sorting time      : ");


        //ClrZBuffer();                                                     // Clear Screen and Z buffer
        Clear(BLACK);

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

