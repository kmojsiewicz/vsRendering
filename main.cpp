#define _USE_MATH_DEFINES
#include <vector>
#include <list>
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

        frontTexture.LoadFromBitmap("high.bmp");
        mesh.LoadFromObjectFile("obj\\spyro2.obj", &frontTexture, WHITE, 1.0);
        //mesh.MakeQube(nullptr, WHITE, 1.0);
        //frontTexture.LoadFromBitmap("negz.bmp");   mesh.triangles[0].SetTexture(&frontTexture);   mesh.triangles[1].SetTexture(&frontTexture);      // front
        //rightTexture.LoadFromBitmap("posx.bmp");   mesh.triangles[2].SetTexture(&rightTexture);   mesh.triangles[3].SetTexture(&rightTexture);      // right
        //backTexture.LoadFromBitmap("posz.bmp");    mesh.triangles[4].SetTexture(&backTexture);    mesh.triangles[5].SetTexture(&backTexture);       // back
        //leftTexture.LoadFromBitmap("negx.bmp");    mesh.triangles[6].SetTexture(&leftTexture);    mesh.triangles[7].SetTexture(&leftTexture);       // left
        //topTexture.LoadFromBitmap("posy.bmp");     mesh.triangles[8].SetTexture(&topTexture);     mesh.triangles[9].SetTexture(&topTexture);        // top
        //bottomTexture.LoadFromBitmap("negy.bmp");  mesh.triangles[10].SetTexture(&bottomTexture); mesh.triangles[11].SetTexture(&bottomTexture);    // bottom

        return true;
    }

    int Triangle_ClipAgainstPlane(TVec3d plane_p, TVec3d plane_n, TTriangle& in_tri, TTriangle& out_tri1, TTriangle& out_tri2)
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

        float d0 = dist(in_tri.V1.p);                                       // Get signed distance of each point in triangle to plane
        float d1 = dist(in_tri.V2.p);
        float d2 = dist(in_tri.V3.p);
        if (d0 >= 0) { 
            inside_points[nInsidePointCount++] = &in_tri.V1.p; 
            inside_tex[nInsideTexCount++] = &in_tri.V1.t;
        }
        else { 
            outside_points[nOutsidePointCount++] = &in_tri.V1.p; 
            outside_tex[nOutsideTexCount++] = &in_tri.V1.t;
        }
        if (d1 >= 0) { 
            inside_points[nInsidePointCount++] = &in_tri.V2.p; 
            inside_tex[nInsideTexCount++] = &in_tri.V2.t;
        }
        else { 
            outside_points[nOutsidePointCount++] = &in_tri.V2.p; 
            outside_tex[nOutsideTexCount++] = &in_tri.V2.t;
        }
        if (d2 >= 0) { 
            inside_points[nInsidePointCount++] = &in_tri.V3.p; 
            inside_tex[nInsideTexCount++] = &in_tri.V3.t;
        }
        else { 
            outside_points[nOutsidePointCount++] = &in_tri.V3.p; 
            outside_tex[nOutsideTexCount++] = &in_tri.V3.t;
        }

        // Now classify triangle points, and break the input triangle into 
        // smaller output triangles if required. There are four possible outcomes...
        if (nInsidePointCount == 0) {                                       // All points lie on the outside of plane, so clip whole triangle
            return 0;                                                       // No returned triangles are valid
        }

        if (nInsidePointCount == 3) {                                       // All points lie on the inside of plane, so do nothing
            out_tri1 = in_tri;                                              // and allow the triangle to simply pass through
            return 1;                                                       // Just the one returned original triangle is valid
        }

        if (nInsidePointCount == 1 && nOutsidePointCount == 2) {            // Triangle should be clipped. As two points lie outside the plane, the triangle simply becomes a smaller triangle
            out_tri1.light = in_tri.light;                                  // Copy appearance info to new triangle
            out_tri1.texture = in_tri.texture;
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

            out_tri1.V1.pixel = in_tri.V1.pixel;                            // out_tri1.V1.pixel = RED;
            out_tri1.V2.pixel = in_tri.V2.pixel;                            // out_tri1.V2.pixel = RED;
            out_tri1.V3.pixel = in_tri.V3.pixel;                            // out_tri1.V3.pixel = RED;
            return 1;                                                       // Return the newly formed single triangle
        }

        if (nInsidePointCount == 2 && nOutsidePointCount == 1) {            // Triangle should be clipped. As two points lie inside the plane, the clipped triangle becomes a "quad". Fortunately, we can
            out_tri1.light = in_tri.light;                                  // represent a quad with two new triangles
            out_tri2.light = in_tri.light;                                  // Copy appearance info to new triangles
            out_tri1.texture = in_tri.texture;
            out_tri2.texture = in_tri.texture;
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

            out_tri1.V1.pixel = in_tri.V1.pixel;                            // out_tri1.V1.pixel = GREEN;
            out_tri1.V2.pixel = in_tri.V2.pixel;                            // out_tri1.V2.pixel = GREEN;
            out_tri1.V3.pixel = in_tri.V3.pixel;                            // out_tri1.V3.pixel = GREEN;
            out_tri2.V1.pixel = in_tri.V1.pixel;                            // out_tri2.V1.pixel = BLUE;
            out_tri2.V2.pixel = in_tri.V2.pixel;                            // out_tri2.V2.pixel = BLUE;
            out_tri2.V3.pixel = in_tri.V3.pixel;                            // out_tri2.V3.pixel = BLUE;
            return 2;                                                       // Return two newly formed triangles which form a quad
        }

        return 0;
    }

    bool OnUserUpdate(float fElapsedTime) override
    {
        if (GetKey(VK_UP).bHeld)    vCamera.y += (GetKey(VK_SHIFT).bHeld) ? 16.0f * fElapsedTime : 8.0f * fElapsedTime;
        if (GetKey(VK_DOWN).bHeld)  vCamera.y -= (GetKey(VK_SHIFT).bHeld) ? 16.0f * fElapsedTime : 8.0f * fElapsedTime;
        if (GetKey(VK_LEFT).bHeld)  vCamera.x -= 8.0f * fElapsedTime;
        if (GetKey(VK_RIGHT).bHeld) vCamera.x += 8.0f * fElapsedTime;
        
        TVec3d vForward = Vector_Mul(vLookDir, 8.0f * fElapsedTime);
        TVec3d vForward2 = Vector_Mul(vForward, 2.0f);
        if (GetKey(L'W').bHeld) {
            vCamera = (GetKey(VK_SHIFT).bHeld) ? Vector_Add(vCamera, vForward2) : Vector_Add(vCamera, vForward);
            cout << vCamera.z << endl;
        }
        if (GetKey(L'S').bHeld) vCamera = (GetKey(VK_SHIFT).bHeld) ? Vector_Sub(vCamera, vForward2) : Vector_Sub(vCamera, vForward);
        if (GetKey(L'A').bHeld) fYaw -= 2.0f * fElapsedTime;
        if (GetKey(L'D').bHeld) fYaw += 2.0f * fElapsedTime;
        if (GetKey(L' ').bHeld) fTheta += 1.0f * fElapsedTime;
        if (GetKey(L'1').bHeld) frontTexture.LoadFromBitmap("high.bmp");

        //vCamera.z = 131;

        get_timepoint();
            TMat4x4 matRotZ = Matrix_MakeRotationZ(fTheta);                 // Rotation Z
            TMat4x4 matRotX = Matrix_MakeRotationX(fTheta * 0.5f);          // Rotation X
            TMat4x4 matTrans = Matrix_MakeTranslation(0.0f, 0.0f, 5.0f);
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
            TTriangle triTransformed, triProjected, triViewed = tri;        // get u and v data into projected triangles
            
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
                triViewed.light = max(0.5f, Vector_DotProduct(light_direction, normal));

                // Convert World Space --> View Space
                triViewed.V1.p = Matrix_MultiplyVector(matView, triTransformed.V1.p);
                triViewed.V2.p = Matrix_MultiplyVector(matView, triTransformed.V2.p);
                triViewed.V3.p = Matrix_MultiplyVector(matView, triTransformed.V3.p);

                // Clip Viewed Triangle against near plane, this could form two additional
                // additional triangles. 
                int nClippedTriangles = 0;
                TTriangle clipped[2];
                nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);

                for (int n = 0; n < nClippedTriangles; n++)                 // We may end up with multiple triangles form the clip, so project as required
                {
                    triProjected = clipped[n];

                    // Project triangles from 3D --> 2D with normalising into cartesian space
                    triProjected.V1.p = Matrix_MultiplyVector(matProj, clipped[n].V1.p);
                    triProjected.V2.p = Matrix_MultiplyVector(matProj, clipped[n].V2.p);
                    triProjected.V3.p = Matrix_MultiplyVector(matProj, clipped[n].V3.p);

                    triProjected.V1.t.u = triProjected.V1.t.u / triProjected.V1.p.w;
                    triProjected.V2.t.u = triProjected.V2.t.u / triProjected.V2.p.w;
                    triProjected.V3.t.u = triProjected.V3.t.u / triProjected.V3.p.w;
                    triProjected.V1.t.v = triProjected.V1.t.v / triProjected.V1.p.w;
                    triProjected.V2.t.v = triProjected.V2.t.v / triProjected.V2.p.w;
                    triProjected.V3.t.v = triProjected.V3.t.v / triProjected.V3.p.w;
                    triProjected.V1.t.w = 1.0f / triProjected.V1.p.w;
                    triProjected.V2.t.w = 1.0f / triProjected.V2.p.w;
                    triProjected.V3.t.w = 1.0f / triProjected.V3.p.w;
                    triProjected.V1.p = Vector_Div(triProjected.V1.p, triProjected.V1.p.w);
                    triProjected.V2.p = Vector_Div(triProjected.V2.p, triProjected.V2.p.w);
                    triProjected.V3.p = Vector_Div(triProjected.V3.p, triProjected.V3.p.w);

                    triProjected.V1.p.x *= -1.0f;                           // X/Y are inverted so put them back
                    triProjected.V2.p.x *= -1.0f;
                    triProjected.V3.p.x *= -1.0f;
                    triProjected.V1.p.y *= -1.0f;
                    triProjected.V2.p.y *= -1.0f;
                    triProjected.V3.p.y *= -1.0f;

                    triProjected.Translate({ 1.0f, 1.0f, 0.0f });           // Offset verts into visible normalised space
                    triProjected.Scale(0.5f * (float)ScreenWidth(), 0.5f * (float)ScreenHeight(), 0.5f * 100.0);

                    vTrianglesToRaster.push_back(triProjected);
                }
            }
        }
        print_diff_timepoint("Matrix projection : ");


        get_timepoint();                                                    // sort triangles from front to back
            sort(vTrianglesToRaster.begin(), vTrianglesToRaster.end(), [](TTriangle& t1, TTriangle t2) {
                float z1 = (t1.V1.p.z + t1.V2.p.z + t1.V3.p.z) / 3.0f;
                float z2 = (t2.V1.p.z + t2.V2.p.z + t2.V3.p.z) / 3.0f;
                return z1 < z2;
             });
        print_diff_timepoint("Sorting time      : ");


        ClrZBuffer();                                                       // Clear Screen and Z buffer
        Clear(BLUE);

        get_timepoint();
        for (auto& triToRaster : vTrianglesToRaster) {                      // Loop through all transformed, viewed, projected, and sorted triangles
            // Clip triangles against all four screen edges, this could yield
            // a bunch of triangles, so create a queue that we traverse to 
            //  ensure we only test new triangles generated against planes
            TTriangle clipped[2];
            list<TTriangle> listTriangles;
            listTriangles.push_back(triToRaster);                           // Add initial triangle
            int nNewTriangles = 1;
            for (int p=0; p<4; p++) {
                int nTrisToAdd = 0;
                while (nNewTriangles > 0) {
                    TTriangle test = listTriangles.front();                 // Take triangle from front of queue
                    listTriangles.pop_front();
                    nNewTriangles--;
                    // Clip it against a plane. We only need to test each 
                    // subsequent plane, against subsequent new triangles
                    // as all triangles after a plane clip are guaranteed
                    // to lie on the inside of the plane. I like how this
                    // comment is almost completely and utterly justified
                    switch (p) {
                    case 0:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    }
                    // Clipping may yield a variable number of triangles, so
                    // add these new ones to the back of the queue for subsequent
                    // clipping against next planes
                    for (int w = 0; w < nTrisToAdd; w++) listTriangles.push_back(clipped[w]);
                }
                nNewTriangles = (int)listTriangles.size();
            }
            for (auto& t : listTriangles) {                                 // Draw the transformed, viewed, clipped, projected, sorted, clipped triangles
                //FillTriangle(t);                                            // Rasterize triangle
                //DrawTriangle(t, t.V1.pixel);
                //TexturedTriangle(t.V1.p.x, t.V1.p.y, t.V1.t.u, t.V1.t.v, t.V1.t.w,
                //    t.V2.p.x, t.V2.p.y, t.V2.t.u, t.V2.t.v, t.V2.t.w,
                //    t.V3.p.x, t.V3.p.y, t.V3.t.u, t.V3.t.v, t.V3.t.w, t.texture);
                TexturedTriangle(&t);
            }
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

