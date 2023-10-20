#define _USE_MATH_DEFINES
#include <vector>
#include <math.h>
#include <string>
#include <fstream>
#include <strstream>
#include <algorithm>

#include "platform.h"
#include "bmploader.h"
#include "texture.h"

#define WND_WIDTH   800
#define WND_HEIGHT  600

using namespace std;

static float z_buffer[WND_WIDTH][WND_HEIGHT];                  // z-buffer uses heap data


float q_rsqrt(float number)
{
    long i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = number * 0.5F;
    y = number;
    i = *(long*)&y;                       // evil floating point bit level hacking
    i = 0x5f3759df - (i >> 1);               // what the fuck?
    y = *(float*)&i;
    y = y * (threehalfs - (x2 * y * y));   // 1st iteration
    //y = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

    return y;
}


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

    void Translate(float x, float y, float z) noexcept
    {
        V1.x += x; V2.x += x; V3.x += x;
        V1.y += y; V2.y += y; V3.y += y;
        V1.z += z; V2.z += z; V3.z += z;
    }

    void Scale(float x, float y, float z) noexcept
    {
        V1.x *= x; V2.x *= x; V3.x *= x;
        V1.y *= y; V2.y *= y; V3.y *= y;
        V1.z *= z; V2.z *= z; V3.z *= z;
    }

    TVec3d NormalVector() noexcept
    {
        TVec3d normal, line1, line2;

        line1.x = V2.x - V1.x;
        line1.y = V2.y - V1.y;
        line1.z = V2.z - V1.z;

        line2.x = V3.x - V1.x;
        line2.y = V3.y - V1.y;
        line2.z = V3.z - V1.z;

        normal.x = line1.y * line2.z - line1.z * line2.y;
        normal.y = line1.z * line2.x - line1.x * line2.z;
        normal.z = line1.x * line2.y - line1.y * line2.x;

        float inv_sqrt_nl = q_rsqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        normal.x *= inv_sqrt_nl;
        normal.y *= inv_sqrt_nl;
        normal.z *= inv_sqrt_nl;

        return normal;
    }
};

Texture defaultTexture;

struct TMesh
{
    vector<TTriangle> triangles;

    bool LoadFromObjectFile(std::string sfilename) {

        defaultTexture.set(128, 128);
        defaultTexture.clear(TRGBColor(255, 255, 255, 255));

        vector<TVertex> verts;
        char c;
        std::string line, key;
        std::ifstream ifs(sfilename.c_str(), std::ifstream::in);
        while (ifs.good() && !ifs.eof() && std::getline(ifs, line)) {
            key = "";
            strstream s;
            s << line;
            s >> key >> std::ws;

            if (key == "v") {                                               // vertex
                TVertex v;
                s >> v.x >> v.y >> v.z;
                v.u = 0.0;
                v.v = 0.0;
                verts.push_back(v);
            }
            else if (key == "vp") {                                         // parameter
            }
            else if (key == "vt") {                                         // texture coordinate
            }
            else if (key == "f") {                                          // face
                int f[3], v[3], t[3], n[3];
                while (!s.eof()) {
                    s >> f[0] >> std::ws;
                    if (s.peek() == '/') {
                        s.get();
                        if (s.peek() == '/') {
                            s.get();
                            s >> n[0] >> std::ws;
                        }
                        else {
                            s >> t[0] >> std::ws;
                            if (s.peek() == '/') {
                                s.get();
                                s >> n[0] >> std::ws;
                            }
                        }
                    }
                    while (!s.eof()) {
                        s >> f[1] >> std::ws;
                        if (s.peek() == '/') {
                            s.get();
                            if (s.peek() == '/') {
                                s.get();
                                s >> n[1] >> std::ws;
                            }
                            else {
                                s >> t[1] >> std::ws;
                                if (s.peek() == '/') {
                                    s.get();
                                    s >> n[1] >> std::ws;
                                }
                            }
                        }
                        while (!s.eof()) {
                            s >> f[2] >> std::ws;
                            if (s.peek() == '/') {
                                s.get();
                                if (s.peek() == '/') {
                                    s.get();
                                    s >> n[2] >> std::ws;
                                }
                                else {
                                    s >> t[2] >> std::ws;
                                    if (s.peek() == '/') {
                                        s.get();
                                        s >> n[2] >> std::ws;
                                    }
                                }
                            }
                            triangles.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1], &defaultTexture , 1.0 });
                        }
                    }

                }
            }
        }
        ifs.close();

        return true;
    }
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
    bool printDurTimeOnce = true;
    std::chrono::high_resolution_clock::time_point tp_start, tp_end;
    std::chrono::duration<float> elapsed_sec;
#define get_timepoint()           { if (printDurTimeOnce) { tp_start = std::chrono::high_resolution_clock::now(); } }
#define print_diff_timepoint(...) { if (printDurTimeOnce) { tp_end = std::chrono::high_resolution_clock::now(); \
                                      elapsed_sec = tp_end - tp_start; \
                                      cout << __VA_ARGS__ << std::chrono::duration_cast<std::chrono::microseconds>(elapsed_sec).count() << " us" << endl; } }



    Texture leftTexture, topTexture, rightTexture, bottomTexture, frontTexture, backTexture;

    TMesh mesh;
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

    inline void putPixel(int x, int y, float z, TRGBColor col, float light)
    {
        if (x < 0) x = 0;
        if (x > (WND_WIDTH - 1)) x = WND_WIDTH - 1;
        if (y < 0) y = 0;
        if (y > (WND_HEIGHT - 1)) y = WND_HEIGHT - 1;

        //Pixel p(col.R * light, col.G * light, col.B * light);
        Pixel p(255 * light, 255 * light, 255 * light);
        CLayer* pDrawTarget = GetDrawTarget();

        //if (z_buffer[x][y] > z) {
        //    z_buffer[x][y] = z;
        pDrawTarget->SetPixel(x, y, p);
        //}
    }
    bool Draw(int32_t x, int32_t y, Pixel p)
    {
        CLayer* pDrawTarget = GetDrawTarget();
        pDrawTarget->SetPixel(x, y, p);
        return true;
    }

    void DrawLine(int32_t x1, int32_t y1, int32_t x2, int32_t y2, Pixel p = WHITE, uint32_t pattern = 0xFFFFFFFF)
    {
        int x, y, dx, dy, dx1, dy1, px, py, xe, ye, i;
        dx = x2 - x1; dy = y2 - y1;

        auto rol = [&](void) { pattern = (pattern << 1) | (pattern >> 31); return pattern & 1; };

        //olc::vi2d p1(x1, y1), p2(x2, y2);
        //if (!ClipLineToScreen(p1, p2))
        //	return;
        //x1 = p1.x; y1 = p1.y;
        //x2 = p2.x; y2 = p2.y;

        // straight lines idea by gurkanctn
        if (dx == 0) // Line is vertical
        {
            if (y2 < y1) std::swap(y1, y2);
            for (y = y1; y <= y2; y++) if (rol()) Draw(x1, y, p);
            return;
        }

        if (dy == 0) // Line is horizontal
        {
            if (x2 < x1) std::swap(x1, x2);
            for (x = x1; x <= x2; x++) if (rol()) Draw(x, y1, p);
            return;
        }

        // Line is Funk-aye
        dx1 = abs(dx); dy1 = abs(dy);
        px = 2 * dy1 - dx1;	py = 2 * dx1 - dy1;
        if (dy1 <= dx1)
        {
            if (dx >= 0)
            {
                x = x1; y = y1; xe = x2;
            }
            else
            {
                x = x2; y = y2; xe = x1;
            }

            if (rol()) Draw(x, y, p);

            for (i = 0; x < xe; i++)
            {
                x = x + 1;
                if (px < 0)
                    px = px + 2 * dy1;
                else
                {
                    if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0)) y = y + 1; else y = y - 1;
                    px = px + 2 * (dy1 - dx1);
                }
                if (rol()) Draw(x, y, p);
            }
        }
        else
        {
            if (dy >= 0)
            {
                x = x1; y = y1; ye = y2;
            }
            else
            {
                x = x2; y = y2; ye = y1;
            }

            if (rol()) Draw(x, y, p);

            for (i = 0; y < ye; i++)
            {
                y = y + 1;
                if (py <= 0)
                    py = py + 2 * dx1;
                else
                {
                    if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0)) x = x + 1; else x = x - 1;
                    py = py + 2 * (dx1 - dy1);
                }
                if (rol()) Draw(x, y, p);
            }
        }
    }

    void drawTriangle(TTriangle t)
    {
        DrawLine(t.V1.x, t.V1.y, t.V2.x, t.V2.y);
        DrawLine(t.V2.x, t.V2.y, t.V3.x, t.V3.y);
        DrawLine(t.V1.x, t.V1.y, t.V3.x, t.V3.y);
    }

    void fillTriangle(TTriangle t, Pixel p)
    {
        auto drawline = [&](int sx, int ex, int ny) { for (int i = sx; i <= ex; i++) Draw(i, ny, p); };

        int t1x, t2x, y, minx, maxx, t1xp, t2xp;
        bool changed1 = false;
        bool changed2 = false;
        int signx1, signx2, dx1, dy1, dx2, dy2;
        int e1, e2;
        // Sort vertices
        if (t.V1.y > t.V2.y) { std::swap(t.V1, t.V2); }
        if (t.V1.y > t.V3.y) { std::swap(t.V1, t.V3); }
        if (t.V2.y > t.V3.y) { std::swap(t.V2, t.V3); }

        int x1, y1, x2, y2, x3, y3;
        x1 = t.V1.x; y1 = t.V1.y;
        x2 = t.V2.x; y2 = t.V2.y;
        x3 = t.V3.x; y3 = t.V3.y;

        t1x = t2x = x1; y = y1;   // Starting points
        dx1 = (int)(x2 - x1);
        if (dx1 < 0) { dx1 = -dx1; signx1 = -1; }
        else signx1 = 1;
        dy1 = (int)(y2 - y1);

        dx2 = (int)(x3 - x1);
        if (dx2 < 0) { dx2 = -dx2; signx2 = -1; }
        else signx2 = 1;
        dy2 = (int)(y3 - y1);

        if (dy1 > dx1) { std::swap(dx1, dy1); changed1 = true; }
        if (dy2 > dx2) { std::swap(dy2, dx2); changed2 = true; }

        e2 = (int)(dx2 >> 1);
        // Flat top, just process the second half
        if (y1 == y2) goto next;
        e1 = (int)(dx1 >> 1);

        for (int i = 0; i < dx1;) {
            t1xp = 0; t2xp = 0;
            if (t1x < t2x) { minx = t1x; maxx = t2x; }
            else { minx = t2x; maxx = t1x; }
            // process first line until y value is about to change
            while (i < dx1) {
                i++;
                e1 += dy1;
                while (e1 >= dx1) {
                    e1 -= dx1;
                    if (changed1) t1xp = signx1;//t1x += signx1;
                    else          goto next1;
                }
                if (changed1) break;
                else t1x += signx1;
            }
            // Move line
        next1:
            // process second line until y value is about to change
            while (1) {
                e2 += dy2;
                while (e2 >= dx2) {
                    e2 -= dx2;
                    if (changed2) t2xp = signx2;//t2x += signx2;
                    else          goto next2;
                }
                if (changed2)     break;
                else              t2x += signx2;
            }
        next2:
            if (minx > t1x) minx = t1x;
            if (minx > t2x) minx = t2x;
            if (maxx < t1x) maxx = t1x;
            if (maxx < t2x) maxx = t2x;
            drawline(minx, maxx, y);    // Draw line from min to max points found on the y
            // Now increase y
            if (!changed1) t1x += signx1;
            t1x += t1xp;
            if (!changed2) t2x += signx2;
            t2x += t2xp;
            y += 1;
            if (y == y2) break;
        }
    next:
        // Second half
        dx1 = (int)(x3 - x2); if (dx1 < 0) { dx1 = -dx1; signx1 = -1; }
        else signx1 = 1;
        dy1 = (int)(y3 - y2);
        t1x = x2;

        if (dy1 > dx1) {   // swap values
            std::swap(dy1, dx1);
            changed1 = true;
        }
        else changed1 = false;

        e1 = (int)(dx1 >> 1);

        for (int i = 0; i <= dx1; i++) {
            t1xp = 0; t2xp = 0;
            if (t1x < t2x) { minx = t1x; maxx = t2x; }
            else { minx = t2x; maxx = t1x; }
            // process first line until y value is about to change
            while (i < dx1) {
                e1 += dy1;
                while (e1 >= dx1) {
                    e1 -= dx1;
                    if (changed1) { t1xp = signx1; break; }//t1x += signx1;
                    else          goto next3;
                }
                if (changed1) break;
                else   	   	  t1x += signx1;
                if (i < dx1) i++;
            }
        next3:
            // process second line until y value is about to change
            while (t2x != x3) {
                e2 += dy2;
                while (e2 >= dx2) {
                    e2 -= dx2;
                    if (changed2) t2xp = signx2;
                    else          goto next4;
                }
                if (changed2)     break;
                else              t2x += signx2;
            }
        next4:

            if (minx > t1x) minx = t1x;
            if (minx > t2x) minx = t2x;
            if (maxx < t1x) maxx = t1x;
            if (maxx < t2x) maxx = t2x;
            drawline(minx, maxx, y);
            if (!changed1) t1x += signx1;
            t1x += t1xp;
            if (!changed2) t2x += signx2;
            t2x += t2xp;
            y += 1;
            if (y > y3) return;
        }
    }

#define SUB_PIX(a) (ceil(a)-a)
    void fillTriangleT(TTriangle t)
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

        mesh.LoadFromObjectFile("test.obj");

        //mesh.triangles = {

        //    // FRONT
        //    { { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 0.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 0.0f },  &frontTexture },
        //    { { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 0.0f, 1.0f, 1.0f },  &frontTexture },

        //    // RIGHT
        //    { { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 0.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },  &rightTexture },
        //    { { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 1.0f, 1.0f, 1.0f },  &rightTexture },

        //    // BACK
        //    { { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 1.0f, 1.0f, 1.0f, 0.0f },  &backTexture },
        //    { { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 0.0f, 0.0f, 1.0f, 1.0f, 1.0f },  &backTexture },

        //    // LEFT
        //    { { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 1.0f, 0.0f, 1.0f, 0.0f },  &leftTexture },
        //    { { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 0.0f, 1.0f, 0.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f },  &leftTexture },

        //    // TOP
        //    { { 0.0f, 1.0f, 0.0f, 0.0f, 1.0f },    { 0.0f, 1.0f, 1.0f, 0.0f, 0.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },  &topTexture },
        //    { { 0.0f, 1.0f, 0.0f, 0.0f, 1.0f },    { 1.0f, 1.0f, 1.0f, 1.0f, 0.0f },    { 1.0f, 1.0f, 0.0f, 1.0f, 1.0f },  &topTexture },

        //    // BOTTOM
        //    { { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 0.0f, 1.0f, 0.0f, 0.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 0.0f },  &bottomTexture },
        //    { { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f },    { 0.0f, 0.0f, 0.0f, 1.0f, 0.0f },    { 1.0f, 0.0f, 0.0f, 1.0f, 1.0f },  &bottomTexture },
        //};

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
        std::chrono::milliseconds elapsedMs;
        vector<TTriangle> vTrianglesToRaster;
        TMat4x4 matRotX, matRotZ, matRotZX;                                 // Set up rotation matrices

        Clear(BLACK);
        clrZBuffer();                                                       // Clear Screen and Z buffer

        fTheta += 1.0f * fElapsedTime;

        get_timepoint();

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

        MultiplyMatrixByMatrix(matRotZX, matRotZ, matRotX);                 // Rotation in Z-Axis and X-Axis

        print_diff_timepoint("Rotation time     : ");

        get_timepoint();

        for (auto tri : mesh.triangles)                                     // Draw Triangles
        {
            TVec3d normal;
            TTriangle triProjected, triTranslated, triRotatedZX;

            triProjected = tri;                                             // get u and v data into projected triangles

            MultiplyMatrixVector(tri.V1, triRotatedZX.V1, matRotZX);
            MultiplyMatrixVector(tri.V2, triRotatedZX.V2, matRotZX);
            MultiplyMatrixVector(tri.V3, triRotatedZX.V3, matRotZX);

            // Offset into the screen
            triTranslated = triRotatedZX;
            triTranslated.V1.z = triRotatedZX.V1.z + 8.0f;
            triTranslated.V2.z = triRotatedZX.V2.z + 8.0f;
            triTranslated.V3.z = triRotatedZX.V3.z + 8.0f;

            normal = triTranslated.NormalVector();

            if (normal.x * (triTranslated.V1.x - vCamera.x) +
                normal.y * (triTranslated.V1.y - vCamera.y) +
                normal.z * (triTranslated.V1.z - vCamera.z) < 0.0)
            {
                // Illumination
                TVec3d light_direction = { 0.0f, 0.0f, -1.0f };

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

        get_timepoint();
        for (auto triProjected : vTrianglesToRaster)                        // Draw Triangles sorted
        {
            fillTriangle(triProjected, WHITE * triProjected.light);         // Rasterize triangle
            //fillTriangleT(triProjected);
        }
        print_diff_timepoint("Rendering time 1  : ");

        get_timepoint();
        for (auto triProjected : vTrianglesToRaster)                        // Draw Triangles sorted
        {
            //fillTriangle(triProjected, WHITE * triProjected.light);         // Rasterize triangle
            fillTriangleT(triProjected);
        }
        print_diff_timepoint("Rendering time 2  : ");

        printDurTimeOnce = false;
        return true;
    }
};

int WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR CmdLine, int nShowCmd)
{
    Engine3D demo;

    if (demo.Construct(WND_WIDTH, WND_HEIGHT, false, false, true)) {
        //if (demo.Construct(WND_WIDTH, WND_HEIGHT, 1, 1, false, false)) {
        demo.Start();
    }

    return 0;
}

