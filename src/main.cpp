#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>

#define MAX_SIZE 1000
#define PI 3.141592
#define X_AXIS { 1.0, 0.0, 0.0 }
#define Y_AXIS { 0.0, 1.0, 0.0 }
#define Z_AXIS { 0.0, 0.0, 1.0 }

char BUF[MAX_SIZE][MAX_SIZE];
double Z_BUF[MAX_SIZE][MAX_SIZE];
const char TEXTURE[13] = { ' ', '.', ',', '-', '~', ':', ';', '=', '!', '*', '#', '$', '@' }; 

struct Vertex {
    double x, y, z;
};

Vertex operator+(const Vertex& lhs, const Vertex& rhs) {
    return { lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
}

Vertex operator-(const Vertex& lhs, const Vertex& rhs) {
    return { lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
}

Vertex operator*(const Vertex& lhs, const Vertex& rhs) {
    return { lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z };
}

Vertex operator/(const Vertex& lhs, const Vertex& rhs) {
    return { lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z };
}

Vertex scale(const Vertex& v, const double s) {
    return { v.x * s, v.y * s, v.z * s };
}

double dot(const Vertex& lhs, const Vertex& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vertex cross(const Vertex& lhs, const Vertex& rhs) {
    return { 
        lhs.y * rhs.z - rhs.y * lhs.z, 
        -(lhs.x * rhs.z - lhs.z * rhs.x), 
        lhs.x * rhs.y - lhs.y * rhs.x 
    };
}

Vertex normalize(const Vertex& v) {
    double size = std::sqrt(dot(v, v));
    return { v.x / size, v.y / size, v.z / size };
}

struct Matrix {
    Vertex u;
    Vertex v;
    Vertex n;
};

Matrix transpose(const Matrix& m) {
    return {
        { m.u.x, m.v.x, m.n.x },
        { m.u.y, m.v.y, m.n.y },
        { m.u.z, m.v.z, m.n.z }
    };
}

Vertex operator*(const Matrix& m, const Vertex& v) {
    Matrix m_T = transpose(m);
    return { dot(m_T.u, v), dot(m_T.v, v), dot(m_T.n, v) };
}

Matrix operator*(const Matrix& m1, const Matrix& m2) {
    Matrix m1_T = transpose(m1);
    return {
        { dot(m1_T.u, m2.u), dot(m1_T.v, m2.u), dot(m1_T.n, m2.u) },
        { dot(m1_T.u, m2.v), dot(m1_T.v, m2.v), dot(m1_T.n, m2.v) },
        { dot(m1_T.u, m2.n), dot(m1_T.v, m2.n), dot(m1_T.n, m2.n) }
    };
}

int main(int argc, char* argv[]) {
    double minor_radius = 1.0;
    double major_radius = 2.0;
    double sample_rate = 0.005;

    double rotation_period = 2.0;
    Vertex rotation_axis = { 1.0, 1.0, 1.0 };

    Vertex camera_pos = { 0.0, 0.0, 5.0 };
    double near_plane = -1.0;

    Vertex lighting = { 1.0, 0.0, 1.0 };

    int viewport_width = 50;
    int viewport_height = 50;
    while (true) {
        memset(BUF, 0, sizeof(BUF));
        memset(Z_BUF, -1, sizeof(Z_BUF));

        double theta = 2 * PI * ((double)std::clock() / CLOCKS_PER_SEC) / rotation_period;
        double sint = std::sin(theta);
        double cost = std::cos(theta);

        Vertex n = normalize(rotation_axis);
        Vertex u = normalize(cross(Y_AXIS, n));
        Vertex v = normalize(cross(n, u));

        Matrix R = { u, v, n };
        Matrix R_z = { 
            {  cost, sint, 0.0 },
            { -sint, cost, 0.0 },
            {   0.0,  0.0, 1.0 } 
        };
        Matrix tmp = R_z * transpose(R);
        Matrix world_matrix = R * R_z * transpose(R);
        
        for (double i = 0; i < 2 * PI; i += 2 * PI * sample_rate) {
            for (double j = 0; j < 2 * PI; j += 2 * PI * sample_rate) {
                double sini = std::sin(i);
                double cosi = std::cos(i);
                double sinj = std::sin(j);
                double cosj = std::cos(j);

                Vertex center = { sini, 0.0, cosi };
                Vertex normal = normalize(scale(center, cosj) + Vertex{ 0.0, sinj, 0.0 });
                Vertex vertex = scale(center, major_radius) + scale(normal, minor_radius);

                vertex = world_matrix * vertex - camera_pos;
                normal = normalize(world_matrix * normal);

                double proj_x = (vertex.x * -near_plane) / -vertex.z;
                double proj_y = (vertex.y * -near_plane) / -vertex.z;

                int pos_x = proj_x * (viewport_width / 2.0) + (viewport_width / 2.0);
                int pos_y = proj_y * (viewport_height / 2.0) + (viewport_height / 2.0);

                lighting = normalize(lighting);
                double diff_term = dot(lighting, normal);

                if (pos_x < 0 || pos_x >= viewport_width || pos_y < 0 || pos_y >= viewport_height) continue;
                if (Z_BUF[pos_x][pos_y] > vertex.z) continue;

                BUF[pos_x][pos_y] = TEXTURE[std::max(0, int(std::round(diff_term * 12.0)))];
                Z_BUF[pos_x][pos_y] = vertex.z;
            }
        }
        std::system("clear");
        for (int i = viewport_height - 1; i >= 0; --i) {
            for (int j = 0; j < viewport_width; ++j) {
                std::putchar(BUF[j][i] ? BUF[j][i] : ' ');
                std::putchar(' ');
            }
            std::putchar('\n');
        }
    }
}