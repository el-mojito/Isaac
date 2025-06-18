#include <GLFW/glfw3.h>

void start_viewer();
void draw_points(double* x, double* y, double* vx, double* vy, int* n);

static GLFWwindow* window;

// Viridis RGB values
static const float viridis[6][3] = {
    {0.267f, 0.005f, 0.329f},
    {0.283f, 0.141f, 0.458f},
    {0.254f, 0.265f, 0.530f},
    {0.207f, 0.372f, 0.553f},
    {0.164f, 0.471f, 0.558f},
    {0.993f, 0.906f, 0.144f}
};

void viridis_colormap(float t, float* r, float* g, float* b) {
    t = fminf(fmaxf(t, 0.0f), 1.0f);
    float pos = t * 5.0f;
    int idx = (int)pos;
    if (idx >= 5) idx = 4;
    float frac = pos - idx;

    *r = (1 - frac) * viridis[idx][0] + frac * viridis[idx + 1][0];
    *g = (1 - frac) * viridis[idx][1] + frac * viridis[idx + 1][1];
    *b = (1 - frac) * viridis[idx][2] + frac * viridis[idx + 1][2];
}

void start_viewer() { 
    if (!glfwInit()) return;
    window = glfwCreateWindow(800, 800, "N-Body Viewer", NULL, NULL);
    if (!window) return;
    glfwMakeContextCurrent(window);
    glPointSize(2.0f);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1e16, 1e16, -1e16, 1e16, -1.0, 1.0); // match coordinate space
    glMatrixMode(GL_MODELVIEW);
}

void draw_points(double* x, double* y, double* vx, double* vy, int* n) {
    if (!window || glfwWindowShouldClose(window)) return;

    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    glBegin(GL_POINTS);

    for (int i = 0; i < *n; ++i) {
        float v = sqrtf((float)(vx[i] * vx[i] + vy[i] * vy[i]));
        float norm_v = fminf(fmaxf(v / 200000.0f, 0.0f), 1.0f);

        float r, g, b;
        viridis_colormap(norm_v, &r, &g, &b);
        glColor3f(r, g, b);
        glVertex2f((float)x[i], (float)y[i]);
    }

    glEnd();
    glFlush();
    glfwSwapBuffers(window);
    glfwPollEvents();
}
