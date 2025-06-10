#include <GLFW/glfw3.h>

void start_viewer();
void draw_points(double* x, double* y, int* n);

static GLFWwindow* window;

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

void draw_points(double* x, double* y, int* n) {
    if (!window || glfwWindowShouldClose(window)) return;

    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    glBegin(GL_POINTS);
    for (int i = 0; i < *n; ++i) {
        glVertex2f((float)x[i], (float)y[i]);
    }
    glEnd();
    glfwSwapBuffers(window);
    glfwPollEvents();
}
