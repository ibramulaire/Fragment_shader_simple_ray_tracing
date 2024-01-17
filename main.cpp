#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>

#include <GL/glew.h>

#include <GL/gl.h>
#include <GLFW/glfw3.h>

#include <math.h>
#include "shader.hpp"
#include <string.h>
#include "../glm/glm.hpp"
#include "../glm/gtc/matrix_transform.hpp"
using namespace glm;

void Resize_Window(GLFWwindow *window, int largeur, int hauteur);
void Input(GLFWwindow *window);
int Initialisation();

int largeur = 600;
int hauteur = 800;
const char *titre = "tuto1";
GLFWwindow *window;

GLuint MyShader;

using namespace std;
unsigned int VBO, VAO, EBO;

float vertices[] = {
    1.f, 1.f, 0.0f,   // top right
    1.f, -1.f, 0.0f,  // bottom right
    -1.f, -1.f, 0.0f, // bottom left
    -1.f, 1.f, 0.0f   // top left
};
unsigned int indices[] = {
    // note that we start from 0!
    0, 1, 3, // first triangle
    1, 2, 3  // second triangle
};

GLuint locCameraPosition;
vec3 cameraPosition(0., 0., 3.);

bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraPosX;
float cameraPosY;
float cameraDistance = 10.0;
float vitesse = 0.08;
float Scroolvitesse = 1;

double deltaTime = 0.0;
double previousTime = 0.0;

GLuint MatrixIDMVP, MatrixIDView, MatrixIDModel, MatrixIDPerspective;
glm::mat4 MVP;                     // justement la voilà
glm::mat4 Model, View, Projection; // Matrices constituant MVP
void Genere_gpu_data()
{
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);

    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3* sizeof(float)));
    // glEnableVertexAttribArray(1);

    // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

int main(int argc, char const *argv[])
{
    if (!Initialisation())
    {
        cout << "creation de fenetre echoué" << endl;
        return -1;
    }

    MyShader = LoadShaders("VS.vert", "FS.frag");

    MatrixIDModel = glGetUniformLocation(MyShader, "MODEL");

    MatrixIDMVP = glGetUniformLocation(MyShader, "MVP");
    MatrixIDView = glGetUniformLocation(MyShader, "VIEW");
    MatrixIDModel = glGetUniformLocation(MyShader, "MODEL");
    MatrixIDPerspective = glGetUniformLocation(MyShader, "PERSPECTIVE");

    /* on recupere l'ID */
    locCameraPosition = glGetUniformLocation(MyShader, "cameraPosition");

    Genere_gpu_data();

    glUseProgram(MyShader);

    int Large = glGetUniformLocation(MyShader, "Largeur");
    int longe = glGetUniformLocation(MyShader, "Hauteur");

    while (!glfwWindowShouldClose(window))
    {
        deltaTime = glfwGetTime() - previousTime;
        previousTime = glfwGetTime();

        Input(window);

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        Projection = glm::perspective(glm::radians(60.f), 1.0f, 1.0f, 1000.0f);

        // Model = glm::translate(Model,glm::vec3(cameraPosX,cameraPosY,0));

        // Model = glm::scale(Model,glm::vec3(.5, .8, .8));
/*
        View = glm::lookAt(vec3( vec4(cameraPosition, 1.0)), // Camera is at (0,0,3), in World Space
                           glm::vec3(0, 0, 0),                      // and looks at the origin
                           glm::vec3(0, 1, 0)                       // Head is up (set to 0,-1,0 to look upside-down)
        );
*/
        Model = glm::mat4(1.0f);
        Model = glm::mat4(1.0f);
        Model = glm::translate(Model, glm::vec3(0,0, -cameraDistance));
       
        Model = glm::rotate(Model, glm::radians(cameraAngleX), glm::vec3(1, 0, 0));
        Model = glm::rotate(Model, glm::radians(cameraAngleY), glm::vec3(0, 1, 0));
        
     
     //  std::cout<<cameraAngleX<<std::endl;

        MVP = Projection * View * Model;
        glUniform1f(Large, float(largeur));
        glUniform1f(longe, float(hauteur));

        glUniformMatrix4fv(MatrixIDMVP, 1, GL_FALSE, &MVP[0][0]);
        glUniformMatrix4fv(MatrixIDView, 1, GL_FALSE, &View[0][0]);
        glUniformMatrix4fv(MatrixIDModel, 1, GL_FALSE, &Model[0][0]);
        glUniformMatrix4fv(MatrixIDPerspective, 1, GL_FALSE, &Projection[0][0]);

        glUniform3f(locCameraPosition, cameraPosition.x, cameraPosition.y, cameraPosition.z);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glDeleteProgram(MyShader);

    glfwTerminate();

    return 0;
}

void Resize_Window(GLFWwindow *window, int larg, int haut)
{
    glViewport(0, 0, larg, haut);
    largeur = larg;
    hauteur = haut;
}

void mouse(GLFWwindow *window, int button, int state, int mods)
{
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    mouseX = x;
    mouseY = y;

    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if (state == GLFW_PRESS)
        {
            mouseLeftDown = true;
        }
        else if (state == GLFW_RELEASE)
            mouseLeftDown = false;
    }

    else if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
        if (state == GLFW_PRESS)
        {
            mouseRightDown = true;
        }
        else if (state == GLFW_RELEASE)
            mouseRightDown = false;
    }

    else if (button == GLFW_MOUSE_BUTTON_MIDDLE)
    {
        if (state == GLFW_PRESS)
        {
            mouseMiddleDown = true;
        }
        else if (state == GLFW_RELEASE)
            mouseMiddleDown = false;
    }
}
void mouseScroll(GLFWwindow *window, double x, double y)
{

    cameraDistance += y * Scroolvitesse;
}

void mouseMotion(GLFWwindow *window, double x, double y)
{

    if (mouseLeftDown)
    {
        cameraAngleY += (x - mouseX) * vitesse;
        cameraAngleX += (y - mouseY) * vitesse;
        mouseX = x;
        mouseY = y;
    }
    if (mouseRightDown)
    {
    }

    if (mouseMiddleDown)
    {
        cameraPosX += (x - mouseX) * vitesse;
        cameraPosY -= (y - mouseY) * vitesse;

        mouseX = x;
        mouseY = y;
    }
}

void Input(GLFWwindow *window)
{

    glfwSetCursorPosCallback(window, mouseMotion);
    // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetMouseButtonCallback(window, mouse);
    glfwSetScrollCallback(window, mouseScroll);

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

int Initialisation() // initialise une fenetrecontex et glew
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    window = glfwCreateWindow(largeur, hauteur, titre, NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return 0;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, Resize_Window);

    if (glewInit() != GLEW_OK)
    {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return 0;
    }

    // info version GLSL
    std::cout << "***** Info GPU *****" << std::endl;
    std::cout << "Fabricant : " << glGetString(GL_VENDOR) << std::endl;
    std::cout << "Carte graphique: " << glGetString(GL_RENDERER) << std::endl;
    std::cout << "Version : " << glGetString(GL_VERSION) << std::endl;
    std::cout << "Version GLSL : " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

    return 1;
}
