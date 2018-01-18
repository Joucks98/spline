#include <conio.h>
#include "myWindow.h"
#include "NurbsBase.h"

//#pragma comment( linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

static GLenum doubleBuffer = GL_TRUE;
myWindow* pWindow = nullptr;

void displayFunc()
{
    pWindow->myDisplay();
}
void reshapeFunc(int w, int h)
{
    pWindow->reshape(w, h);
}

void motionFunc(int x, int y)
{
    pWindow->myMotion(x, y);
}
void passiveMotionFunc(int x, int y)
{
    pWindow->passiveMotion(x, y);
}

void entryFunc(int state)
{
    pWindow->entryFunc(state);
}

void mouseFunc(int button, int state, int x, int y)
{
    pWindow->myMouse(button, state, x, y);
}
void keyFunc(unsigned char key, int x, int y)
{
    pWindow->myKey(key, x, y);
}
void specialKeyFunc(GLint key, GLint x, GLint y)
{
    pWindow->mySpecialKey(key, x, y);
}
void idleFunc()
{
    pWindow->idle();
}

void Args(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-sb") == 0) doubleBuffer = GL_FALSE;
        if (strcmp(argv[i], "-db") == 0) doubleBuffer = GL_TRUE;
    }
}
/////////////////main///////////////
void main(int argc, char **argv)
{
    Args(argc, argv);
    glutInit(&argc, argv);
    GLenum type = GLUT_RGBA;
    type |= (doubleBuffer ? GLUT_DOUBLE : GLUT_SINGLE);
    //type |= GLUT_DEPTH;
    glutInitDisplayMode(type);
    glutInitWindowSize(500, 500/*mw.getImageWidth(), mw.getImageHeight()*/);
    glutInitWindowPosition(150, 50);
    glutCreateWindow("B Spline interpolation");


    //int menuID = glutCreateMenu(ProcessMenu);
    //glutAddMenuEntry("Save Image", 1);
    //glutAddMenuEntry("Flip", 2);
    //glutAddMenuEntry("zoom pixel fill window", 3);
    //glutAddMenuEntry("Just Red", 4);
    //glutAddMenuEntry("Just Green", 5);
    //glutAddMenuEntry("Just Blue", 6);
    //glutAddMenuEntry("black & white", 7);
    //glutAddMenuEntry("invert map", 8);
    //glutAttachMenu(GLUT_RIGHT_BUTTON);

    myWindow mw;
    pWindow = &mw;

    glutDisplayFunc(displayFunc);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);
    glutReshapeFunc(reshapeFunc);
    glutKeyboardFunc(keyFunc);
    glutSpecialFunc(specialKeyFunc);
    glutPassiveMotionFunc(passiveMotionFunc);
    glutEntryFunc(entryFunc);
    glutIdleFunc(idleFunc);
    glutMainLoop();

    
    _getch();
}