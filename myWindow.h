#pragma once
#include <gl/glut.h>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "InterpolationCurve.h"

using namespace std;
using namespace Eigen;

class myWindow
{
    typedef enum {
        MOVE,
        DELPT,
        INSERTPT
    }EDITFUNC;

public:
    myWindow();
    ~myWindow();

    void myInit();
    void myDisplay();
    void reshape(int w, int h);

    // motion function will not run if the mouse button is in up state
    void myMotion(int x, int y);
    void passiveMotion(int x, int y);

    void entryFunc(int state);

    void myMouse(int button, int state, int x, int y);
    void myKey(unsigned char key, int x, int y);
    void mySpecialKey(GLint key, GLint x, GLint y);
    void idle();


    GLint getImageWidth() const { return imagewidth; }
    GLint getImageHeight() const { return imageheight; }

private:

    void selectFont(int size, int charset, const char * face);
    void drawCNString(const char* str);
    void XPrintString(const char *s);
    void coordsToStr(int x, int y, string * str);
    pair<double, double> strToCoords(const string& str);

    void showMinDist();

    void readImageFile();
    /* 函数grab
    * 抓取窗口中的像素
    */
    void grab(void);
    pair<int, int> findMoveIndex(const vector<InterpolationCurve>& crvVec, double x, double y);
    void offsetCurve(double lamda);


    //static myWindow* pWindow;

    GLuint TextFont;

    vector<InterpolationCurve> m_crvVec;
    GLint imagewidth;
    GLint imageheight;
    GLint pixellength;
    GLubyte* pixeldata;


    
    double modelMat[16], projMat[16];
    GLfloat yOffset = 0, xOffset = 0, zoomRatio = 1.f, oOffset = 0.f;

    bool toShowCp = false, toShowInter = true, toShowDer = false, toShowHull = false, toBiHe = false, toOffset = false;
    bool toEdit = false, bCtrl = false;

    pair<int, int> crv_pt_idxs = make_pair(-1, -1); // first indicate curveid, second indicate pointid.

    std::vector<double> chosenPt, hoverPt;
    string coordString;

    std::vector<double> m_minDistPt;
    InterpolationCurve m_shadowCrv;
};