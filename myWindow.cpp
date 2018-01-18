#include "myWindow.h"

#define BMP_Header_Length 54
#define MAX_CHAR    128

//myWindow* myWindow::pWindow = nullptr;

myWindow::myWindow():chosenPt(2)
{
    //pWindow = this;
    readImageFile();
    myInit();
    
}

myWindow::~myWindow()
{
    if (pixeldata != nullptr)
        free(pixeldata);
}

void myWindow::myInit()
{
    glClearColor(0.0, 0.5, 0.5, 1.0);
    //glDepthFunc(GL_LEQUAL);
    //glEnable(GL_DEPTH_TEST);
    //glEnable(GL_AUTO_NORMAL);
    //glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    //glShadeModel(GL_SMOOTH);    

    //����MAX_CHAR����������ʾ�б����
    TextFont = glGenLists(MAX_CHAR);
    //��ÿ���ַ��Ļ������װ����Ӧ����ʾ�б���
    wglUseFontBitmaps(wglGetCurrentDC(), 0, MAX_CHAR, TextFont);
}

void myWindow::myDisplay()
{
    GLint viewPort[4];
    glGetIntegerv(GL_VIEWPORT, viewPort);
    //glViewport(0, 0, (GLsizei)viewPort[2] * zoomRatio, (GLsizei)viewPort[3] * zoomRatio);

    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0, yOffset, 0);
    glTranslatef(xOffset, 0, 0);
    //glTranslatef(0, 0, -10);

    glScalef(zoomRatio, zoomRatio, zoomRatio);
    //glDrawPixels(imagewidth, imageheight, GL_BGR_EXT, GL_UNSIGNED_BYTE, pixeldata


    for (int i = 0; i < curveVec.size(); ++i)
    {
        auto &iCurve = curveVec[i];
        if (iCurve.getReadyFlag()) // iCurve has get all control points
        {
            iCurve.display(toEdit);
        }
        if (toEdit || toShowInter)
            iCurve.showInterPoints(toEdit);
        // check if in edit mode
        if (toEdit || toShowCp)
            iCurve.showControlPoints(toEdit);
        if (toShowDer)
            iCurve.showDerivates();
        if (toShowHull)
            iCurve.showHull();
        
    }

    if (toEdit)
    {
        selectFont(48, ANSI_CHARSET, "Comic Sans MS");
        glColor3f(1.0f, 1.0f, 1.0f);

        if (viewPort[2] <= viewPort[3])
        {
            glRasterPos2f(-15.0f, -15.0f*viewPort[3] / viewPort[2]);
        }
        else
        {
            glRasterPos2f(-15.0f*viewPort[2] / viewPort[3], -15.0f);
        }

        drawCNString("In Editing...");

        if (bCtrl)
        {
            auto tmp = strToCoords(coordString);
            if (tmp.first != -1)
            {
                glRasterPos3f((GLfloat)tmp.first, (GLfloat)tmp.second, (GLfloat)0.0);
                XPrintString(coordString.data());
            }

            if (!minDistPt.empty())
            {
                showMinDist();
            }
        }

        // highlight picked point
        if (crv_pt_idxs.first >= 0)
        {
            InterpolationCurve & chosenCrv = curveVec[crv_pt_idxs.first];
            glColor3f(0.7f, 0.6f, 1.f);
            glPointSize(20.0);
            glDisable(GL_POINT_SMOOTH);
            glBegin(GL_POINTS);
            glVertex2d(chosenCrv.getInterPointCoords()[2 * crv_pt_idxs.second], chosenCrv.getInterPointCoords()[2 * crv_pt_idxs.second + 1]);
            glEnd();
        }
        if (!hoverPt.empty())
        {
            glColor3f(1.0f, 1.f, 1.f);
            glPointSize(5.0);
            glDisable(GL_POINT_SMOOTH);
            glBegin(GL_POINTS);
            glVertex2d(hoverPt[0], hoverPt[1]);
            glEnd();
        }

        if (shadowCrv.getReadyFlag()) 
        {
            shadowCrv.display(toEdit);
            shadowCrv.showInterPoints(toEdit);
            //shadowCrv.showControlPoints();
        }
        
    }
    if (1)
        glutSwapBuffers();
    else
        glFlush();
}

void myWindow::reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    glPixelZoom((GLfloat)w / imagewidth, (GLfloat)h / imageheight);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (w <= h)
        //gluOrtho2D(-5.0, 5.0, -5.0 * (GLfloat)h / (GLfloat)w, 5.0 * (GLfloat)h / (GLfloat)w);
        glOrtho(-15.0, 15.0, -15.0 * (GLfloat)h / (GLfloat)w, 15.0 * (GLfloat)h / (GLfloat)w, -10.0, 10.0);
    else
        //gluOrtho2D(-5.0*(GLfloat)w / (GLfloat)h, 5.0*(GLfloat)w / (GLfloat)h, -5.0, 5.0);
        glOrtho(-15.0*(GLfloat)w / (GLfloat)h, 15.0*(GLfloat)w / (GLfloat)h, -15.0, 15.0, -10.0, 10.0);
    //gluPerspective(60.0, (GLfloat)w / (GLfloat)h, 1.0, 20.0);
}

void myWindow::myMotion(int x, int y)
{
    GLint viewPort[4];
    glGetIntegerv(GL_VIEWPORT, viewPort);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMat);
    glGetDoublev(GL_PROJECTION_MATRIX, projMat);
    double qCoord[3];
    gluUnProject(x, viewPort[3] - y, 0, modelMat, projMat, viewPort,
        &qCoord[0], &qCoord[1], &qCoord[2]);

    if (curveVec.empty())
        return;

    if (!toEdit)
    {
        InterpolationCurve &iCurve = curveVec.back();
        if (!iCurve.appendable())
            return;
        std::vector<double> copied(iCurve.getInterPointCoords());
        iCurve.update(iCurve.addInterPointCoords(qCoord, iCurve.getDimension()));
        myDisplay(); // this line should be retain for display instantly
        iCurve.update(iCurve.setInterPointCoords(std::move(copied)));  // should restore for the true insert point.
    }
    else
    {
        coordsToStr(x, y, &coordString); // update coordstring in Motion in Edit
        if (crv_pt_idxs.second >= 0 && crv_pt_idxs.second < curveVec[crv_pt_idxs.first].getInterPointNum() 
            && curveVec[crv_pt_idxs.first].isInFocus()) // only inFocus can move and remove for to show in stripple line.
        {
            InterpolationCurve &eCurve = curveVec[crv_pt_idxs.first];
            //glutSetCursor(GLUT_CURSOR_CROSSHAIR); // ��ʾ��׽���в����Ķ�

            
            if ((eCurve.getInterPointNum() > 1) && (Interpolation::twoPointDist(&chosenPt[0], qCoord, 2) < 20.0 /*/ (eCurve.getInterPointNum() - 1)*/))
            {
                eCurve.update(eCurve.modifyInerPointCoords(crv_pt_idxs.second, qCoord));
            }
            else
            {
                eCurve.update(eCurve.removeInterPointCoords(crv_pt_idxs.second));
                crv_pt_idxs = make_pair(-1, -1);
            }
            myDisplay();
        }
        //glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    }
}

void myWindow::passiveMotion(int x, int y)
{
    coordsToStr(x, y, &coordString);

    auto coords = strToCoords(coordString);
    crv_pt_idxs = findMoveIndex(curveVec, coords.first, coords.second);

    if (toEdit && crv_pt_idxs.second >= 0)
    {
        glutSetCursor(GLUT_CURSOR_FULL_CROSSHAIR/*GLUT_CURSOR_CROSSHAIR*/); // ��ʾ��׽�������
    }
    else
    {
        glutSetCursor(GLUT_CURSOR_RIGHT_ARROW);
    }

    glutPostRedisplay();
}

void myWindow::entryFunc(int state)
{
    if (state == GLUT_ENTERED)
    {
    }
    else
    {
        coordString.clear();
    }
}

void myWindow::myMouse(int button, int state, int x, int y)
{
    GLint viewPort[4];
    glGetIntegerv(GL_VIEWPORT, viewPort);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMat);
    glGetDoublev(GL_PROJECTION_MATRIX, projMat);
    double qCoord[3];
    gluUnProject(x, viewPort[3] - y, 0, modelMat, projMat, viewPort,
        &qCoord[0], &qCoord[1], &qCoord[2]);


    if (curveVec.empty() || !curveVec.back().appendable())
    {
        // create a new curve
        curveVec.push_back(InterpolationCurve());
    }

    InterpolationCurve& iCurve = curveVec.back();
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        if (!toEdit)
        {
            glutSetCursor(GLUT_CURSOR_FULL_CROSSHAIR/*GLUT_CURSOR_CROSSHAIR*/);
        }
        else
        {
            //crv_pt_idxs= findMoveIndex(curveVec, qCoord[0], qCoord[1]);

            if (crv_pt_idxs.first != -1)
            {
                curveVec[crv_pt_idxs.first].setFocus(!curveVec[crv_pt_idxs.first].isInFocus());
                
                shadowCrv = curveVec[crv_pt_idxs.first];
                shadowCrv.setFocus(false);
               

                chosenPt[0] = curveVec[crv_pt_idxs.first].getInterPointCoords()[2 * crv_pt_idxs.second];
                chosenPt[1] = curveVec[crv_pt_idxs.first].getInterPointCoords()[2 * crv_pt_idxs.second + 1];
            }
            else
            {
                
            }


            int mod = glutGetModifiers();
            if (mod == GLUT_ACTIVE_CTRL && crv_pt_idxs.second > 0) // can't insert to be end.
            {
                curveVec[crv_pt_idxs.first].update(curveVec[crv_pt_idxs.first].insertInterPointCoords(qCoord, crv_pt_idxs.second));
            }

            if (bCtrl)
            {
                minDistPt.clear();
                double minDist = MAXLONG;
                for (auto & m : curveVec)
                {
                    std::vector<double> tmp;
                    if (m.FindNearestCurvePoint(qCoord, &tmp) != 0)
                        continue;
                    double tmpDist = Interpolation::twoPointDist(qCoord, &tmp[0], 2);
                    if ( tmpDist < minDist)
                    {
                        minDist = tmpDist;
                        minDistPt.swap(tmp);
                    }
                }                
                chosenPt[0] = qCoord[0];
                chosenPt[1] = qCoord[1];
            }
            
        }
    }



    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    {
        glutSetCursor(GLUT_CURSOR_RIGHT_ARROW);
        if (!toEdit)
        {
            // get click point coordinates
            iCurve.update(iCurve.addInterPointCoords(qCoord, iCurve.getDimension()));

        }
        else
        {
            if (crv_pt_idxs.second != -1 && curveVec[crv_pt_idxs.first].isInFocus())
            {
                curveVec[crv_pt_idxs.first].update(curveVec[crv_pt_idxs.first].modifyInerPointCoords(crv_pt_idxs.second, qCoord));
                crv_pt_idxs = make_pair(-1, -1);

                shadowCrv.clear();

            }

        }

    }

    if (button == GLUT_RIGHT_BUTTON)
    {
        if (!toEdit)
        {
            system("cls");
            iCurve.setAppendFlag(false);
        }
        glutSetCursor(GLUT_CURSOR_RIGHT_ARROW);
    }

    //if (!toEdit && button == GLUT_MIDDLE_BUTTON)
    //{        
    //    iCurve.setAppendFlag(false);
    //    // create a new curve
    //    // avoid the case : middle button click many times.
    //    // assure only one void curve append curveVec
    //    //if (!iCurve.appendable()) 
    //    //curveVec.push_back(InterpolationCurve());
    //}

    myDisplay();

    return;
}

void myWindow::myKey(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 's':
    case 'S':
        // swap buffer, or grab may have a wrong screen shot.
        glutSwapBuffers();
        grab();
        // swap buffer back to show the correct pic
        glutSwapBuffers();
        break;
    case 'b':
    case 'B':
        toBiHe = !toBiHe;
        for (auto iter = curveVec.begin(); iter != curveVec.end(); ++iter)
        {
            if (iter->isInFocus())
            {
                iter->update(iter->SetClose(toBiHe));
            }
        }        
        glutPostRedisplay();
        break;
    case 'c':
    case 'C':
        toShowCp = !toShowCp;
        glutPostRedisplay();
        break;
    case 'd':
    case 'D':
        toShowDer = !toShowDer;
        /*for (auto iter = curveVec.begin(); iter != curveVec.end(); ++iter)
        {
            if (iter->isInFocus())
            {
                std::vector<double> norTmp;
                iter->getDerNorEndPts(iter->getUparam())
                iter->update(iter->SetClose(toBiHe));
            }
        }*/
        glutPostRedisplay();
        break;
    case 'h':
    case 'H':
        toShowHull = !toShowHull;
        glutPostRedisplay();
        break;
    case 'i':
    case 'I':
        toShowInter = !toShowInter;
        glutPostRedisplay();
        break;
    case 'e':
    case 'E':
        toEdit = !toEdit;
        glutPostRedisplay();
    case 'o':
    case 'O':
        toOffset = !toOffset;
        
        /*if (toOffset)
        {
            oOffset *= -1.0f;
        }*/
        
        glutPostRedisplay();
        break;
    case 'r':
    case 'R':
        /*if (!curveVec.empty())
        curveVec.pop_back();*/
        for (auto iter = curveVec.begin(); iter != curveVec.end(); )
        {
            if (iter->isInFocus())
                iter = curveVec.erase(iter);
            else
                ++iter;
        }
        crv_pt_idxs = make_pair(-1, -1); // for the index is not correct
        glutPostRedisplay();
        break;

    case 27: // Esc
        exit(0);
    default:
        break;
    }
}

void myWindow::mySpecialKey(GLint key, GLint x, GLint y)
{
    if (key == GLUT_KEY_UP)
    {
        yOffset += .1f;
    }
    if (key == GLUT_KEY_DOWN)
    {
        yOffset -= .1f;
    }
    if (key == GLUT_KEY_LEFT)
    {
        xOffset -= .1f;
    }
    if (key == GLUT_KEY_RIGHT)
    {
        xOffset += .1f;
    }
    if (key == GLUT_KEY_PAGE_UP)
    {
        zoomRatio *= 1.2f;
    }
    if (key == GLUT_KEY_PAGE_DOWN)
    {
        zoomRatio /= 1.2f;
    }
    if (key == GLUT_KEY_HOME)
    {
        oOffset += (GLfloat)0.1;        
        offsetCurve(oOffset);
    }
    if (key == GLUT_KEY_END)
    {
        oOffset -= (GLfloat)0.1;
        offsetCurve(oOffset);
    }
    myDisplay();
    //glutPostRedisplay();
}

void myWindow::idle()
{
    if (GetAsyncKeyState(VK_CONTROL)/*glutGetModifiers() == GLUT_ACTIVE_CTRL*/)
    {
        bCtrl = true;
    }
    else
    {
        bCtrl = false;
    }
    glutPostRedisplay();
}

void myWindow::selectFont(int size, int charset, const char * face)
{
    HFONT hFont = CreateFontA(size, 0, 0, 0, FW_MEDIUM, 0, 0, 0,
        charset, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
        DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, face);

    HFONT hOldFont = (HFONT)SelectObject(wglGetCurrentDC(), hFont);
    DeleteObject(hOldFont);
}

void myWindow::drawCNString(const char * str)
{
    int len, i;
    wchar_t* wstring;
    HDC hDC = wglGetCurrentDC();
    GLuint list = glGenLists(1);
    len = 0;
    for (i = 0; str[i] != '\0'; ++i)
    {
        if (IsDBCSLeadByte(str[i]))
            ++i;
        ++len;
    }
    wstring = (wchar_t*)malloc((len + 1) * sizeof(wchar_t));
    MultiByteToWideChar(CP_ACP, MB_PRECOMPOSED, str, -1, wstring, len);
    wstring[len] = L'\0';
    for (i = 0; i<len; ++i)
    {
        wglUseFontBitmapsW(hDC, wstring[i], 1, list);
        glCallList(list);
    }
    free(wstring);
    glDeleteLists(list, 1);
}

void myWindow::XPrintString(const char * s)
{
    glPushAttrib(GL_LIST_BIT);
    //����ÿ���ַ���Ӧ����ʾ�б�������ÿ���ַ�
    for (; *s != '\0'; ++s)
        glCallList(TextFont + *s);
    glPopAttrib();
}

void myWindow::coordsToStr(int x, int y, string * str)
{
    assert(str != nullptr);
    GLint viewPort[4];
    GLdouble modelMat[16], projMat[16];
    glGetIntegerv(GL_VIEWPORT, viewPort);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMat);
    glGetDoublev(GL_PROJECTION_MATRIX, projMat);
    double qCoord[3];
    gluUnProject(x, viewPort[3] - y, 0, modelMat, projMat, viewPort,
        &qCoord[0], &qCoord[1], &qCoord[2]);
    char ch[MAX_CHAR];
    sprintf_s(ch, sizeof(ch), "%.5f, %.5f", qCoord[0], qCoord[1]);
    str->assign(ch);
}

pair<double, double> myWindow::strToCoords(const string & str)
{
    if (str.length() != 0)
    {
        char tmpCh[MAX_CHAR];
        strcpy_s(tmpCh, MAX_CHAR, str.data());
        const char *split = ",";
        char *p, *buf;
        double x = 0, y = 0;
        if (p = strtok_s(tmpCh, split, &buf))
        {
            x = atof(p);
            if (p = strtok_s(NULL, split, &buf))
                y = atof(p); //strtod(p)
        }
        return make_pair(x, y);
    }
    return make_pair(-1, -1);
}

void myWindow::showMinDist()
{
    glBegin(GL_LINES);
    glVertex2d(minDistPt[0], minDistPt[1]);
    glVertex2d(chosenPt[0], chosenPt[1]);
    glEnd();
    glBegin(GL_POINTS);
    glVertex2d(minDistPt[0], minDistPt[1]);
    glVertex2d(chosenPt[0], chosenPt[1]);
    glEnd();

    char ch[MAX_CHAR];
    /*glRasterPos2f((GLfloat)minDistPt[0], (GLfloat)minDistPt[1]);
    sprintf_s(ch, sizeof(ch), "(%.2f, %.2f)", minDistPt[0], minDistPt[1]);
    XPrintString(ch);
    glRasterPos2f((GLfloat)chosenPt[0], (GLfloat)chosenPt[1]);
    sprintf_s(ch, sizeof(ch), "(%.2f, %.2f)", chosenPt[0], chosenPt[1]);
    XPrintString(ch);*/
    glRasterPos2f((GLfloat)((minDistPt[0] + chosenPt[0])*.5f), (GLfloat)((minDistPt[1] + chosenPt[1])*.5f));
    sprintf_s(ch, sizeof(ch), "%.5f", Interpolation::twoPointDist(&minDistPt[0], &chosenPt[0], 2));
    XPrintString(ch);
}

void myWindow::readImageFile()
{
    //���ļ�
    FILE* pfile;
    int err = fopen_s(&pfile, "460.bmp", "rb");
    if (err != 0) exit(0);
    //��ȡͼ���С
    fseek(pfile, 0x0012, SEEK_SET);
    fread(&imagewidth, sizeof(imagewidth), 1, pfile);
    fread(&imageheight, sizeof(imageheight), 1, pfile);
    //�����������ݳ���
    pixellength = imagewidth * 3;
    //while (pixellength % 4 != 0)pixellength++;
    pixellength += (4 - pixellength % 4);
    pixellength *= imageheight;
    //��ȡ��������
    pixeldata = (GLubyte*)malloc(pixellength);
    if (pixeldata == 0) exit(0);
    fseek(pfile, 54, SEEK_SET);
    fread(pixeldata, pixellength, 1, pfile);

    //�ر��ļ�
    fclose(pfile);
}

void myWindow::grab(void)
{
    FILE*    pDummyFile = nullptr;
    FILE*    pWritingFile = nullptr;
    GLubyte* pPixelData = nullptr;
    GLubyte  BMP_Header[BMP_Header_Length];

    GLint    PixelDataLength;

    GLint viewPort[4];
    glGetIntegerv(GL_VIEWPORT, viewPort);
    GLsizei ColorChannel = 3;
    PixelDataLength = viewPort[2] * ColorChannel;
    // �����������ݵ�ʵ�ʳ���
    PixelDataLength += (4 - PixelDataLength % 4);
    PixelDataLength *= viewPort[3];

    // �����ڴ�ʹ��ļ�
    pPixelData = (GLubyte*)malloc(PixelDataLength);
    if (pPixelData == nullptr)
        exit(0);
    int err = fopen_s(&pDummyFile, "460.bmp", "rb");
    if (err != 0 || pDummyFile == nullptr)
        exit(0);
    char strBmp[256] = { '\0' };
    static int ii = 0;
    sprintf_s(strBmp, "ScreenShot%d.bmp", ii++);
    err = fopen_s(&pWritingFile, strBmp, "wb");
    if (err != 0 || pWritingFile == nullptr)
        exit(0);

    // ��ȡ����
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glReadPixels(viewPort[0], viewPort[1], viewPort[2], viewPort[3],
        GL_BGR_EXT, GL_UNSIGNED_BYTE, pPixelData);

    // ��dummy.bmp���ļ�ͷ����Ϊ���ļ����ļ�ͷ
    fread(BMP_Header, sizeof(BMP_Header), 1, pDummyFile);
    fwrite(BMP_Header, sizeof(BMP_Header), 1, pWritingFile);
    fseek(pWritingFile, 0x0012, SEEK_SET);
    GLint i = viewPort[2];
    GLint j = viewPort[3];
    fwrite(&i, sizeof(i), 1, pWritingFile);
    fwrite(&j, sizeof(j), 1, pWritingFile);

    // д����������
    fseek(pWritingFile, 0, SEEK_END);
    if (strlen((char*)pPixelData) != 0)
    {
        fwrite(pPixelData, PixelDataLength, 1, pWritingFile);
    }

    // �ͷ��ڴ�͹ر��ļ�
    fclose(pDummyFile);
    fclose(pWritingFile);
    free(pPixelData);
}

pair<int, int> myWindow::findMoveIndex(const vector<InterpolationCurve>& crvVec, double x, double y)
{
    int crvId = -1, ptId = -1;
    double distMin = LONG_MAX;
    std::vector<double> Q = { x, y };
    for (int k = 0; k < crvVec.size(); ++k)
    {
        const InterpolationCurve& crv = crvVec[k];
        if (crv.getInterPointCoords().empty())
            continue;

        const std::vector<double> &iVec = crv.getInterPointCoords();

        for (int i = 0; i < crv.getInterPointNum(); ++i)
        {
            std::vector<double> m{ iVec[i * 2], iVec[i * 2 + 1] };
            double tmp = Interpolation::twoPointDist(&Q[0], &m[0], 2);
            if (tmp < distMin && tmp < 0.4)
            {
                distMin = tmp;
                crvId = k;
                ptId = i;
            }
        }
        //std::vector<double> crvPt;
        //crv.FindNearestCurvePoint(&Q[0], &crvPt);
        //if (!crvPt.empty() && Interpolation::twoPointDist(&Q[0], &crvPt[0], 2) < .3)
        //{
        //    //const_cast<InterpolationCurve&>(crv).setFocus(1);
        //    hoverPt.swap(crvPt);
        //}
        //else
        //{
        //    hoverPt.clear();
        //}
    }

    return make_pair(crvId, ptId);
}

void myWindow::offsetCurve(double lamda)
{
    //const InterpolationCurve* crv = nullptr;
    for (auto& crv : curveVec)
    {
        if (!crv.isInFocus())
            continue;
        
        if (!shadowCrv.getReadyFlag())
        {
            shadowCrv = crv;
            shadowCrv.setFocus(false);
        }
        crv.setOffsetLength(lamda);

    }
}

