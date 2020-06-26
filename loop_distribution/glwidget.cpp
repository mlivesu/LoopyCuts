/****************************************************************************
**
** Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include "glwidget.h"
#include <wrap/qt/trackball.h>
#include <wrap/gl/picking.h>
#include <wrap/qt/anttweakbarMapper.h>
#include <wrap/gl/gl_field.h>
#include <tracing_field/graph/anisotropic_graph.h>
#include <tracing_field/loop_finder.h>
#include <tracing_field/sharp_feature_sampler.h>
#include <loop_stats.h>
#include <loop_mesher.h>
#include <loop_splitter.h>

TwBar *bar;
char * filename;/// filename of the mesh to load
CMesh mesh;     /// the active mesh instance
vcg::GlTrimesh<CMesh> glWrap;    /// the active mesh opengl wrapper
vcg::Trackball track;     /// the active manipulator
GLW::DrawMode drawmode=GLW::DMFlatWire;     /// the current drawmode

std::string pathM,pathF,pathS;

bool to_pick=false;
int xMouse,yMouse;

StatsCollector<typename CMesh::ScalarType> Stats;
AnisotropicGraph<CMesh> aniGraph(mesh,Stats);
SharpFeaturesManager<CMesh> SharpF(mesh);
LoopSplitter<CMesh> LoopSplit(aniGraph,SharpF);

LoopFinder<CMesh> LFind(aniGraph,SharpF);
typedef Path<CMesh::ScalarType> PathType;

LoopFinder<CMesh>::SampleParams SParam;

int GraphRes=3;
bool drawSharpF=true;
bool drawSing=true;
bool drawField=true;
bool drawPaths=false;
bool splitted=false;
bool batch_process=false;
//bool has_field=false;
bool has_features=false;
bool delete_unref=true;
bool add_sing_nodes=false;

MyEMesh EdgeM;

std::vector<PathType> ChoosenPaths;
std::vector<bool> IsLoop;

PathStats<CMesh>::Stats S;


vcg::GridStaticPtr<CFace,double> Gr;
bool select_loops=false;

void InitializeGraph()
{
    aniGraph.InitGraph(GraphRes,SParam.MaxDev,add_sing_nodes,delete_unref);
}


void TW_CALL InitGraph(void *)
{
    //if(!init_color)
    SharpF.CheckConnectivity();
        //mesh.CheckConnectivity();

    InitializeGraph();

    //   CTable.Init(aniGraph);
    //   SharpSampl.TestConvexClose();

    //    SharpSampl.SampleLoopTest();
    //    SharpSampl.InitCandidatesLoopTest(SParam.PenaltyDrift,SParam.AvoidSelfIntersections);
    //    SharpSampl.SampleLoopTest(SParam.PenaltyDrift,SParam.AvoidSelfIntersections);
}

void CopyPathOnQuality()
{
    for (size_t i=0;i<EdgeM.edge.size();i++)
        EdgeM.edge[i].Q()=EdgeM.edge[i].PathIndex;
}


void SaveSetup()
{
    std::string SetupPath=pathM;
    SetupPath.erase(SetupPath.find_last_of("."));
    SetupPath.append("_setup.txt");
    FILE *F=NULL;
    F=fopen(SetupPath.c_str(),"wt");
    assert(F!=NULL);
    fprintf(F,"Graph Resolution %d \n",GraphRes);
    fprintf(F,"Max Angle %5.5f \n",SParam.MaxDev);
    fprintf(F,"Deviance Penalty %5.5f \n",SParam.PenaltyDrift);
    fprintf(F,"OverSample Rate %d \n",SParam.OversampleRate);
    fprintf(F,"Ortho Loop Samples %d \n",SParam.OrthoLoopsSample);

    if (SParam.FindConcaveLoops)
        fprintf(F,"Complete Concave 1 \n");
    else
        fprintf(F,"Complete Concave 0 \n");

    if (SParam.AvoidSelfIntersections)
        fprintf(F,"Avoid Self 1 \n");
    else
        fprintf(F,"Avoid Self 0 \n");

    if (SParam.OnlyOrthogonalClose)
        fprintf(F,"Only Ortho Close 1 \n");
    else
        fprintf(F,"Only Ortho Close 0 \n");

    if (SParam.ResampleAdditionalLoops)
        fprintf(F,"Resample Loops 1 \n");
    else
        fprintf(F,"Resample Loops 0 \n");

    if (SParam.UpdateLoops)
        fprintf(F,"Update Loops 1 \n");
    else
        fprintf(F,"Update Loops 0 \n");

    if (SParam.CloseConvexEndPoints)
        fprintf(F,"Close Convex 1 \n");
    else
        fprintf(F,"Close Convex 0 \n");

    if (SParam.FixInterCross)
        fprintf(F,"Fix Cross 1 \n");
    else
        fprintf(F,"Fix Cross 0 \n");

    fprintf(F,"Additional %d \n",SParam.AdditionalLoops);
    fclose(F);
}

void LoadSetup()
{
    std::string SetupPath=pathM;
    SetupPath.erase(SetupPath.find_last_of("."));
    SetupPath.append("_setup.txt");
    FILE *F=NULL;
    F=fopen(SetupPath.c_str(),"rt");
    if(F==NULL){std::cout<<"Setup File non found"<<std::endl;return;}
    std::cout<<"Setup File Loaded"<<std::endl;
    fscanf(F,"Graph Resolution %d \n",&GraphRes);
    float data;

    fscanf(F,"Max Angle %f \n",&data);
    SParam.MaxDev=(CMesh::ScalarType)data;

    fscanf(F,"Deviance Penalty %f \n",&data);
    SParam.PenaltyDrift=(CMesh::ScalarType)data;

    int DataI;
    fscanf(F,"OverSample Rate %d \n",&DataI);
    SParam.OversampleRate=DataI;

    fscanf(F,"Ortho Loop Samples %d \n",&DataI);
    SParam.OrthoLoopsSample=DataI;

    fscanf(F,"Complete Concave %d \n",&DataI);
    if (DataI==1)
        SParam.FindConcaveLoops=true;
    else
        SParam.FindConcaveLoops=false;

    fscanf(F,"Avoid Self %d \n",&DataI);
    if (DataI==1)
        SParam.AvoidSelfIntersections=true;
    else
        SParam.AvoidSelfIntersections=false;

    fscanf(F,"Only Ortho Close %d \n",&DataI);
    if (DataI==1)
        SParam.OnlyOrthogonalClose=true;
    else
        SParam.OnlyOrthogonalClose=false;

    fscanf(F,"Resample Loops %d \n",&DataI);
    if (DataI==1)
        SParam.ResampleAdditionalLoops=true;
    else
        SParam.ResampleAdditionalLoops=false;

    fscanf(F,"Update Loops %d \n",&DataI);
    if (DataI==1)
        SParam.UpdateLoops=true;
    else
        SParam.UpdateLoops=false;

    fscanf(F,"Close Convex %d \n",&DataI);
    if (DataI==1)
        SParam.CloseConvexEndPoints=true;
    else
        SParam.CloseConvexEndPoints=false;

    fscanf(F,"Fix Cross %d \n",&DataI);
    if (DataI==1)
        SParam.FixInterCross=true;
    else
        SParam.FixInterCross=false;

    fscanf(F,"Additional %d \n",&DataI);
    SParam.AdditionalLoops=DataI;

     fclose(F);

}

void LoadAll()
{
    printf("Loading the mesh \n");
    bool loadedMesh=mesh.LoadFromFile(pathM);
    mesh.UpdateAttributes();
    if (!loadedMesh)
    {
        std::cout<<"*** ERROR LOADING MESH ***"<<std::endl;
        exit(0);
    }
    std::cout<<"Loaded "<<mesh.fn<<" faces and "<<mesh.vn<<" edges"<<std::endl;

     //FIELD LOAD
    bool loadedField=mesh.LoadField(pathF);
    if (!loadedField){
        std::cout<<"*** ERROR LOADING FIELD ***"<<std::endl;
        exit(0);
    }

    if (!has_features)return;

    bool loadedFeatures=SharpF.LoadSharpFeatures(pathS);
    if (!loadedFeatures){
        std::cout<<"*** ERROR LOADING FEATURES ***"<<std::endl;
        exit(0);
    }
}


void TW_CALL SaveSetupButton(void *)
{SaveSetup();}

void Reload()
{
    drawSing=true;
    drawField=true;
    drawPaths=false;
    //drawSmoothPath=true;
    splitted=false;
    //colorByKind=false;

    EdgeM.Clear();

    ChoosenPaths.clear();
    IsLoop.clear();

    //EdgePos.clear();
    //FaceEdgeLoops.clear();

//    CurrLoop=-1;
//    CurrGroup=-1;
//    CurrIndex=-1;

    S.Clear();

//    mesh.LoadFromFile(pathM);
//    mesh.LoadField(pathF);
    LoadAll();

    Gr.Set(mesh.face.begin(),mesh.face.end());
}

void TW_CALL ReloadAll(void *)
{
    Reload();
}


void SaveAllData()
{
    std::string SplittedPath=pathM;
    SplittedPath.erase(SplittedPath.find_last_of("."));
    std::string LoopPath=SplittedPath;
    SplittedPath.append("_splitted.obj");

    LoopPath.append("_loop.txt");
    //vcg::tri::io::ExporterPLY<CMesh>::Save(mesh,SplittedPath.c_str());
    vcg::tri::io::ExporterOBJ<CMesh>::Save(mesh,SplittedPath.c_str(),0);

    LoopSplit.SaveLoopInfo(LoopPath);

    SaveSetup();

}


void TW_CALL SaveData(void *)
{
   SaveAllData();
}

void SplitStep()
{
    std::vector<bool> Essential,CrossOK;
    LFind.GetFinalPaths(ChoosenPaths,IsLoop,Essential,CrossOK);
    LoopSplit.Split(ChoosenPaths,IsLoop,CrossOK);
}

void TW_CALL SplitLoop(void *)
{
  SplitStep();
}


void SampleLoop()
{
    LFind.SampleLoops(SParam);
    drawField=false;
    drawPaths=true;
    drawmode=GLW::DMFlat;
}

void TW_CALL InitLoop(void *)
{
    drawSharpF=false;
    SampleLoop();
    //mesh.InitFaceEdgeSelFromFeatureGeo();
    //SharpF.InitFaceEdgeSelFromFeatureGeo();
    drawPaths=true;
}

void TW_CALL SnapLoop(void *)
{
    LFind.SnapLoop();
}

void BatchLoops()
{
    SharpF.CheckConnectivity();
    InitializeGraph();

    drawSharpF=false;
    SampleLoop();
    drawPaths=true;

    LFind.SnapLoop();

    SplitStep();
}

void TW_CALL BatchProcess(void *)
{
//    InitFeatures();
//    //CleanFeatureStep();
//    //SmoothFieldStep();
//    Reload();
    BatchLoops();
//    //mesh.InitFaceEdgeSelFromFeatureGeo();
//    SharpF.InitFaceEdgeSelFromFeatureGeo();
//    LFind.SnapLoop();
//    SplitStep();
//    colorByKind=true;
//    drawSharpF=false;
//    DrawCandidates=false;
}

//void TW_CALL DeleteSelectdLoop(void *)
//{
//    LFind.DeleteSelectedLoop();
//}

void InitLoopBar(QWidget *w)
{
    (void) w;
    bar = TwNewBar("LoopTool");
    TwDefine("LoopTool size='700 1100' ");
    TwDefine("LoopTool position='40 40' ");

    TwCopyCDStringToClientFunc (CopyCDStringToClient);


    TwAddButton(bar,"Reload All",ReloadAll,0," label='Reload All' ");

     TwEnumVal drawmodes[4] = { {GLW::DMSmooth, "Smooth"}, {GLW::DMPoints, "Per Points"},{GLW::DMFlatWire, "FlatWire"},{GLW::DMFlat, "Flat"}};
    // Create a type for the enum shapeEV
    TwType drawMode = TwDefineEnum("DrawMode", drawmodes, 4);
    TwAddVarRW(bar, "Draw Mode", drawMode, &drawmode, " keyIncr='<' keyDecr='>' help='Change draw mode.' ");

    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"DrawSharpF",TW_TYPE_BOOLCPP,&drawSharpF,"label='Draw Sharp Features'");

    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"DrawSing",TW_TYPE_BOOLCPP,&drawSing,"label='Draw Singularities'");
    TwAddVarRW(bar,"DrawField",TW_TYPE_BOOLCPP,&drawField,"label='Draw Field'");

    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"Resolution",TW_TYPE_INT32,&GraphRes,"label='Graph Resolution'");
    TwAddButton(bar,"Init Graph",InitGraph,0,	" label='Initialize Graph' ");

    TwAddSeparator(bar,NULL,NULL);

    TwAddVarRW(bar,"MaxAngle",TW_TYPE_DOUBLE,&SParam.MaxDev,"label='Max Angle'");
    TwAddVarRW(bar,"DeviancePenalty",TW_TYPE_DOUBLE,&SParam.PenaltyDrift,"label='Deviance Penalty'");
    TwAddSeparator(bar,NULL,NULL);
    TwAddVarRW(bar,"OversampleRate",TW_TYPE_INT32,&SParam.OversampleRate,"label='Oversample Rate'");
    TwAddVarRW(bar,"OrthosampleRate",TW_TYPE_INT32,&SParam.OrthoLoopsSample,"label='Ortho Loop Samples'");
    TwAddVarRW(bar,"ConcaveLoops",TW_TYPE_BOOLCPP,&SParam.FindConcaveLoops,"label='Complete Concave Loops'");
    TwAddVarRW(bar,"SelfIntersection",TW_TYPE_BOOLCPP,&SParam.AvoidSelfIntersections,"label='Avoid Self Intersections'");
    TwAddVarRW(bar,"OnlyOrtho",TW_TYPE_BOOLCPP,&SParam.OnlyOrthogonalClose,"label='Only Ortho Close'");
    TwAddVarRW(bar,"ResLoops",TW_TYPE_BOOLCPP,&SParam.ResampleAdditionalLoops,"label='Resample Loops'");
    TwAddVarRW(bar,"UpdateLoops",TW_TYPE_BOOLCPP,&SParam.UpdateLoops,"label='Update Loops'");
    TwAddVarRW(bar,"CloseConvex",TW_TYPE_BOOLCPP,&SParam.CloseConvexEndPoints,"label='Close Convex'");
    TwAddVarRW(bar,"InterCross",TW_TYPE_BOOLCPP,&SParam.FixInterCross,"label='Final Fix Inter-crossing'");
    TwAddSeparator(bar,NULL,NULL);
    TwAddVarRW(bar,"AdditionalLoop",TW_TYPE_INT32,&SParam.AdditionalLoops,"label='Additional Loops'");
    TwAddButton(bar,"Init Loop",InitLoop,0,	" label='Initialize Loops' ");
    //TwAddButton(bar,"Delete Loop",DeleteSelectdLoop,0,	" label='Delete Sel Loop' ");
    //TwAddVarRW(bar,"DrawPath",TW_TYPE_BOOLCPP,&drawPaths,"label='Draw Paths'");

    TwAddSeparator(bar,NULL,NULL);

    TwAddButton(bar,"Snap Loop",SnapLoop,0,	" label='Snap Loops' ");
    TwAddButton(bar,"Split Loop",SplitLoop,0,	" label='Split Loops' ");

    TwAddSeparator(bar,NULL,NULL);

    TwAddButton(bar,"Batch Process",BatchProcess,0,	" label='BatchProcess' ");
    TwAddButton(bar,"Save Data",SaveData,0,	" label='Save Data' ");
    TwAddButton(bar,"Save Setup",SaveSetupButton,0,	" label='Save Setup' ");
}

void InitBar(QWidget *w)
{
    InitLoopBar(w);
}


GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    filename=0;
    hasToPick=false;
    InitBar(this);
    LoadSetup();
    LoadAll();
    Gr.Set(mesh.face.begin(),mesh.face.end());
    if (batch_process)
    {

        std::cout<<"*** BATCH PROCESSING START ***"<<std::endl;
        int t0=clock();
        std::string pathConsole=pathM;
        pathConsole.erase(pathConsole.find_last_of("."));
        pathConsole.append("_console.txt");

        std::ofstream out(pathConsole.c_str());

        std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
        std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

        BatchLoops();
        int t1=clock();
        std::cout<<"*** TIME :"<<(float)(t1-t0)/(float)CLOCKS_PER_SEC<<" ***"<<std::endl;
        std::cout.rdbuf(coutbuf); //reset to standard output again

        SaveAllData();
        batch_process=false;
        std::cout<<"*** BATCH PROCESSING COMPLETE ***"<<std::endl;
        exit(0);
    }
    //mesh.FlatDegree=60;
    //SharpF.FlatDegree=60;
}

void GLWidget::initializeGL ()
{
//    FieldParam.alpha_curv=0.2;
//    FieldParam.curvRing=4;
//    glewInit();
    glClearColor(0, 0, 0, 0);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
}

void GLWidget::resizeGL (int w, int h)
{
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    TwWindowSize(w, h);
    initializeGL();
}

void GLWidget::paintGL ()
{

    glClearColor(255,255,255,255);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, GLWidget::width()/(float)GLWidget::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();
    glPushMatrix();
    track.Apply();
    //glPushMatrix();
    glWrap.m = &mesh;
    if(mesh.vert.size()>0)
    {
        vcg::glScale(2.0f/mesh.bbox.Diag());
        glTranslate(-mesh.bbox.Center());

        vcg::glColor(vcg::Color4b(220,220,220,255));
        //glWrap.Draw(GLW::DrawMode(drawmode),GLW::CMNone,GLW::TMNone);

        glWrap.Draw(GLW::DrawMode(drawmode),GLW::CMPerFace,GLW::TMNone);

        if (drawField)//&&(has_field))
            vcg::GLField<CMesh>::GLDrawFaceField(mesh,false,false,0.007);

        if (drawSing)//&&(has_field))
            vcg::GLField<CMesh>::GLDrawSingularity(mesh);

        if (drawPaths && (!splitted))
        {
            LFind.DrawLoops();
        }

        LoopSplit.GLDrawLoopInfo();

//        if (drawSmoothPath)
//            GLDrawLoops();


        if (drawSharpF)
            SharpF.GLDrawSharpFeatures();

        aniGraph.GlDrawWrongTris();

        //std::cout<<"8"<<std::endl<<std::flush;
//        if (DrawCandidates)
//            LFind.DrawCandidates();

        //LFind.DrawGenerativeNodes();
        //LFind.DrawSingNodes();
        //LFind.GlDrawDisabledLinks();
        //LFind.DrawIntersections();
        //LFind.DrawFeatures();
        //std::cout<<"9"<<std::endl<<std::flush;

        //SharpSampl.GLDrawSharpDebug();
    }

    //glPopMatrix();
    //track.DrawPostApply();
    if(hasToPick)
    {
        hasToPick=false;
        typename CFace::CoordType pp;
        if(Pick<CFace::CoordType>(pointToPick[0],pointToPick[1],pp))
        {
            typename CFace::CoordType closPt;
            typename CFace::ScalarType minD;
            //std::cout<<pp.X()<<" "<<pp.Y()<<" "<<pp.Z()<<std::endl;
            CFace *f=vcg::tri::GetClosestFaceBase(mesh,Gr,pp,mesh.bbox.Diag(),minD,closPt);
            if ((f!=NULL)&&(drawPaths))
            {
                if (select_loops)
                    LFind.SelectClosestLoop(closPt);
                else
                {
                    if (LFind.selected_loop==-1)
                        LFind.SelectClosestLoopNode(closPt);
                    else
                    {
                        LFind.MoveSelectedToNodeGraph(closPt);
                        LFind.selected_loop=-1;
                    }
                }
            }
            if ((f!=NULL)&&(!drawPaths))
            {
                std::cout<<"flip"<<std::endl;
                int IndexE=0;
                vcg::Segment3d S0(f->P(0),f->P(1));
                vcg::Segment3d S1(f->P(1),f->P(2));
                vcg::Segment3d S2(f->P(2),f->P(0));
                double D0,D1,D2;
                typename CFace::CoordType closS;
                vcg::SegmentPointDistance(S0,pp,closS,D0);
                vcg::SegmentPointDistance(S1,pp,closS,D1);
                vcg::SegmentPointDistance(S2,pp,closS,D2);
                if ((D1<D0)&&(D1<D2))IndexE=1;
                if ((D2<D0)&&(D2<D1))IndexE=2;
                vcg::face::FlipEdge((*f),IndexE);
                mesh.UpdateAttributes();
                vcg::tri::io::ExporterOBJ<CMesh>::Save(mesh,"remeshed_flip.obj",0);

                //f->C()=vcg::Color4b(255,0,0,255);
            }
        }
    }

    glPopMatrix();

    TwDraw();
}

void GLWidget::keyReleaseEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt) track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
    updateGL ();
}


void GLWidget::keyPressEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control) track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));

    TwKeyPressQt(e);
    updateGL ();
}

void GLWidget::mousePressEvent (QMouseEvent * e)
{
    if(!TwMousePressQt(this,e))
    {
        e->accept ();
        setFocus ();
        track.MouseDown(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG (e->button (), e->modifiers ()));
    }
    updateGL ();
}

void GLWidget::mouseMoveEvent (QMouseEvent * e)
{
    if (e->buttons ()) {
        track.MouseMove(QT2VCG_X(this, e), QT2VCG_Y(this, e));
        updateGL ();
    }
    TwMouseMotion(QTLogicalToDevice(this, e->x()), QTLogicalToDevice(this, e->y()));
}

void GLWidget::mouseDoubleClickEvent (QMouseEvent * e)
{
    if (e->buttons ())
    {
        xMouse=QT2VCG_X(this, e);
        yMouse=QT2VCG_Y(this, e);
        //pointToPick=Point2i(e->x(),height()-e->y());
        pointToPick=Point2i(xMouse,yMouse);
        hasToPick=true;
        updateGL ();
    }
    updateGL();
}


void GLWidget::mouseReleaseEvent (QMouseEvent * e)
{
    track.MouseUp(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG(e->button (), e->modifiers ()));
    TwMouseReleaseQt(this,e);
    updateGL ();
}

void GLWidget::wheelEvent (QWheelEvent * e)
{
    const int WHEEL_STEP = 120;
    track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
    updateGL ();
}
