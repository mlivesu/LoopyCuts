/***************************************************************************/
/* Copyright(C) 2020

 Marco Livesu
 Italian National Research Council

 and

 Nico Pietroni
 University Of Technology Sydney

 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/


#include <QApplication>
#include <QDesktopWidget>
#include <wrap/qt/anttweakbarMapper.h>
#include "glwidget.h"
#include <QWindow>
#include <QFileInfo>

extern CMesh mesh;
extern std::string pathM;
extern std::string pathF;
extern std::string pathS;
extern bool has_features;
extern bool batch_process;
extern bool delete_unref;
extern bool add_sing_nodes;

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QWindow dummy;
    QString def_string = QString("GLOBAL fontscaling=%1").arg((int)dummy.devicePixelRatio());
    TwDefine(def_string.toStdString().c_str());
    printf("%s\n",qPrintable(def_string));
    fflush(stdout);

    // Set functions to handle string copy
    TwCopyCDStringToClientFunc(CopyCDStringToClient);
    TwCopyStdStringToClientFunc(CopyStdStringToClient);

    if( !TwInit(TW_OPENGL, NULL) )
    {
        fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
        return 1;
    }

    //PARAMETERS CHECK
    if(argc<2)
    {
        printf("error: pass mesh name as parameter \n");
        fflush(stdout);
        exit(0);
    }

    //MESH LOAD
    pathM=std::string(argv[1]);
    QString pathMQ=QString(pathM.c_str());
    QFileInfo f_infoM(pathMQ);
    if (!f_infoM.exists())
    {
        std::cout<<"error: mesh fileneme wrong"<<std::endl;
        fflush(stdout);
        exit(0);
    }
    else
        std::cout<<"Mesh file correct"<<std::endl;


    //FIELD LOAD
    pathF=pathM;
    pathF.erase(pathF.find_last_of("."));
    pathF.append(".rosy");

    QString pathFQ=QString(pathF.c_str());
    QFileInfo f_infoF(pathFQ);
    if (!f_infoF.exists())
    {
        printf("error: field fileneme wrong\n");
        fflush(stdout);
        exit(0);
    }
    else
        std::cout<<"Field file correct"<<std::endl;

    pathS=pathM;
    pathS.erase(pathS.find_last_of("."));
    pathS.append(".sharp");
    QString pathSQ=QString(pathS.c_str());
    QFileInfo f_infoS(pathSQ);
    if (!f_infoS.exists())
    {
        printf("no feature line \n");
        has_features=false;
        fflush(stdout);
        //exit(0);
    }
    else
    {
       has_features=true;
       std::cout<<"Sharp file correct"<<std::endl;
    }


    batch_process=false;
    if (argc>=3)
    {
        //then check if it mush batch process
        std::string pathComm;
        pathComm=std::string(argv[2]);
        if (pathComm==std::string("batch"))batch_process=true;
        if (batch_process)
        {
            std::cout<<"*** BATCH PROCESSING ***"<<std::endl;
            std::cout<<"* DATASET "<<pathM.c_str()<<" *"<<std::endl;
        }
    }
    if (argc>=4)
    {
       delete_unref=(bool)atoi(argv[3]);
       if (delete_unref)
          std::cout<<"must delete unref nodes"<<std::endl;
       else
          std::cout<<"keep unref nodes"<<std::endl;
    }
    if (argc>=5)
    {
       add_sing_nodes=(bool)atoi(argv[4]);
       if (add_sing_nodes)
          std::cout<<"must add sing nodes"<<std::endl;
       else
          std::cout<<"no sing nodes"<<std::endl;
    }
    GLWidget window;
    window.show();
    return app.exec();
}
