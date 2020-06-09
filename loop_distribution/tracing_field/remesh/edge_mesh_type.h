#ifndef EDGE_MESH_TYPE
#define EDGE_MESH_TYPE

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/polygon_support.h>
#include <vcg/simplex/face/component_polygon.h>
#include <vcg/complex/algorithms/update/normal.h>
#ifndef NO_TRACING_OPENGL
#include <wrap/gl/trimesh.h>
#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/gl/gl_field.h>
//#include <para_text.h>
#include <wrap/io_trimesh/export_field.h>
#endif
#include <vcg/complex/algorithms/polygonal_algorithms.h>
//#include <igl/harmonic.h>
//#include <wrap/igl/lscm_parametrization.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <vcg/space/distance3.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/simplex/edge/topology.h>
#include <vcg/complex/algorithms/closest.h>
#include <tracing_field/remesh/edge_splitter.h>
#ifndef NO_TRACING_OPENGL
#include <wrap/gl/addons.h>
#endif
class MyEEdge;
class MyEVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyEVertex>   ::AsVertexType,
        vcg::Use<MyEEdge>     ::AsEdgeType>{};

class MyEVertex  : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3d,
        vcg::vertex::Normal3d, vcg::vertex::VEAdj,
        vcg::vertex::Qualityd,vcg::vertex::BitFlags  >
{
public:
    bool IsSing;


    MyEVertex(){IsSing=false;}

    void ImportData(const MyEVertex  & left )
    {
        vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3d,
                vcg::vertex::Normal3d, vcg::vertex::VEAdj,
                vcg::vertex::Qualityd,vcg::vertex::BitFlags  >::ImportData(left);

        IsSing=left.IsSing;
    }
};

class MyEEdge    : public vcg::Edge<MyUsedTypes,vcg::edge::VertexRef,
        vcg::edge::VEAdj, vcg::edge::EEAdj,
        vcg::edge::BitFlags, vcg::edge::Qualityd,
        vcg::edge::Color4b>
{
public:

    int Face;
    int M4Dir;
    int PathIndex;

    void ImportData(const MyEEdge  & left )
    {
        vcg::Edge<MyUsedTypes,vcg::edge::VertexRef,
                vcg::edge::VEAdj, vcg::edge::EEAdj,
                vcg::edge::BitFlags, vcg::edge::Qualityd,
                vcg::edge::Color4b>::ImportData(left);

        Face=left.Face;
        M4Dir=left.M4Dir;
        PathIndex=left.PathIndex;
    }
};

class MyEMesh    : public vcg::tri::TriMesh< vector<MyEVertex> , vector<MyEEdge>  >
{

    std::vector<CoordType> Elevation;

    template <class Trimesh>
    void EvaluateHLevel(Trimesh &base_mesh)
    {
        Elevation.clear();
        Elevation.resize(vert.size(),CoordType(0,0,0));

        std::vector<CoordType> vertNorm;
        vertNorm.resize(vert.size(),CoordType(0,0,0));


        //sort the vertex quality as minimum path
        for (size_t i=0;i<vert.size();i++)
            vert[i].Q()=std::numeric_limits<ScalarType>::max();

        for (size_t i=0;i<edge.size();i++)
        {
            VertexType *v0=edge[i].V(0);
            VertexType *v1=edge[i].V(1);
            ScalarType pathI=(ScalarType)edge[i].PathIndex;
            vert[i].Q()=std::max(vert[i].Q(),pathI);
            int indexV0=vcg::tri::Index(*this,v0);
            int indexV1=vcg::tri::Index(*this,v1);
            int FIndex=edge[i].Face;
            assert(FIndex>=0);
            assert(FIndex<base_mesh.face.size());
            typename Trimesh::FaceType *f=&base_mesh.face[FIndex];
            vertNorm[indexV0]+=f->N();
            vertNorm[indexV1]+=f->N();
        }

        //then sort vertices by considering quality
        std::vector<std::pair<ScalarType,int> > VertQ;
        for (size_t i=0;i<vert.size();i++)
        {
            VertQ.push_back(std::pair<ScalarType,int>(vert[i].Q(),i));
            vertNorm[i].Normalize();
        }
        std::sort(VertQ.begin(),VertQ.end());


        //save the magnitudo in the quality
        std::map<CoordType,int> CurrL;
        ScalarType distN=base_mesh.bbox.Diag()*0.001;
        for (size_t i=0;i<VertQ.size();++i)
        {
            int currV=VertQ[i].second;
            CoordType pos=vert[currV].P();
            if (CurrL.count(pos)==0){CurrL[pos]=0;vert[currV].Q()=0;continue;}
            CurrL[pos]++;

            if (vert[currV].IsSing)
                vert[currV].Q()=0;
            else
                vert[currV].Q()=((ScalarType)CurrL[pos]*distN);
        }
        //then to some smooth step
        for (int s=0;s<3;s++)
        {
            std::vector<int> Num(vert.size(),0);
            std::vector<ScalarType> SwapQ(vert.size(),0);
            std::vector<CoordType> SwapN(vert.size(),CoordType(0,0,0));

            for (size_t i=0;i<edge.size();i++)
            {
                VertexType *v0=edge[i].V(0);
                VertexType *v1=edge[i].V(1);
                ScalarType Q0=v0->Q();
                ScalarType Q1=v1->Q();

                int indexV0=vcg::tri::Index(*this,v0);
                int indexV1=vcg::tri::Index(*this,v1);
                CoordType N0=vertNorm[indexV0];
                CoordType N1=vertNorm[indexV1];
                Num[indexV0]++;
                Num[indexV1]++;
                SwapQ[indexV0]+=Q1;
                SwapQ[indexV1]+=Q0;
                SwapN[indexV0]+=N1;
                SwapN[indexV1]+=N0;
            }
            for (size_t i=0;i<vert.size();i++)
            {
                vertNorm[i]=vertNorm[i]*0.5+(SwapN[i]/(ScalarType)Num[i])*0.5;
                vertNorm[i].Normalize();
                if (vert[i].IsSing) continue;
                vert[i].Q()=vert[i].Q()*0.5+(SwapQ[i]/(ScalarType)Num[i])*0.5;

            }
        }
        for (size_t i=0;i<vert.size();i++)
            Elevation[i]=vertNorm[i]*vert[i].Q();
    }


public:

    void InitColorByPath()
    {
        //first find the maximum Q
        int MaxP=0;
        for (size_t i=0;i<edge.size();i++)
            if (edge[i].PathIndex>MaxP)MaxP=edge[i].PathIndex;

        //then colorize by path
        for (size_t i=0;i<edge.size();i++)
        {
            int currP=edge[i].PathIndex;
            edge[i].C()=vcg::Color4b::Scatter(MaxP+1,currP);
        }
    }

#ifndef NO_TRACING_OPENGL

    void GLDrawSolid(ScalarType Width)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99995);
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        for (size_t i=0;i<edge.size();++i)
        {
            CoordType pos0=edge[i].P(0);
            CoordType pos1=edge[i].P(1);
            vcg::Point3f pos0f,pos1f;
            pos0f.Import(pos0);
            pos1f.Import(pos1);
            vcg::Color4b col=edge[i].C();
            col[3]=100;
            vcg::glColor(col);
            vcg::Add_Ons::glCylinder<vcg::Add_Ons::DMSolid>(pos0f,pos1f,Width);

            vcg::Add_Ons::glPoint<vcg::Add_Ons::DMSolid>(pos0f,Width);
            vcg::Add_Ons::glPoint<vcg::Add_Ons::DMSolid>(pos1f,Width);
        }
        glPopAttrib();
    }

    void GLDraw()
    {
        const ScalarType Line_Width=6;

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9997);

        //        glDepthRange(0,0.99995);
        //        glEnable (GL_BLEND);
        //        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glLineWidth(Line_Width);
        for (size_t i=0;i<edge.size();++i)
        {
            //CoordType bary=edge[i].P(0)*0.5+edge[i].P(1)*0.5;
            vcg::Color4b col=edge[i].C();
            //col[3]=200;
            vcg::glColor(col);
            glBegin(GL_LINES);
            vcg::glVertex(edge[i].P(0));
            vcg::glVertex(edge[i].P(1));
            //                        vcg::glVertex(bary*0.1+edge[i].P(0)*0.9);
            //                        vcg::glVertex(bary*0.1+edge[i].P(1)*0.9);
            glEnd();
        }

        glDepthRange(0,0.9998);
        glLineWidth(Line_Width+2);
        for (size_t i=0;i<edge.size();++i)
        {
            vcg::glColor(vcg::Color4b(0,0,0,255));
            glBegin(GL_LINES);
            vcg::glVertex(edge[i].P(0));
            vcg::glVertex(edge[i].P(1));
            //                        vcg::glVertex(bary*0.1+edge[i].P(0)*0.9);
            //                        vcg::glVertex(bary*0.1+edge[i].P(1)*0.9);
            glEnd();
        }

        //        glDepthRange(0,0.9995);
        //        glPointSize(Line_Width*2);
        //        vcg::glColor(vcg::Color4b(255,0,0,255));
        //        for (size_t i=0;i<vert.size();++i)
        //        {
        //            if (!vert[i].IsSing)continue;
        //            glBegin(GL_POINTS);
        //            vcg::glVertex(vert[i].P());
        //            glEnd();
        //        }
        //        glDepthRange(0,0.9996);
        //        glPointSize(Line_Width*2+2);
        //        vcg::glColor(vcg::Color4b(0,0,0,255));
        //        for (size_t i=0;i<vert.size();++i)
        //        {
        //            if (!vert[i].IsSing)continue;
        //            glBegin(GL_POINTS);
        //            vcg::glVertex(vert[i].P());
        //            glEnd();
        //        }
        glPopAttrib();
    }

#endif

    void InitFromCoords(const std::vector<std::pair<CoordType,CoordType> > &edge3D,
                        const std::vector<std::pair<size_t,size_t> > &FaceDir,
                        const std::vector<int> &PathIndex)
    {
        Clear();
        //add vertices
        vcg::tri::Allocator<MyEMesh>::AddVertices(*this,edge3D.size()*2);
        //add edges
        vcg::tri::Allocator<MyEMesh>::AddEdges(*this,edge3D.size());
        //initialize the mesh with vertices
        for (size_t i=0;i<edge3D.size();i++)
        {
            vert[i*2].P()=edge3D[i].first;
            vert[(i*2)+1].P()=edge3D[i].second;
            edge[i].V(0)=&vert[i*2];
            edge[i].V(1)=&vert[(i*2)+1];
            edge[i].Face=FaceDir[i].first;
            edge[i].M4Dir=FaceDir[i].second;
            edge[i].PathIndex=PathIndex[i];
        }
        vcg::tri::Clean<MyEMesh>::MergeCloseVertex(*this,0);
        vcg::tri::Clean<MyEMesh>::RemoveUnreferencedVertex(*this);
        //vcg::tri::Clean<MyEMesh>::RemoveDuplicateEdge(*this);
        vcg::tri::Allocator<MyEMesh>::CompactVertexVector(*this);
        vcg::tri::Allocator<MyEMesh>::CompactEdgeVector(*this);
    }

    void SelectNonValence2Vert()
    {
        vcg::tri::UpdateFlags<MyEMesh>::VertexClearS(*this);

        std::vector<size_t> Num(vert.size(),0);
        for(size_t i=0;i<edge.size();i++)
        {
            int Index0=vcg::tri::Index(*this,edge[i].V(0));
            int Index1=vcg::tri::Index(*this,edge[i].V(1));
            Num[Index0]++;
            Num[Index1]++;
        }

        for(size_t i=0;i<Num.size();i++)
        {
            if (Num[i]==2)continue;
            vert[i].SetS();
        }
    }

    void SelectSingularVert()
    {
        vcg::tri::UpdateFlags<MyEMesh>::VertexClearS(*this);

        std::vector<size_t> Num(vert.size(),0);
        std::vector<int> Path(vert.size(),-1);
        for(size_t i=0;i<edge.size();i++)
        {
            int Index0=vcg::tri::Index(*this,edge[i].V(0));
            int Index1=vcg::tri::Index(*this,edge[i].V(1));
            Num[Index0]++;
            Num[Index1]++;
            int IndexPath=edge[i].PathIndex;

            if (Path[Index0]==-1)
                Path[Index0]=IndexPath;
            else
            {
                if (Path[Index0]!=IndexPath)vert[Index0].SetS();
            }
            if (Path[Index1]==-1)
                Path[Index1]=IndexPath;
            else
            {
                if (Path[Index1]!=IndexPath)vert[Index1].SetS();
            }

        }

        for(size_t i=0;i<Num.size();i++)
        {
            if (Num[i]==2)continue;
            vert[i].SetS();
        }

        //finally select singularities
        for(size_t i=0;i<vert.size();i++)
        {
            if (!vert[i].IsSing)continue;
            vert[i].SetS();
        }


    }

    void Smooth(int step=3,
                bool FixSelectedV=true)
    {
        for(int i=0;i<step;++i)
        {
            std::vector<CoordType> Pos(vert.size(),CoordType(0,0,0));
            std::vector<size_t> Num(vert.size(),0);
            for(size_t i=0;i<edge.size();i++)
            {
                int Index0=vcg::tri::Index(*this,edge[i].V(0));
                int Index1=vcg::tri::Index(*this,edge[i].V(1));
                CoordType Pos0=edge[i].V(0)->P();
                CoordType Pos1=edge[i].V(1)->P();
                Pos[Index0]+=Pos1;
                Pos[Index1]+=Pos0;
                Num[Index0]++;
                Num[Index1]++;
            }
            for(size_t i=0;i<vert.size();i++)
            {
                if ((FixSelectedV)&&(vert[i].IsS()))continue;
                vert[i].P()=vert[i].P()*0.5+
                        (Pos[i]/(ScalarType)Num[i])*0.5;
            }
        }
    }

    void SmoothVertNorm(int step=3)
    {
        for(int i=0;i<step;++i)
        {
            std::vector<CoordType> Norm(vert.size(),CoordType(0,0,0));
            std::vector<size_t> Num(vert.size(),0);
            for(size_t i=0;i<edge.size();i++)
            {
                int Index0=vcg::tri::Index(*this,edge[i].V(0));
                int Index1=vcg::tri::Index(*this,edge[i].V(1));
                CoordType Pos0=edge[i].V(0)->N();
                CoordType Pos1=edge[i].V(1)->N();
                Norm[Index0]+=Pos1;
                Norm[Index1]+=Pos0;
                Num[Index0]++;
                Num[Index1]++;
            }
            for(size_t i=0;i<vert.size();i++)
            {
                vert[i].N()=vert[i].N()*0.5+
                        (Norm[i]/(ScalarType)Num[i])*0.5;
                vert[i].N().Normalize();
            }
        }
    }

    template <class TrimeshType>
    void CollectSelfSplitOperation(TrimeshType &baseMesh,
                                   std::vector<std::pair<int,int> > &SplitOP,
                                   std::vector<CoordType> &SplitPos)
    {
        typedef typename TrimeshType::FaceType FaceType;

        std::map<int,std::vector<int> > FaceEdges;

        //initialize per face edges
        for (size_t i=0;i<edge.size();i++)
            FaceEdges[edge[i].Face].push_back(i);

        //clear all the rest
        SplitOP.clear();
        SplitPos.clear();
        vcg::tri::UpdateFlags<MyEMesh>::EdgeClear(*this);

        //for each edge
        for (size_t i=0;i<edge.size();i++)
        {
            int IndexE0=i;
            EdgeType *e0=&edge[IndexE0];

            //check with others if there's a common face
            int IndexF=edge[IndexE0].Face;

            //get the first direction
            int M4Dir0=edge[IndexE0].M4Dir;

            for (size_t j=0;j<FaceEdges[IndexF].size();j++)
            {
                int IndexE1=FaceEdges[IndexF][j];
                EdgeType *e1=&edge[IndexE1];
                if (IndexE0==IndexE1)continue;
                //safety check
                assert(e0->Face==e1->Face);

                //check if that edge has already been selected for split
                if (e1->IsS())continue;
                //check if that edge has already been selected for split
                if (e0->IsS())continue;

                //get the second direction
                int M4Dir1=edge[IndexE1].M4Dir;

                //                //check direction
                if ((M4Dir0%2)==(M4Dir1%2))continue;

                //see if they are adjacent, in such case no split is needed
                if ((e0->V(0)==e1->V(0))||
                        (e0->V(0)==e1->V(1)))
                    continue;

                if ((e0->V(1)==e1->V(0))||
                        (e0->V(1)==e1->V(1)))
                    continue;

                vcg::Segment3<ScalarType> s0(e0->P(0),e0->P(1));
                vcg::Segment3<ScalarType> s1(e1->P(0),e1->P(1));
                FaceType *f=&baseMesh.face[IndexF];

                CoordType IntPos;
                if (!SegSegTangSpaceIntersect<CoordType>(s0,s1,f->N(),IntPos))continue;

                e0->SetS();
                e1->SetS();

                //CoordType AvPos=(closest0+closest1)/2;
                SplitOP.push_back(std::pair<int,int>(IndexE0,IndexE1));
                SplitPos.push_back(IntPos);

                //SplitPos.push_back(AvPos);
            }
        }
    }

    void PerformSplitOperations(std::vector<std::pair<int,int> > &SplitOP,
                                std::vector<CoordType> &SplitPos)
    {
        for (size_t i=0;i<SplitOP.size();i++)
        {
            int IndexE0=SplitOP[i].first;
            int IndexE1=SplitOP[i].second;

            //add the new vertex
            CoordType NewP=SplitPos[i];
            vcg::tri::Allocator<MeshType>::AddVertex(*this,NewP);
            vcg::tri::Allocator<MeshType>::AddEdges(*this,2);

            //get the 4 edges
            EdgeType *e0=&edge[IndexE0];
            EdgeType *e1=&edge[IndexE1];
            EdgeType *e2=&edge[edge.size()-2];
            EdgeType *e3=&edge[edge.size()-1];


            int Face0=e0->Face;
            int Face1=e1->Face;

            //safety check
            assert(Face0==Face1);
            int M4Dir0=e0->M4Dir;
            int M4Dir1=e1->M4Dir;
            assert((M4Dir0%2)!=(M4Dir1%2));

            int PathIndex0=e0->PathIndex;
            int PathIndex1=e1->PathIndex;


            //get the 4 vertices
            VertexType *v0=e0->V(0);
            VertexType *v1=e0->V(1);
            VertexType *v2=e1->V(0);
            VertexType *v3=e1->V(1);
            VertexType *v5=&vert.back();

            e0->V(0)=v0;
            e0->V(1)=v5;
            e0->PathIndex=PathIndex0;
            e0->M4Dir=M4Dir0;
            e0->Face=Face0;

            e1->V(0)=v1;
            e1->V(1)=v5;
            e1->PathIndex=PathIndex0;
            e1->M4Dir=M4Dir0;
            e1->Face=Face0;

            e2->V(0)=v2;
            e2->V(1)=v5;
            e2->PathIndex=PathIndex1;
            e2->M4Dir=M4Dir1;
            e2->Face=Face1;

            e3->V(0)=v3;
            e3->V(1)=v5;
            e3->PathIndex=PathIndex1;
            e3->M4Dir=M4Dir1;
            e3->Face=Face1;
        }
    }

    template <class TrimeshType>
    void SplitSelfIntersections(TrimeshType &baseMesh)
    {

        std::vector<std::pair<int,int> > SplitOP;
        std::vector<CoordType> SplitPos;

        bool splitted=false;

        do
        {

            CollectSelfSplitOperation(baseMesh,SplitOP,SplitPos);

            splitted=(SplitOP.size()>0);

            PerformSplitOperations(SplitOP,SplitPos);

            if (SplitOP.size()>0)
                std::cout << "SPLIT ITERATION " << SplitOP.size() << std::endl;

        }while (splitted);

    }

    void SplitEdge(size_t IndexE,CoordType NewP)
    {
        assert(IndexE>=0);
        assert(IndexE<edge.size());
        EdgeType *curr_e=&edge[IndexE];

        VertexType *v0=curr_e->V(0);
        VertexType *v1=curr_e->V(1);

        int indexv0=vcg::tri::Index(*this,v0);
        int indexv1=vcg::tri::Index(*this,v1);

        //add a new vertex
        vcg::tri::Allocator<MyEMesh>::AddVertices(*this,1);
        //set the poistion as the average
        vert.back().P()=NewP;
        //and update edges
        vcg::tri::Allocator<MyEMesh>::AddEdges(*this,1);
        edge[IndexE].V(0)=&vert[indexv0];
        edge[IndexE].V(1)=&vert.back();
        edge.back().V(0)=&vert.back();
        edge.back().V(1)=&vert[indexv1];

        edge.back().PathIndex=edge[IndexE].PathIndex;
        edge.back().M4Dir=edge[IndexE].M4Dir;
        edge.back().Face=edge[IndexE].Face;
        //edge.back().Q()=edge[IndexE].Q();
    }

    void SplitEdge(size_t IndexE)
    {
        assert(IndexE>=0);
        assert(IndexE<edge.size());
        EdgeType *curr_e=&edge[IndexE];

        VertexType *v0=curr_e->V(0);
        VertexType *v1=curr_e->V(1);

        CoordType NewP=(v0->P()+v1->P())/2;

        SplitEdge(IndexE,NewP);
    }

    struct InterpInfo
    {
        VertexType *v0;
        VertexType *v1;
        vcg::Segment3<ScalarType> ReprojSeg;
        vcg::Matrix33<ScalarType> RotMatrix;

        InterpInfo()
        {
            v0=NULL;v1=NULL;
        }
    };

    template <class BaseTriMesh>
    void SmoothLocally(BaseTriMesh &basemesh,
                       int steps=10,
                       bool checklimits=true)
    {
        typedef typename BaseTriMesh::FaceType FaceType;
        //update topology
        vcg::tri::UpdateTopology<MyEMesh>::VertexEdge(*this);

        //first select all singular vertices
        SelectSingularVert();

        //create the vector of interpolating values
        std::vector<InterpInfo> InterpOp;
        InterpOp.resize(vert.size());

        //then for each vertex non singular
        //get the two neighbors
        for (size_t i=0;i<vert.size();i++)
        {
            //avoid the singularities
            if (vert[i].IsS())continue;

            //get the star of vertices
            std::vector<VertexType *> starVec;
            vcg::edge::VVStarVE<VertexType>(&vert[i],starVec);

            //get the star of edges
            std::vector<EdgeType *> starEdge;
            vcg::edge::VEStarVE<EdgeType>(&vert[i],starEdge);

            //get the shared edge
            assert(starVec.size()<=2);

            if (starVec.size()==1)continue;

            int IndexF0=starEdge[0]->Face;
            int IndexF1=starEdge[1]->Face;

            FaceType *f0=&basemesh.face[IndexF0];
            FaceType *f1=&basemesh.face[IndexF1];
            if (f0==f1)
            {
                vert[i].SetS();
                continue;
            }
            //assert(f0!=f1);

            //get the shared edge
            int i0,i1;
            bool share=vcg::face::ShareEdgeFF(f0,f1,&i0,&i1);
            if(!share)//due to split of internal nodes in faces may happens
            {
                vert[i].SetS();
                continue;
            }
            vcg::Segment3<ScalarType> Seg3(f0->P0(i0),f0->P1(i0));

            InterpOp[i].v0=starVec[0];
            InterpOp[i].v1=starVec[1];
            InterpOp[i].ReprojSeg=Seg3;
            InterpOp[i].RotMatrix=vcg::RotationMatrix(f1->N(),f0->N());
            assert(InterpOp[i].v0!=InterpOp[i].v1);
        }

        //then smooth n steps
        for (int s=0;s<steps;s++)
        {
            std::vector<CoordType> TargetPos(vert.size(),CoordType(0,0,0));

            for (size_t i=0;i<vert.size();i++)
            {
                if (vert[i].IsS())continue;
                assert(InterpOp[i].v0!=NULL);
                assert(InterpOp[i].v1!=NULL);
                VertexType *v0=InterpOp[i].v0;
                VertexType *v1=InterpOp[i].v1;

                vcg::Segment3<ScalarType> SegMesh=InterpOp[i].ReprojSeg;
                //set the Origin
                CoordType Origin=SegMesh.P0();
                SegMesh.P0()=CoordType(0,0,0);
                SegMesh.P1()-=Origin;

                //get the edge mesh
                vcg::Segment3<ScalarType> SegEdge(v0->P(),v1->P());
                SegEdge.P0()-=Origin;
                SegEdge.P1()-=Origin;
                //then rotate
                SegEdge.P1()=InterpOp[i].RotMatrix*SegEdge.P1();

                //then get the closest
                CoordType closest0,closest1;
                ScalarType dist;
                bool parallel;
                vcg::SegmentSegmentDistance(SegEdge,SegMesh,dist,parallel,closest0,closest1);


                //check distance wrt segment extremes
                ScalarType ratio0=(SegMesh.P0()-closest1).Norm()/SegMesh.Length();
                if (checklimits)
                {
                    if (ratio0>=0.95)
                        closest1=SegMesh.P0()*0.05+SegMesh.P1()*0.95;
                    if (ratio0<=0.05)
                        closest1=SegMesh.P0()*0.95+SegMesh.P1()*0.05;
                }
                //sum back the origin
                closest1+=Origin;

                TargetPos[i]=closest1;
            }

            for (size_t i=0;i<vert.size();i++)
            {
                if (vert[i].IsS())continue;
                vert[i].P()=vert[i].P()*0.5+TargetPos[i]*0.5;
            }
        }
        //finally select all singular vertices
        SelectSingularVert();
    }

    template <class BaseTriMesh>
    void InitFaceEdgeMap(BaseTriMesh &basemesh)
    {
        typedef typename BaseTriMesh::FaceType FaceType;
        typedef vcg::GridStaticPtr<FaceType, ScalarType> TriMeshGrid;
        TriMeshGrid grid;

        //set maximum spatial search distance
        ScalarType MaxD=basemesh.bbox.Diag();

        //initialize grid
        grid.Set(basemesh.face.begin(),basemesh.face.end());

        //iterate for each edge
        for (size_t i=0;i<edge.size();i++)
        {
            //find the barycenter
            CoordType AvEdge=(edge[i].P(0)+edge[i].P(1))/2;

            //get closest point
            CoordType closestPt;
            ScalarType minDist;
            FaceType *f=NULL;
            f=vcg::tri::GetClosestFaceBase(basemesh,grid,AvEdge,MaxD,minDist,closestPt);
            if (minDist>0.00001)
            {
                std::cout << "WARNING.. vertex on edge mesh missing along edge dist " << minDist << std::endl;
                edge[i].V(0)->SetS();
                edge[i].V(1)->SetS();
            }
            else
                edge[i].ClearS();

            assert(f!=NULL);
            //get the index
            int IndexF=vcg::tri::Index(basemesh,f);

            //set face to edge relation
            edge[i].Face=IndexF;
        }
    }

    void RefineEdgeIfSmaller(ScalarType min_size)
    {
        std::cout << "*** Number of edges " << edge.size () << std::endl;
        bool refined=false;
        do
        {
            refined=false;
            std::vector<size_t> to_refine;
            for (size_t i=0;i<edge.size();i++)
            {
                if ((edge[i].P(0)-edge[i].P(1)).Norm()<=min_size)continue;
                to_refine.push_back(i);
            }
            if (to_refine.size()>0)refined=true;

            //then split the edges the new edges
            for (size_t i=0;i<to_refine.size();i++)
            {
                int IndexE=to_refine[i];
                SplitEdge(IndexE);
            }
        }while (refined);
        std::cout << "*** Number of edges after refinement " << edge.size () << std::endl;
    }

    template <class Trimesh>
    void ComputeElevation(Trimesh &base_mesh)
    {
        EvaluateHLevel(base_mesh);
    }

    void ApplyElevation()
    {
        assert(Elevation.size()==vert.size());
        for (size_t i=0;i<vert.size();i++)
            vert[i].P()+=Elevation[i];
    }

    template <class MeshType>
    void RemoveEdgeFromMesh(MeshType &mesh)
    {
        //get Sharp features
        //then remove from edge mesh the ones that are incident on original vertices
        std::set<std::pair<CoordType,CoordType> > OrigEdges;
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                CoordType pos0=mesh.face[i].P0(j);
                CoordType pos1=mesh.face[i].P1(j);
                std::pair<CoordType,CoordType> entry(std::min(pos0,pos1),std::max(pos0,pos1));
                OrigEdges.insert(entry);
            }
        for (size_t i=0;i<edge.size();i++)
        {
            typename MeshType::CoordType pos0=edge[i].P(0);
            typename MeshType::CoordType pos1=edge[i].P(1);
            std::pair<typename MeshType::CoordType,typename MeshType::CoordType> entry(std::min(pos0,pos1),std::max(pos0,pos1));
            if (OrigEdges.count(entry)==0)continue;

             vcg::tri::Allocator<MyEMesh>::DeleteEdge(*this,edge[i]);
        }
        vcg::tri::Clean<MyEMesh>::RemoveUnreferencedVertex(*this);
        vcg::tri::Allocator<MyEMesh>::CompactEveryVector(*this);
    }

};



#endif
