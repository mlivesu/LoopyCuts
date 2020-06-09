#ifndef EDGE_SPLITTER
#define EDGE_SPLITTER

#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/space/intersection2.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <vcg/simplex/face/pos.h>
#include <tracing_field/triangle_mesh_functions.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>


template <class CoordType>
void TrasformMatrices(CoordType DirX,
                      CoordType Norm,
                      vcg::Matrix33<typename CoordType::ScalarType> &RotM,
                      vcg::Matrix33<typename CoordType::ScalarType> &InvRotM)
{
    DirX.Normalize();
    Norm.Normalize();
    typename CoordType::ScalarType dotN=DirX*Norm;
    if (!(fabs(dotN)<0.001))
    {
        DirX-=Norm*dotN;
        DirX.Normalize();
        assert(fabs(DirX*Norm)<0.001);
        //        std::cout<<"dirX "<<DirX.X()<<" "<<DirX.Y()<<" "<<DirX.Z()<<std::endl;
        //        std::cout<<"Norm "<<Norm.X()<<" "<<Norm.Y()<<" "<<Norm.Z()<<std::endl;
        //        std::cout<<"Dot "<<fabs(DirX*Norm)<<std::endl;
        //        assert(0);
    }
    //assert(fabs(DirX*Norm)<0.1);//should be orthogonal

    CoordType DirZ=Norm;
    CoordType DirY=DirX^DirZ;
    DirY.Normalize();

    RotM.SetColumn(0,DirX);
    RotM.SetColumn(1,DirY);
    RotM.SetColumn(2,DirZ);
    RotM.transposeInPlace();

    InvRotM=vcg::Inverse<typename CoordType::ScalarType>(RotM);

}

template <class CoordType>
bool SegSegTangSpaceIntersect(const vcg::Segment3<typename CoordType::ScalarType> &seg0,
                              const vcg::Segment3<typename CoordType::ScalarType> &seg1,
                              CoordType Norm,
                              CoordType &IntPoint)
{
    typedef typename CoordType::ScalarType ScalarType;

    ScalarType delta=std::min(seg0.Length(),seg1.Length())*0.00001;

    //compute the rotation matrix
    CoordType DirX=seg0.P1()-seg0.P0();

    vcg::Matrix33<ScalarType> RotM,InvRotM;
    TrasformMatrices(DirX,Norm,RotM,InvRotM);

    //express the segment in the local reference frame
    CoordType Origin=seg0.P0();

    //the two segments in 2D
    vcg::Segment2<ScalarType> seg0_2D,seg1_2D;
    seg0_2D.P0()=vcg::Point2<ScalarType>(0,0);
    CoordType Swap=seg0.P1()-Origin;
    Swap=RotM*Swap;
    //this because axis X is along this segment
    assert(fabs(Swap.Z())<=0.0001);
    assert(fabs(Swap.Y())<=0.0001);
    seg0_2D.P1()=vcg::Point2<ScalarType>(Swap.X(),Swap.Y());

    Swap=seg1.P0()-Origin;
    Swap=RotM*Swap;
    assert(fabs(Swap.Z())<=0.0001);
    seg1_2D.P0()=vcg::Point2<ScalarType>(Swap.X(),Swap.Y());

    Swap=seg1.P1()-Origin;
    Swap=RotM*Swap;
    assert(fabs(Swap.Z())<=0.0001);
    seg1_2D.P1()=vcg::Point2<ScalarType>(Swap.X(),Swap.Y());

    vcg::Point2<ScalarType> IntPoint2D;
    bool Intersected=vcg::SegmentSegmentIntersection(seg0_2D,seg1_2D,IntPoint2D);
    if (!Intersected)return false;

    //then see if the intersection point is one of the extremes
    if ((IntPoint2D-seg0_2D.P0()).Norm()<delta)return false;
    if ((IntPoint2D-seg0_2D.P1()).Norm()<delta)return false;
    if ((IntPoint2D-seg1_2D.P0()).Norm()<delta)return false;
    if ((IntPoint2D-seg1_2D.P1()).Norm()<delta)return false;

    vcg::Point2<ScalarType> Dir2D_0=seg0_2D.P1()-seg0_2D.P0();
    Dir2D_0.Normalize();

    vcg::Point2<ScalarType> Dir2D_1=seg1_2D.P1()-seg1_2D.P0();
    Dir2D_1.Normalize();

    //see if they are aligned
    if (fabs(Dir2D_0*Dir2D_1)>0.99999) return false;

    //transform back to 3D
    IntPoint.X()=IntPoint2D.X();
    IntPoint.Y()=IntPoint2D.Y();
    IntPoint.Z()=0;

    IntPoint=InvRotM*IntPoint;
    IntPoint+=Origin;
    return true;
}

template <class MESH_TYPE>
class SelEdgePred
{
    typedef MESH_TYPE TriMeshType;
    typedef typename TriMeshType::CoordType    CoordType;
    typedef typename TriMeshType::VertexType   VertexType;
    typedef typename TriMeshType::FaceType     FaceType;
    typedef typename TriMeshType::ScalarType   ScalarType;
public:

    std::map<std::pair<CoordType,CoordType>, CoordType> *SplitMap;

    bool operator()(vcg::face::Pos<typename MESH_TYPE::FaceType> ep)
    {
        CoordType Pos0=ep.f->P0(ep.z);
        CoordType Pos1=ep.f->P1(ep.z);
        std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        return (SplitMap->count(Key)>0);
    }
};


template<class MESH_TYPE>
class SelMidPointFunctor : public std::unary_function<vcg::face::Pos<typename MESH_TYPE::FaceType> , typename MESH_TYPE::CoordType>
{
    typedef MESH_TYPE TriMeshType;
    typedef typename TriMeshType::CoordType    CoordType;
    typedef typename TriMeshType::VertexType   VertexType;
    typedef typename TriMeshType::FaceType     FaceType;
    typedef typename TriMeshType::ScalarType   ScalarType;
public:
    TriMeshType *m;
    std::map<std::pair<CoordType,CoordType>, CoordType> *SplitMap;

    void operator()(typename MESH_TYPE::VertexType &nv,
                    const vcg::face::Pos<typename MESH_TYPE::FaceType> &ep)
    {
        VertexType *v0=ep.V();
        VertexType *v1=ep.VFlip();
        ScalarType Q0=v0->Q();
        ScalarType Q1=v1->Q();

        CoordType Pos0=v0->P();
        CoordType Pos1=v1->P();

        vcg::Point2<ScalarType> t0=v0->T().P();
        vcg::Point2<ScalarType> t1=v1->T().P();
        std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        assert(SplitMap->count(Key)>0);
        nv.P()=(*SplitMap)[Key];

        //        nv.OriginalV=false;
        OriginalV(*m,nv)=false;

        //        IsSing<MESH_TYPE>(*m,nv)=false;
        //nv.IsSing=false;
        ScalarType alpha=(Pos0-nv.P()).Norm()/(Pos0-Pos1).Norm();
        nv.T().P()=t0*(1-alpha)+t1*alpha;
        nv.Q()=Q0*(1-alpha)+Q1*alpha;
        nv.SetV();

    }

    template<class FL_TYPE>
    vcg::TexCoord2<FL_TYPE,1> WedgeInterp(vcg::TexCoord2<FL_TYPE,1> &t0,
                                          vcg::TexCoord2<FL_TYPE,1> &t1)
    {
        vcg::TexCoord2<FL_TYPE,1> tmp;
        assert(t0.n()== t1.n());
        tmp.n()=t0.n();
        tmp.t()=(t0.t()+t1.t())/2;

        return tmp;
    }

    SelMidPointFunctor(TriMeshType *_m){m=_m;}
};

template < class TriMeshType,class EdgeMeshType >
class EdgeSplitter
{
    typedef typename TriMeshType::CoordType    CoordType;
    typedef typename TriMeshType::VertexType   VertexType;
    typedef typename TriMeshType::FaceType     FaceType;
    typedef typename TriMeshType::ScalarType   ScalarType;

    typedef vcg::GridStaticPtr<FaceType, ScalarType> TriMeshGrid;
    TriMeshGrid grid;

    TriMeshType &trimesh;
    EdgeMeshType &emesh;

    bool SplitTriStep()
    {
        //then Split the edges
        SelEdgePred<TriMeshType> SelEP;
        SelMidPointFunctor<TriMeshType> MidEP(&trimesh);
        SelEP.SplitMap=&SplitMap;
        MidEP.SplitMap=&SplitMap;

        bool refined=vcg::tri::RefineE<TriMeshType,SelMidPointFunctor<TriMeshType>,SelEdgePred<TriMeshType> >(trimesh,MidEP,SelEP);

        return refined;
    }

    bool SplitEdgeIfNeeded(int IndexE,
                           const int FIndex,
                           int &IndexE0,
                           int &IndexE1)
    {
        IndexE0=-1;
        IndexE1=-1;

        //initialize first segment
        vcg::Segment3<ScalarType> seg0(emesh.edge[IndexE].P(0),
                                       emesh.edge[IndexE].P(1));

        for (size_t i=0;i<3;i++)
        {
            //get face's segment on edge
            CoordType p0=trimesh.face[FIndex].P0(i);
            CoordType p1=trimesh.face[FIndex].P1(i);
            vcg::Segment3<ScalarType> seg1(p0,p1);

            CoordType Norm=trimesh.face[FIndex].N();

            //check if there's intersection
            CoordType IntPoint;
            if (!SegSegTangSpaceIntersect<CoordType>(seg0,seg1,Norm,IntPoint))continue;

            //in that case split the edge
            emesh.SplitEdge(IndexE,IntPoint);
            IndexE0=IndexE;
            IndexE1=emesh.edge.size()-1;
            return true;
        }
        return false;
    }

    void SplitEdgeIfNeeded(int IndexE,
                           std::vector<int> Faces)
    {
        //set stack of edges that need to be checked
        std::vector<int> To_check_edges;
        To_check_edges.push_back(IndexE);

        do
        {
            //get the one that need to be tested
            int CurrEdge=To_check_edges.back();
            //erase from the stack
            To_check_edges.pop_back();
            for (size_t i=0;i<Faces.size();i++)
            {
                int IndexE0,IndexE1;
                bool splitted=SplitEdgeIfNeeded(CurrEdge,Faces[i],IndexE0,IndexE1);
                if (!splitted)continue;

                //safety check
                assert((IndexE0>=0)&&(IndexE0<(int)emesh.edge.size()));
                assert((IndexE1>=0)&&(IndexE1<(int)emesh.edge.size()));

                //otherwise add new possible split operations
                To_check_edges.push_back(IndexE0);
                To_check_edges.push_back(IndexE1);
            }
        }while (!To_check_edges.empty());
    }

    void SplitEdgeStep()
    {
        //bool splitted=false;
        size_t edge_size=emesh.edge.size();
        for (size_t i=0;i<edge_size;i++)
        {
            //get the index of the old faces
            int IndexF=EdgeFace[i];

            //then find all new faces with such index
            std::vector<int> NewFIndex;
            NewFIndex=std::vector<int>(FaceOriginalMap[IndexF].begin(),
                                       FaceOriginalMap[IndexF].end());

            SplitEdgeIfNeeded(i,NewFIndex);
        }
    }

    std::map<std::pair<CoordType,CoordType>, CoordType> SplitMap;

    void AddPossibleSplittedEges(const FaceType &CurrF,
                                 const size_t face_edge,
                                 const std::vector<size_t> &edgesI)
    {
        assert((face_edge>=0)&&(face_edge<3));
        //create the segment
        CoordType PosF0=std::min(CurrF.cP0(face_edge),CurrF.cP1(face_edge));
        CoordType PosF1=std::max(CurrF.cP0(face_edge),CurrF.cP1(face_edge));

        //create first segment
        vcg::Segment3<ScalarType> seg0(PosF0,PosF1);

        for (size_t i=0;i<edgesI.size();i++)
        {
            size_t IndexE=edgesI[i];
            CoordType p0=emesh.edge[IndexE].P(0);
            CoordType p1=emesh.edge[IndexE].P(1);

            //get second segment
            //enlarge it for numerical reasons
            CoordType bary=(p0+p1)/2;
            CoordType p0_ext=bary+(p0-bary)*1.25;
            CoordType p1_ext=bary+(p1-bary)*1.25;
            vcg::Segment3<ScalarType> seg1(p0_ext,p1_ext);

            CoordType IntPoint;
            if (!SegSegTangSpaceIntersect<CoordType>(seg0,seg1,CurrF.cN(),IntPoint))continue;

            std::pair<CoordType,CoordType> Key(std::min(PosF0,PosF1),
                                               std::max(PosF0,PosF1));
            SplitMap[Key]=IntPoint;
        }
    }



    void AddPossibleSplittedEges(FaceType &CurrF)
    {
        int IndexF=vcg::tri::Index(trimesh,CurrF);

        //see if it has trasversed by one path
        if (FaceEdges.count(IndexF)==0)return;

        //otherwise check real intersection
        //with all possible segments
        for (size_t j=0;j<3;j++)
            AddPossibleSplittedEges(CurrF,j,FaceEdges[IndexF]);
    }

    void AddPossibleSplittedEges()
    {
        //clear the map
        SplitMap.clear();

        //then add the intersections
        for (size_t i=0;i<trimesh.face.size();i++)
            AddPossibleSplittedEges(trimesh.face[i]);

        std::cout << "Added splitting operations " << SplitMap.size() << std::endl;
    }


    //for each face collect all the edges that transverse it
    std::map<size_t,std::vector<size_t> > FaceEdges;
    //for each face collect the face it pass over
    std::vector<size_t> EdgeFace;

    void InitFaceEdgeMap()
    {
        //allocate space
        EdgeFace.resize(emesh.edge.size(),-1);
        //clear the map
        FaceEdges.clear();
        //set maximum spatial search distance
        ScalarType MaxD=trimesh.bbox.Diag();

        //initialize grid
        grid.Set(trimesh.face.begin(),trimesh.face.end());

        //iterate for each edge
        for (size_t i=0;i<emesh.edge.size();i++)
        {
            //find the barycenter
            CoordType AvEdge=(emesh.edge[i].P(0)+emesh.edge[i].P(1))/2;

            //get closest point
            CoordType closestPt;
            ScalarType minDist;
            FaceType *f=NULL;
            f=vcg::tri::GetClosestFaceBase(trimesh,grid,AvEdge,MaxD,minDist,closestPt);
            if (minDist>0.00001)
            {
                std::cout << "WARNING.. vertex on edge mesh missing along edge dist " << minDist << std::endl;
                emesh.edge[i].V(0)->SetS();
                emesh.edge[i].V(1)->SetS();
            }
            else
                emesh.edge[i].ClearS();

            assert(f!=NULL);
            //get the index
            int IndexF=vcg::tri::Index(trimesh,f);

            //set face to edge relation
            FaceEdges[IndexF].push_back(i);

            //set edge to face
            EdgeFace[i]=IndexF;
        }
        std::cout << "Setted crossing faces " << FaceEdges.size() << std::endl;
    }

    //set per face quality as index of the face
    void InitFaceIndexOnQ()
    {
        for (size_t i=0;i<trimesh.face.size();i++)
            trimesh.face[i].Q()=i;
    }

    std::map<int,std::vector<int> > FaceOriginalMap;

    //set for each index stored into quality set all
    //the faces that have such index
    void InitFaceIndexMap()
    {
        FaceOriginalMap.clear();
        for (size_t i=0;i<trimesh.face.size();i++)
        {
            int Index=(int)trimesh.face[i].Q();
            FaceOriginalMap[Index].push_back(i);
        }
    }

    void SetVertClosestFace(std::vector<FaceType*> &closestF)
    {
        ScalarType MaxD=trimesh.bbox.Diag();
        closestF.clear();
        for (size_t i=0;i<emesh.vert.size();i++)
        {
            CoordType p=emesh.vert[i].P();
            //get closest point
            CoordType closestPt;
            ScalarType minDist;
            FaceType *f=NULL;
            f=vcg::tri::GetClosestFaceBase(trimesh,grid,p,MaxD,minDist,closestPt);
            assert(f!=NULL);
            closestF.push_back(f);
        }
    }

    ScalarType AVEdgeSize()
    {
        ScalarType sum=0;
        int num=0;
        for (size_t i=0;i<trimesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                sum+=(trimesh.face[i].P0(j)-trimesh.face[i].P1(j)).Norm();
                num++;
            }
        return (sum/(ScalarType)num);
    }

    bool SplitInternalStep()
    {
        bool splitted=false;
        vcg::tri::UpdateFlags<EdgeMeshType>::VertexClearV(emesh);

        ScalarType bary_delta=0.0001;
        ScalarType delta=AVEdgeSize()*0.0001;

        //initialize grid
        grid.Set(trimesh.face.begin(),trimesh.face.end());

        std::vector<FaceType*> closestF;
        std::vector<int> closestFIndex;
        SetVertClosestFace(closestF);

        //transform to index
        for (size_t i=0;i<closestF.size();i++)
            closestFIndex.push_back(vcg::tri::Index(trimesh,closestF[i]));

        //then each non regular vertex see if fall inside or not
        //emesh.SelectNonValence2Vert();
        emesh.SelectSingularVert();

        for (size_t i=0;i<emesh.vert.size();i++)
        {
            if (!emesh.vert[i].IsS())continue;
            if (emesh.vert[i].IsV())continue;

            CoordType currP=emesh.vert[i].P();
            FaceType *closF=&trimesh.face[closestFIndex[i]];
            if (closF->IsV())continue;


            //find distances to extremes
            bool found_extr=false;
            for (size_t j=0;j<3;j++)
            {

                //if the vertex is on edge then
                if ((closF->P(j)-currP).Norm()<delta)
                {
                    emesh.vert[i].P()=closF->P(j);
                    emesh.vert[i].SetV();//mark as solved
                    found_extr=true;
                    continue;
                }
            }
            if (found_extr)continue;

            //check on edge
            CoordType bary;
            vcg::InterpolationParameters(*closF,currP,bary);

            for (int j=0;j<3;j++)
            {
                if (fabs(bary.V(j))<=bary_delta)
                {
                    emesh.vert[i].SetV();//mark as solved
                    bary.V(j)=0;
                    bary.V((j+1)%3)=1-bary.V((j+2)%3);
                    emesh.vert[i].P()=closF->P(0)*bary.X()+
                            closF->P(1)*bary.Y()+
                            closF->P(2)*bary.Z();
                    //put on the edge
                    found_extr=true;
                    continue;
                }
            }
            if (found_extr)continue;
            //add two new faces
            vcg::tri::Allocator<TriMeshType>::AddFaces(trimesh,2);
            vcg::tri::Allocator<TriMeshType>::AddVertices(trimesh,1);

            //reupdate the pointer in case it has been changed
            closF=&trimesh.face[closestFIndex[i]];

            //otherwise split the face
            int sizeF=trimesh.face.size();
            FaceType *f0=&trimesh.face[sizeF-1];
            FaceType *f1=&trimesh.face[sizeF-2];
            f0->ImportData(*closF);
            f1->ImportData(*closF);
            VertexType *nv=&trimesh.vert.back();
            vcg::tri::CenterPointBarycenter<TriMeshType> CPoint;
            vcg::tri::TriSplit<TriMeshType>::Apply(closF,f0,f1,nv,CPoint);
            nv->P()=currP;
            closF->SetV();
            emesh.vert[i].SetV();//mark as solved
            splitted=true;
        }
        //        vcg::tri::UpdateNormal<TriMeshType>::PerFaceNormalized(trimesh);
        //        vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(trimesh);
        return splitted;
    }


    std::vector<std::pair<FaceType*,int> > ClosestEdge;

    ScalarType EvaluateDist(const FaceType &f,int edgeIndex,
                            const typename EdgeMeshType::EdgeType &e)
    {
        ScalarType delta_norm=0.99;

        CoordType Pe0=e.cP(0);
        CoordType Pe1=e.cP(1);
        CoordType Pf0=f.cP0(edgeIndex);
        CoordType Pf1=f.cP1(edgeIndex);
        ScalarType dist0=(Pe0-Pf0).Norm()+(Pe1-Pf1).Norm();
        ScalarType dist1=(Pe0-Pf1).Norm()+(Pe1-Pf0).Norm();

        CoordType dir0=Pe1-Pe0;
        CoordType dir1=Pf1-Pf0;
        dir0.Normalize();
        dir1.Normalize();
        if ((fabs(dir0*dir1))<delta_norm)return std::numeric_limits<ScalarType>::max();
        return (std::min(dist0,dist1));
    }

    void InitClosestEdge()
    {
        //ScalarType AvE=trimesh.AVEdgeSize()*0.5;

        ScalarType AvE=AVEdgeSize();
        std::set<std::pair<FaceType*,int> > InsertedEdges;

        //initialize the grid
        grid.Set(trimesh.face.begin(),trimesh.face.end());

        //allocate space
        ClosestEdge.resize(emesh.edge.size(),std::pair<FaceType*,int>(NULL,-1));

        //then find closest edge for each face
        for (size_t i=0;i<emesh.edge.size();i++)
        {
            //get the bounding box and inflate it
            vcg::Box3<ScalarType> bbox;
            bbox.Add(emesh.edge[i].P(0));
            bbox.Add(emesh.edge[i].P(1));
            bbox.Offset(AvE);

            std::vector<FaceType*> closeF;
            //get the faces into the bounding box
            vcg::tri::GetInBoxFace(trimesh,grid,bbox,closeF);

            //safety check
            assert(closeF.size()>0);

            //then find the closest edge
            ScalarType minD=std::numeric_limits<ScalarType>::max();
            for (size_t j=0;j<closeF.size();j++)
                for (size_t k=0;k<3;k++)
                {
                    ScalarType testD=EvaluateDist(*closeF[j],k,emesh.edge[i]);
                    if (testD<minD)
                    {
                        minD=testD;
                        ClosestEdge[i]=std::pair<FaceType*,int>(closeF[j],k);
                    }
                }

            //finally check that is not null
            assert(ClosestEdge[i].first!=NULL);

            //and then check it has not been already inserted
            if(InsertedEdges.count(ClosestEdge[i])>0)
            {
                std::cout<<"WARNING an Edge has been selected twice"<<std::endl;
            }

            //then add to the set of already inserted ones
            InsertedEdges.insert(ClosestEdge[i]);
        }
    }

    void MarkPathEdges()
    {
        //#ifndef TRACING_EXTERNAL
        for (size_t i=0;i<trimesh.face.size();i++)
            for (size_t j=0;j<3;j++)
                //trimesh.face[i].PathIndex[j]=-1;
                PathIndex(trimesh,trimesh.face[i],j)=-1;
        //#endif

        //vcg::tri::UpdateFlags<TriMeshType>::FaceClearCreases(trimesh);
        vcg::tri::UpdateFlags<TriMeshType>::FaceClearFaceEdgeS(trimesh);
        vcg::tri::UpdateTopology<TriMeshType>::FaceFace(trimesh);

        for (size_t i=0;i<ClosestEdge.size();i++)
        {
            FaceType *f=ClosestEdge[i].first;
            int IndexE=ClosestEdge[i].second;
            //f->SetCrease(IndexE);
            f->SetFaceEdgeS(IndexE);
            FaceType *f1=f->FFp(IndexE);
            int IndexE1=f->FFi(IndexE);
            if (f1==f)continue;
            //f1->SetCrease(IndexE1);
            f1->SetFaceEdgeS(IndexE1);

            //save the original path
            int IndexF0=vcg::tri::Index(trimesh,f);
            int IndexF1=vcg::tri::Index(trimesh,f1);

            //#ifndef TRACING_EXTERNAL
            //trimesh.face[IndexF0].PathIndex[IndexE]=emesh.edge[i].PathIndex;
            //trimesh.face[IndexF1].PathIndex[IndexE1]=emesh.edge[i].PathIndex;

            PathIndex<TriMeshType>(trimesh,trimesh.face[IndexF0],IndexE)=emesh.edge[i].PathIndex;
            PathIndex<TriMeshType>(trimesh,trimesh.face[IndexF1],IndexE1)=emesh.edge[i].PathIndex;
            //#endif
        }
    }

public:

    void SplitFace(FaceType *to_split,
                   CoordType newPos,
                   VertexType *&newV)
    {
        std::map<std::pair<CoordType,CoordType>,int> PathMap;
        //to_split->GetPathMap(PathMap);
        GetPathMap<TriMeshType>(trimesh,*to_split,PathMap);

        //add two new faces and vertices
        vcg::tri::Allocator<TriMeshType>::AddFaces(trimesh,2);
        vcg::tri::Allocator<TriMeshType>::AddVertices(trimesh,1);

        //then split the face
        int sizeF=trimesh.face.size();
        FaceType *f0=to_split;
        FaceType *f1=&trimesh.face[sizeF-1];
        FaceType *f2=&trimesh.face[sizeF-2];

        //import the values
        f1->ImportData(*f0);
        f2->ImportData(*f0);

        //get the new vertex
        newV=&trimesh.vert.back();

        vcg::tri::CenterPointBarycenter<TriMeshType> CPoint;
        vcg::tri::TriSplit<TriMeshType>::Apply(f0,f1,f2,newV,CPoint);
        newV->P()=newPos;

        //set the path Index
        //        f0->SetFromPathMap(PathMap);
        //        f1->SetFromPathMap(PathMap);
        //        f2->SetFromPathMap(PathMap);
        SetFromPathMap<TriMeshType>(trimesh,*f0,PathMap);
        SetFromPathMap<TriMeshType>(trimesh,*f1,PathMap);
        SetFromPathMap<TriMeshType>(trimesh,*f2,PathMap);

        vcg::tri::UpdateTopology<TriMeshType>::FaceFace(trimesh);
        vcg::tri::UpdateTopology<TriMeshType>::VertexFace(trimesh);

        vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(trimesh);
    }

    void SplitEdges(std::vector<std::pair<CoordType,CoordType> > &Edges,
                    std::vector<CoordType> &MiddlePos)
    {
        std::map<std::pair<CoordType,CoordType>,int> PathMap;

        //get path map
        GetPathMap(trimesh,PathMap);

        //add the edges to be splitted
        for (size_t i=0;i<Edges.size();i++)
        {
            //add the splitting action
            SplitMap[ Edges[i] ]=MiddlePos[i];

            //erase old edge
            assert(PathMap.count(Edges[i])>0);
            int IndexPath=PathMap[ Edges[i] ];
            PathMap.erase(Edges[i]);

            //add the now ones
            CoordType pos0=Edges[i].first;
            CoordType pos1=Edges[i].second;
            std::pair<CoordType,CoordType> Key0(std::min(pos0,MiddlePos[i]),
                                                std::max(pos0,MiddlePos[i]));

            std::pair<CoordType,CoordType> Key1(std::min(pos1,MiddlePos[i]),
                                                std::max(pos1,MiddlePos[i]));

            PathMap[Key0]=IndexPath;
            PathMap[Key1]=IndexPath;
        }

        //then split the edge mesh
        SelEdgePred<TriMeshType> SelEP;
        SelMidPointFunctor<TriMeshType> MidEP(&trimesh);
        SelEP.SplitMap=&SplitMap;
        MidEP.SplitMap=&SplitMap;

        vcg::tri::RefineE<TriMeshType,SelMidPointFunctor<TriMeshType>,SelEdgePred<TriMeshType> >(trimesh,MidEP,SelEP);

        //finally reupdate the paths
        SetFromPathMap(trimesh,PathMap);

        vcg::tri::UpdateTopology<TriMeshType>::FaceFace(trimesh);
        vcg::tri::UpdateTopology<TriMeshType>::VertexFace(trimesh);

        vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(trimesh);
    }

    void SplitInternal()
    {
        bool splitted=false;
        do
        {
            vcg::tri::UpdateFlags<TriMeshType>::FaceClearV(trimesh);
            //initialize original face index
            InitFaceIndexOnQ();
            //initialize face to edge map
            InitFaceEdgeMap();
            //split interior faces
            splitted=SplitInternalStep();
            //initialize face to edge map
            InitFaceIndexMap();
            //then split the edge mesh
            SplitEdgeStep();
        }while(splitted);
        emesh.InitColorByPath();

        vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(trimesh);
    }

    //    void SplitAlignedEdges()
    //    {
    //       //first select only
    //    }

    void Split(bool MarkPathE=true)
    {
        //split internal ones for safety
        //SplitInternal();

        //save the singularities
        bool refined=true;
        do{
            std::cout << "PERFORMIN SPLIT STEP " << std::endl;
            //save into quality the index of the face
            InitFaceIndexOnQ();
            //initialize face to edge map
            InitFaceEdgeMap();
            //initialize all possible face split operations
            AddPossibleSplittedEges();
            //split the tri mesh
            refined=SplitTriStep();
            //then set all the faces that have such index
            InitFaceIndexMap();
            //then split the edge mesh considering all the new generated faces
            SplitEdgeStep();
            //refined|=SplitEdgeStep();
        }while (refined);

        std::cout<<"SPLITTING OPS TERMINATED"<<std::endl;
        ////        //finally collapse and erase duplicated vertices
        //        ScalarType delta=trimesh.AVEdgeSize()*0.001;
        //        vcg::tri::Clean<TriMeshType>::MergeCloseVertex(trimesh,delta);
        //        vcg::tri::Clean<TriMeshType>::RemoveDegenerateFace(trimesh);
        //        vcg::tri::Allocator<TriMeshType>::CompactEveryVector(trimesh);

        vcg::tri::UpdateTopology<TriMeshType>::FaceFace(trimesh);
        vcg::tri::UpdateTopology<TriMeshType>::VertexFace(trimesh);

        vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(trimesh);

        vcg::tri::CrossField<TriMeshType>::UpdateSingularByCross(trimesh);


        InitClosestEdge();

        if (MarkPathE)
            MarkPathEdges();

    }

    void GetPathEdges(std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgePaths)
    {
        int NumPaths=-1;
        FaceEdgePaths.clear();
        for (size_t i=0;i<trimesh.face.size();i++)
            for (size_t j=0;j<3;j++)
                NumPaths=std::max(NumPaths,PathIndex(trimesh,trimesh.face[i],j));

        if (NumPaths<0)return;
        FaceEdgePaths.resize(NumPaths+1);
        for (size_t i=0;i<trimesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                int IndexPath=PathIndex(trimesh,trimesh.face[i],j);
                if (IndexPath<0)continue;
                FaceEdgePaths[IndexPath].push_back(std::pair<size_t,size_t>(i,j));
            }
    }

    EdgeSplitter( TriMeshType &_trimesh,
                  EdgeMeshType &_emesh):trimesh(_trimesh),emesh(_emesh)
    {}
};

#endif //SEPARATRIX_OPTIMIZER
