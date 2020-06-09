#ifndef MY_TRI_MESH_FUNCTIONS
#define MY_TRI_MESH_FUNCTIONS
#include <vcg/complex/complex.h>
#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/voronoi_processing.h>

enum VertType{VRegular,VPath,VSingular,VCross,VTJunction};


template <class MeshType>
void ValidateAdditionalParameters(MeshType m)
{
    //typename MeshType::template PerVertexAttributeHandle<bool> IsSing = vcg::tri::Allocator<MeshType>:: template GetPerVertexAttribute<bool> (m);
    typename MeshType::template PerVertexAttributeHandle<bool> IsSing =  vcg::tri::Allocator<MeshType>:: template FindPerVertexAttribute<bool>(m,std::string("IsSing"));
    if(!MeshType::template PerVertexAttributeHandle<bool>::IsValidHandle(m,IsSing))assert(0);
}

struct EdgePath
{
  int PerEdgePath[3];

  EdgePath()
  {
      PerEdgePath[0]=-1;
      PerEdgePath[1]=-1;
      PerEdgePath[2]=-1;
  }
};

template <class MeshType>
void AllocateAdditionalParameters(MeshType &m)
{
    if (!vcg::tri::HasPerVertexAttribute(m,std::string("Singular")))
        vcg::tri::Allocator<MeshType>:: template AddPerVertexAttribute<bool> (m,std::string("Singular"));
    if (!vcg::tri::HasPerVertexAttribute(m,std::string("VertexKind")))
        vcg::tri::Allocator<MeshType>:: template AddPerVertexAttribute<VertType> (m,std::string("VertexKind"));
    if (!vcg::tri::HasPerVertexAttribute(m,std::string("IsConstrained")))
        vcg::tri::Allocator<MeshType>:: template AddPerVertexAttribute<bool> (m,std::string("IsConstrained"));
    if (!vcg::tri::HasPerVertexAttribute(m,std::string("OriginalV")))
        vcg::tri::Allocator<MeshType>:: template AddPerVertexAttribute<bool> (m,std::string("OriginalV"));
    if (!vcg::tri::HasPerFaceAttribute(m,std::string("Region")))
        vcg::tri::Allocator<MeshType>:: template AddPerFaceAttribute<int> (m,std::string("Region"));
    if (!vcg::tri::HasPerFaceAttribute(m,std::string("PathIndex")))
        vcg::tri::Allocator<MeshType>:: template AddPerFaceAttribute<EdgePath> (m,std::string("PathIndex"));
    //typename MeshType::template PerVertexAttributeHandle<bool> IsSing = vcg::tri::Allocator<MeshType>:: template GetPerVertexAttribute<bool> (m,std::string("IsSing"));
}


template <class MeshType>
bool &OriginalV(MeshType &m, typename MeshType::VertexType &v)
{
    assert(vcg::tri::HasPerVertexAttribute(m,std::string("OriginalV")));
    typename MeshType::template PerVertexAttributeHandle<bool> OrigV =  vcg::tri::Allocator<MeshType>:: template FindPerVertexAttribute<bool>(m,std::string("OriginalV"));
    return OrigV[v];
}

template <class MeshType>
bool &IsSing(MeshType &m, typename MeshType::VertexType &v)
{
    assert(vcg::tri::HasPerVertexAttribute(m,std::string("Singular")));
    typename MeshType::template PerVertexAttributeHandle<bool> IsSing =  vcg::tri::Allocator<MeshType>:: template FindPerVertexAttribute<bool>(m,std::string("Singular"));
    return IsSing[v];
}

template <class MeshType>
VertType &VertexKind(MeshType &m, typename MeshType::VertexType &v)
{
    assert(vcg::tri::HasPerVertexAttribute(m,std::string("VertexKind")));
    typename MeshType::template PerVertexAttributeHandle<VertType> VKind =  vcg::tri::Allocator<MeshType>:: template FindPerVertexAttribute<VertType>(m,std::string("VertexKind"));
    return VKind[v];
}

template <class MeshType>
bool &IsConstrained(MeshType &m, typename MeshType::VertexType &v)
{
    assert(vcg::tri::HasPerVertexAttribute(m,std::string("IsConstrained")));
    typename MeshType::template PerVertexAttributeHandle<bool> VConstr =  vcg::tri::Allocator<MeshType>:: template FindPerVertexAttribute<bool>(m,std::string("IsConstrained"));
    return VConstr[v];
}

template <class MeshType>
int &Region(MeshType &m, typename MeshType::FaceType &f)
{
    assert(vcg::tri::HasPerFaceAttribute(m,std::string("Region")));
    typename MeshType::template PerFaceAttributeHandle<int> FRegion =  vcg::tri::Allocator<MeshType>:: template FindPerFaceAttribute<int>(m,std::string("Region"));
    return FRegion[f];
}

template <class MeshType>
int &PathIndex(MeshType &m, typename MeshType::FaceType &f,int IndexV)
{
    assert(IndexV>=0);
    assert(IndexV<3);
    assert(vcg::tri::HasPerFaceAttribute(m,std::string("PathIndex")));
    typename MeshType::template PerFaceAttributeHandle<EdgePath> PIndex =  vcg::tri::Allocator<MeshType>:: template FindPerFaceAttribute<EdgePath>(m,std::string("PathIndex"));
    return (PIndex[f].PerEdgePath[IndexV]);
}


template <class TriMeshType>
void ClearPathIndexes(TriMeshType &m,typename TriMeshType::FaceType &f)
{

    for (int j=0;j<3;j++)
    {
        //f.ClearCrease(j);
        f.SetFaceEdgeS(j);
        //PathIndex[j]=-1;
        PathIndex(m,f,j)=-1;
    }
}

template <class TriMeshType>
void GetPathMap(TriMeshType &m,typename TriMeshType::FaceType &f,
                std::map<std::pair<typename TriMeshType::CoordType,typename TriMeshType::CoordType>,int> &PathMap)
{
    typedef typename TriMeshType::CoordType CoordType;
    PathMap.clear();
    for (int j=0;j<3;j++)
    {
        //if (PathIndex[j]==-1)continue;
        if (PathIndex(m,f,j)==-1)continue;

        CoordType pos0,pos1;
        pos0=f.P0(j);
        pos1=f.P1(j);
        std::pair<CoordType,CoordType> key(std::min(pos0,pos1),
                                           std::max(pos0,pos1));
        //PathMap[key]=PathIndex[j];
        PathMap[key]=PathIndex(m,f,j);
    }
}

template <class TriMeshType>
void SetFromPathMap(TriMeshType &m,typename TriMeshType::FaceType &f,
                    std::map<std::pair<typename TriMeshType::CoordType,typename TriMeshType::CoordType>,int> &PathMap)
{
    typedef typename TriMeshType::CoordType CoordType;
    //clear all creases
    ClearPathIndexes(m,f);
    //then uptate it
    for (int j=0;j<3;j++)
    {
        CoordType Pos0=f.P0(j);
        CoordType Pos1=f.P1(j);
        std::pair<CoordType,CoordType> key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        if (PathMap.count(key)==0)continue;
        //f.SetCrease(j);
        f.SetFaceEdgeS(j);
        //PathIndex[j]=PathMap[key];
        PathIndex(m,f,j)=PathMap[key];
    }
}

template <class TriMeshType>
size_t numRegions(TriMeshType &m)
{
    std::vector<size_t> VertRegion;
    for (size_t i=0;i<m.face.size();i++)
    {
        typename TriMeshType::FaceType* currF=&m.face[i];
        VertRegion.push_back(Region(m,*currF));
    }
    std::sort(VertRegion.begin(),VertRegion.end());
    std::vector<size_t>::iterator last = std::unique(VertRegion.begin(),VertRegion.end());
    VertRegion.erase(last, VertRegion.end());
    return (VertRegion.size());
}

template <class TriMeshType>
void ColorFaceRegions(TriMeshType &m)
{
    int NumR=numRegions(m);
    for (size_t i=0;i<m.face.size();i++)
        m.face[i].C()=vcg::Color4b::Scatter(NumR,(size_t)Region(m,m.face[i]));
}

template <class TriMeshType>
void UpdateAttributes(TriMeshType &m)
{
    vcg::tri::UpdateNormal<TriMeshType>::PerFaceNormalized(m);
    vcg::tri::UpdateNormal<TriMeshType>::PerVertexNormalized(m);
    vcg::tri::UpdateBounding<TriMeshType>::Box(m);
    vcg::tri::UpdateTopology<TriMeshType>::FaceFace(m);
    vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);
    vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromFF(m);
    vcg::tri::UpdateFlags<TriMeshType>::VertexBorderFromFaceBorder(m);
}

template <class TriMeshType>
typename TriMeshType::ScalarType Area(TriMeshType &m)
{
    typename TriMeshType::ScalarType AreaVal=0;
    for (size_t i=0;i<m.face.size();i++)
        AreaVal+=vcg::DoubleArea(m.face[i]);
    return AreaVal;
}

//template <class TriMeshType>
//void Parameterize(TriMeshType &m)
//{
//    vcg::tri::OptimizeUV_LSCM(m,TriMeshType::VertexType::SELECTED);
//    vcg::tri::UV_Utils<TriMeshType>::CopyVertUVWedge(m);
//}

template <class TriMeshType>
void AssignFaceRegion(TriMeshType &m)
{
    size_t numRegions=0;
    vcg::tri::UpdateFlags<TriMeshType>::FaceClearV(m);
    for (size_t i=0;i<m.face.size();i++)
    {
        if (m.face[i].IsV())continue;
        std::vector<typename TriMeshType::FaceType*> explFaces;
        explFaces.push_back(&m.face[i]);
        do{
            typename TriMeshType::FaceType* currF=explFaces.back();
            explFaces.pop_back();
            //currF->Q()=numRegions;
            //currF->Region=numRegions;
            Region(m,*currF)=numRegions;

            currF->SetV();
            for (size_t j=0;j<3;j++)
            {
                //if (currF->IsCrease(j))continue;
                if (currF->IsFaceEdgeS(j))continue;
                typename TriMeshType::FaceType* nextF=currF->FFp(j);
                if (nextF->IsV())continue;
                explFaces.push_back(nextF);
            }
        }while (!explFaces.empty());
        numRegions++;
    }
    std::cout << "Found " << numRegions <<" regions" << std::endl;
}

template <class TriMeshType>
void GetNeighRegions(TriMeshType &m,int indexV,std::vector<size_t> &VertReg)
{
    std::vector<typename TriMeshType::FaceType *> faceVec;
    std::vector<int> indexes;
    //get the faces sharing that vertex
    vcg::face::VFStarVF<typename TriMeshType::FaceType>(&m.vert[indexV],faceVec,indexes);

    //then get the number of regions in the neighborhood
    VertReg.clear();
    for (size_t j=0;j<faceVec.size();j++)
        //VertReg.push_back(faceVec[j]->Region);
        VertReg.push_back(Region(m,*faceVec[j]));


    std::sort(VertReg.begin(),VertReg.end());
    std::vector<size_t>::iterator it;
    it = std::unique (VertReg.begin(), VertReg.end());
    VertReg.resize( std::distance(VertReg.begin(),it) );
}

template <class TriMeshType>
void SetVertexPathOnQ(TriMeshType &m)
{
    vcg::tri::UpdateQuality<TriMeshType>::VertexConstant(m,-1);

    //set on paths
    for (size_t i=0;i<m.face.size();i++)
    {
        for (size_t j=0;j<3;j++)
        {
            //if (face[i].PathIndex[j]!=-1)

            if (PathIndex<TriMeshType>(m,m.face[i],(int)j)!=-1)
            {
                //int IndexP=face[i].PathIndex[j];
                int IndexP=PathIndex<TriMeshType>(m,m.face[i],(int)j);
                //if (face[i].V0(j)->VertexKind==VPath)
                if (VertexKind(m,*m.face[i].V0(j))==VPath)
                    m.face[i].V0(j)->Q()=IndexP;
                //if (face[i].V1(j)->VertexKind==VPath)
                if (VertexKind(m,*m.face[i].V1(j))==VPath)
                    m.face[i].V1(j)->Q()=IndexP;
            }
        }
    }
}

template <class TriMeshType>
void SelectCorners(TriMeshType &m)
{
    std::vector<typename TriMeshType::CoordType> sampleVec;
    vcg::tri::TrivialSampler<TriMeshType> mps(sampleVec);
    vcg::tri::SurfaceSampling<TriMeshType,vcg::tri::TrivialSampler<TriMeshType> >::VertexBorderCorner(m,mps,vcg::math::ToRad(100.f));

    std::vector<typename TriMeshType::VertexType*> vertSamples;
    vcg::tri::VoronoiProcessing<TriMeshType>::SeedToVertexConversion(m,sampleVec,vertSamples);
    for (size_t i=0;i<vertSamples.size();i++)
    {
        vertSamples[i]->SetS();
        //vertSamples[i]->VertexKind=VCross;
        VertexKind(m,*vertSamples[i])=VCross;
    }
}

template <class TriMeshType>
void SetVertexKind(TriMeshType &m)
{
    vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromFF(m);
    vcg::tri::UpdateFlags<TriMeshType>::VertexBorderFromFaceBorder(m);
    vcg::tri::UpdateFlags<TriMeshType>::VertexClearS(m);


    //set the regular ones
    for (size_t i=0;i<m.vert.size();i++)
        //vert[i].VertexKind=VRegular;
        VertexKind(m,m.vert[i])=VRegular;

    //then set the path ones
    for (size_t i=0;i<m.face.size();i++)
    {
        for (size_t j=0;j<3;j++)
            //if (face[i].PathIndex[j]!=-1)
            if (PathIndex<TriMeshType>(m,m.face[i],(int)j)!=-1)
            {
                // face[i].V0(j)->VertexKind=VPath;
                // face[i].V1(j)->VertexKind=VPath;
                VertexKind<TriMeshType>(m,*m.face[i].V0(j))=VPath;
                VertexKind<TriMeshType>(m,*m.face[i].V1(j))=VPath;
            }
    }

    //selet the crease ones on boder
    SelectCorners(m);

    //set the Cross or T
    for (size_t i=0;i<m.vert.size();i++)
    {
        //then get the number of regions in the neighborhood
        std::vector<size_t> VertReg;
        GetNeighRegions(m,i,VertReg);

        //if the number of neighbors is bigger than 2
        //it could be a t-vertex or a cross vertex
        if (VertReg.size()==4)
        {
            //then select it
            m.vert[i].SetS();
            //vert[i].VertexKind=VCross;
            VertexKind(m,m.vert[i])=VCross;
        }

        if (VertReg.size()==3)
        {
            m.vert[i].SetS();
            //vert[i].VertexKind=VTJunction;
            VertexKind(m,m.vert[i])=VTJunction;
        }

        if ((VertReg.size()==2)&&(m.vert[i].IsB()))
        {
            m.vert[i].SetS();
            //vert[i].VertexKind=VCross;
            VertexKind(m,m.vert[i])=VCross;
        }
    }

    //set the singularities
    for (size_t i=0;i<m.vert.size();i++)
    {
        //if it is not a singular then go on
        //if (!IsSingularVert(vert[i]))continue;
        if (!IsSing(m,m.vert[i]))continue;
        //if (vcg::tri::CrossField<MeshType>::IsSingular(*this,vert[i]))continue;

        //set the singularities
        //vert[i].VertexKind=VSingular;
        VertexKind(m,m.vert[i])=VSingular;
        //and select them
        m.vert[i].SetS();
    }
}

template <class TriMeshType>
typename TriMeshType::ScalarType AVEdgeSize(TriMeshType &m)
{
    typedef typename TriMeshType::ScalarType ScalarType;
    ScalarType sum=0;
    int num=0;
    for (size_t i=0;i<m.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            sum+=(m.face[i].P0(j)-m.face[i].P1(j)).Norm();
            num++;
        }
    return (sum/(ScalarType)num);
}

template <class TriMeshType>
bool OnTVert(TriMeshType &m, vcg::face::Pos<typename TriMeshType::FaceType> CurrP)
{
    assert(CurrP.V()->IsS());//should be on a vertex

    //should not be on a sing or cross
    //int IndexV=vcg::tri::Index(splittedMesh,CurrP.V());

    //if (VertexKind[IndexV]!=VTJunction)return false;

    //if (CurrP.V()->VertexKind!=VTJunction)return false;
    if (VertexKind(m,*CurrP.V())!=VTJunction)return false;

    //then check if from current face is a T-junction or not

    //then get the Path index
    size_t indexF0=vcg::tri::Index(m,CurrP.F());
    size_t edgeI0=CurrP.E();

    //go to the next on border
    CurrP.NextEdgeS();
    //CurrP.Next();
    size_t indexF1=vcg::tri::Index(m,CurrP.F());
    size_t edgeI1=CurrP.E();

    //        int currPath0=EdgeSep[indexF0][edgeI0];
    //        int currPath1=EdgeSep[indexF1][edgeI1];

//        int currPath0=face[indexF0].PathIndex[edgeI0];
//        int currPath1=face[indexF1].PathIndex[edgeI1];
    int currPath0=PathIndex<TriMeshType>(m,m.face[indexF0],(int)edgeI0);
    int currPath1=PathIndex<TriMeshType>(m,m.face[indexF1],(int)edgeI1);

    return (currPath0==currPath1);
    return true;
}

template <class TriMeshType>
bool OnPatchVert(TriMeshType &m,vcg::face::Pos<typename TriMeshType::FaceType> CurrP)
{
    if (!CurrP.V()->IsS())return false;
    //if (!CurrP.IsCrease())return false;
    if (!CurrP.IsEdgeS())return false;

    return (!OnTVert(m,CurrP));
}

template <class TriMeshType>
void InitOriginalV(TriMeshType &m)
{
    for (size_t i=0;i<m.vert.size();i++)
        //vert[i].OriginalV=true;
        OriginalV(m,m.vert[i])=true;
}

template <class TriMeshType>
void EraseSinglePaths(TriMeshType &m)
{
    int erased_edges=0;
    bool corrected=false;
    do{
        corrected=false;
        std::vector<int> CreaseCount(m.vert.size(),0);
        for (size_t i=0;i<m.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!m.face[i].IsCrease(j))continue;
                int IndexV0=vcg::tri::Index(m,m.face[i].V0(j));
                int IndexV1=vcg::tri::Index(m,m.face[i].V1(j));
                if (IndexV0>IndexV1)continue;

                CreaseCount[IndexV0]++;
                CreaseCount[IndexV1]++;
            }
        //then set the vertices to be unCReased
        std::set<typename TriMeshType::VertexType*> to_clean;
        for (size_t i=0;i<m.vert.size();i++)
        {
            if (m.vert[i].IsB())continue;
            if (CreaseCount[i]!=1)continue;
            to_clean.insert(&m.vert[i]);
        }

        for (size_t i=0;i<m.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!m.face[i].IsCrease(j))continue;

                typename TriMeshType::VertexType *v0=m.face[i].V0(j);
                typename TriMeshType::VertexType *v1=m.face[i].V1(j);

                if ((to_clean.count(v0)==0)&&
                        (to_clean.count(v1)==0))continue;

                m.face[i].ClearCrease(j);
                //face[i].PathIndex[j]=-1;
                PathIndex<TriMeshType>(m,m.face[i],(int)j)=-1;

                //go to the other side
                typename TriMeshType::FaceType *f1=m.face[i].FFp(j);
                int OppE=m.face[i].FFi(j);
                f1->ClearCrease(OppE);
                //f1->PathIndex[OppE]=-1;
                PathIndex<TriMeshType>(m,(*f1),OppE)=-1;
                corrected=true;
                erased_edges++;
            }
    }while (corrected);

    std::cout << "*** ERASED "<< erased_edges <<" SINGLE VERTICES ***" <<  std::endl;
}

template <class TriMeshType>
size_t MaxPath(TriMeshType &m)
{
    int MaxP=0;
    for (size_t i=0;i<m.face.size();i++)
        for (int j=0;j<3;j++)
        {
            //if (face[i].PathIndex[j]<MaxP)continue;
            if (PathIndex<TriMeshType>(m,m.face[i],j)<MaxP)continue;

            //MaxP=face[i].PathIndex[j];
            MaxP=PathIndex<TriMeshType>(m,m.face[i],j);
        }
    return MaxP;
}

template <class TriMeshType>
void ClearPathIndexes(TriMeshType &m)
{
    for (size_t i=0;i<m.face.size();i++)
        for (int j=0;j<3;j++)
        {
            //m.face[i].ClearCrease(j);
            m.face[i].ClearFaceEdgeS(j);
            //face[i].PathIndex[j]=-1;
            PathIndex<TriMeshType>(m,m.face[i],j)=-1;
        }
}

template <class TriMeshType>
void GetPathMap(TriMeshType &m,
                std::map<std::pair<typename TriMeshType::CoordType,typename TriMeshType::CoordType>,int> &PathMap)
{
    typedef typename TriMeshType::CoordType CoordType;
    PathMap.clear();
    for (size_t i=0;i<m.face.size();i++)
        for (int j=0;j<3;j++)
        {
            //if (face[i].PathIndex[j]==-1)continue;
            if (PathIndex<TriMeshType>(m,m.face[i],j)==-1)continue;

            CoordType pos0,pos1;
            pos0=m.face[i].P0(j);
            pos1=m.face[i].P1(j);
            std::pair<CoordType,CoordType> key(std::min(pos0,pos1),
                                               std::max(pos0,pos1));
           //PathMap[key]=face[i].PathIndex[j];
            PathMap[key]=PathIndex<TriMeshType>(m,m.face[i],j);
        }
}

template <class TriMeshType>
void SetFromPathMap(TriMeshType &m,
                    std::map<std::pair<typename TriMeshType::CoordType,typename TriMeshType::CoordType>,int> &PathMap)
{
   typedef typename TriMeshType::CoordType CoordType;

    //clear all creases
    ClearPathIndexes(m);
    //then uptate it
    for (size_t i=0;i<m.face.size();i++)
        for (int j=0;j<3;j++)
        {
            CoordType Pos0=m.face[i].P0(j);
            CoordType Pos1=m.face[i].P1(j);
            std::pair<CoordType,CoordType> key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
            if (PathMap.count(key)==0)continue;
            //m.face[i].SetCrease(j);
            m.face[i].SetFaceEdgeS(j);
            //face[i].PathIndex[j]=PathMap[key];
            PathIndex<TriMeshType>(m,m.face[i],j)=PathMap[key];
        }
}

template <class TriMeshType>
void SetBorderAsPath(TriMeshType &m)
{
    vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromFF(m);
    //then update it
    size_t MPath=MaxPath(m);
    for (size_t i=0;i<m.face.size();i++)
        for (int j=0;j<3;j++)
        {
            if (!m.face[i].IsB(j))continue;
            //m.face[i].SetCrease(j);
            m.face[i].SetFaceEdgeS(j);
            //face[i].PathIndex[j]=MPath+1;
            PathIndex<TriMeshType>(m,m.face[i],j)=MPath+1;
        }
}

template <class TriMeshType>
void SelectFacesOfPatch(TriMeshType &m,int PatchIndex)
{
    vcg::tri::UpdateFlags<TriMeshType>::FaceClearS(m);
    for (size_t i=0;i<m.face.size();i++)
        //if (face[i].Region==PatchIndex)face[i].SetS();
        if (Region(m,m.face[i])==PatchIndex)m.face[i].SetS();
}

template <class TriMeshType>
void SelectPathVertices(TriMeshType &m)
{
    vcg::tri::UpdateFlags<TriMeshType>::VertexClearS(m);
    for (size_t i=0;i<m.vert.size();i++)
    {
        //if (vert[i].VertexKind==VRegular)continue;
        if (VertexKind(m,m.vert[i])==VRegular)continue;
        m.vert[i].SetS();
    }
}

#endif
