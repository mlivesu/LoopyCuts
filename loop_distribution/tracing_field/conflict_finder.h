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

#ifndef CONFLICT_FINDER
#define CONFLICT_FINDER

//#include <anisotropic_geodesic.h>
#include <tracing_field/graph/anisotropic_graph.h>

/**
 * @brief The class that performs parametrization by considering the connections within the graph
 */
template < class MeshType >
class ConflictFinder {
public:
    /**
     * @brief AnisotropicGraph Constructor.
     * @param _anigraph
     */
    ConflictFinder(AnisotropicGraph<MeshType> &_anigraph) : anigraph(_anigraph){ }

    /**
     * @brief ~AnisotropicGraph Destructor.
     */
    ~ConflictFinder() {}

private:
    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;

    AnisotropicGraph<MeshType> &anigraph;               //the mesh on which the separatrix graph is calculated

    typedef typename AnisotropicGraph<MeshType>::Path PathType;

public:

    struct ConflictInfo
    {
        std::vector<size_t> Paths;

        ConflictInfo(const std::vector<size_t> &_Paths)
        {
            Paths=std::vector<size_t>(_Paths.begin(),_Paths.end());
            std::sort(Paths.begin(),Paths.end());
        }

        ConflictInfo(const ConflictInfo &CInfo)
        {
            Paths=std::vector<size_t>(CInfo.Paths.begin(),CInfo.Paths.end());
        }

        ConflictInfo(){}

        inline bool operator <(const ConflictInfo &CInfo)const
        {
            return (Paths<CInfo.Paths);
        }

        inline bool operator ==(const ConflictInfo &CInfo)const
        {
            return (Paths==CInfo.Paths);
        }
    };

private:

    struct PathStep
    {
        size_t IndexPathNode0;
        size_t IndexPathNode1;
        size_t Node0;
        size_t Node1;
        CoordType Pos0;
        CoordType Pos1;
        size_t FaceIndex;
        size_t M4Dir;
        size_t PathIndex;

        bool ConflictWith(PathStep &PathS1,
                          bool checkRealInt,
                          int NDir=2)
        {
            if (Node0==PathS1.Node0)return true;
            if (Node0==PathS1.Node1)return true;
            if (Node1==PathS1.Node0)return true;
            if (Node1==PathS1.Node1)return true;

            if (FaceIndex!=PathS1.FaceIndex)return false;

            if (NDir>0)
            {
                if ((M4Dir % NDir)!=(PathS1.M4Dir % NDir))return false;
            }
            if (!checkRealInt)return true;

            CoordType dir0=Pos1-Pos0;
            CoordType dir1=PathS1.Pos1-PathS1.Pos0;
            dir0.Normalize();
            dir1.Normalize();

            assert(Pos0!=Pos1);
            assert(PathS1.Pos0!=PathS1.Pos1);
            if (Pos0==PathS1.Pos0)return true;
            if (Pos0==PathS1.Pos1)return true;
            if (Pos1==PathS1.Pos0)return true;
            if (Pos1==PathS1.Pos1)return true;

            //check real intersection
            CoordType Dir0=(Pos0-PathS1.Pos0)^(Pos0-PathS1.Pos1);
            CoordType Dir1=(Pos1-PathS1.Pos0)^(Pos1-PathS1.Pos1);
            Dir0.Normalize();
            Dir1.Normalize();

            //the two direction should be different
            //assert(Dir0*Dir1)
            if ((Dir0*Dir1)<0)return true;

            return false;
        }

        PathStep(size_t _Node0,
                 size_t _Node1,
                 CoordType _Pos0,
                 CoordType _Pos1,
                 size_t _FaceIndex,
                 size_t _M4Dir,
                 size_t _PathIndex,
                 size_t _IndexPathNode0,
                 size_t _IndexPathNode1)
        {
            Node0=_Node0;
            Node1=_Node1;
            Pos0=_Pos0;
            Pos1=_Pos1;
            FaceIndex=_FaceIndex;
            M4Dir=_M4Dir;
            PathIndex=_PathIndex;
            IndexPathNode0=_IndexPathNode0;
            IndexPathNode1=_IndexPathNode1;
        }

        PathStep(const PathStep &PStep)
        {
            Node0=PStep.Node0;
            Node1=PStep.Node1;
            Pos0=PStep.Pos0;
            Pos1=PStep.Pos1;
            FaceIndex=PStep.FaceIndex;
            M4Dir=PStep.M4Dir;
            PathIndex=PStep.PathIndex;
            IndexPathNode0=PStep.IndexPathNode0;
            IndexPathNode1=PStep.IndexPathNode1;
        }
    };

    void AddPossibleConflicts(std::vector<PathStep> &FacePath,
                              std::vector<ConflictInfo> &Conflicts,
                              bool checkRealInt,
                              bool colorize=false,
                              int NDir=2)
    {
        for (int i=0;i<(int)FacePath.size();i++)
        {
            std::vector<size_t> conflPath;
            PathStep P0=FacePath[i];
            for (int j=0;j<(int)FacePath.size();j++)
            {
                if (i==j)continue;
                PathStep P1=FacePath[j];
                if (P0.ConflictWith(P1,checkRealInt,NDir))
                    conflPath.push_back(P1.PathIndex);
            }
            if (conflPath.size()>0)
            {
                if (colorize)anigraph.Mesh().face[P0.FaceIndex].SetS();

                conflPath.push_back(P0.PathIndex);
                ConflictInfo CInfo(conflPath);
                Conflicts.push_back(CInfo);
            }
        }
    }

public:

    struct IntPoint
    {
      CoordType Pos;
      int Path0;
      int Path1;
      std::pair<int,int> Interval0;
      std::pair<int,int> Interval1;

      void CheckValid()const
      {
         assert(Path0>=0);
         assert(Path1>=0);
         assert(Interval0.first>=0);
         assert(Interval0.second>=0);
         assert(Interval1.first>=0);
         assert(Interval1.second>=0);
      }

      inline bool operator ==(const IntPoint &IPoint)const
      {
          CheckValid();
          IPoint.CheckValid();
          return (Pos==IPoint.Pos);
      }

      inline bool operator <(const IntPoint &IPoint)const
      {
          CheckValid();
          IPoint.CheckValid();
          return (Pos<IPoint.Pos);
      }

      IntPoint(CoordType _Pos,
              size_t _Path0,
              size_t _Path1,
              std::pair<size_t,size_t> _Interval0,
              std::pair<size_t,size_t> _Interval1)
      {
          Pos=_Pos;
          Path0=_Path0;
          Path1=_Path1;
          Interval0=_Interval0;
          Interval1=_Interval1;
          CheckValid();
      }

      IntPoint(const IntPoint &IPoint)
      {
          Pos=IPoint.Pos;
          Path0=IPoint.Path0;
          Path1=IPoint.Path1;
          Interval0=IPoint.Interval0;
          Interval1=IPoint.Interval1;
          CheckValid();
          IPoint.CheckValid();
      }

      IntPoint()
      {
          Pos=CoordType(0,0,0);
          Path0=-1;
          Path1=-1;
          Interval0=std::pair<int,int>(-1,-1);
          Interval1=std::pair<int,int>(-1,-1);
      }
    };

private:

    void FindGeometricIntersections(std::vector<PathStep> &FacePath,
                                    std::vector<IntPoint> &IPoint)
    {
        if (FacePath.size()<2)return;
        for (int i=0;i<(int)FacePath.size()-1;i++)
        {
            //std::vector<size_t> conflPath;
            PathStep P0=FacePath[i];
            size_t PIndex0=P0.PathIndex;
            for (int j=(i+1);j<(int)FacePath.size();j++)
            {
                PathStep P1=FacePath[j];
                size_t PIndex1=P1.PathIndex;
                //size_t IndexPathNode1=P1.IndexPathNode;
                //check point to point equivalence
                if (P0.Pos0==P1.Pos0)
                {
                    std::pair<size_t,size_t> Interval0(P0.IndexPathNode0,P0.IndexPathNode0);
                    std::pair<size_t,size_t> Interval1(P1.IndexPathNode0,P1.IndexPathNode0);
                    IntPoint IP(P0.Pos0,PIndex0,PIndex1,Interval0,Interval1);
                    IPoint.push_back(IP);
//                    IntPos.push_back(P0.Pos0);

//                    Path0.push_back(PIndex0);
//                    Path1.push_back(PIndex1);

//                    //put the same one
//                    NodesIndex0.push_back(std::pair<size_t,size_t>(P0.IndexPathNode0,P0.IndexPathNode0));
//                    NodesIndex1.push_back(std::pair<size_t,size_t>(P1.IndexPathNode0,P1.IndexPathNode0));
                    continue;
                }
                if (P0.Pos0==P1.Pos1)
                {
                    std::pair<size_t,size_t> Interval0(P0.IndexPathNode0,P0.IndexPathNode0);
                    std::pair<size_t,size_t> Interval1(P1.IndexPathNode1,P1.IndexPathNode1);
                    IntPoint IP(P0.Pos0,PIndex0,PIndex1,Interval0,Interval1);
                    IPoint.push_back(IP);

//                    IntPos.push_back(P0.Pos0);
//                    Path0.push_back(PIndex0);
//                    Path1.push_back(PIndex1);

//                    //put the same one
//                    NodesIndex0.push_back(std::pair<size_t,size_t>(P0.IndexPathNode0,P0.IndexPathNode0));
//                    NodesIndex1.push_back(std::pair<size_t,size_t>(P1.IndexPathNode1,P1.IndexPathNode1));
                    continue;
                }
                if (P0.Pos1==P1.Pos0)
                {
                    std::pair<size_t,size_t> Interval0(P0.IndexPathNode1,P0.IndexPathNode1);
                    std::pair<size_t,size_t> Interval1(P1.IndexPathNode0,P1.IndexPathNode0);
                    IntPoint IP(P0.Pos1,PIndex0,PIndex1,Interval0,Interval1);
                    IPoint.push_back(IP);

//                    IntPos.push_back(P0.Pos1);
//                    Path0.push_back(PIndex0);
//                    Path1.push_back(PIndex1);

//                    //put the same one
//                    NodesIndex0.push_back(std::pair<size_t,size_t>(P0.IndexPathNode1,P0.IndexPathNode1));
//                    NodesIndex1.push_back(std::pair<size_t,size_t>(P1.IndexPathNode0,P1.IndexPathNode0));
                    continue;
                }
                if (P0.Pos1==P1.Pos1)
                {
                    std::pair<size_t,size_t> Interval0(P0.IndexPathNode1,P0.IndexPathNode1);
                    std::pair<size_t,size_t> Interval1(P1.IndexPathNode1,P1.IndexPathNode1);
                    IntPoint IP(P0.Pos1,PIndex0,PIndex1,Interval0,Interval1);
                    IPoint.push_back(IP);
//                    IntPos.push_back(P0.Pos1);
//                    Path0.push_back(PIndex0);
//                    Path1.push_back(PIndex1);

//                    //put the same one
//                    NodesIndex0.push_back(std::pair<size_t,size_t>(P0.IndexPathNode1,P0.IndexPathNode1));
//                    NodesIndex1.push_back(std::pair<size_t,size_t>(P1.IndexPathNode1,P1.IndexPathNode1));
                    continue;
                }

                if (P0.ConflictWith(P1,true,0))
                {

                    vcg::Segment3<ScalarType> S0(P0.Pos0,P0.Pos1);
                    vcg::Segment3<ScalarType> S1(P1.Pos0,P1.Pos1);

                    CoordType clos0,clos1;
                    ScalarType dist;
                    bool parallel;
                    vcg::SegmentSegmentDistance(S0,S1,dist,parallel,clos0,clos1);

                    std::pair<size_t,size_t> Interval0(P0.IndexPathNode0,P0.IndexPathNode1);
                    std::pair<size_t,size_t> Interval1(P1.IndexPathNode0,P1.IndexPathNode1);
                    IntPoint IP((clos0+clos1)/2,PIndex0,PIndex1,Interval0,Interval1);
                    IPoint.push_back(IP);

//                    IntPos.push_back((clos0+clos1)/2);

//                    Path0.push_back(PIndex0);
//                    Path1.push_back(PIndex1);
//                    NodesIndex0.push_back(std::pair<size_t,size_t>(P0.IndexPathNode0,P0.IndexPathNode1));
//                    NodesIndex1.push_back(std::pair<size_t,size_t>(P1.IndexPathNode0,P1.IndexPathNode1));

                    //NodeIndex.push_back();
                }
            }
        }

    }

    std::vector<std::vector<size_t> > perPathConfl;


    void AddFaceDir(const size_t &indexPath,
                    const PathType &Path,
                    std::vector<std::vector<PathStep> > &PerFacePath,
                    bool IsLoop=false)
    {
        std::vector<std::pair<size_t,size_t> > FaceDir;

        anigraph.GetFacesDir(Path,FaceDir,IsLoop);

        if (!IsLoop)
            assert(FaceDir.size()==Path.nodes.size()-1);
        else
            assert(FaceDir.size()==Path.nodes.size());

        int limit=Path.nodes.size();
        if (!IsLoop)
            limit=Path.nodes.size()-1;

        for (int j=0;j<limit;j++)
        {
            size_t Index0=j;
            size_t Index1=(j+1) % Path.nodes.size();
            //get the two nodes
            size_t N0=Path.nodes[Index0];
            size_t N1=Path.nodes[Index1];

            if(!IsLoop)
            {
                assert(!anigraph.IsSink(N0));
                //this to check the conflicts correctly when they start from the same singularity
                if (anigraph.IsSink(N1))
                    N1=anigraph.OppositeSink(N1);
            }

            CoordType pos0=anigraph.NodePos(N0);
            CoordType pos1=anigraph.NodePos(N1);

            //and the face it passes throught
            int FaceI=FaceDir[j].first;
            int M4Dir=FaceDir[j].second;

            //add the entry in the vector
            assert(FaceI>=0);
            assert(FaceI<(int)PerFacePath.size());
            PerFacePath[FaceI].push_back(PathStep(N0,N1,pos0,pos1,FaceI,M4Dir,indexPath,Index0,Index1));
        }
    }

public:

    bool SelfConflict(const PathType &Path,
                      size_t N_dir=2,
                      bool Loop=true)
    {
        //std::map<size_t,size_t> PerFacePath;
        std::map<CoordType,std::vector<std::pair<size_t,size_t> > > PerNodePath;

        //check if multiple pass throught the same vertex
        size_t Limit=Path.nodes.size();
        int start=0;
        if (!Loop)
        {
            Limit-=2;
            start=1;
        }
        for (size_t i=start;i<Limit;i++)
        {
            int IndexN=Path.nodes[i];
            CoordType currP=anigraph.NodePos(IndexN);
            size_t currD=anigraph.NodeDir(IndexN)%N_dir;
            size_t seqI=(i+1)%Path.nodes.size();
            if (PerNodePath.count(currP)>0)
            {
                for (size_t j=0;j<PerNodePath[currP].size();j++)
                {
                    size_t dir0=PerNodePath[currP][j].first;
                    size_t indexN=PerNodePath[currP][j].second;
                    if ((dir0==currD)&&(indexN!=seqI))return true;
                }
            }
            PerNodePath[currP].push_back(std::pair<size_t,size_t>(currD,i));
        }

        //then check paths
        std::vector<std::vector<PathStep> > PerFacePath;
        PerFacePath.resize(anigraph.NumFaces());

        std::vector<ConflictInfo>  SelfConflicts;
        AddFaceDir(0,Path,PerFacePath,Loop);

        for (size_t i=0;i<PerFacePath.size();i++)//for each face
            AddPossibleConflicts(PerFacePath[i],SelfConflicts,true,false,N_dir);//false);

        return (SelfConflicts.size()>0);
    }


    void SetGeodesicBarrier(const std::vector<PathType> &Paths,
                            bool checkRealInt,
                            bool IsLoop=false,
                            size_t NDir=2)
    {
        //first disable all the nodes that belongs to the Paths
        for (size_t i=0;i<Paths.size();i++)
            for (size_t j=0;j<Paths[i].nodes.size();j++)
            {
                size_t IndexN=Paths[i].nodes[j];
                anigraph.SetValid(IndexN,false);
                //if is a singularity or a sink invalidate also the opposite
                if (anigraph.IsSink(IndexN))
                {
                    size_t OppNode=anigraph.OppositeSink(IndexN);
                    anigraph.SetValid(OppNode,false);
                }
                if (anigraph.IsSingularity(IndexN))
                {
                    size_t OppNode=anigraph.OppositeSing(IndexN);
                    anigraph.SetValid(OppNode,false);
                }
            }

        //then disable all arcs that goes or exit from non valid nodes
        anigraph.InvalidateArcsOnNonValidNodes();

        //initialize per face list of conflicts
        std::vector<std::vector<PathStep> > PerFacePath;
        PerFacePath.resize(anigraph.NumFaces());

        //then run over all paths and add for each face the path passing through
        for (size_t i=0;i<Paths.size();i++)
            AddFaceDir(i,Paths[i],PerFacePath,IsLoop);

        //then check for each path passing from the inserted faces
        //disable all the arcs that cross one of the inserted face path
        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            std::vector<size_t> NeighNodes;
            std::vector<size_t> FacesId;
            std::vector<size_t> M4Dirs;
            std::vector<bool> ValidNeigh;

            //get neighbors
            anigraph.GetNeighInfo(i,NeighNodes,FacesId,M4Dirs,ValidNeigh);

            size_t IndexN0=i;
            CoordType pos0=anigraph.NodePos(IndexN0);

            //then for each neighbor
            for (size_t j=0;j<NeighNodes.size();j++)
            {
                if (!ValidNeigh[j])continue;//is already nn valid
                size_t M4Dir0=M4Dirs[j];
                size_t IndexN1=NeighNodes[j];
                size_t IndexF=FacesId[j];
                CoordType pos1=anigraph.NodePos(IndexN1);
                //create a temporary path step
                size_t LoopI=0;
                PathStep currStep(IndexN0,IndexN1,pos0,pos1,IndexF,M4Dir0,LoopI);
                //then check versus othr path step within the same face
                for (size_t k=0;k<PerFacePath[IndexF].size();k++)
                {
                    if (!currStep.ConflictWith(PerFacePath[IndexF][k],checkRealInt,NDir))continue;
                    ValidNeigh[j]=false;//just found a conflict no need to check anymore
                    break;
                }
            }
            //finally set validity of neighbors for the node
            anigraph.SetNeighValidity(IndexN0,ValidNeigh);
        }
    }

    //THIS WORKS IN GENERAL BUT IS NOT THE MOST CLEAN IMPLEMENTATION
    //MUST CHECK REAL INTERSECTION TO BE SAFE
    bool CrossIntersect(PathType &Path,bool IsLoop)
    {
        return (SelfConflict(Path,1,IsLoop));
        //        std::vector<std::pair<size_t,size_t> > FaceDir;
        //        anigraph.GetFacesDir(Path,FaceDir,IsLoop);
        //        std::set<size_t> FaceCount;

        //        //then run over all paths
        //        for (size_t i=0;i<FaceDir.size();i++)
        //        {
        //            int FaceIndex=FaceDir[i].first;
        //            if (FaceCount.count(FaceIndex)==0)
        //                FaceCount.insert(FaceIndex);
        //            else
        //               return true;
        //        }

        //        return false;
    }

//    void FindCrossIntersections(const PathType &P0,const PathType &P1,
//                                bool IsLoop0,bool IsLoop1,
//                                std::vector<CoordType> &IntersPos)
//    {
//        IntersPos.clear();

//        std::vector<std::pair<size_t,size_t> > CrossPairs,CrossNodes;

//        //initialize per face list of conflicts
//        std::vector<std::vector<PathStep> > PerFacePath;
//        PerFacePath.resize(anigraph.NumFaces());

//        //then add per face paths
//        AddFaceDir(0,P0,PerFacePath,IsLoop0);
//        AddFaceDir(1,P1,PerFacePath,IsLoop1);

//        std::vector<size_t> Path0,Path1;
//        std::vector<std::pair<size_t,size_t> > NodesIndex0,NodesIndex1;

//        for (size_t i=0;i<PerFacePath.size();i++)//for each face
//            FindGeometricIntersections(PerFacePath[i],IntersPos,Path0,Path1,NodesIndex0,NodesIndex1);

//        std::sort(IntersPos.begin(),IntersPos.end());
//        typename std::vector<CoordType>::iterator it;
//        it = std::unique (IntersPos.begin(), IntersPos.end());
//        IntersPos.resize( std::distance(IntersPos.begin(),it) );
//    }

    void FindCrossIntersections(const PathType &P0,const PathType &P1,
                                bool IsLoop0,bool IsLoop1,
                                std::vector<CoordType> &IntersPos)
    {
        IntersPos.clear();

        //std::vector<std::pair<size_t,size_t> > CrossPairs,CrossNodes;

        //initialize per face list of conflicts
        std::vector<std::vector<PathStep> > PerFacePath;
        PerFacePath.resize(anigraph.NumFaces());

        //then add per face paths
        AddFaceDir(0,P0,PerFacePath,IsLoop0);
        AddFaceDir(1,P1,PerFacePath,IsLoop1);

//        std::vector<size_t> Path0,Path1;
//        std::vector<std::pair<size_t,size_t> > NodesIndex0,NodesIndex1;

        std::vector<IntPoint> IPoint;
        for (size_t i=0;i<PerFacePath.size();i++)//for each face
            FindGeometricIntersections(PerFacePath[i],IPoint);

        std::sort(IPoint.begin(),IPoint.end());
        typename std::vector<IntPoint>::iterator it;
        it = std::unique (IPoint.begin(), IPoint.end());
        IPoint.resize( std::distance(IPoint.begin(),it) );

        for (size_t i=0;i<IPoint.size();i++)
            IntersPos.push_back(IPoint[i].Pos);

        size_t size0=IntersPos.size();
        std::sort(IntersPos.begin(),IntersPos.end());
        typename std::vector<CoordType>::iterator it1;
        it1 = std::unique (IntersPos.begin(), IntersPos.end());
        IntersPos.resize( std::distance(IntersPos.begin(),it1) );
        size_t size1=IntersPos.size();
        assert(size0==size1);
    }



    void FindCrossIntersections(const std::vector<PathType> &Paths,
                                const std::vector<bool> &IsLoop,
                                std::vector<IntPoint> &IPoint)
    {

        //initialize per face list of conflicts
        std::vector<std::vector<PathStep> > PerFacePath;
        PerFacePath.resize(anigraph.NumFaces());
        for (size_t i=0;i<Paths.size();i++)
            AddFaceDir(i,Paths[i],PerFacePath,IsLoop[i]);

        IPoint.clear();
        for (size_t i=0;i<PerFacePath.size();i++)
        {
            FindGeometricIntersections(PerFacePath[i],IPoint);
        }

        //then remove duplicates from adjacent faces
        std::sort(IPoint.begin(),IPoint.end());
        typename std::vector<IntPoint>::iterator it;
        it = std::unique (IPoint.begin(), IPoint.end());
        IPoint.resize( std::distance(IPoint.begin(),it) );

    }

    void FindConflicts(const std::vector<PathType> &Paths,
                       std::vector<ConflictInfo> &Conflicts,
                       bool checkRealInt,
                       bool IsLoop=false,
                       bool colorize=false,
                       bool vert_sing=false,
                       bool write_msg=false)
    {
        perPathConfl.clear();
        perPathConfl.resize(Paths.size());

        int t0 = clock();
        if (write_msg)
        std::cout << "Initializing Conflicts" << std::endl;

        Conflicts.clear();

        //initialize per face list of conflicts
        std::vector<std::vector<PathStep> > PerFacePath;
        PerFacePath.resize(anigraph.NumFaces());

        //then run over all paths
        for (size_t i=0;i<Paths.size();i++)
        {
            AddFaceDir(i,Paths[i],PerFacePath,IsLoop);
        }

        //then find conflicts
        for (size_t i=0;i<PerFacePath.size();i++)//for each face
            AddPossibleConflicts(PerFacePath[i],Conflicts,checkRealInt,colorize);


        //add singularity conflicts
        if ((!IsLoop)&&(vert_sing))
        {
            std::map<size_t,std::vector<size_t> > PerExtremeNodePath;
            for (size_t i=0;i<Paths.size();i++)
            {
                size_t N_init=Paths[i].nodes.front();
                size_t N_end=Paths[i].nodes.back();

                assert(anigraph.IsSingularity(N_init));
                assert(anigraph.IsSink(N_end));

                N_end=anigraph.OppositeSink(N_end);

                assert(anigraph.IsSingularity(N_end));

                PerExtremeNodePath[N_init].push_back(i);
                PerExtremeNodePath[N_end].push_back(i);
            }

            typename std::map<size_t,std::vector<size_t> >::iterator ItExtr;
            for (ItExtr=PerExtremeNodePath.begin();ItExtr!=PerExtremeNodePath.end();ItExtr++)
            {
                if ((*ItExtr).second.size()==1)continue;//no conflict.. only one path arrive there!

                ConflictInfo CInfo((*ItExtr).second);
                Conflicts.push_back(CInfo);
            }

        }
        //sort and make them unique
        std::sort(Conflicts.begin(),Conflicts.end());
        typename std::vector<ConflictInfo>::iterator it;
        it = std::unique (Conflicts.begin(), Conflicts.end());
        Conflicts.resize( std::distance(Conflicts.begin(),it) );
        int t1 = clock();
        if (write_msg)
        {
            std::cout << "done. Time " << (float)(t1 - t0)/CLOCKS_PER_SEC << "sec" << std::endl;
            std::cout << "Found " << Conflicts.size() << " Conflicts" << std::endl;
        }
    }

};

#endif //CONFLICT_FINDER
