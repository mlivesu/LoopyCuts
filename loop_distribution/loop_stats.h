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

#ifndef LOOP_STATS_H
#define LOOP_STATS_H

template <class MeshType>
class PathStats
{

    static void IndexFromPos(const MeshType & mesh,
                             const std::vector<std::pair<size_t,size_t> > &FaceEdgePath,
                             std::vector<std::pair<size_t,size_t> > &EdgePair)
    {
        EdgePair.clear();
        for (size_t i=0;i<FaceEdgePath.size();i++)
        {
            size_t IndexF=FaceEdgePath[i].first;
            size_t IndexE=FaceEdgePath[i].second;

            int IndexV0=vcg::tri::Index(mesh,mesh.face[IndexF].cV0(IndexE));
            int IndexV1=vcg::tri::Index(mesh,mesh.face[IndexF].cV1(IndexE));
            std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));

            EdgePair.push_back(Key);
        }
    }

    static void IndexFromPos(const MeshType & mesh,
                             const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgePath,
                             std::vector<std::vector<std::pair<size_t,size_t> > > &EdgePair)
    {
        EdgePair.clear();
        EdgePair.resize(FaceEdgePath.size());
        for (size_t i=0;i<FaceEdgePath.size();i++)
            IndexFromPos(mesh,FaceEdgePath[i],EdgePair[i]);
    }

    static void ComputeOverlapPaths(const MeshType & mesh,
                                    const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgePath,
                                    std::vector< std::pair<size_t,size_t> > &OverlapPath)
    {
        std::cout<<"*** CHECK MUTUAL OVERLAPS ***"<<std::endl;

        OverlapPath.clear();

        std::vector<std::vector<std::pair<size_t,size_t> > > EdgePair;
        IndexFromPos(mesh,FaceEdgePath,EdgePair);

        for (size_t i=0;i<EdgePair.size();i++)
            std::sort(EdgePair[i].begin(),EdgePair[i].end());

        //do intersection
        for (size_t i=0;i<EdgePair.size()-1;i++)
            for (size_t j=(i+1);j<EdgePair.size();j++)
            {
                std::vector<std::pair<size_t,size_t> > OverlapEdges;
                std::set_intersection(EdgePair[i].begin(),EdgePair[i].end(),
                                      EdgePair[j].begin(),EdgePair[j].end(),
                                      std::back_inserter(OverlapEdges));
                if (OverlapEdges.size()==0)continue;

                OverlapPath.push_back(std::pair<size_t,size_t>(i,j));
            }
    }

    static void ComputeSelfOverlapPaths(const MeshType & mesh,
                                        const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgePath,
                                        std::vector<size_t > &SelfOverlapPath)
    {
        std::cout<<"*** CHECK SELF OVERLAPS ***"<<std::endl;
        SelfOverlapPath.clear();

        std::vector<std::vector<std::pair<size_t,size_t> > > EdgePair;
        IndexFromPos(mesh,FaceEdgePath,EdgePair);

        //sort them
        for (size_t i=0;i<EdgePair.size();i++)
        {
            size_t size0=EdgePair[i].size();
            std::sort(EdgePair[i].begin(),EdgePair[i].end());
            auto last = std::unique(EdgePair[i].begin(),EdgePair[i].end());
            EdgePair[i].erase(last, EdgePair[i].end());
            size_t size1=EdgePair[i].size();
            if (size1==size0)continue;
            SelfOverlapPath.push_back(i);
        }
    }

//    static void ComputeSelfOverlapPaths(const MeshType & mesh,
//                                        const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgePath,
//                                        std::vector<size_t > &SelfOverlapLoop)
//    {
//        OverlapLoop.clear();

//        std::vector<std::vector<std::pair<size_t,size_t> > > EdgePair;
//        IndexFromPos(mesh,FaceEdgePath,EdgePair);

//        //sort them
//        for (size_t i=0;i<EdgePair.size();i++)
//        {
//            size_t size0=EdgePair[i].size();
//            std::sort(EdgePair[i].begin(),EdgePair[i].end());
//            auto last = std::unique(EdgePair[i].begin(),EdgePair[i].end());
//            EdgePair[i].erase(last, v.end());
//            size_t size1=EdgePair[i].size();
//            if (size1==size0)continue;
//            SelfOverlapLoop.push_back(i);
//        }
//    }

    static bool SelfIntersect(std::vector<std::pair<size_t,size_t> > &EdgePair)
    {
        //remove overlapping ones
        std::sort(EdgePair.begin(),EdgePair.end());
        auto last = std::unique(EdgePair.begin(),EdgePair.end());
        EdgePair.erase(last, EdgePair.end());

        //then set the count for each vertex the amount
        std::map<size_t,size_t> Count;
        for (size_t i=0;i<EdgePair.size();i++)
        {
            int IndexV0=EdgePair[i].first;
            int IndexV1=EdgePair[i].second;
            if (Count.count(IndexV0)==0)
                Count[IndexV0]=1;
            else
                Count[IndexV0]++;

            if (Count.count(IndexV1)==0)
                Count[IndexV1]=1;
            else
                Count[IndexV1]++;
        }

        std::map<size_t,size_t>::iterator IteInters;
        for (IteInters=Count.begin();IteInters!=Count.end();IteInters++)
        {
            if ((*IteInters).second>2)return true;
        }
        return false;
    }


    static void ComputeSelfIntersectingPaths(const MeshType & mesh,
                                             const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgePath,
                                             std::vector<size_t > &SelfIntersectLoop)
    {
        std::cout<<"*** CHECK SELF INTERSECTING PATHS ***"<<std::endl;

        SelfIntersectLoop.clear();

        std::vector<std::vector<std::pair<size_t,size_t> > > EdgePair;
        IndexFromPos(mesh,FaceEdgePath,EdgePair);

        for (size_t i=0;i<EdgePair.size();i++)
        {
            if (!SelfIntersect(EdgePair[i]))continue;
            SelfIntersectLoop.push_back(i);
        }
    }

    static void RemoveEndPoints(std::vector<size_t> &VertSeq)
    {
        std::map<size_t,size_t> VertCount;
        for (size_t i=0;i<VertSeq.size();i++)
        {
            size_t IndexV=VertSeq[i];
            if (VertCount.count(IndexV)==0)
                VertCount[IndexV]=1;
            else
                VertCount[IndexV]++;
        }

        std::vector<size_t> NewVertSeq;
        for (size_t i=0;i<VertSeq.size();i++)
        {
            size_t IndexV=VertSeq[i];
            if (VertCount[IndexV]==1)continue;
            NewVertSeq.push_back(IndexV);
        }
        VertSeq=NewVertSeq;
    }

    static size_t ComputeIntersectingCardinality(const std::vector<std::pair<size_t,size_t> >  &EdgePair0,
                                                 const std::vector<std::pair<size_t,size_t> >  &EdgePair1)
    {
        std::vector<std::pair<size_t,size_t> > EdgePair0T=EdgePair0;
        std::vector<std::pair<size_t,size_t> > EdgePair1T=EdgePair1;

        //remove self overlaps
        std::sort(EdgePair0T.begin(),EdgePair0T.end());
        auto last = std::unique(EdgePair0T.begin(),EdgePair0T.end());
        EdgePair0T.erase(last, EdgePair0T.end());

        //remove self overlaps
        std::sort(EdgePair1T.begin(),EdgePair1T.end());
        last = std::unique(EdgePair1T.begin(),EdgePair1T.end());
        EdgePair1T.erase(last, EdgePair1T.end());

        std::vector<std::pair<size_t,size_t> > NonOverlapEdges;
        std::set_difference(EdgePair1T.begin(),EdgePair1T.end(),
                            EdgePair0T.begin(),EdgePair0T.end(),
                            std::back_inserter(NonOverlapEdges));

        EdgePair1T=NonOverlapEdges;

        //then get vertices
        std::vector<size_t> Vert0,Vert1;
        for (size_t i=0;i<EdgePair0T.size();i++)
        {
            Vert0.push_back(EdgePair0T[i].first);
            Vert0.push_back(EdgePair0T[i].second);
        }

        RemoveEndPoints(Vert0);

        //erase self intersecrtions
        std::sort(Vert0.begin(),Vert0.end());
        auto last1 = std::unique(Vert0.begin(),Vert0.end());
        Vert0.erase(last1, Vert0.end());

        for (size_t i=0;i<EdgePair1T.size();i++)
        {
            Vert1.push_back(EdgePair1T[i].first);
            Vert1.push_back(EdgePair1T[i].second);
        }
        RemoveEndPoints(Vert1);

        //erase self intersecrtions
        std::sort(Vert1.begin(),Vert1.end());
        last1 = std::unique(Vert1.begin(),Vert1.end());
        Vert1.erase(last1, Vert1.end());

        //finally count!
        Vert0.insert(Vert0.end(),Vert1.begin(),Vert1.end());
        size_t size0=Vert0.size();
        std::sort(Vert0.begin(),Vert0.end());
        last1 = std::unique(Vert0.begin(),Vert0.end());
        Vert0.erase(last1, Vert0.end());
        size_t size1=Vert0.size();
        assert(size1<=size0);
        return (size0-size1);
    }

    static void ComputeIntersectingByCardinality(const MeshType & mesh,
                                                    const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgePath,
                                                    std::vector< std::vector<std::pair<size_t,size_t> > > &IntersectCardinality)
    {
        std::cout<<"*** CHECK CARDINALITY INTERSECTING PATHS ***"<<std::endl;

        IntersectCardinality.clear();

        std::vector<std::vector<std::pair<size_t,size_t> > > EdgePair;
        IndexFromPos(mesh,FaceEdgePath,EdgePair);
        for (size_t i=0;i<EdgePair.size()-1;i++)
            for (size_t j=(i+1);j<EdgePair.size();j++)
            {
               int NumInt=ComputeIntersectingCardinality(EdgePair[i],EdgePair[j]);
               if (IntersectCardinality.size()<(NumInt+1))
                   IntersectCardinality.resize(NumInt+1);
               IntersectCardinality[NumInt].push_back(std::pair<size_t,size_t>(i,j));
            }
    }

public:

    struct Stats
    {
        //the set of paths that self overlao
        std::vector< size_t> SelfOverlapPaths;

        //the set of paths that self intersect
        std::vector< size_t> SelfIntersectPaths;

        //the set of paths that overlaps each other
        std::vector< std::pair<size_t,size_t> > OverlapPaths;

        //the set of paths that intersect by number of cardinality
        std::vector< std::vector<std::pair<size_t,size_t> > > IntersectCardinality;

        void Clear()
        {
            SelfOverlapPaths.clear();
            SelfIntersectPaths.clear();
            OverlapPaths.clear();
            IntersectCardinality.clear();
        }
    };


    static void ComputeStats(const MeshType & mesh,
                             const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgePath,
                             Stats &PStats)
    {

        ComputeSelfOverlapPaths(mesh,FaceEdgePath,PStats.SelfOverlapPaths);

        ComputeOverlapPaths(mesh,FaceEdgePath,PStats.OverlapPaths);

        ComputeSelfIntersectingPaths(mesh,FaceEdgePath,PStats.SelfIntersectPaths);

        ComputeIntersectingByCardinality(mesh,FaceEdgePath,PStats.IntersectCardinality);

    }
};
#endif
