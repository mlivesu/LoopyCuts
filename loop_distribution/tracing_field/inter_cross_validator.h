#ifndef INTER_CROSS_VALIDATOR
#define INTER_CROSS_VALIDATOR

#include <tracing_field/conflict_finder.h>

//template < class MeshType >
//class InterCrossValidator {
//    typedef typename MeshType::CoordType    CoordType;
//    typedef typename MeshType::VertexType   VertexType;
//    typedef typename MeshType::FaceType     FaceType;
//    typedef typename MeshType::ScalarType   ScalarType;

//    AnisotropicGraph<MeshType> &anigraph;

//    typedef typename AnisotropicGraph<MeshType>::Path PathType;


//public:

//    //    struct UnsolvedInfo
//    //    {
//    //        size_t PathIndex;
//    //        size_t Start0,End0;
//    //        size_t Start1,End1;

//    //        UnsolvedInfo(size_t _PathIndex0,
//    //                     size_t _PathIndex1,
//    //                     size_t _Start0,
//    //                     size_t _End0,
//    //                     size_t _Start1,
//    //                     size_t _End1)
//    //        {
//    //            PathIndex0=_PathIndex0;
//    //            PathIndex1=_PathIndex1;
//    //            Start0=_Start0;
//    //            Start1=_Start1;
//    //            End0=_End0;
//    //            End1=_End1;
//    //        }
//    //    };

//    struct UnsolvedInfo
//    {
//        size_t PathIndex;
//        size_t Start,End;

//        UnsolvedInfo(size_t _PathIndex,
//                     size_t _Start,
//                     size_t _End)
//        {
//            PathIndex=_PathIndex;
//            Start=_Start;
//            End=_End;
//        }
//    };

//    struct IntersectionInterval
//    {
//        int PathIndex;
//        int NodeIndex0;
//        int NodeIndex1;

//        IntersectionInterval(size_t &_PathIndex,
//                             size_t &_NodeIndex0,
//                             size_t &_NodeIndex1)
//        {
//            PathIndex=_PathIndex;
//            NodeIndex0=_NodeIndex0;
//            NodeIndex1=_NodeIndex1;
//        }

//        inline bool operator <(const IntersectionInterval &IInt)const
//        {
//            if (NodeIndex0==IInt.NodeIndex0)
//                return (NodeIndex1<IInt.NodeIndex1);
//            return (NodeIndex0<IInt.NodeIndex0);
//        }

//        void Check()
//        {
//            assert(PathIndex>=0);
//            assert(NodeIndex0>=0);
//            assert(NodeIndex1>=0);
//        }

//        IntersectionInterval()
//        {
//            PathIndex=-1;
//            NodeIndex0=-1;
//            NodeIndex1=-1;
//        }
//    };

//    void GetCrossUnsolved(std::vector<PathType> &TestPath,
//                          std::vector<UnsolvedInfo> &CrossUnsolved,
//                          bool PrintDebug=true)const
//    {
//        CrossUnsolved.clear();

//        ConflictFinder<MeshType> CFinder(anigraph);
//        std::vector<bool> IsLoop(TestPath.size(),true);

//        if (PrintDebug)
//        {
//            std::cout<<"There are "<<TestPath.size()<<" loops to test"<<std::endl;
//        }
//        std::vector<typename ConflictFinder<MeshType>::IntPoint> IPoint;
//        CFinder.FindCrossIntersections(TestPath,IsLoop,IPoint);

//        if (PrintDebug)
//        {
//            std::cout<<"There are "<<IPoint.size()<<" intersections"<<std::endl;
//        }

//        //then sort the intersections per loop, per loop pairs, node other loops
//        std::vector<std::vector<IntersectionInterval> > SortedCross(TestPath.size());
//        for (size_t i=0;i<IPoint.size();i++)
//        {
//            size_t IndexP0=IPoint[i].Path0;
//            size_t IndexP1=IPoint[i].Path1;

//            assert(IndexP0<SortedCross.size());
//            assert(IndexP1<SortedCross.size());
//            assert(IndexP0<TestPath.size());
//            assert(IndexP1<TestPath.size());

//            size_t IndexP0N0=IPoint[i].Interval0.first;
//            size_t IndexP0N1=IPoint[i].Interval0.second;

//            size_t IndexP1N0=IPoint[i].Interval1.first;
//            size_t IndexP1N1=IPoint[i].Interval1.second;

//            assert(IndexP0N0<TestPath[IndexP0].nodes.size());
//            assert(IndexP0N1<TestPath[IndexP0].nodes.size());
//            assert(IndexP1N0<TestPath[IndexP1].nodes.size());
//            assert(IndexP1N1<TestPath[IndexP1].nodes.size());

//            SortedCross[IndexP0].push_back(IntersectionInterval(IndexP1,IndexP0N0,IndexP0N1));
//            SortedCross[IndexP1].push_back(IntersectionInterval(IndexP0,IndexP1N0,IndexP1N1));
//        }
//        for (size_t i=0;i<SortedCross.size();i++)
//            std::sort(SortedCross[i].begin(),SortedCross[i].end());

//        for (size_t i=0;i<SortedCross.size();i++)
//        {
//            //no crossing, must be inserted
//            if (SortedCross[i].size()<=2)
//            {
//                CrossUnsolved.push_back(UnsolvedInfo(i,0,TestPath[i].nodes.size()-1));
//                if (PrintDebug)
//                    std::cout<<"Smaller or Equal to 2 Intersections"<<std::endl;
//                continue;
//            }
//            //check the crossing sequence
//            for (size_t j=0;j<SortedCross[i].size();j++)
//            {
//                size_t CurrI=j;
//                size_t NextI=(j+1)%SortedCross[i].size();
//                size_t IndexPath0=SortedCross[i][CurrI].PathIndex;
//                size_t IndexPath1=SortedCross[i][NextI].PathIndex;
//                if (IndexPath0!=IndexPath1)continue;

//                if (PrintDebug)
//                {
//                    std::cout<<"Current Path "<<i<<std::endl;
//                    std::cout<<"Index0 "<<IndexPath0<<std::endl;
//                    std::cout<<"Index1 "<<IndexPath1<<std::endl;
//                }
//                //then check the interval for the node
//                size_t NodeI0=SortedCross[i][CurrI].NodeIndex1;
//                size_t NodeI1=SortedCross[i][NextI].NodeIndex0;
//                CrossUnsolved.push_back(UnsolvedInfo(i,NodeI0,NodeI1));
//            }
//        }
//        if (!PrintDebug)return;

//        //debug test
//        std::cout<<"*** UNSOLVED CROSS ***"<<std::endl;
//        std::cout<<"There are "<<CrossUnsolved.size()<<" unsolved crossed loops"<<std::endl;
//        for (size_t i=0;i<CrossUnsolved.size();i++)
//        {
//            std::cout<<"* Loop "<<CrossUnsolved[i].PathIndex<<std::endl;
//            std::cout<<"-Lenght "<<TestPath[CrossUnsolved[i].PathIndex].nodes.size()<<std::endl;
//            std::cout<<"    - between "<<CrossUnsolved[i].Start<<" and "<<CrossUnsolved[i].End<<std::endl;
//        }
//    }

//    InterCrossValidator(AnisotropicGraph<MeshType> &_anigraph) : anigraph(_anigraph)
//    {}

//};


#include <tracing_field/conflict_finder.h>

template < class MeshType >
class InterCrossValidator {
    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;

    AnisotropicGraph<MeshType> &anigraph;

    typedef typename AnisotropicGraph<MeshType>::Path PathType;


public:

    struct UnsolvedInfo
    {
        size_t PathIndex;
        size_t Start,End;

        UnsolvedInfo(size_t _PathIndex,
                     size_t _Start,
                     size_t _End)
        {
            PathIndex=_PathIndex;
            Start=_Start;
            End=_End;
        }
    };

    struct IntersectionInterval
    {
        int PathIndex;
        int NodeIndex0;
        int NodeIndex1;

        IntersectionInterval(int _PathIndex,
                             int _NodeIndex0,
                             int _NodeIndex1)
        {
            PathIndex=_PathIndex;
            NodeIndex0=_NodeIndex0;
            NodeIndex1=_NodeIndex1;
        }

        inline bool operator <(const IntersectionInterval &IInt)const
        {
            if (NodeIndex0==IInt.NodeIndex0)
                return (NodeIndex1<IInt.NodeIndex1);
            return (NodeIndex0<IInt.NodeIndex0);
        }

        void Check()
        {
            assert(PathIndex>=0);
            assert(NodeIndex0>=0);
            assert(NodeIndex1>=0);
        }

        IntersectionInterval()
        {
            PathIndex=-1;
            NodeIndex0=-1;
            NodeIndex1=-1;
        }
    };

    //    void GetCrossUnsolved(std::vector<PathType> &TestPath,
    //                          std::vector<UnsolvedInfo> &CrossUnsolved,
    //                          bool PrintDebug=true)const
    void GetCrossUnsolved(std::vector<std::vector<vcg::face::Pos<FaceType> > > &SharpFeatures,
                          std::vector<PathType> &TestPath,
                          std::vector<UnsolvedInfo> &CrossUnsolved,
                          bool PrintDebug=true)const
    {
        CrossUnsolved.clear();

        ConflictFinder<MeshType> CFinder(anigraph);
        std::vector<bool> IsLoop(TestPath.size(),true);

        if (PrintDebug)
        {
            std::cout<<"There are "<<TestPath.size()<<" loops to test"<<std::endl;
        }

        std::vector<typename ConflictFinder<MeshType>::IntPoint> IPoint;
        CFinder.FindCrossIntersections(TestPath,IsLoop,IPoint);

        if (PrintDebug)
        {
            std::cout<<"There are "<<IPoint.size()<<" intersections"<<std::endl;
        }

        //also get the intersections with sharp features
        std::map<std::pair<int,int>, int > FeatureFaceEdges;
        for (size_t i=0;i<SharpFeatures.size();i++)
            for (size_t j=0;j<SharpFeatures[i].size();j++)
            {
                vcg::face::Pos<FaceType> SharpPos=SharpFeatures[i][j];
                FaceType *f=SharpPos.F();
                int FaceI=vcg::tri::Index(anigraph.Mesh(),f);
                int EdgeI=SharpPos.E();
                SharpPos.FlipF();
                f=SharpPos.F();
                int FaceIOpp=vcg::tri::Index(anigraph.Mesh(),f);
                int EdgeIOpp=SharpPos.E();
                assert(FaceI!=FaceIOpp);
                assert(EdgeI>=0);
                assert(EdgeIOpp>=0);
                std::pair<int,int> key0(FaceI,EdgeI);
                std::pair<int,int> key1(FaceIOpp,EdgeIOpp);
                //only one edge per feature
                assert(FeatureFaceEdges.count(key0)==0);
                assert(FeatureFaceEdges.count(key1)==0);
                FeatureFaceEdges[key0]=i;
                FeatureFaceEdges[key1]=i;
            }

        std::vector<std::vector<IntersectionInterval> > SortedCross(TestPath.size());
        for (size_t i=0;i<TestPath.size();i++)
            for (size_t j=0;j<TestPath[i].nodes.size();j++)
            {
                size_t currN=TestPath[i].nodes[j];
                size_t FaceI,EdgeI,M4Dir;
                anigraph.GetFaceEdgeDir(currN,FaceI,EdgeI,M4Dir);
                std::pair<int,int> key((int)FaceI,(int)EdgeI);
                if (FeatureFaceEdges.count(key)==0)continue;
                int FeatureIndex=-FeatureFaceEdges[key];
                //negative for the sharp features
                SortedCross[i].push_back(IntersectionInterval(FeatureIndex,(int)j,(int)j));
            }

        //then sort the intersections per loop, per loop pairs, node other loops
        for (size_t i=0;i<IPoint.size();i++)
        {

            size_t IndexP0=IPoint[i].Path0;
            size_t IndexP1=IPoint[i].Path1;

            assert(IndexP0<SortedCross.size());
            assert(IndexP1<SortedCross.size());
            assert(IndexP0<TestPath.size());
            assert(IndexP1<TestPath.size());

            size_t IndexP0N0=IPoint[i].Interval0.first;
            size_t IndexP0N1=IPoint[i].Interval0.second;

            size_t IndexP1N0=IPoint[i].Interval1.first;
            size_t IndexP1N1=IPoint[i].Interval1.second;

            assert(IndexP0N0<TestPath[IndexP0].nodes.size());
            assert(IndexP0N1<TestPath[IndexP0].nodes.size());
            assert(IndexP1N0<TestPath[IndexP1].nodes.size());
            assert(IndexP1N1<TestPath[IndexP1].nodes.size());

            SortedCross[IndexP0].push_back(IntersectionInterval((int)IndexP1,(int)IndexP0N0,(int)IndexP0N1));
            SortedCross[IndexP1].push_back(IntersectionInterval((int)IndexP0,(int)IndexP1N0,(int)IndexP1N1));
        }
        for (size_t i=0;i<SortedCross.size();i++)
            std::sort(SortedCross[i].begin(),SortedCross[i].end());

        for (size_t i=0;i<SortedCross.size();i++)
        {
            //no crossing, must be inserted
            if (SortedCross[i].size()<=2)
            {
                CrossUnsolved.push_back(UnsolvedInfo(i,0,TestPath[i].nodes.size()-1));
                if (PrintDebug)
                    std::cout<<"Smaller or Equal to 2 Intersections"<<std::endl;
                continue;
            }
            //check the crossing sequence
            for (size_t j=0;j<SortedCross[i].size();j++)
            {
                size_t CurrI=j;
                size_t NextI=(j+1)%SortedCross[i].size();
                size_t IndexPath0=SortedCross[i][CurrI].PathIndex;
                size_t IndexPath1=SortedCross[i][NextI].PathIndex;
                if (IndexPath0!=IndexPath1)continue;

                if (PrintDebug)
                {
                    std::cout<<"Current Path "<<i<<std::endl;
                    std::cout<<"Index0 "<<IndexPath0<<std::endl;
                    std::cout<<"Index1 "<<IndexPath1<<std::endl;
                }
                //then check the interval for the node
                size_t NodeI0=SortedCross[i][CurrI].NodeIndex1;
                size_t NodeI1=SortedCross[i][NextI].NodeIndex0;
//                assert(NodeI0>=0);
//                assert(NodeI1>=0);
//                assert(NodeI0!=NodeI1);
                if (NodeI0==NodeI1)//single intersection
                {
                    assert(NodeI0>=0);
                    assert(NodeI1>=0);
                    size_t NumNodes=TestPath[i].nodes.size();
                    size_t Step1_4=std::max((size_t)1,NumNodes/4);
                    size_t Step1_2=std::max((size_t)1,NumNodes/2);
                    NodeI0=(NodeI0+Step1_4)%NumNodes;
                    NodeI1=(NodeI0+Step1_2)%NumNodes;
                    assert(NodeI0!=NodeI1);
                }
                CrossUnsolved.push_back(UnsolvedInfo(i,NodeI0,NodeI1));
            }
        }
        if (!PrintDebug)return;

        //debug test
        std::cout<<"*** UNSOLVED CROSS ***"<<std::endl;
        std::cout<<"There are "<<CrossUnsolved.size()<<" unsolved crossed loops"<<std::endl;
        for (size_t i=0;i<CrossUnsolved.size();i++)
        {
            std::cout<<"* Loop "<<CrossUnsolved[i].PathIndex<<std::endl;
            std::cout<<"-Lenght "<<TestPath[CrossUnsolved[i].PathIndex].nodes.size()<<std::endl;
            std::cout<<"    - between "<<CrossUnsolved[i].Start<<" and "<<CrossUnsolved[i].End<<std::endl;
        }
    }

    InterCrossValidator(AnisotropicGraph<MeshType> &_anigraph) : anigraph(_anigraph)
    {}

};
#endif //LOOP_OPTIMIZER
