#ifndef STATS_COLLECTOR
#define STATS_COLLECTOR

#include <iostream>

/**
 * @brief The class that removes small degenerate artifacts from the mesh
 */
template < class ScalarType >
class StatsCollector
{
public:
    /// general info on the field
    int NumSingularities;

    /// this first part is relative to the graph setup
    int NumNodes;
    int EdgeNodes;
    int BoundaryNodes;
    int SingularityNodes;
    int SinkNodes;
    ScalarType timeGraph;

    /// some info on the graph elabouration
    int WrongOpposite;
    int CorrectedSeparatrix;

    /// this part is relative to the separatrix elabouration
    int NumIterations;        ///the number of solving iterations
    int NumNeighbors;         ///the number of per sink neighbors to look at
    int UnresolvedDirections0;/// this are the unresolved directions before the t-junction
    int UnresolvedDirections1;/// this are the unresolved after the t-junction are added
    ScalarType timeSeparatrix;

    ///this part is relative to number of patches per size
    int NumPatches;
    std::vector<size_t> ValencePatches;

    ///this part is relative to parametrization
    int NumVariables;
    int NumUnresolvedT;
    int subDiv;
    ScalarType Percentile;
    ScalarType timeParam;

    ///eventual subdivision steps
    int NumSubPatches;

    /// finally the quadrangulation step
    int NumQuads;
    int NumVert;
    int NumIrregular;

    //print the stats values
    void PrintValues(std::string filepath="",
                     bool printdebug=true)
    {
        std::streambuf *coutbuf;
        std::ofstream out;
        if (filepath!="")
        {
            out=std::ofstream(filepath.c_str());
            coutbuf = std::cout.rdbuf(); //save old buf
            std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
        }

        std::cout << "_______________________"<< std::endl;
        std::cout << "*** FIELD STATS ***" <<  std::endl;
        std::cout << "Num Singularities " << NumSingularities << std::endl;

        std::cout << "_______________________"<< std::endl;
        std::cout << "*** GRAPH STATS ***" <<  std::endl;
        std::cout << "Num Nodes " << NumNodes << std::endl;
        std::cout << "Num Edge Nodes " << EdgeNodes << std::endl;
        std::cout << "Num Boundary Nodes " << BoundaryNodes << std::endl;
        std::cout << "Num Singularity Nodes " << SingularityNodes << std::endl;
        std::cout << "Num Sink Nodes " << SinkNodes << std::endl;
        std::cout << "Time " << timeGraph << std::endl;

        std::cout << "_______________________ " << std::endl;
        std::cout << "*** SEPARATRIX STATS ***" <<  std::endl;
        std::cout << "Num Iterations  " << NumIterations << std::endl;
        std::cout << "Num Neighbors  " << NumNeighbors << std::endl;
        std::cout << "Num Resolved  " << SinkNodes-UnresolvedDirections0 << std::endl;
        std::cout << "Num TJunctions  " << UnresolvedDirections0 << std::endl;
        std::cout << "Time " << timeSeparatrix << std::endl;

        std::cout << "_______________________ "  << std::endl;
        std::cout << "*** PATCH STATS ***" <<  std::endl;
        std::cout << "Num Patches " << NumPatches << std::endl;
        for (size_t i=0;i<ValencePatches.size();i++)
        {
            if (ValencePatches[i]==0)continue;
            std::cout << "Patches with Valence " << i << " :  " << ValencePatches[i] << std::endl;
        }

        std::cout << "_______________________ "  << std::endl;
        std::cout << "*** PARAM STATS ***" <<  std::endl;
        std::cout << "Num Variables " << NumVariables << std::endl;
        std::cout << "Percentile " << Percentile << std::endl;
        std::cout << "Subdivision Factor " << subDiv << std::endl;
        std::cout << "Num Unresolved T " << NumUnresolvedT << std::endl;
        std::cout << "Time " << timeParam << std::endl;

        std::cout << "_______________________ "  << std::endl;
        std::cout << "*** QUAD STATS ***" <<  std::endl;
        std::cout << "Num Quads " << NumQuads << std::endl;
        std::cout << "Num Vertices " << NumVert << std::endl;
        std::cout << "Num Irregular " << NumIrregular << std::endl;

        std::cout << "_______________________ "  << std::endl;
        std::cout << "*** SUB STATS ***" <<  std::endl;
        std::cout << "Num Subdivided Faces " << NumSubPatches << std::endl;
        std::cout << "_______________________ " << std::endl;
        std::cout << "_______________________ " << std::endl;
        std::cout << "_______________________ " << std::endl;
        std::cout << "_______________________ " << std::endl;

        if (printdebug)
        {
            std::cout << "debug graph" <<  std::endl;
            std::cout << "Num Wrong Opposites " << WrongOpposite << std::endl;
            std::cout << "Num Corrected Separatrix " << CorrectedSeparatrix << std::endl;
            std::cout << "_______________________ " << std::endl;
            std::cout << "debug separatrices" <<  std::endl;
            std::cout << "Num Unresolved after TJunctions " << UnresolvedDirections1 << std::endl;
        }
        std::cout.flush();

        if (filepath!="")
        {
            out.close();
            std::cout.rdbuf(coutbuf); //reset to standard output again
        }
    }

//    //print the stats values
//    void PrintValues(std::string filepath="",
//                     bool printdebug=true)
//    {
//        std::streambuf *coutbuf;
//        if (filepath!="")
//        {
//            std::ofstream out(filepath.c_str());
//            coutbuf = std::cout.rdbuf(); //save old buf
//            std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
//        }

//        std::cout << "_______________________"<< std::endl;
//        std::cout << "*** FIELD STATS ***" <<  std::endl;
//        std::cout << "Num Singularities " << NumSingularities << std::endl;

//        std::cout << "_______________________"<< std::endl;
//        std::cout << "*** GRAPH STATS ***" <<  std::endl;
//        std::cout << "Num Nodes " << NumNodes << std::endl;
//        std::cout << "Num Edge Nodes " << EdgeNodes << std::endl;
//        std::cout << "Num Boundary Nodes " << BoundaryNodes << std::endl;
//        std::cout << "Num Singularity Nodes " << SingularityNodes << std::endl;
//        std::cout << "Num Sink Nodes " << SinkNodes << std::endl;
//        std::cout << "Time " << timeGraph << std::endl;

//        std::cout << "_______________________ " << std::endl;
//        std::cout << "*** SEPARATRIX STATS ***" <<  std::endl;
//        std::cout << "Num Iterations  " << NumIterations << std::endl;
//        std::cout << "Num Neighbors  " << NumNeighbors << std::endl;
//        std::cout << "Num Resolved  " << SinkNodes-UnresolvedDirections0 << std::endl;
//        std::cout << "Num TJunctions  " << UnresolvedDirections0 << std::endl;
//        std::cout << "Time " << timeSeparatrix << std::endl;

//        std::cout << "_______________________ "  << std::endl;
//        std::cout << "*** PATCH STATS ***" <<  std::endl;
//        std::cout << "Num Patches " << NumPatches << std::endl;
//        for (size_t i=0;i<ValencePatches.size();i++)
//        {
//            if (ValencePatches[i]==0)continue;
//            std::cout << "Patches with Valence " << i << " :  " << ValencePatches[i] << std::endl;
//        }

//        std::cout << "_______________________ "  << std::endl;
//        std::cout << "*** PARAM STATS ***" <<  std::endl;
//        std::cout << "Num Variables " << NumVariables << std::endl;
//        std::cout << "Percentile " << Percentile << std::endl;
//        std::cout << "Subdivision Factor " << subDiv << std::endl;
//        std::cout << "Num Unresolved T " << NumUnresolvedT << std::endl;
//        std::cout << "Time " << timeParam << std::endl;

//        std::cout << "_______________________ "  << std::endl;
//        std::cout << "*** QUAD STATS ***" <<  std::endl;
//        std::cout << "Num Quads " << NumQuads << std::endl;
//        std::cout << "Num Vertices " << NumVert << std::endl;
//        std::cout << "Num Irregular " << NumIrregular << std::endl;

//        std::cout << "_______________________ "  << std::endl;
//        std::cout << "*** SUB STATS ***" <<  std::endl;
//        std::cout << "Num Subdivided Faces " << NumSubPatches << std::endl;
//        std::cout << "_______________________ " << std::endl;
//        std::cout << "_______________________ " << std::endl;
//        std::cout << "_______________________ " << std::endl;
//        std::cout << "_______________________ " << std::endl;

//        if (printdebug)
//        {
//            std::cout << "debug graph" <<  std::endl;
//            std::cout << "Num Wrong Opposites " << WrongOpposite << std::endl;
//            std::cout << "Num Corrected Separatrix " << CorrectedSeparatrix << std::endl;
//            std::cout << "_______________________ " << std::endl;
//            std::cout << "debug separatrices" <<  std::endl;
//            std::cout << "Num Unresolved after TJunctions " << UnresolvedDirections1 << std::endl;
//        }
//        std::cout.flush();

//        if (filepath!="") std::cout.rdbuf(coutbuf); //reset to standard output again
//    }

    StatsCollector()
    {
        ValencePatches.resize(10,0);

        NumSingularities=0;
        NumNodes=0;
        EdgeNodes=0;
        BoundaryNodes=0;
        SingularityNodes=0;
        SinkNodes=0;
        timeGraph=0;
        WrongOpposite=0;
        CorrectedSeparatrix=0;
        NumIterations=0;
        NumNeighbors=0;
        UnresolvedDirections0=0;
        UnresolvedDirections1=0;
        timeSeparatrix=0;

        NumPatches=0;

        NumVariables=0;
        NumUnresolvedT=0;
        subDiv=0;
        Percentile=0;
        timeParam=0;

        NumSubPatches=0;

        NumQuads=0;
        NumVert=0;
        NumIrregular=0;

    }
};
#endif //STATS_COLLECTOR
