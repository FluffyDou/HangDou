/*******************************************************************/
// This is the main function to test Spectral Matching.
// Assumptions:
// (a) We read in a standard PDB file. Data truncation is allowed but
// the header must be complete.
// (b) ChainID starts from 'A' and increase alphabetically

// Hang 10/06/2014
/*******************************************************************/

// This header contains all the typedef, other necessary headers
// and some helper functions
#include "SpectralMatching.h"
//#include <iomanip>
//#include "Eigen\Eigenvalues"

/********* Global variables **************/

// Setter up the timer
double   elapsedTime = 0;
CPUTimer cpuTimer;

// Export module and data file name
char exportModelFileName[100]  = "..\\Release\\ModuleDataHelixAndSequence.csv";
char exportResultFileName[100] = "..\\Release\\MatchResult.csv";

// If we switch the input module and data.
bool switchFlag = false;

PDBObject     pdbReader;
HelixGeometry moduleHelix;
HelixGeometry dataHelix;
// Real input handles to score matrix generation
HelixGeometry *inputModuleHelix = 0;
HelixGeometry *inputDataHelix   = 0;

// Double line direction
ScoreMatrixDoubleLineDir scoreMatrixDLD;
SMSolverDLD smSolverDLD;

// Single line direction
ScoreMatrixSingleLineDir scoreMatrixSLD;
SMSolverSLD smSolverSLD;

// 2D Points
ScoreMatrix2DPointCloud  scoreMatrix2DP;
SMSolver2DP smSolver2DP;

// Bounding box diagonal
double moduleBBD;
double dataBBD;

// For visualization in Mathematica
Vector3f moduleOffset(50,50,50);
Vector3f dataOffset(110,110,50);
// S.t. in Mathematica, 1 unit means 100 angstrom
float    scale = 0.01f;
// Parameters for score matrix construction
std::vector<double> thresholds(7,0);
// Gaussian by default
int scoreFunFlag = 0;

/******************************* Basic Functions **********************************/
void ReadInTheFile(const char* fileName, PDBObject & pdbReader);
// Retrieve geometry for one chain
void RetrieveHelixGeometry(HelixGeometry& helixGeo, PDBObject& pdbReader, char chainID, Vector3f modelOffset);
// Retrieve geometry for specific chains
void RetrieveHelixGeometry(HelixGeometry& helixGeo, PDBObject& pdbReader, Vector3f modelOffset, std::string chainIDs = "All");

 // Make sure that the module set contains more helices
bool SwtichModuleData(HelixGeometry* & inputModuleHelix, HelixGeometry* & inputDataHelix, HelixGeometry* moduleHelix, HelixGeometry* dataHelix);

// Read in new parameter list to compute score matrix
void SetParameterForScoreMatrix(std::vector<double> & thresholds, std::string parameters, int & scoreFunFlag, int sfFlag);

////// For double line direction
void ConstructScoreMatrixDLD(ScoreMatrixDoubleLineDir & scoreMatrixDLD, HelixGeometry *moduleHelix, HelixGeometry *dataHelix, const std::vector<double> & thresholds, int & scoreFunFlag);
// Solve for matching result
void SolveForResultDLD(SMSolverDLD & smSolverDLD, const ScoreMatrixDoubleLineDir& matrix, int methodFlag = 0, double eigenRatio = 1.0);

////// For single line direction
void ConstructScoreMatrixSLD(ScoreMatrixSingleLineDir & scoreMatrixSLD, HelixGeometry *moduleHelix, HelixGeometry *dataHelix, const std::vector<double> & thresholds, int & scoreFunFlag);
// Construct the score matrix for orient single line direction
void ConstructScoreMatrixSLDOrient(ScoreMatrixSingleLineDir & scoreMatrixSLD, HelixGeometry *moduleHelix, HelixGeometry *dataHelix, const std::vector<double> & thresholds, int & scoreFunFlag);
// Solve for matching result
void SolveForResultSLD(SMSolverSLD & smSolverSLD, const ScoreMatrixSingleLineDir& matrix, int methodFlag = 0, double eigenRatio = 1.0);

///// For 2D point cloud
void ConstructScoreMatrix2DP(ScoreMatrix2DPointCloud & scoreMatrix2DP, PointCloud2D* modulePoints, PointCloud2D* dataPoints, const std::vector<double> & thresholds, int & scoreFunFlag);
// Solve for matching result
void SolveForResult2DP(SMSolver2DP & smSolver2DP, const ScoreMatrix2DPointCloud& matrix, int methodFlag = 0);

// Export module and data geometry for Mathematica
void ExportModuleAndDataGeometry(HelixGeometry *inputModuleHelix, HelixGeometry *inputDataHelix, bool switchFlag);
// Export 2D point cloud. ToDo.
void ExportModuleAndDataGeometry(PointCloud2D *module, PointCloud2D *data);
// Export match result for Mathematica
void ExportMatchResult();


/********************************************************** Scripts *****************************************************************/

// Run test for new stop TPS oriented line model
void RunForTPSOrientLineNewStop();

// Gather accuracy, function score, and time for TPS double line direction
void RunForTPSLineNewStop();

// Run test for TPS orient double line direction
void RunForTPSLine();
// Run test for TPS orient single line direction
void RunForTPSOrientLine();

// Test on 2D point cloud
void RunForTest2DPoint();
// Test on TPS 2D point cloud
void RunForTest2DPointTPS();

// Gather accuracy, function score, and time for single line direction
void RunForTestSingleDir();

// Record accuracy, score sum and time for given method with given parameters. (Iterative Pick, Simple Greedy)
void RunForAccuracyAndScoreSum();
// Record accuracy, score sum and time for given method with given parameters. (Iterative Pick, Iterative pick new stop)
void RunForAccuracyAndScoreSumNewStop();

// Run the batch to obtain accuracy volume (angle, midpoint distance and helix length)
void RunForAccuracyVolumeADL();
// Compute the accuracy
float ComputeAccuracyDLD(const std::vector<int> & matchTruth, const std::vector<int> & matchResult, const std::vector<int> & m_matchFlag);
float ComputeAccuracySLD(const std::vector<int> & matchTruth, const std::vector<int> & matchResult);

// Run the batch to obtain accuracy volume (distance weight, midpoint distance and helix length)
// Fix angle. Absolute distance weight
void RunForAccuracyVolumeDWDL();
// Fix angle. BBD ratio distance weight

void RunForAccuracyVolumeDWDLBBD();
// Compute the bounding box diagonal
void ComputeBBD(HelixGeometry* helixGeo, double & diagonal);

// Fix distance weight threshold.
void RunForAccuracyVolumeADLDW();

// Compute the accuracy volume for simple greedy
void RunForAccuracyVolumeSGBADL();

using Eigen::ArrayXd;
using Eigen::SelfAdjointEigenSolver;

//#define TEMP_NUM 10. Flag

// This program runs in a loop. It reads in PDB files and run SM to do helices
// (line segments) matching. It exports following information for Mathematica:
// ModuleHelix, module helix sequence, DataHelix, data helix sequence,
// matching result and matching result's picked order
// So far, it does not export intermediate results of Iterative Picking

// Pass in command: ToDo.
int main( int args, char * argv[] )
{       
/*
    std::string   testFile1 = "../Data/1n0u_A.pdb";
    std::string   testFile2 = "../Data/1n0v_C.pdb";
    std::string   testPara  = "1 170.0 170.0, 5.0 1.308996941 1.308996941 50.0"; // "0.0 50.0 50.0 10.0 2.0943951057 2.0943951057 100.0"  "0.0 50.0 50.0 5.0 0.3 0.3 5.0" 1, 1.0472, 1.0472, 16.3333""

    // Read in module set
    ReadInTheFile(testFile1.c_str(), pdbReader);
    //std::cout<<"chainID offset: "<<pdbReader.GetChainIdOffset()<<"\n";
    // Retrieve Helix model --- update so that we can choose all chains or specific chains. ToDo.
    //RetrieveHelixGeometry(moduleHelix, pdbReader, 'A', moduleOffset);
    RetrieveHelixGeometry(moduleHelix, pdbReader, moduleOffset, "All");
    // Read in data set
    ReadInTheFile(testFile2.c_str(), pdbReader);
    // Retrieve Helix model 
    //RetrieveHelixGeometry(dataHelix, pdbReader, 'A', dataOffset);
    RetrieveHelixGeometry(dataHelix, pdbReader, dataOffset, "All");
    // Make sure that the module set contains more helices --- Not necessary. Remove it. ToDo.
    SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

    // Read in new parameter list to compute score matrix
    SetParameterForScoreMatrix( thresholds, testPara, scoreFunFlag, 0 );    

    // Construct the score matrix
//    ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);
    // Solve for result
//    SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 0);

    ConstructScoreMatrixSLD(scoreMatrixSLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);

    SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 2);

    // Compute the score sum
    smSolverSLD.ComputeFunScore(scoreMatrixSLD);
    double funScore = smSolverSLD.GetFunScore();

    //// Compute indicator vector based on matching truth
    //VectorXd indicatorVec = VectorXd::Zero( scoreMatrixSLD.GetDim() );
    const int mLen = scoreMatrixSLD.GetModuleNum();
    const int dLen = scoreMatrixSLD.GetDataNum();
    //// Generate indicator vector --- input to IPFP
    //for(int i = 0; i < dLen; ++i)
    //{
    //    if(0 !=smSolverSLD.GetMatchResult()[i])
    //        indicatorVec[ mLen*i + smSolverSLD.GetMatchResult()[i] - 1 ] = 1.0;
    //}

    //MatrixXd fullMatrix = scoreMatrixSLD.GetMatrix();
    //fullMatrix .triangularView<Eigen::Upper>() = fullMatrix.transpose();

    VectorXd indicatorVec = smSolverSLD.GetIndicatorVec(scoreMatrixSLD);
    MatrixXd fullMatrix   = smSolverSLD.GetFullMatrix(scoreMatrixSLD);
    IPFP ipfp(fullMatrix, indicatorVec, funScore, 10, mLen, dLen);

    // Time for matrix (graph) construction and solve
    cpuTimer.Start();

    ipfp.DoIPFP();

    elapsedTime = cpuTimer.Tick();
    std::cout<<"Time to do IPFP: "<<elapsedTime*(1e-3)<<"s.\n";

    std::cout<<"Initial score:  "<<funScore<<"\n";
    std::cout<<"Refined score:  "<<ipfp.GetRefinedScore()<<"\n";
    std::cout<<"Refined result:\n"<<ipfp.GetRefinedResult().transpose()<<"\n";

    //std::vector<int> tempMatchTruth(9, 0);
    //tempMatchTruth[0] = 1;
    //tempMatchTruth[1] = 2;
    //tempMatchTruth[2] = 3;
    //tempMatchTruth[3] = 4;
    //tempMatchTruth[4] = 6;
    //tempMatchTruth[5] = 7;
    //tempMatchTruth[6] = 8;
    //tempMatchTruth[7] = 0;
    //tempMatchTruth[8] = 9;

    //// Update the match truth based on different parameter setting. 
    //for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
    //{
    //    if( 0 != tempMatchTruth[h] )
    //    {
    //        // The truth has the same format as in Mathematica --- index starts from 1
    //        float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
    //        float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

    //        if( fabs(moduleLen - datalen) > thresholds[6] && thresholds[6] > 0 )
    //            tempMatchTruth[h] = 0;
    //    }
    //}
   
    //std::cout<<"Matrix Sum: "<<scoreMatrixDLD.GetMatrixSum()<<"\n";
    //std::cout<<"Match Truth:\n";

    //for(auto it = tempMatchTruth.begin(); it != tempMatchTruth.end(); ++it)
    //{
    //    std::cout<<*it<<", ";
    //}
    //std::cout<<"\n";

    //std::cout<<"Parameter list:\n";
    //for(auto it = scoreMatrixSLD.GetParaList().begin(); it != scoreMatrixSLD.GetParaList().end(); ++it)
    //{
    //    std::cout<<*it<<", ";
    //}
    //std::cout<<"\n";

    std::cout<<"Match result:\n";
    for(auto it = smSolverSLD.GetMatchResult().begin(); it != smSolverSLD.GetMatchResult().end(); ++it)
    {
        std::cout<<*it<<", ";
    }
    std::cout<<"\n";

    //std::cout<<"Match flag:\n";
    //for(auto it = smSolverDLD.GetMatchFlag().begin(); it != smSolverDLD.GetMatchFlag().end(); ++it)
    //{
    //    std::cout<<*it<<", ";
    //}
    //std::cout<<"\n";

    //float accuracyValue = ComputeAccuracy(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
    //std::cout<<"accuracy value: "<<accuracyValue<<"\n";

    //MatrixXd tm(3,3);
    //tm << 1,2,3,4,5,6,7,8,9;

    //VectorXd tv(3);
    //tv << 1,2,3;

    //double tr = tv.transpose() * tm * tv;

    //MatrixXd testM(3,3);
    //testM<< 1, 2, 3,
    //4, 5, 6,
    //7, 8, 9;

    //std::cout<<"testM:\n"<<testM<<"\n";

    //testM.triangularView<Eigen::Upper>() = testM.transpose();

    //std::cout<<"testM:\n"<<testM<<"\n";

    //MatrixXd F = MatrixXd::Random(2,2);
    //MatrixXd G = MatrixXd::Random(2,2);
    //MatrixXd H(1,1);

    //int* tempM = new int[TEMP_NUM*TEMP_NUM];

    //int nThreads = 0;
    //int i = 0; int j = 0;
    //int offset = TEMP_NUM;
    ////#pragma omp parallel for
    //#pragma omp parallel default(shared) private(i,j)
    //{
    //#pragma omp master
    //nThreads = omp_get_num_threads();

    //    // Enter the loop to compute the accuracy volume --- #pragma omp parallel for
    //#pragma omp for
    //for (i = 0; i < TEMP_NUM; ++i)
    //    for(j = 0; j < TEMP_NUM; ++j)
    //    {
    //        // Construct the score matrix
    //        //ConstructScoreMatrix(scoreMatrixDLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);

    //        //SolveForResult(smSolverDLD, scoreMatrixDLD);
    //        //std::cout <<"i,j: "<<i<<", "<<j<< std::endl;
    //        tempM[offset*i+j] = 10*i+j;
    //        //H.resize(i+1,i+1);
    //        //for(int h = 0; h < i; ++h)
    //        //    for(int k = 0; k < i; ++k)
    //        //        H(h,k) = (double)(h+k);

    //    }

    //}

    //for (i = 0; i < TEMP_NUM; ++i)
    //{
    //    for(j = 0; j < TEMP_NUM; ++j)
    //    {
    //        std::cout <<tempM[offset*i+j]<<", ";
    //    }

    //    std::cout<<"\n";
    //}
    //std::cout<<nThreads<<" OpenMP threads were used.\n";
    //std::cout << "Done!" << std::endl;

    //delete[] tempM;
    ////Eigen::initParallel();


    std::cout<<"Try Point Cloud.\n";
    
    // Create Module
    PointCloud2D modulePC(15);

    // Create Data
    PointCloud2D dataPC;
    dataPC.SetPointCloud2D(modulePC);

    dataPC.AddWhiteNoise(0.0f, 2.0f);

    dataPC.AddOutLiersUniform(5);

    dataPC.TransForm(Vector2f(10.5, 20.0), 0.25*PI);
    //dataPC.TransForm(Vector2f(10.5, 20.0), 0.0);

    // Add outliers for Module
    modulePC.AddOutLiersUniform(5);

    std::cout<<"Module:\n";

    for(auto it = modulePC.GetPointCloud().begin(); it != modulePC.GetPointCloud().end(); ++it)
    {
        std::cout<<"("<<it->transpose()<<"), ";
    }
    std::cout<<"\n";

    std::cout<<"Data:\n";

    for(auto it = dataPC.GetPointCloud().begin(); it != dataPC.GetPointCloud().end(); ++it)
    {
        std::cout<<"("<<it->transpose()<<"), ";
    }
    std::cout<<"\n";

    std::string   testPara  = "0.0 0.0 0.0 5.0 0.0 0.0 0.0"; 
    // Read in new parameter list to compute score matrix --- Quadratic (2) for 2D point cloud
    SetParameterForScoreMatrix( thresholds, testPara, scoreFunFlag, 2 );    

    // Construct the score matrix
    ConstructScoreMatrix2DP(scoreMatrix2DP, &modulePC, &dataPC, thresholds, scoreFunFlag);

    std::cout<<"Point number test: "<<scoreMatrix2DP.GetModuleNum()<<","<<scoreMatrix2DP.GetDataNum()<<"\n";

    SolveForResult2DP(smSolver2DP, scoreMatrix2DP, 2);

    std::cout<<"Match result:\n";
    for(auto it = smSolver2DP.GetMatchResult().begin(); it != smSolver2DP.GetMatchResult().end(); ++it)
    {
        std::cout<<*it<<", ";
    }
    std::cout<<"\n";
*/

    //MatrixXd testM(3,3);
    //testM<< 1, 2, 3,
    //        4, 5, 6,
    //        7, 8, 9;

    //// Initialize the matlab library
    //if( !LeadingEigenvecInitialize() )
    //{
    //    std::cerr << "Could not initialize the library properly."<< std::endl;
    //    return -1;
    //}

    //// Input matrix
    //mwArray matrix_mat(3,3, mxDOUBLE_CLASS, mxREAL);

    //// Output
    //mwArray out;

    //// Pass in data to matlab    
    //matrix_mat.SetData(testM.data(), testM.cols()*testM.cols());

    //// Matlab kernel computation --- 1 means one output
    //LeadingEigenVec(1, out, matrix_mat);

    //VectorXd prunedVec(testM.cols());
    ////double* prunedVec = new double[];
    //out.GetData(prunedVec.data(), testM.cols());

    //std::cout<<"result:\n"<<prunedVec<<"\n";

    //LeadingEigenvecTerminate();

    std::cout<<"Begin process. Press q or Q to exit.\n";

    // Take all the external input commands (from Mathematica)
    //std::string   testFile = "importModuleAndData, ../Data/1oel_A.pdb, A, ../Data/2c7c_A.pdb, A";
    //std::string   testPara = "0.0 50.0 50.0 10.0 2.0943951057 2.0943951057 100.0"; // "0.0 50.0 50.0 10.0 2.0943951057 2.0943951057 100.0"  "0.0 50.0 50.0 5.0 0.3 0.3 5.0"
    // constructAndSolve, 0.0 50.0 50.0 10.0 2.0943951057 2.0943951057 100.0
    // constructAndSolve, 0.0 50.0 50.0 5.0 0.3 0.3 5.0
    // constructAndSolve, 0. 0. 0. 10.0 1.91986 1.91986 100.0

    std::string   inputString = "";

    // Main loop of the program
    while(true)
    {
        // Wait to read in the next command line
        std::getline( std::cin, inputString );

        // Parse the current line by comma
        //std::istringstream buf(inputString);
        //std::istream_iterator<std::string> beg(buf), end;
        //std::vector<std::string> tokens(beg, end);

        // Protect from empty input
        if("" == inputString)
            continue;

        // Parse the input command by ','
        std::vector<std::string> tokens;
        ParseByCharater(inputString, tokens, ',');

        // Exit when capturing q or Q
        if( "q" ==  tokens[0] || "Q" ==  tokens[0] )
            break;

        // Run for accuracy volume (angle, midpoint distance, helix length)
        if( "RunBatch" == tokens[0] )
        {
            // Run for adl volume through Iterative Picking
//            RunForAccuracyVolumeADL();

            // Run for adl volume through Simple Greedy
//            RunForAccuracyVolumeSGBADL();

            // Run for dwdl volume through Iterative Picking
//            RunForAccuracyVolumeDWDL();

            // Fix distance weight
//            RunForAccuracyVolumeADLDW();

            // Record accuracy, score sum and time for double line dir
//            RunForAccuracyAndScoreSum();

            // Record accuracy, function score and time for single line dir
//            RunForTestSingleDir();

            // Run test for 2D point cloud
//            RunForTest2DPoint();

            // Run test for TPS 2D point cloud
//            RunForTest2DPointTPS();

            // Run test for TPS orient single line direction
//            RunForTPSOrientLine();

            // Run test for new stop real protein data
//            RunForAccuracyAndScoreSumNewStop();

            // Run test for new stop TPS model
            RunForTPSLineNewStop();

            // Run test for new stop TPS oriented line model
//            RunForTPSOrientLineNewStop();


            break;
        }

        // Read in module and data sets
        if( "importModuleAndData" == tokens[0] )
        {
            // Read in module set
            ReadInTheFile(tokens[1].c_str(), pdbReader);
            // Retrieve Helix model
            RetrieveHelixGeometry(moduleHelix, pdbReader, moduleOffset, trimString(tokens[2]));
            //RetrieveHelixGeometry(moduleHelix, pdbReader, 'A', moduleOffset);
            // Read in data set
            ReadInTheFile(tokens[3].c_str(), pdbReader);
            // Retrieve Helix model 
            RetrieveHelixGeometry(dataHelix, pdbReader, dataOffset, trimString(tokens[4]));
            //RetrieveHelixGeometry(dataHelix, pdbReader, 'A', dataOffset);
            // Make sure that the module set contains more helices --- Not necessary. Remove it. ToDo.
            switchFlag = SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

            std::cout<<"Finish import models.\n";
        }

        // Export module and data geometry for Mathematica (include: transform the models to world space before dump the geometry).
        if( "exportModuleAndData" == tokens[0] )
        {
            // Export the model geometry
            ExportModuleAndDataGeometry(inputModuleHelix, inputDataHelix, switchFlag);

            // Inform the console we have exported model geometry
            std::cout<<"Finish export models.\n";
        }

        // Construct and solve the score matrix
        if( "constructAndSolve" == tokens[0] )
        {
            // Read in new parameter list to compute score matrix
            SetParameterForScoreMatrix( thresholds, tokens[1], scoreFunFlag, 0 );

            // Construct the score matrix
            ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);

            // Iterative Picking
            SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 0);

            // Export the result (including line direction) and picked order for Mathematica
            ExportMatchResult();

            // Inform the console we have exported matching result
            std::cout<<"Finish export result.\n";
        }

    } // end while

    std::cout<<"End process.\n";

    return 0;

}

/**************************************************************************************************************************/
/**************************************************** Scripts *************************************************************/
/**************************************************************************************************************************/

void RunForTPSOrientLineNewStop()
{
    /******************************** Read in all the data ********************************/
    std::cout<<"Begin to test on TPS Orient Line...\n";

    std::string inputListFile  = "../Data/TPSLine/fileList.txt";
    std::string dataDir        = "../Data/TPSLine/3tgl/";

    // Source model
    std::vector< std::string > inputFileList;
    // Deformed models by bending energy (each energy level has 50 bended models)
    std::vector< std::string > outputFileList;

    // Only one input source
    inputFileList.push_back("inputLineSegs.txt");

    std::ifstream file(inputListFile);
    std::string   line;

    // Read the deformed model list by energy level
    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputListFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            // Attention: no protection for space character
            outputFileList.push_back(line);
        }

        std::cout<<"Finish reading in all the file names!\n";
    }

    file.close();

    /////////  Process all the files /////////

    // Read in input line source
    std::vector< HelixGeometry > inputSetList;
    std::vector<lineSeg> lineSegList;

    std::string inputFile = dataDir + inputFileList[0];
    file.open(inputFile, std::ifstream::in);

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, '\t');

            // Each line has 6 coordinates (2 points)
            if("e" != tokens[0])
            {
                // Read in the two points
                lineSeg ls = {Vector3f(std::stof(tokens[0]),std::stof(tokens[1]),std::stof(tokens[2])), Vector3f(std::stof(tokens[3]),std::stof(tokens[4]),std::stof(tokens[5]))};
                lineSegList.push_back( ls );
            }

        }

        std::cout<<"Finish reading in the inputPoints file!\n";
    }

    file.close();

    // Build the helix geometry
    HelixGeometry tempGeo(lineSegList.size());
    //tempGeo.Reset(lineSegList.size());

    int dummyIndex = 0;
    for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
    {
        tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
    }

    //// Well, not really necessary.
    //tempGeo.ComputeModelCenter();
    ////tempGeo.SetModelOffset(modelOffset);
    //tempGeo.TransformModel();

    inputSetList.push_back(tempGeo);
    lineSegList.clear();

    // Process all the samples by energy level
    unsigned int energyIndex = 0;
    while(energyIndex < outputFileList.size())
    {
        // Read in the data and model from the file
        // ATTENTION!! "inputPoints" are source and "outputPoints" are deformed points
        std::string outputFile = dataDir+outputFileList[energyIndex];

        // The list holding all the deformed models
        std::vector< HelixGeometry > outputSetList;

        ////////////////////// Read output points //////////////////////
        file.open(outputFile, std::ifstream::in);

        if( !file.is_open() )
            std::cout<<"Unable to open "<<outputFile<<"...\n";
        else
        {
            std::cout<<"Read in the file...\n";

            bool initFlag = true;

            std::vector<lineSeg> lineSegList;

            while( std::getline(file, line) )
            {
                // Protect from empty lines
                if( "" == line )
                    continue;

                std::vector<std::string> tokens;
                ParseByCharater(line, tokens, '\t');

                if("e" == tokens[0])
                {
                    // We create new input and outputs
                    if( !initFlag )
                    {
                        // Build the helix geometry
                        HelixGeometry tempGeo;
                        tempGeo.Reset(lineSegList.size());

                        int dummyIndex = 0;
                        for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
                        {
                            tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
                        }

                        outputSetList.push_back(tempGeo);
                        lineSegList.clear();
                    }
                    // Avoid the first e
                    initFlag = false;
                }
                else
                {
                    // Read in the two points
                    lineSeg ls = {Vector3f(std::stof(tokens[0]),std::stof(tokens[1]),std::stof(tokens[2])), Vector3f(std::stof(tokens[3]),std::stof(tokens[4]),std::stof(tokens[5]))};
                    lineSegList.push_back( ls );
                }

            }

            // Record the last set in the file
            HelixGeometry tempGeo;
            tempGeo.Reset(lineSegList.size());

            int dummyIndex = 0;
            for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
            {
                tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
                //std::cout<<(*it)[0].transpose()<<", "<<(*it)[1].transpose()<<"\n";
                //std::cout<<tempGeo.GetLineSegment(dummyIndex).GetP0().transpose()<<", "<<tempGeo.GetLineSegment(dummyIndex).GetP1().transpose()<<"\n\n";
            }

            outputSetList.push_back(tempGeo);
            lineSegList.clear();

            std::cout<<"Finish reading in the deformed models file!\n";
        }

        file.close();  

        std::cout<<"Processing "<<energyIndex+1<<"th bending energy.\n";

        // Parameter settings
        double distanceWeight = 120.0;
        scoreFunFlag  = 0;    
        thresholds[0] = 1.0; // wf
        thresholds[1] = distanceWeight;
        thresholds[2] = distanceWeight;
        thresholds[3] = 5.0;   // td
        thresholds[4] = 90.0*PI / 180.0;       // ta
        thresholds[5] = 90.0*PI / 180.0;       // ta
        thresholds[6] = 40.0; // tl

        // Initialize data we want to record
        float accuracy_IPB = 0.0f;
        float accuracy_SGB = 0.0f;
        float accuracy_HUB = 0.0f;

        float accuracy_IPB_IPFP = 0.0f;
        float accuracy_SGB_IPFP = 0.0f;
        float accuracy_HUB_IPFP = 0.0f;

        double scoreRatio_IPBOverSGB = 0.0;
        double scoreRatio_HUBOverSGB = 0.0;

        double scoreRatio_IPB_IPFP_Over_SGB = 0.0;
        double scoreRatio_HUB_IPFP_Over_SGB = 0.0;
        double scoreRatio_SGB_IPFP_Over_SGB = 0.0;

        double scoreRatio_IPBOverSGB_IPFP = 0.0f;
        double scoreRatio_HUBOverSGB_IPFP = 0.0f;

        int index = 0; // We have 50 samples for each energy level.
        while( index < (int)outputSetList.size() )
        {
            //// Bind module and data
            SwtichModuleData(inputModuleHelix, inputDataHelix, &(inputSetList[0]), &(outputSetList[index]));

            //// Create match result
            int lineSegNum = 10;
            std::vector<int> tempMatchTruth(lineSegNum);

            // Fix match truth (start from 1, not 0)
            for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
            {
                tempMatchTruth[h] = h+1;
            }

            // Update the match truth based on different parameter setting. 
            for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
            {
                if( 0 != tempMatchTruth[h] )
                {
                    // The truth has the same format as in Mathematica --- index starts from 1
                    float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                    float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                    if( fabs(moduleLen - datalen) > thresholds[6] && thresholds[6] > 0 )
                        tempMatchTruth[h] = 0;
                }
            }

            //// Construct the score matrix
            ConstructScoreMatrixSLDOrient(scoreMatrixSLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);
            MatrixXd fullMatrix = smSolverSLD.GetFullMatrix(scoreMatrixSLD);

            VectorXd indicatorVec;
            IPFP ipfp;

            VectorXi tempResultVecXd;
            std::vector<int> tempResult;

            const int mLen = scoreMatrixSLD.GetModuleNum();
            const int dLen = scoreMatrixSLD.GetDataNum();

            ///***************************** Solve by IPB New Stop *****************************/
            double eigenRatio = 0.0;
            SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 3, eigenRatio);

            //// Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_IPB += ComputeAccuracySLD( tempMatchTruth, smSolverSLD.GetMatchResult() );

            //// Retrieve function score
            smSolverSLD.ComputeFunScore(scoreMatrixSLD);
            double score_IPB = smSolverSLD.GetFunScore();

            //// Refine through IPFP
            //// Prepare indicator vector and score matrix
            indicatorVec = smSolverSLD.GetIndicatorVec(scoreMatrixSLD);

            //// Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolverSLD.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            //// Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_IPB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_IPB_IPFP = ipfp.GetRefinedScore();

            ///***************************** Solve by SGB *****************************/
            SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 1);

            //// Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_SGB += ComputeAccuracySLD( tempMatchTruth, smSolverSLD.GetMatchResult() );

            //// Retrieve function score
            smSolverSLD.ComputeFunScore(scoreMatrixSLD);
            double score_SGB = smSolverSLD.GetFunScore();

            //// Refine through IPFP
            //// Prepare indicator vector and score matrix
            indicatorVec = smSolverSLD.GetIndicatorVec(scoreMatrixSLD);

            //// Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolverSLD.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            //// Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_SGB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_SGB_IPFP = ipfp.GetRefinedScore();

            ///***************************** Solve by HUB *****************************/
            SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 2);

            //// Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_HUB += ComputeAccuracySLD( tempMatchTruth, smSolverSLD.GetMatchResult() );

            //// Retrieve function score
            smSolverSLD.ComputeFunScore(scoreMatrixSLD);
            double score_HUB = smSolverSLD.GetFunScore();

            //// Refine through IPFP
            //// Prepare indicator vector and score matrix
            indicatorVec = smSolverSLD.GetIndicatorVec(scoreMatrixSLD);

            //// Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolverSLD.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            //// Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_HUB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_HUB_IPFP = ipfp.GetRefinedScore();

            ///************ Update information **************/
            scoreRatio_IPBOverSGB += score_IPB/score_SGB;
            scoreRatio_HUBOverSGB += score_HUB/score_SGB;

            scoreRatio_IPBOverSGB_IPFP += score_IPB_IPFP/score_SGB_IPFP;
            scoreRatio_HUBOverSGB_IPFP += score_HUB_IPFP/score_SGB_IPFP;

            scoreRatio_IPB_IPFP_Over_SGB += score_IPB_IPFP/score_SGB;
            scoreRatio_HUB_IPFP_Over_SGB += score_HUB_IPFP/score_SGB;
            scoreRatio_SGB_IPFP_Over_SGB += score_SGB_IPFP/score_SGB;

            //// Move to next round
            ++index;
            ////std::cout<<"Done with "<<index<<"th set.\n";
        }
        //std::cout<<"We have "<<index<<" sets for this energy level.\n";

        // Compute average
        scoreRatio_IPBOverSGB /= (double)outputSetList.size();
        scoreRatio_HUBOverSGB /= (double)outputSetList.size();

        scoreRatio_IPBOverSGB_IPFP /= (double)outputSetList.size();
        scoreRatio_HUBOverSGB_IPFP /= (double)outputSetList.size();

        accuracy_IPB /= (float)outputSetList.size();
        accuracy_SGB /= (float)outputSetList.size();
        accuracy_HUB /= (float)outputSetList.size();

        accuracy_IPB_IPFP /= (float)outputSetList.size();
        accuracy_SGB_IPFP /= (float)outputSetList.size();
        accuracy_HUB_IPFP /= (float)outputSetList.size();

        scoreRatio_IPB_IPFP_Over_SGB /= (double)outputSetList.size();
        scoreRatio_HUB_IPFP_Over_SGB /= (double)outputSetList.size();
        scoreRatio_SGB_IPFP_Over_SGB /= (double)outputSetList.size();

        std::cout<<"Function Score:\n";
        std::cout<<"IPB over SGB:     "<<scoreRatio_IPBOverSGB<<"\n";
        std::cout<<"HUB over SGB:     "<<scoreRatio_HUBOverSGB<<"\n";
        std::cout<<"IPB IPFP over SGB:"<<scoreRatio_IPB_IPFP_Over_SGB<<"\n";
        std::cout<<"HUB IPFP over SGB:"<<scoreRatio_HUB_IPFP_Over_SGB<<"\n";
        std::cout<<"SGB IPFP over SGB:"<<scoreRatio_SGB_IPFP_Over_SGB<<"\n";

        std::cout<<"IPB over SGB IPFP:"<<scoreRatio_IPBOverSGB_IPFP<<"\n";
        std::cout<<"HUB over SGB IPFP:"<<scoreRatio_HUBOverSGB_IPFP<<"\n";

        std::cout<<"IPB average accuracy:      "<<accuracy_IPB<<"\n";
        std::cout<<"IPB average accuracy IPFP: "<<accuracy_IPB_IPFP<<"\n";
        std::cout<<"HUB average accuracy:      "<<accuracy_HUB<<"\n";
        std::cout<<"HUB average accuracy IPFP: "<<accuracy_HUB_IPFP<<"\n";
        std::cout<<"SGB average accuracy:      "<<accuracy_SGB<<"\n";
        std::cout<<"SGB average accuracy IPFP: "<<accuracy_SGB_IPFP<<"\n";

        std::cout<<"Finish testing this batch!\n";

        // Move to the next pair of files
        ++energyIndex;
    }
}

// Gather accuracy, function score, and time for TPS double line direction
void RunForTPSLineNewStop()
{
    /******************************** Read in all the data ********************************/
    std::cout<<"Begin to test on TPS Line New Stop...\n";

    std::string inputListFile  = "../Data/TPSLine/fileList.txt";
    std::string dataDir        = "../Data/TPSLine/3tgl/";

    // Source model
    std::vector< std::string > inputFileList;
    // Deformed models by bending energy (each energy level has 50 bended models)
    std::vector< std::string > outputFileList;

    // Only one input source
    inputFileList.push_back("inputLineSegs.txt");

    std::ifstream file(inputListFile);
    std::string   line;

    // Read the deformed model list by energy level
    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputListFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            // Attention: no protection for space character
            outputFileList.push_back(line);
        }

        std::cout<<"Finish reading in all the file names!\n";
    }

    file.close();

    /////////  Process all the files /////////

    // Read in input line source
    std::vector< HelixGeometry > inputSetList;
    std::vector<lineSeg> lineSegList;

    std::string inputFile = dataDir + inputFileList[0];
    file.open(inputFile, std::ifstream::in);

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, '\t');

            // Each line has 6 coordinates (2 points)
            if("e" != tokens[0])
            {
                // Read in the two points
                lineSeg ls = {Vector3f(std::stof(tokens[0]),std::stof(tokens[1]),std::stof(tokens[2])), Vector3f(std::stof(tokens[3]),std::stof(tokens[4]),std::stof(tokens[5]))};
                lineSegList.push_back( ls );
            }

        }

        std::cout<<"Finish reading in the inputPoints file!\n";
    }

    file.close();

    // Build the helix geometry
    HelixGeometry tempGeo(lineSegList.size());
    //tempGeo.Reset(lineSegList.size());

    int dummyIndex = 0;
    for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
    {
        tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
    }

    //// Well, not really necessary.
    //tempGeo.ComputeModelCenter();
    ////tempGeo.SetModelOffset(modelOffset);
    //tempGeo.TransformModel();

    inputSetList.push_back(tempGeo);
    lineSegList.clear();

    // Process all the samples by energy level
    unsigned int energyIndex = 0;
    while(energyIndex < outputFileList.size())
    {
        // Read in the data and model from the file
        // ATTENTION!! "inputPoints" are source and "outputPoints" are deformed points
        std::string outputFile = dataDir+outputFileList[energyIndex];

        // The list holding all the deformed models
        std::vector< HelixGeometry > outputSetList;

        ////////////////////// Read output points //////////////////////
        file.open(outputFile, std::ifstream::in);

        if( !file.is_open() )
            std::cout<<"Unable to open "<<outputFile<<"...\n";
        else
        {
            std::cout<<"Read in the file...\n";

            bool initFlag = true;

            std::vector<lineSeg> lineSegList;

            while( std::getline(file, line) )
            {
                // Protect from empty lines
                if( "" == line )
                    continue;

                std::vector<std::string> tokens;
                ParseByCharater(line, tokens, '\t');

                if("e" == tokens[0])
                {
                    // We create new input and outputs
                    if( !initFlag )
                    {
                        // Build the helix geometry
                        HelixGeometry tempGeo;
                        tempGeo.Reset(lineSegList.size());

                        int dummyIndex = 0;
                        for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
                        {
                            tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
                        }

                        outputSetList.push_back(tempGeo);
                        lineSegList.clear();
                    }
                    // Avoid the first e
                    initFlag = false;
                }
                else
                {
                    // Read in the two points
                    lineSeg ls = {Vector3f(std::stof(tokens[0]),std::stof(tokens[1]),std::stof(tokens[2])), Vector3f(std::stof(tokens[3]),std::stof(tokens[4]),std::stof(tokens[5]))};
                    lineSegList.push_back( ls );
                }

            }

            // Record the last set in the file
            HelixGeometry tempGeo;
            tempGeo.Reset(lineSegList.size());

            int dummyIndex = 0;
            for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
            {
                tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
                //std::cout<<(*it)[0].transpose()<<", "<<(*it)[1].transpose()<<"\n";
                //std::cout<<tempGeo.GetLineSegment(dummyIndex).GetP0().transpose()<<", "<<tempGeo.GetLineSegment(dummyIndex).GetP1().transpose()<<"\n\n";
            }

            outputSetList.push_back(tempGeo);
            lineSegList.clear();

            std::cout<<"Finish reading in the deformed models file!\n";
        }

        file.close();  

        std::cout<<"Processing "<<energyIndex+1<<"th bending energy.\n";

        // Parameter settings
        double distanceWeight = 120.0;
        scoreFunFlag  = 0;    
        thresholds[0] = 1.0; // wf
        thresholds[1] = distanceWeight;
        thresholds[2] = distanceWeight;
        thresholds[3] = 5.0;   // td
        thresholds[4] = 90.0*PI / 180.0;       // ta
        thresholds[5] = 90.0*PI / 180.0;       // ta
        thresholds[6] = 40.0; // tl

        // Initialize data we want to record
        float accuracy_IPB = 0.0f;
        float accuracy_SGB = 0.0f;

        double scoreRatio_IPBOverSGB = 0.0;

        int index = 0; // We have 50 samples for each energy level.
        while( index < (int)outputSetList.size() )
        {
            //// Bind module and data
            SwtichModuleData(inputModuleHelix, inputDataHelix, &(inputSetList[0]), &(outputSetList[index]));

            //// Create match result
            int lineSegNum = 10;
            std::vector<int> tempMatchTruth(lineSegNum);

            // Fix match truth (start from 1, not 0)
            for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
            {
                tempMatchTruth[h] = h+1;
            }

            // Update the match truth based on different parameter setting. 
            for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
            {
                if( 0 != tempMatchTruth[h] )
                {
                    // The truth has the same format as in Mathematica --- index starts from 1
                    float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                    float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                    if( fabs(moduleLen - datalen) > thresholds[6] && thresholds[6] > 0 )
                        tempMatchTruth[h] = 0;
                }
            }

            //// Construct the score matrix
            ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);
            MatrixXd fullMatrix = smSolverDLD.GetFullMatrix(scoreMatrixDLD); 

            const int mLen = scoreMatrixSLD.GetModuleNum();
            const int dLen = scoreMatrixSLD.GetDataNum();

            ///***************************** Solve by IPB New Stop *****************************/
            double eigenRatio = 0.3;
            SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 2, eigenRatio);

            //// Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_IPB += ComputeAccuracyDLD( tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag() );

            //// Retrieve function score
            smSolverDLD.ComputeFunScore(scoreMatrixDLD);
            double score_IPB = smSolverDLD.GetFunScore();

            ///***************************** Solve by SGB *****************************/
            SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 1);

            //// Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_SGB += ComputeAccuracyDLD( tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag() );

            //// Retrieve function score
            smSolverDLD.ComputeFunScore(scoreMatrixDLD);
            double score_SGB = smSolverDLD.GetFunScore();

            ///************ Update information **************/
            scoreRatio_IPBOverSGB += score_IPB/score_SGB;

            //// Move to next round
            ++index;
            ////std::cout<<"Done with "<<index<<"th set.\n";
        }
        //std::cout<<"We have "<<index<<" sets for this energy level.\n";

        // Compute average
        scoreRatio_IPBOverSGB /= (double)outputSetList.size();

        accuracy_IPB /= (float)outputSetList.size();
        accuracy_SGB /= (float)outputSetList.size();

        std::cout<<"Function Score:\n";
        std::cout<<"IPB over SGB:         "<<scoreRatio_IPBOverSGB<<"\n";
        std::cout<<"IPB average accuracy: "<<accuracy_IPB<<"\n";
        std::cout<<"SGB average accuracy: "<<accuracy_SGB<<"\n";

        std::cout<<"Finish testing this batch!\n";

        // Move to the next pair of files
        ++energyIndex;
    }
}

// Gather accuracy, function score, and time for TPS single line direction
void RunForTPSOrientLine()
{
    // Wrap all this information into a class. ToDo. Flag Ha.
    std::cout<<"Begin to test on TPS Orient Line...\n";

    std::string inputListFile  = "../Data/TPSLine/fileList.txt";
    std::string dataDir        = "../Data/TPSLine/3tgl/";

    // Source model
    std::vector< std::string > inputFileList;
    // Deformed models by bending energy (each energy level has 50 bended models)
    std::vector< std::string > outputFileList;

    // Only one input source
    inputFileList.push_back("inputLineSegs.txt");

    std::ifstream file(inputListFile);
    std::string   line;

    // Read the deformed model list by energy level
    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputListFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            // Attention: no protection for space character
            outputFileList.push_back(line);
        }

        std::cout<<"Finish reading in all the file names!\n";
    }

    file.close();

    /////////  Process all the files /////////

    // Read in input line source
    std::vector< HelixGeometry > inputSetList;
    std::vector<lineSeg> lineSegList;

    std::string inputFile = dataDir + inputFileList[0];
    file.open(inputFile, std::ifstream::in);

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, '\t');

            // Each line has 6 coordinates (2 points)
            if("e" != tokens[0])
            {
                // Read in the two points
                lineSeg ls = {Vector3f(std::stof(tokens[0]),std::stof(tokens[1]),std::stof(tokens[2])), Vector3f(std::stof(tokens[3]),std::stof(tokens[4]),std::stof(tokens[5]))};
                lineSegList.push_back( ls );
            }

        }

        std::cout<<"Finish reading in the inputPoints file!\n";
    }

    file.close();

    // Build the helix geometry
    HelixGeometry tempGeo(lineSegList.size());
    //tempGeo.Reset(lineSegList.size());

    int dummyIndex = 0;
    for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
    {
        tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
    }

    //// Well, not really necessary.
    //tempGeo.ComputeModelCenter();
    //tempGeo.SetModelOffset(modelOffset);
    //tempGeo.TransformModel();

    inputSetList.push_back(tempGeo);
    lineSegList.clear();

    // Process all the samples by energy level
    unsigned int energyIndex = 0;
    while(energyIndex < outputFileList.size())
    {
        // Read in the data and model from the file
        // ATTENTION!! "inputPoints" are source and "outputPoints" are deformed points
        std::string outputFile = dataDir+outputFileList[energyIndex];

        // The list holding all the deformed models
        std::vector< HelixGeometry > outputSetList;
        
        ////////////////////// Read output points //////////////////////
        file.open(outputFile, std::ifstream::in);

        if( !file.is_open() )
            std::cout<<"Unable to open "<<outputFile<<"...\n";
        else
        {
            std::cout<<"Read in the file...\n";

            bool initFlag = true;

            std::vector<lineSeg> lineSegList;

            while( std::getline(file, line) )
            {
                // Protect from empty lines
                if( "" == line )
                    continue;

                std::vector<std::string> tokens;
                ParseByCharater(line, tokens, '\t');

                if("e" == tokens[0])
                {
                    // We create new input and outputs
                    if( !initFlag )
                    {
                        // Build the helix geometry
                        HelixGeometry tempGeo;
                        tempGeo.Reset(lineSegList.size());

                        int dummyIndex = 0;
                        for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
                        {
                            tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
                        }

                        outputSetList.push_back(tempGeo);
                        lineSegList.clear();
                    }
                    // Avoid the first e
                    initFlag = false;
                }
                else
                {
                    // Read in the two points
                    lineSeg ls = {Vector3f(std::stof(tokens[0]),std::stof(tokens[1]),std::stof(tokens[2])), Vector3f(std::stof(tokens[3]),std::stof(tokens[4]),std::stof(tokens[5]))};
                    lineSegList.push_back( ls );
                }

            }

            // Record the last set in the file
            HelixGeometry tempGeo;
            tempGeo.Reset(lineSegList.size());

            int dummyIndex = 0;
            for(auto it = lineSegList.begin(); it != lineSegList.end(); ++it, ++dummyIndex)
            {
                tempGeo.AddLineSeg( dummyIndex, LineSegment(it->at(0), it->at(1)) );
                //std::cout<<(*it)[0].transpose()<<", "<<(*it)[1].transpose()<<"\n";
                //std::cout<<tempGeo.GetLineSegment(dummyIndex).GetP0().transpose()<<", "<<tempGeo.GetLineSegment(dummyIndex).GetP1().transpose()<<"\n\n";
            }

            outputSetList.push_back(tempGeo);
            lineSegList.clear();

            std::cout<<"Finish reading in the deformed models file!\n";
        }

        file.close();  

        std::cout<<"Processing "<<energyIndex+1<<"th bending energy.\n";

        // Parameter settings
        double distanceWeight = 120.0;
        scoreFunFlag  = 0;    
        thresholds[0] = 1.0; // wf
        thresholds[1] = distanceWeight;
        thresholds[2] = distanceWeight;
        thresholds[3] = 5.0;   // td
        thresholds[4] = 90.0*PI / 180.0;       // ta
        thresholds[5] = 90.0*PI / 180.0;       // ta
        thresholds[6] = 40.0; // tl

        // Initialize data we want to record
        float accuracy_IPB = 0.0f;
        float accuracy_SGB = 0.0f;
        float accuracy_HUB = 0.0f;

        float accuracy_IPB_IPFP = 0.0f;
        float accuracy_SGB_IPFP = 0.0f;
        float accuracy_HUB_IPFP = 0.0f;

        double scoreRatio_IPBOverSGB = 0.0;
        double scoreRatio_HUBOverSGB = 0.0;

        double scoreRatio_IPB_IPFP_Over_SGB = 0.0;
        double scoreRatio_HUB_IPFP_Over_SGB = 0.0;
        double scoreRatio_SGB_IPFP_Over_SGB = 0.0;

        double scoreRatio_IPBOverSGB_IPFP = 0.0f;
        double scoreRatio_HUBOverSGB_IPFP = 0.0f;

        int index = 0; // We have 50 samples for each energy level.
        while( index < (int)outputSetList.size() )
        {
            //// Bind module and data
            SwtichModuleData(inputModuleHelix, inputDataHelix, &(inputSetList[0]), &(outputSetList[index]));

            //// Create match result
            int lineSegNum = 10;
            std::vector<int> tempMatchTruth(lineSegNum);

            // Fix match truth (start from 1, not 0)
            for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
            {
                tempMatchTruth[h] = h+1;
            }

            // Update the match truth based on different parameter (threshold) setting. 
            for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
            {
                if( 0 != tempMatchTruth[h] )
                {
                    // The truth has the same format as in Mathematica --- index starts from 1
                    float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                    float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                    if( fabs(moduleLen - datalen) > thresholds[6] && thresholds[6] > 0 )
                        tempMatchTruth[h] = 0;
                }
            }

            //// Construct the score matrix
            ConstructScoreMatrixSLDOrient(scoreMatrixSLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);
            MatrixXd fullMatrix = smSolverSLD.GetFullMatrix(scoreMatrixSLD);

            VectorXd indicatorVec;
            IPFP ipfp;

            VectorXi tempResultVecXd;
            std::vector<int> tempResult;

            const int mLen = scoreMatrixSLD.GetModuleNum();
            const int dLen = scoreMatrixSLD.GetDataNum();

            ///***************************** Solve by IPB *****************************/
            SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 0, 1.0);

            //// Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_IPB += ComputeAccuracySLD( tempMatchTruth, smSolverSLD.GetMatchResult() );

            //// Retrieve function score
            smSolverSLD.ComputeFunScore(scoreMatrixSLD);
            double score_IPB = smSolverSLD.GetFunScore();

            //// Refine through IPFP
            //// Prepare indicator vector and score matrix
            indicatorVec = smSolverSLD.GetIndicatorVec(scoreMatrixSLD);

            //// Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolverSLD.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            //// Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_IPB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_IPB_IPFP = ipfp.GetRefinedScore();

            ///***************************** Solve by SGB *****************************/
            SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 1, 1.0);

            //// Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_SGB += ComputeAccuracySLD( tempMatchTruth, smSolverSLD.GetMatchResult() );

            //// Retrieve function score
            smSolverSLD.ComputeFunScore(scoreMatrixSLD);
            double score_SGB = smSolverSLD.GetFunScore();

            //// Refine through IPFP
            //// Prepare indicator vector and score matrix
            indicatorVec = smSolverSLD.GetIndicatorVec(scoreMatrixSLD);

            //// Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolverSLD.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            //// Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_SGB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_SGB_IPFP = ipfp.GetRefinedScore();

            ///***************************** Solve by HUB *****************************/
            SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 2, 1.0);

            //// Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_HUB += ComputeAccuracySLD( tempMatchTruth, smSolverSLD.GetMatchResult() );

            //// Retrieve function score
            smSolverSLD.ComputeFunScore(scoreMatrixSLD);
            double score_HUB = smSolverSLD.GetFunScore();

            //// Refine through IPFP
            //// Prepare indicator vector and score matrix
            indicatorVec = smSolverSLD.GetIndicatorVec(scoreMatrixSLD);

            //// Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolverSLD.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            //// Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_HUB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_HUB_IPFP = ipfp.GetRefinedScore();

            ///************ Update information **************/
            scoreRatio_IPBOverSGB += score_IPB/score_SGB;
            scoreRatio_HUBOverSGB += score_HUB/score_SGB;

            scoreRatio_IPBOverSGB_IPFP += score_IPB_IPFP/score_SGB_IPFP;
            scoreRatio_HUBOverSGB_IPFP += score_HUB_IPFP/score_SGB_IPFP;

            scoreRatio_IPB_IPFP_Over_SGB += score_IPB_IPFP/score_SGB;
            scoreRatio_HUB_IPFP_Over_SGB += score_HUB_IPFP/score_SGB;
            scoreRatio_SGB_IPFP_Over_SGB += score_SGB_IPFP/score_SGB;

            //// Move to next round
            ++index;
            ////std::cout<<"Done with "<<index<<"th set.\n";
        }
        //std::cout<<"We have "<<index<<" sets for this energy level.\n";

        // Compute average
        scoreRatio_IPBOverSGB /= (double)outputSetList.size();
        scoreRatio_HUBOverSGB /= (double)outputSetList.size();

        scoreRatio_IPBOverSGB_IPFP /= (double)outputSetList.size();
        scoreRatio_HUBOverSGB_IPFP /= (double)outputSetList.size();

        accuracy_IPB /= (float)outputSetList.size();
        accuracy_SGB /= (float)outputSetList.size();
        accuracy_HUB /= (float)outputSetList.size();

        accuracy_IPB_IPFP /= (float)outputSetList.size();
        accuracy_SGB_IPFP /= (float)outputSetList.size();
        accuracy_HUB_IPFP /= (float)outputSetList.size();

        scoreRatio_IPB_IPFP_Over_SGB /= (double)outputSetList.size();
        scoreRatio_HUB_IPFP_Over_SGB /= (double)outputSetList.size();
        scoreRatio_SGB_IPFP_Over_SGB /= (double)outputSetList.size();

        std::cout<<"Function Score:\n";
        std::cout<<"IPB over SGB:     "<<scoreRatio_IPBOverSGB<<"\n";
        std::cout<<"HUB over SGB:     "<<scoreRatio_HUBOverSGB<<"\n";
        std::cout<<"IPB IPFP over SGB:"<<scoreRatio_IPB_IPFP_Over_SGB<<"\n";
        std::cout<<"HUB IPFP over SGB:"<<scoreRatio_HUB_IPFP_Over_SGB<<"\n";
        std::cout<<"SGB IPFP over SGB:"<<scoreRatio_SGB_IPFP_Over_SGB<<"\n";

        std::cout<<"IPB over SGB IPFP:"<<scoreRatio_IPBOverSGB_IPFP<<"\n";
        std::cout<<"HUB over SGB IPFP:"<<scoreRatio_HUBOverSGB_IPFP<<"\n";

        std::cout<<"IPB average accuracy:      "<<accuracy_IPB<<"\n";
        std::cout<<"IPB average accuracy IPFP: "<<accuracy_IPB_IPFP<<"\n";
        std::cout<<"HUB average accuracy:      "<<accuracy_HUB<<"\n";
        std::cout<<"HUB average accuracy IPFP: "<<accuracy_HUB_IPFP<<"\n";
        std::cout<<"SGB average accuracy:      "<<accuracy_SGB<<"\n";
        std::cout<<"SGB average accuracy IPFP: "<<accuracy_SGB_IPFP<<"\n";

        std::cout<<"Finish testing this batch!\n";

        // Move to the next pair of files
        ++energyIndex;
    }
}

// Gather accuracy, function score, and time for 2D point cloud
void RunForTest2DPoint()
{
    // Wrap all this information into a class. ToDo. Flag Ha.

    std::cout<<"Begin to test on 2D point cloud...\n";

    // Parameter settings
    double distanceWeight = 170.0;
    scoreFunFlag = 0;    
    thresholds[0] = 1.0; // wf
    thresholds[1] = distanceWeight;
    thresholds[2] = distanceWeight;
    thresholds[3] = 5.0;   // td
    thresholds[4] = 75.0*PI / 180.0;       // ta
    thresholds[5] = 75.0*PI / 180.0;       // ta
    thresholds[6] = 50.0; // tl

    int numOfRuns     = 50;
    int numOfInliers  = 15;
    int numOfOutliers = 0;

    std::cout<<"Parameter setting: td: "<<thresholds[3]<<".\n";
    std::cout<<"Number of runs:        "<<numOfRuns<<".\n";

    std::vector<float> deformNoiseList(11);
    for(int i = 0; i < (int)deformNoiseList.size(); ++i)
    {
        deformNoiseList[i] = (float)i;
    }

    // Set up random generator
    // Set up uniformly distributed random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> uniformDisTrans(10.0, 300.0);
    std::uniform_real_distribution<float> uniformDisRot(-PI, PI);

    // Loop over one parameter
    for(auto deformNoise = deformNoiseList.begin(); deformNoise != deformNoiseList.end(); ++deformNoise)
    {
        // Initialize accuracy and function score
        std::vector<float> accuracy(6);
        for(int k = 0; k < (int)accuracy.size(); ++k)
            accuracy[k] = 0.0f;

        std::vector<double> funScore(6);
        for(int k = 0; k < (int)funScore.size(); ++k)
            funScore[k] = 0.0;

        std::cout<<"(+) Deform noise: "<<*deformNoise<<"\n";

        // For each parameter, test many times
        for(int i = 0; i < numOfRuns; ++i)
        {
            // Create Module
            PointCloud2D modulePC(numOfInliers);

            // Create Data
            PointCloud2D dataPC;
            dataPC.SetPointCloud2D(modulePC);

            dataPC.AddWhiteNoise(0.0f, *deformNoise);

            dataPC.AddOutLiersUniform(numOfOutliers);

            // Random translate and rotate
            dataPC.TransForm(Vector2f(uniformDisTrans(gen), uniformDisTrans(gen)), uniformDisRot(gen));
            //dataPC.TransForm(Vector2f(10.5, 20.0), 0.0);

            // Add outliers for Module
            modulePC.AddOutLiersUniform(numOfOutliers);  

            // Create match result
            std::vector<int> tempMatchTruth(dataPC.GetInliersNum());
            for(int h = 0; h < (int)tempMatchTruth.size(); ++h)
            {
                tempMatchTruth[h] = h+1;
            }

            // Construct the score matrix
            ConstructScoreMatrix2DP(scoreMatrix2DP, &modulePC, &dataPC, thresholds, scoreFunFlag);
            MatrixXd fullMatrix   = smSolver2DP.GetFullMatrix(scoreMatrix2DP);
            VectorXd indicatorVec;
            IPFP ipfp;

            VectorXi tempResultVecXd;
            std::vector<int> tempResult;

            const int mLen = scoreMatrix2DP.GetModuleNum();
            const int dLen = scoreMatrix2DP.GetDataNum();

            /***************************** Solve by IPB *****************************/
            SolveForResult2DP(smSolver2DP, scoreMatrix2DP, 0);

            // Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy[0] += ComputeAccuracySLD( tempMatchTruth, smSolver2DP.GetMatchResult() );

            // Retrieve function score
            smSolver2DP.ComputeFunScore(scoreMatrix2DP);
            funScore[0] += smSolver2DP.GetFunScore();

            // Refine through IPFP
            // Prepare indicator vector and score matrix
            indicatorVec = smSolver2DP.GetIndicatorVec(scoreMatrix2DP);

            // Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolver2DP.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            // Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy[3] += ComputeAccuracySLD(tempMatchTruth, tempResult);
            funScore[3] += ipfp.GetRefinedScore();

            /***************************** Solve by SGB *****************************/
            SolveForResult2DP(smSolver2DP, scoreMatrix2DP, 1);

            // Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy[1] += ComputeAccuracySLD( tempMatchTruth, smSolver2DP.GetMatchResult() );

            // Retrieve function score
            smSolver2DP.ComputeFunScore(scoreMatrix2DP);
            funScore[1] += smSolver2DP.GetFunScore();

            // Refine through IPFP
            // Prepare indicator vector and score matrix
            indicatorVec = smSolver2DP.GetIndicatorVec(scoreMatrix2DP);

            // Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolver2DP.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            // Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy[4] += ComputeAccuracySLD(tempMatchTruth, tempResult);
            funScore[4] += ipfp.GetRefinedScore();

            /***************************** Solve by HUB *****************************/
            SolveForResult2DP(smSolver2DP, scoreMatrix2DP, 2);

            // Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy[2] += ComputeAccuracySLD( tempMatchTruth, smSolver2DP.GetMatchResult() );

            // Retrieve function score
            smSolver2DP.ComputeFunScore(scoreMatrix2DP);
            funScore[2] += smSolver2DP.GetFunScore();

            // Refine through IPFP
            // Prepare indicator vector and score matrix
            indicatorVec = smSolver2DP.GetIndicatorVec(scoreMatrix2DP);

            // Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolver2DP.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            // Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy[5] += ComputeAccuracySLD(tempMatchTruth, tempResult);
            funScore[5] += ipfp.GetRefinedScore();
        }

        std::cout<<"IPB average accuracy:      "<<accuracy[0]/(float)numOfRuns<<"\n";
        std::cout<<"IPB IPFP average accuracy: "<<accuracy[3]/(float)numOfRuns<<"\n";
        std::cout<<"IPB average funScore:      "<<funScore[0]/(float)numOfRuns<<"\n";
        std::cout<<"IPB IPFP average funScore: "<<funScore[3]/(float)numOfRuns<<"\n";

        std::cout<<"SGB average accuracy:      "<<accuracy[1]/(float)numOfRuns<<"\n";
        std::cout<<"SGB IPFP average accuracy: "<<accuracy[4]/(float)numOfRuns<<"\n";
        std::cout<<"SGB average funScore:      "<<funScore[1]/(float)numOfRuns<<"\n";
        std::cout<<"SGB IPFP average funScore: "<<funScore[4]/(float)numOfRuns<<"\n";

        std::cout<<"HUB average accuracy:      "<<accuracy[2]/(float)numOfRuns<<"\n";
        std::cout<<"HUB IPFP average accuracy: "<<accuracy[5]/(float)numOfRuns<<"\n";
        std::cout<<"HUB average funScore:      "<<funScore[2]/(float)numOfRuns<<"\n";
        std::cout<<"HUB IPFP average funScore: "<<funScore[5]/(float)numOfRuns<<"\n";
        std::cout<<"Done.\n";
    }

    std::cout<<"Finish testing on 2D point cloud!\n";
}

// Test one TPS point cloud
void RunForTest2DPointTPS()
{
    // Wrap all this information into a class. ToDo. Flag Ha.
    std::cout<<"Begin to test on TPS 2D point cloud...\n";

    std::string inputListFile  = "../Data/TPSPoint/fileList.txt";
    std::string dataDir        = "../Data/TPSPoint/withOutliers/30inliers_10outliers_30trails/";
    std::vector< std::string > inputFileList;
    std::vector< std::string > outputFileList;

    std::ifstream file(inputListFile);
    std::string   line;

    // Read the file list
    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputListFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            // Parse each line for file name out by space
            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ' ');

            inputFileList.push_back(tokens[0]);
            outputFileList.push_back(tokens[1]);
        }

        std::cout<<"Finish reading in all the files!\n";
    }

    file.close();

    /////////  Process on all the files /////////
    unsigned int i = 0;
    while(i < inputFileList.size())
    {
        // Read in the data and model from the file
        // ATTENTION!! "inputPoints" are source and "outputPoints" are deformed points
        std::string inputFile  = dataDir+inputFileList[i];
        std::string outputFile = dataDir+outputFileList[i];

        std::vector< PointCloud2D > inputSetList;
        std::vector< PointCloud2D > outputSetList;

        ////////////////////// Read input points //////////////////////
        file.open(inputFile, std::ifstream::in);

        if( !file.is_open() )
            std::cout<<"Unable to open "<<inputFile<<"...\n";
        else
        {
            std::cout<<"Read in the file...\n";

            bool initFlag = true;

            std::vector<Vector2f> tempPointList;

            while( std::getline(file, line) )
            {
                // Protect from empty lines
                if( "" == line )
                    continue;

                std::vector<std::string> tokens;
                ParseByCharater(line, tokens, '\t');

                if("e" == tokens[0])
                {
                    // We create new input and outputs
                    if( !initFlag )
                    {
                        inputSetList.push_back( PointCloud2D(tempPointList) );
                        tempPointList.clear();
                    }
                    // Avoid the first e
                    initFlag = false;
                }
                else
                {
                    // Read in the two points
                    tempPointList.push_back( Vector2f(std::stof(tokens[0]), std::stof(tokens[1])) );
                }

            }
            // Record the last set in the file
            inputSetList.push_back( PointCloud2D(tempPointList) );

            std::cout<<"Finish reading in the inputPoints file!\n";
        }

        file.close();

        ////////////////////// Read output points //////////////////////
        file.open(outputFile, std::ifstream::in);

        if( !file.is_open() )
            std::cout<<"Unable to open "<<outputFile<<"...\n";
        else
        {
            std::cout<<"Read in the file...\n";

            bool initFlag = true;

            std::vector<Vector2f> tempPointList;

            while( std::getline(file, line) )
            {
                // Protect from empty lines
                if( "" == line )
                    continue;

                std::vector<std::string> tokens;
                ParseByCharater(line, tokens, '\t');

                if("e" == tokens[0])
                {
                    // We create new input and outputs
                    if( !initFlag )
                    {
                        outputSetList.push_back( PointCloud2D(tempPointList) );
                        tempPointList.clear();
                    }
                    // Avoid the first e
                    initFlag = false;
                }
                else
                {
                    // Read in the two points
                    tempPointList.push_back( Vector2f(std::stof(tokens[0]), std::stof(tokens[1])) );
                }

            }
            // Record the last set in the file
            outputSetList.push_back( PointCloud2D(tempPointList) );

            std::cout<<"Finish reading in the outputPoints file!\n";
        }

        file.close();  

        std::cout<<i<<"th file. "<<(int)inputSetList.size()<<" sets in all.\n";


        // Parameter settings
        double distanceWeight = 170.0;
        scoreFunFlag  = 2;    
        thresholds[0] = 0.0; // wf
        thresholds[1] = distanceWeight;
        thresholds[2] = distanceWeight;
        thresholds[3] = 20.0;   // td
        thresholds[4] = 90.0*PI / 180.0;       // ta
        thresholds[5] = 90.0*PI / 180.0;       // ta
        thresholds[6] = 50.0; // tl

        // Initialize data we want to record
        float accuracy_IPB = 0.0f;
        float accuracy_SGB = 0.0f;
        float accuracy_HUB = 0.0f;

        float accuracy_IPB_IPFP = 0.0f;
        float accuracy_SGB_IPFP = 0.0f;
        float accuracy_HUB_IPFP = 0.0f;

        double scoreRatio_IPBOverSGB = 0.0;
        double scoreRatio_HUBOverSGB = 0.0;

        double scoreRatio_IPB_IPFP_Over_SGB = 0.0;
        double scoreRatio_HUB_IPFP_Over_SGB = 0.0;
        double scoreRatio_SGB_IPFP_Over_SGB = 0.0;

        double scoreRatio_IPBOverSGB_IPFP = 0.0f;
        double scoreRatio_HUBOverSGB_IPFP = 0.0f;

        int index = 0;
        while( index < (int)inputSetList.size() )
        {

            // Create Module
            PointCloud2D modulePC(inputSetList[index]);
            // Create Data
            PointCloud2D dataPC(outputSetList[index]);

            // Create match result
            std::vector<int> tempMatchTruth(dataPC.GetInliersNum());
            for(int h = 0; h < (int)tempMatchTruth.size(); ++h)
            {
                tempMatchTruth[h] = 0;
            }
            // We have 10 outliers
            for(int h = 0; h < (int)tempMatchTruth.size()-10; ++h)
            {
                tempMatchTruth[h] = h+1;
            }

            // Construct the score matrix
            ConstructScoreMatrix2DP(scoreMatrix2DP, &modulePC, &dataPC, thresholds, scoreFunFlag);
            MatrixXd fullMatrix = smSolver2DP.GetFullMatrix(scoreMatrix2DP);
            VectorXd indicatorVec;
            IPFP ipfp;

            VectorXi tempResultVecXd;
            std::vector<int> tempResult;

            const int mLen = scoreMatrix2DP.GetModuleNum();
            const int dLen = scoreMatrix2DP.GetDataNum();

            /***************************** Solve by IPB *****************************/
            SolveForResult2DP(smSolver2DP, scoreMatrix2DP, 0);

            // Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_IPB += ComputeAccuracySLD( tempMatchTruth, smSolver2DP.GetMatchResult() );

            // Retrieve function score
            smSolver2DP.ComputeFunScore(scoreMatrix2DP);
            double score_IPB = smSolver2DP.GetFunScore();

            // Refine through IPFP
            // Prepare indicator vector and score matrix
            indicatorVec = smSolver2DP.GetIndicatorVec(scoreMatrix2DP);

            // Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolver2DP.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            // Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_IPB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_IPB_IPFP = ipfp.GetRefinedScore();

            /***************************** Solve by SGB *****************************/
            SolveForResult2DP(smSolver2DP, scoreMatrix2DP, 1);

            // Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_SGB += ComputeAccuracySLD( tempMatchTruth, smSolver2DP.GetMatchResult() );

            // Retrieve function score
            smSolver2DP.ComputeFunScore(scoreMatrix2DP);
            double score_SGB = smSolver2DP.GetFunScore();

            // Refine through IPFP
            // Prepare indicator vector and score matrix
            indicatorVec = smSolver2DP.GetIndicatorVec(scoreMatrix2DP);

            // Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolver2DP.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            // Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_SGB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_SGB_IPFP = ipfp.GetRefinedScore();

            /***************************** Solve by HUB *****************************/
            SolveForResult2DP(smSolver2DP, scoreMatrix2DP, 2);

            // Retrieve accuracy --- we can use accuracy function for SLD.
            accuracy_HUB += ComputeAccuracySLD( tempMatchTruth, smSolver2DP.GetMatchResult() );

            // Retrieve function score
            smSolver2DP.ComputeFunScore(scoreMatrix2DP);
            double score_HUB = smSolver2DP.GetFunScore();

            // Refine through IPFP
            // Prepare indicator vector and score matrix
            indicatorVec = smSolver2DP.GetIndicatorVec(scoreMatrix2DP);

            // Do IPFP
            ipfp.Reset(fullMatrix, indicatorVec, smSolver2DP.GetFunScore(), 10, mLen, dLen);
            ipfp.DoIPFP();
            // Obtain refined accuracy
            tempResultVecXd = ipfp.GetRefinedResult();
            tempResult = std::vector<int>(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];

            accuracy_HUB_IPFP    += ComputeAccuracySLD(tempMatchTruth, tempResult);
            double score_HUB_IPFP = ipfp.GetRefinedScore();

            /************ Update information **************/
            scoreRatio_IPBOverSGB += score_IPB/score_SGB;
            scoreRatio_HUBOverSGB += score_HUB/score_SGB;

            scoreRatio_IPBOverSGB_IPFP += score_IPB_IPFP/score_SGB_IPFP;
            scoreRatio_HUBOverSGB_IPFP += score_HUB_IPFP/score_SGB_IPFP;

            scoreRatio_IPB_IPFP_Over_SGB += score_IPB_IPFP/score_SGB;
            scoreRatio_HUB_IPFP_Over_SGB += score_HUB_IPFP/score_SGB;
            scoreRatio_SGB_IPFP_Over_SGB += score_SGB_IPFP/score_SGB;

            // Move to next round
            ++index;
            //std::cout<<"Done with "<<index<<"th set.\n";
        }

        // Compute average
        scoreRatio_IPBOverSGB /= (double)inputSetList.size();
        scoreRatio_HUBOverSGB /= (double)inputSetList.size();

        scoreRatio_IPBOverSGB_IPFP /= (double)inputSetList.size();
        scoreRatio_HUBOverSGB_IPFP /= (double)inputSetList.size();

        accuracy_IPB /= (float)inputSetList.size();
        accuracy_SGB /= (float)inputSetList.size();
        accuracy_HUB /= (float)inputSetList.size();

        accuracy_IPB_IPFP /= (float)inputSetList.size();
        accuracy_SGB_IPFP /= (float)inputSetList.size();
        accuracy_HUB_IPFP /= (float)inputSetList.size();

        scoreRatio_IPB_IPFP_Over_SGB /= (double)inputSetList.size();
        scoreRatio_HUB_IPFP_Over_SGB /= (double)inputSetList.size();
        scoreRatio_SGB_IPFP_Over_SGB /= (double)inputSetList.size();

        std::cout<<"Function Score:\n";
        std::cout<<"IPB over SGB:     "<<scoreRatio_IPBOverSGB<<"\n";
        std::cout<<"HUB over SGB:     "<<scoreRatio_HUBOverSGB<<"\n";
        std::cout<<"IPB IPFP over SGB:"<<scoreRatio_IPB_IPFP_Over_SGB<<"\n";
        std::cout<<"HUB IPFP over SGB:"<<scoreRatio_HUB_IPFP_Over_SGB<<"\n";
        std::cout<<"SGB IPFP over SGB:"<<scoreRatio_SGB_IPFP_Over_SGB<<"\n";

        std::cout<<"IPB over SGB IPFP:"<<scoreRatio_IPBOverSGB_IPFP<<"\n";
        std::cout<<"HUB over SGB IPFP:"<<scoreRatio_HUBOverSGB_IPFP<<"\n";

        std::cout<<"IPB average accuracy:      "<<accuracy_IPB<<"\n";
        std::cout<<"IPB average accuracy IPFP: "<<accuracy_IPB_IPFP<<"\n";
        std::cout<<"HUB average accuracy:      "<<accuracy_HUB<<"\n";
        std::cout<<"HUB average accuracy IPFP: "<<accuracy_HUB_IPFP<<"\n";
        std::cout<<"SGB average accuracy:      "<<accuracy_SGB<<"\n";
        std::cout<<"SGB average accuracy IPFP: "<<accuracy_SGB_IPFP<<"\n";

        std::cout<<"Finish testing on TSP 2D point cloud!\n";

        // Move to the next pair of files
        ++i;
    }
}

// Gather accuracy, function score, and time for single line direction
void RunForTestSingleDir()
{
    // Wrap all this information into a class. ToDo. Flag Ha.

    std::cout<<"Begin to test improvement of score and accuracy with IPFP for IPB, SGB and HUB...\n";

    std::string inputFile = "../Data/AccuracyScoreSumBatch.txt";
    std::string dataDir   = "../Data/";
    std::string fileExt   = ".pdb";

    //std::string outPutType    = "record";
    //std::string outPutDirGau  = "../Data/AccuracyScoreSumBatch/";
    //std::string outPutDirEpan = "../Data/IterativePicking/EpanAccuracy/";
    //std::string outPutExt     = ".txt";

    std::vector< std::vector<std::string> > fileList;
    std::vector< std::vector<int> >         matchTruthList;

    // Read in all the batch data from AccuracyBatch.txt
    // The information in the file has to be line by line, in the order above.
    // The angles in the AccuracyBatch file are measured by degree
    std::ifstream file(inputFile);
    std::string line;

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ':');

            if("fileName" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<std::string> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( content[i] );
                }
                fileList.push_back(temp);
            }
            else if("matchTruth" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<int> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stoi(content[i]) );
                }
                matchTruthList.push_back(temp);
            }

        }

        std::cout<<"Finish reading in the file!\n";
    }

    file.close();        

    // Parameter settings
    double distanceWeight = 170.0;
    scoreFunFlag = 0;    
    thresholds[0] = 1.0; // wf
    thresholds[1] = distanceWeight;
    thresholds[2] = distanceWeight;
    thresholds[3] = 5.0;   // td
    thresholds[4] = 75.0*PI / 180.0;       // ta
    thresholds[5] = 75.0*PI / 180.0;       // ta
    thresholds[6] = 50.0; // tl

    std::cout<<"Parameter setting: wf, dw, td, ta, tl"<<thresholds[0]<<", "<<thresholds[1]<<", "<<thresholds[3]<<", "<<thresholds[4]<<", "<<thresholds[6]<<".\n";

    unsigned int i = 0;
    while(i < fileList.size())
    {
        std::cout<<"(+)"<<i<<" Begin to process "<<fileList[i][0]<<" & "<<fileList[i][2]<<".\n";

        //// Read in file and retrieve geometry
        // Read in module set
        std::string fullPath = dataDir + fileList[i][0] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        char chainID = *(fileList[i][1].c_str());
        RetrieveHelixGeometry(moduleHelix, pdbReader, chainID, moduleOffset);

        // Read in data set
        fullPath = dataDir + fileList[i][2] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        chainID = *(fileList[i][3].c_str());
        RetrieveHelixGeometry(dataHelix, pdbReader, chainID, dataOffset);
        // Make sure module set has more helices
        SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

        std::vector<int> tempMatchTruth = matchTruthList[i];

        // Update the match truth based on different parameter setting. 
        for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
        {
            if( 0 != tempMatchTruth[h] )
            {
                // The truth has the same format as in Mathematica --- index starts from 1
                float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                if( fabs(moduleLen - datalen) > thresholds[6] && thresholds[6] > 0 )
                    tempMatchTruth[h] = 0;
            }
        }

        double accuracy = 0.0;
        double funScore = 0.0;

        // Time for matrix (graph) construction and solve
        //cpuTimer.Start();

        // Compute accuracy
        ConstructScoreMatrixSLD(scoreMatrixSLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);

        const int mLen = scoreMatrixSLD.GetModuleNum();
        const int dLen = scoreMatrixSLD.GetDataNum();

        //elapsedTime = cpuTimer.Tick();
        //std::cout<<"Time to construct matrix: "<<elapsedTime*(1e-3)<<"s.\n";

        /*************** Record and time Iterative Picking *****************/
        std::vector<std::string> methodId(3);
        methodId[0] = "(IPB)"; methodId[1] = "(SGB)"; methodId[2] = "(HUB)";

        for(int index = 0; index < 3; ++index)
        {
            // Solve for match result
            SolveForResultSLD(smSolverSLD, scoreMatrixSLD, index);

            //std::cout<<"Match result:\n";
            //for(auto it = smSolverSLD.GetMatchResult().begin(); it != smSolverSLD.GetMatchResult().end(); ++it)
            //{
            //    std::cout<<*it<<", ";
            //}
            //std::cout<<"\n";

            // Compute the accuracy
            //accuracy = ComputeAccuracySLD(tempMatchTruth, smSolverSLD.GetMatchResult());
            //std::cout<<methodId[index]+" accuracy: "<<accuracy<<"\n";

            // Compute the score sum
            //smSolverSLD.ComputeFunScore(scoreMatrixSLD);
            //funScore = smSolverSLD.GetFunScore();
            //std::cout<<methodId[index]+" function score: "<<funScore<<"\n";

            // Compute the refined score sum (x'Mx)
            smSolverSLD.ComputeFunScore(scoreMatrixSLD);
            funScore = smSolverSLD.GetFunScore();
            // Prepare indicator vector and score matrix
            VectorXd indicatorVec = smSolverSLD.GetIndicatorVec(scoreMatrixSLD);
            MatrixXd fullMatrix   = smSolverSLD.GetFullMatrix(scoreMatrixSLD);
            // Do IPFP
            IPFP ipfp(fullMatrix, indicatorVec, funScore, 10, mLen, dLen);
            ipfp.DoIPFP();
            // Obtain refined accuracy
            VectorXi tempResultVecXd = ipfp.GetRefinedResult();
            std::vector<int> tempResult(tempResultVecXd.size());
            for(int i = 0; i < tempResultVecXd.size(); ++i)
                tempResult[i] = tempResultVecXd[i];
            accuracy = ComputeAccuracySLD(tempMatchTruth, tempResult);

            // Output refined result
            std::cout<<methodId[index] + "Refined score:    "<<ipfp.GetRefinedScore()<<"\n";
            std::cout<<methodId[index] + "Refined accuracy: "<<accuracy<<"\n";
        }

        /*************** Record and Time Simple Greedy *****************/

        //// Solve for match result
        //SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 1);

        //// Compute the accuracy
        //accuracy = ComputeAccuracySLD(tempMatchTruth, smSolverSLD.GetMatchResult());
        //std::cout<<"(SGB) accuracy: "<<accuracy<<"\n";

        //// Compute the score sum
        //smSolverSLD.ComputeFunScore(scoreMatrixSLD);
        //funScore = smSolverSLD.GetFunScore();
        //std::cout<<"(SGB) function score: "<<funScore<<"\n";

        /*************** Record and Time Hungarian *****************/
        //// Time for matrix (graph) construction and solve
        //cpuTimer.Start();

        //// Solve for match result
        //SolveForResultSLD(smSolverSLD, scoreMatrixSLD, 2);

        //// Compute elapsed time
        //elapsedTime = cpuTimer.Tick();
        //std::cout<<"(HUB) Time to solve for result: "<<elapsedTime*(1e-3)<<"s.\n";

        //// Compute the accuracy
        //accuracy = ComputeAccuracySLD(tempMatchTruth, smSolverSLD.GetMatchResult());
        //std::cout<<"(HUB) accuracy: "<<accuracy<<"\n";

        //// Compute the score sum
        //smSolverSLD.ComputeFunScore(scoreMatrixSLD);
        //funScore = smSolverSLD.GetFunScore();
        //std::cout<<"(HUB) function score: "<<funScore<<"\n";

        // Move to the next sample set
        std::cout<<"Done.\n";

        ++i;
    }

    std::cout<<"Finish processing the batch!\n";
}

// Record accuracy, score sum and time for given method with given parameters. (Iterative Pick, Simple Greedy)
void RunForAccuracyAndScoreSum()
{
    // Wrap all this information into a class. ToDo. Flag Ha.

    std::cout<<"Begin to compute accuracy and score sum for IPB and SGB...\n";

    std::string inputFile = "../Data/AccuracyScoreSumBatch.txt";
    std::string dataDir   = "../Data/";
    std::string fileExt   = ".pdb";

    //std::string outPutType    = "record";
    //std::string outPutDirGau  = "../Data/AccuracyScoreSumBatch/";
    //std::string outPutDirEpan = "../Data/IterativePicking/EpanAccuracy/";
    //std::string outPutExt     = ".txt";

    std::vector< std::vector<std::string> > fileList;
    std::vector< std::vector<int> >         matchTruthList;

    // Read in all the batch data from AccuracyBatch.txt
    // The information in the file has to be line by line, in the order above.
    // The angles in the AccuracyBatch file are measured by degree
    std::ifstream file(inputFile);
    std::string line;

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ':');

            if("fileName" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<std::string> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( content[i] );
                }
                fileList.push_back(temp);
            }
            else if("matchTruth" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<int> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stoi(content[i]) );
                }
                matchTruthList.push_back(temp);
            }
         
        }

        std::cout<<"Finish reading in the file!\n";
    }

    file.close();        

    // Parameter settings
    double distanceWeight = 120.0;
    scoreFunFlag = 0;    
    thresholds[0] = 1.0; // wf
    thresholds[1] = distanceWeight;
    thresholds[2] = distanceWeight;
    thresholds[3] = 5.0;   // td
    thresholds[4] = 90.0*PI / 180.0;       // ta
    thresholds[5] = 90.0*PI / 180.0;       // ta
    thresholds[6] = 40.0; // tl
    
    std::cout<<"Parameter setting: wf, dw, td, ta, tl"<<thresholds[0]<<", "<<thresholds[1]<<", "<<thresholds[3]<<", "<<thresholds[4]<<", "<<thresholds[6]<<".\n";

    unsigned int i = 0;
    while(i < fileList.size())
    {
        std::cout<<"(+) Begin to process "<<fileList[i][0]<<" & "<<fileList[i][2]<<".\n";

        //// Read in file and retrieve geometry
        // Read in module set
        std::string fullPath = dataDir + fileList[i][0] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        char chainID = *(fileList[i][1].c_str());
        RetrieveHelixGeometry(moduleHelix, pdbReader, chainID, moduleOffset);

        // Read in data set
        fullPath = dataDir + fileList[i][2] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        chainID = *(fileList[i][3].c_str());
        RetrieveHelixGeometry(dataHelix, pdbReader, chainID, dataOffset);
        // Make sure module set has more helices
        SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

        std::vector<int> tempMatchTruth = matchTruthList[i];

        // Update the match truth based on different parameter setting. 
        for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
        {
            if( 0 != tempMatchTruth[h] )
            {
                // The truth has the same format as in Mathematica --- index starts from 1
                float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                if( fabs(moduleLen - datalen) > thresholds[6] && thresholds[6] > 0 )
                    tempMatchTruth[h] = 0;
            }
        }

        double accuracy = 0.0;
        double funScore = 0.0;

        // Time for matrix (graph) construction and solve
        cpuTimer.Start();

        // Compute accuracy
        ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);

        elapsedTime = cpuTimer.Tick();
        std::cout<<"Time to construct matrix: "<<elapsedTime*(1e-3)<<"s.\n";

        /*************** Time Iterative Picking *****************/
        // Time for matrix (graph) construction and solve
        cpuTimer.Start();

        // Solve for match result
        SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 0);

        // Compute elapsed time
        elapsedTime = cpuTimer.Tick();
        std::cout<<"(IPB) Time to solve for result: "<<elapsedTime*(1e-3)<<"s.\n";

        // Compute the accuracy
        accuracy = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
        std::cout<<"(IPB) accuracy: "<<accuracy<<"\n";

        // Compute the score sum
        smSolverDLD.ComputeFunScore(scoreMatrixDLD);
        funScore = smSolverDLD.GetFunScore();
        std::cout<<"(IPB) function score: "<<funScore<<"\n";


        /*************** Record and Time Simple Greedy *****************/
        // Time for matrix (graph) construction and solve
        cpuTimer.Start();

        // Solve for match result
        SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 1);

        // Compute elapsed time
        elapsedTime = cpuTimer.Tick();
        std::cout<<"(SGB) Time to solve for result: "<<elapsedTime*(1e-3)<<"s.\n";

        // Compute the accuracy
        accuracy = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
        std::cout<<"(SGB) accuracy: "<<accuracy<<"\n";

        // Compute the score sum
        smSolverDLD.ComputeFunScore(scoreMatrixDLD);
        funScore = smSolverDLD.GetFunScore();
        std::cout<<"(SGB) function score: "<<funScore<<"\n";

        // Move to the next sample set
        std::cout<<"Done with "<<fileList[i][0]<<" & "<<fileList[i][2]<<".\n";

        ++i;
    }

    std::cout<<"Finish processing the batch!\n";
}

// Record accuracy, score sum and time for given method with given parameters. (Iterative Pick, Iterative pick new stop)
void RunForAccuracyAndScoreSumNewStop()
{
    // Wrap all this information into a class. ToDo. Flag Ha.

    std::cout<<"Begin to compute accuracy and score sum for IPB and SGB...\n";

    std::string inputFile = "../Data/AccuracyScoreSumBatch.txt";
    std::string dataDir   = "../Data/";
    std::string fileExt   = ".pdb";

    //std::string outPutType    = "record";
    //std::string outPutDirGau  = "../Data/AccuracyScoreSumBatch/";
    //std::string outPutDirEpan = "../Data/IterativePicking/EpanAccuracy/";
    //std::string outPutExt     = ".txt";

    std::vector< std::vector<std::string> > fileList;
    std::vector< std::vector<int> >         matchTruthList;

    // Read in all the batch data from AccuracyBatch.txt
    // The information in the file has to be line by line, in the order above.
    // The angles in the AccuracyBatch file are measured by degree
    std::ifstream file(inputFile);
    std::string line;

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ':');

            if("fileName" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<std::string> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( content[i] );
                }
                fileList.push_back(temp);
            }
            else if("matchTruth" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<int> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stoi(content[i]) );
                }
                matchTruthList.push_back(temp);
            }

        }

        std::cout<<"Finish reading in the file!\n";
    }

    file.close();        

    // Parameter settings
    double distanceWeight = 120.0;
    scoreFunFlag = 0;    
    thresholds[0] = 1.0; // wf
    thresholds[1] = distanceWeight;
    thresholds[2] = distanceWeight;
    thresholds[3] = 5.0;   // td
    thresholds[4] = 90.0*PI / 180.0;       // ta
    thresholds[5] = 90.0*PI / 180.0;       // ta
    thresholds[6] = 40.0; // tl

    std::cout<<"Parameter setting: wf, dw, td, ta, tl"<<thresholds[0]<<", "<<thresholds[1]<<", "<<thresholds[3]<<", "<<thresholds[4]<<", "<<thresholds[6]<<".\n";

    unsigned int i = 0;
    while(i < fileList.size())
    {
        std::cout<<"(+) Begin to process "<<fileList[i][0]<<" & "<<fileList[i][2]<<".\n";

        //// Read in file and retrieve geometry
        // Read in module set
        std::string fullPath = dataDir + fileList[i][0] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        char chainID = *(fileList[i][1].c_str());
        RetrieveHelixGeometry(moduleHelix, pdbReader, chainID, moduleOffset);

        // Read in data set
        fullPath = dataDir + fileList[i][2] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        chainID = *(fileList[i][3].c_str());
        RetrieveHelixGeometry(dataHelix, pdbReader, chainID, dataOffset);
        // Make sure module set has more helices
        SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

        std::vector<int> tempMatchTruth = matchTruthList[i];

        // Update the match truth based on different parameter setting. 
        for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
        {
            if( 0 != tempMatchTruth[h] )
            {
                // The truth has the same format as in Mathematica --- index starts from 1
                float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                if( fabs(moduleLen - datalen) > thresholds[6] && thresholds[6] > 0 )
                    tempMatchTruth[h] = 0;
            }
        }

        double accuracy = 0.0;
        double funScore = 0.0;

        // Compute accuracy
        ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);

        /*************** Record and Time Iterative Picking *****************/
        //// Time for matrix (graph) construction and solve
        //cpuTimer.Start();

        //// Solve for match result
        //SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 0);

        //// Compute elapsed time
        //elapsedTime = cpuTimer.Tick();
        //std::cout<<"(IPB) Time to solve for result: "<<elapsedTime*(1e-3)<<"s.\n";

        //// Compute the accuracy
        //accuracy = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
        //std::cout<<"(IPB) accuracy: "<<accuracy<<"\n";

        //// Compute the score sum
        //smSolverDLD.ComputeFunScore(scoreMatrixDLD);
        //funScore = smSolverDLD.GetFunScore();
        //std::cout<<"(IPB) function score: "<<funScore<<"\n";


        /*************** Record and Time New Stop *****************/

        //// Initialize the matlab library
        //if( !LeadingEigenvecInitialize() )
        //{
        //    std::cerr << "Could not initialize the library properly."<< std::endl;
        //    return;
        //}

        // Time for matrix (graph) construction and solve
        cpuTimer.Start();

        // Solve for match result
        SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 2, 0.0);

        // Compute elapsed time
        elapsedTime = cpuTimer.Tick();
        std::cout<<"(New Stop) Time to solve for result: "<<elapsedTime*(1e-3)<<"s.\n";

        //// End of matlab kernel
        //LeadingEigenvecTerminate();

        // Compute the accuracy
        accuracy = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
        std::cout<<"(New Stop) accuracy: "<<accuracy<<"\n";

        // Compute the score sum
        smSolverDLD.ComputeFunScore(scoreMatrixDLD);
        funScore = smSolverDLD.GetFunScore();
        std::cout<<"(New Stop) function score: "<<funScore<<"\n";

        // Move to the next sample set
        std::cout<<"Done with "<<fileList[i][0]<<" & "<<fileList[i][2]<<".\n";

        ++i;
    }

    std::cout<<"Finish processing the batch!\n";
}

// Compute the accuracy volume. Fix distance weight threshold
void RunForAccuracyVolumeADLDW()
{
    // Wrap all this information into a class. ToDo.

    std::cout<<"Begin to compute acuracy volume for ADL (IPB)...\n";

    std::string inputFile = "../Data/AccuracyBatch_ADL_DW.txt";
    std::string dataDir   = "../Data/";
    std::string fileExt   = ".pdb";

    std::string outPutType    = "volume";
    std::string outPutDirGau  = "../Data/IterativePicking/GaussianAccuracyADLDW/";
    std::string outPutDirEpan = "../Data/IterativePicking/EpanAccuracy/";
    std::string outPutExt     = ".csv";

    std::vector< std::vector<std::string> > fileList;
    std::vector< std::vector<int> >         matchTruthList;
    std::vector< double >                   scoreFunList;
    std::vector< std::vector<double> >      angleList;
    std::vector< std::vector<double> >      midPDisList;
    std::vector< std::vector<double> >      helixLenList;

    // Read in all the batch data from AccuracyBatch.txt
    // The information in the file has to be line by line, in the order above.
    // The angles in the AccuracyBatch file are measured by degree
    std::ifstream file(inputFile);
    std::string line;

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ':');

            if("fileName" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<std::string> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( content[i] );
                }
                fileList.push_back(temp);
            }
            else if("matchTruth" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<int> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stoi(content[i]) );
                }
                matchTruthList.push_back(temp);
            }
            else if("scoreFun" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                // Now the score fun means distance weight threshold in Angstrom
                scoreFunList.push_back(std::stod(content[0]));
            }
            else if("angle" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( DToR( std::stod(content[i]) ) );
                }
                angleList.push_back(temp);
            }
            else if("midPointDistance" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                midPDisList.push_back(temp);
            }
            else if("helixLenght" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                helixLenList.push_back(temp);
            }            
        }

        std::cout<<"Finish reading in the file!\n";
    }

    file.close();        

    unsigned int i = 0;
    while(i < fileList.size())
    {
        std::cout<<"Begin to process the "<<i+1<<"th sample set...\n";

        //// Time it
        cpuTimer.Start();

        //// Read in file and retrieve geometry
        // Read in module set
        std::string fullPath = dataDir + fileList[i][0] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        char chainID = *(fileList[i][1].c_str());
        RetrieveHelixGeometry(moduleHelix, pdbReader, chainID, moduleOffset);

        // Read in data set
        fullPath = dataDir + fileList[i][2] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        chainID = *(fileList[i][3].c_str());
        RetrieveHelixGeometry(dataHelix, pdbReader, chainID, dataOffset);
        // Make sure module set has more helices
        SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

        // Initialize the accuracy volume --- at least 1x1x1 volume        
        int angleNum    = static_cast<int>( (angleList[i][1]    - angleList[i][0])    / angleList[i][2]    + indexEpsilon ) + 1;
        int midPDisNum  = static_cast<int>( (midPDisList[i][1]  - midPDisList[i][0])  / midPDisList[i][2]  + indexEpsilon ) + 1;
        int helixLenNum = static_cast<int>( (helixLenList[i][1] - helixLenList[i][0]) / helixLenList[i][2] + indexEpsilon ) + 1;

        scoreFunFlag = 0;
        double distanceWeight = scoreFunList[i];
        //std::cout<<angleNum<<", "<<midPDisNum<<", "<<helixLenNum<<"\n";

        // Allocate memory to hold the volume
        float *accuracyVolume = new float[angleNum*midPDisNum*helixLenNum];

        int angleOffset   = midPDisNum*helixLenNum;
        int midPDisOffset = helixLenNum;

        int nThreads = 0;
        //int askedThreads = 2;
        //omp_set_num_threads(askedThreads);

        int a = 0; int d = 0; int l = 0;
        ScoreMatrixDoubleLineDir scoreMatrixDLD;
        SMSolverDLD                 smSolverDLD;
        //unsigned int h = 0;
        //#pragma omp parallel default(shared) private(a,d,l,thresholds,scoreMatrixDLD,smSolverDLD)
        #pragma omp parallel default(shared) private(a,d,l,scoreMatrixDLD,smSolverDLD)
        {
            #pragma omp master
            nThreads = omp_get_num_threads();

            // Enter the loop to compute the accuracy volume --- #pragma omp parallel for
            //#pragma omp for
            //#pragma omp parallel for private(scoreMatrixDLD,smSolverDLD)       
            #pragma omp for
            for(a = 0; a < angleNum; ++a)
                for(d = 0; d < midPDisNum; ++d)
                    for(l = 0; l < helixLenNum; ++l)
                    {
                        // Set up new parameters --- No distance weighting so far
                        //thresholds[0] = 0; // wf
                        //thresholds[3] = midPDisList[i][0] + d*midPDisList[i][2];   // td
                        //thresholds[4] = angleList[i][0] + a*angleList[i][2];       // ta
                        //thresholds[5] = thresholds[4]; // ta
                        //thresholds[6] = helixLenList[i][0] + l*helixLenList[i][2]; // tl

                        std::vector<double> tempThresholds(7,0.0);
                        tempThresholds[0] = 1.0; // wf
                        tempThresholds[1] = distanceWeight;
                        tempThresholds[2] = distanceWeight;
                        tempThresholds[3] = midPDisList[i][0] + d*midPDisList[i][2];   // td
                        tempThresholds[4] = angleList[i][0] + a*angleList[i][2];       // ta
                        tempThresholds[5] = angleList[i][0] + a*angleList[i][2];       // ta
                        tempThresholds[6] = helixLenList[i][0] + l*helixLenList[i][2]; // tl

                        std::vector<int> tempMatchTruth = matchTruthList[i];

                        // Update the match truth based on different parameter setting. 
                        for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
                        {
                            if( 0 != tempMatchTruth[h] )
                            {
                                // The truth has the same format as in Mathematica --- index starts from 1
                                float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                                float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                                if( fabs(moduleLen - datalen) > tempThresholds[6] && tempThresholds[6] > 0 )
                                    tempMatchTruth[h] = 0;
                            }
                        }

                        //std::cout<<"Thresholds: "<<thresholds[0]<<", "<<thresholds[3]<<", "<<thresholds[4]<<", "<<thresholds[5]<<", "<<thresholds[6]<<"\n";
                        // Compute accuracy
                        ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, tempThresholds, scoreFunFlag);
                        // Solve for match result
                        SolveForResultDLD(smSolverDLD, scoreMatrixDLD);
                        // Compute and record accuracy
                        accuracyVolume[a*angleOffset + d*midPDisOffset + l] = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());


                    }// End for loop

            #pragma omp barrier
        } // End Parallel Block

        //// Export the result to the disk.
        // Make output name and choose output folder
        std::string exportFile = "";
        // Include the distance value in file name
        std::ostringstream strs;
        strs << scoreFunList[i];
        exportFile = outPutDirGau + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + "_" + strs.str() + "a" + outPutExt;

        std::ofstream outputFile(exportFile, std::ios::trunc);

        // Use to convert float to string
        std::stringstream stream;
        // Initialize
        stream.str(std::string());

        if( !outputFile.is_open() )
            std::cout<<"Unable to open "<<exportFile<<"...\n";
        else
        {
            // Output volume
            for(int a = 0; a < angleNum; ++a)
            {
                // angle dimension
                stream<<"\"{";
                for(int d = 0; d < midPDisNum; ++d)
                {
                    // midpoint distance dimension
                    stream<<"{";
                    for(int l = 0; l < helixLenNum; ++l)
                    {
                        // helix length dimension
                        stream<<accuracyVolume[a*angleOffset + d*midPDisOffset + l];
                        if(l < helixLenNum - 1)
                            stream<<", ";
                    }

                    stream<<"}";
                    if(d < midPDisNum - 1)
                        stream<<", ";
                }

                stream<<"}\"";
                if(a < angleNum - 1)
                    stream<<",";
            }

            stream<<"\n";

            // Output Parameters --- a, d, l
            // To match for unit in Mathematica, scale by 0.01
            float tempScale = 0.01f;
            stream<<"\"{"<<angleList[i][0]<<", "<<angleList[i][1]<<", "<<angleList[i][2]<<"}\",";
            stream<<"\"{"<<midPDisList[i][0]*tempScale<<", "<<midPDisList[i][1]*tempScale<<", "<<midPDisList[i][2]*tempScale<<"}\",";
            stream<<"\"{"<<helixLenList[i][0]*tempScale<<", "<<helixLenList[i][1]*tempScale<<", "<<helixLenList[i][2]*tempScale<<"}\"";
        }

        // Actually write into the file
        outputFile<<stream.str();
        outputFile.close();

        // Free the memory
        delete[] accuracyVolume;

        // Compute elapsed time
        elapsedTime = cpuTimer.Tick();
        std::cout<<"Done with "<<i+1<<"th sample set. It takes "<<elapsedTime*(1e-3)<<"s. "<<nThreads<<" OpenMP threads were used.\n";

        // Move to the next sample set
        ++i;
    }

    std::cout<<"Finish processing the batch!\n";
}

// Compute the accuracy volume for simple greedy
void RunForAccuracyVolumeSGBADL()
{
    // Wrap all this information into a class. ToDo.

    std::cout<<"Begin to compute acuracy volume for ADL (SGB)...\n";

    std::string inputFile = "../Data/AccuracyBatch.txt";
    std::string dataDir   = "../Data/";
    std::string fileExt   = ".pdb";

    std::string outPutType    = "volume";
    std::string outPutDirGau  = "../Data/SimpleGreedy/GaussianAccuracy/";
    std::string outPutDirEpan = "../Data/SimpleGreedy/EpanAccuracy/";
    std::string outPutExt     = ".csv";

    std::vector< std::vector<std::string> > fileList;
    std::vector< std::vector<int> >         matchTruthList;
    std::vector< int >                      scoreFunList;
    std::vector< std::vector<double> >      angleList;
    std::vector< std::vector<double> >      midPDisList;
    std::vector< std::vector<double> >      helixLenList;

    // Read in all the batch data from AccuracyBatch.txt
    // The information in the file has to be line by line, in the order above.
    // The angles in the AccuracyBatch file are measured by degree
    std::ifstream file(inputFile);
    std::string line;

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ':');

            if("fileName" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<std::string> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( content[i] );
                }
                fileList.push_back(temp);
            }
            else if("matchTruth" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<int> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stoi(content[i]) );
                }
                matchTruthList.push_back(temp);
            }
            else if("scoreFun" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                scoreFunList.push_back(std::stoi(content[0]));
            }
            else if("angle" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( DToR( std::stod(content[i]) ) );
                }
                angleList.push_back(temp);
            }
            else if("midPointDistance" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                midPDisList.push_back(temp);
            }
            else if("helixLenght" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                helixLenList.push_back(temp);
            }            
        }

        std::cout<<"Finish reading in the file!\n";
    }

    file.close();

    unsigned int i = 0;
    while(i < fileList.size())
    {
        std::cout<<"Begin to process the "<<i+1<<"th sample set...\n";

        //// Time it
        cpuTimer.Start();

        //// Read in file and retrieve geometry
        // Read in module set
        std::string fullPath = dataDir + fileList[i][0] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        char chainID = *(fileList[i][1].c_str());
        RetrieveHelixGeometry(moduleHelix, pdbReader, chainID, moduleOffset);

        // Read in data set
        fullPath = dataDir + fileList[i][2] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        chainID = *(fileList[i][3].c_str());
        RetrieveHelixGeometry(dataHelix, pdbReader, chainID, dataOffset);
        // Make sure module set has more helices
        SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

        // Initialize the accuracy volume --- at least 1x1x1 volume        
        int angleNum    = static_cast<int>( (angleList[i][1]    - angleList[i][0])    / angleList[i][2]    + indexEpsilon ) + 1;
        int midPDisNum  = static_cast<int>( (midPDisList[i][1]  - midPDisList[i][0])  / midPDisList[i][2]  + indexEpsilon ) + 1;
        int helixLenNum = static_cast<int>( (helixLenList[i][1] - helixLenList[i][0]) / helixLenList[i][2] + indexEpsilon ) + 1;

        scoreFunFlag = scoreFunList[i];
        //std::cout<<angleNum<<", "<<midPDisNum<<", "<<helixLenNum<<"\n";

        // Allocate memory to hold the volume
        float *accuracyVolume = new float[angleNum*midPDisNum*helixLenNum];

        int angleOffset   = midPDisNum*helixLenNum;
        int midPDisOffset = helixLenNum;

        int nThreads = 0;
        int a = 0; int d = 0; int l = 0;
        //unsigned int h = 0;
        //#pragma omp parallel default(shared) private(a,d,l,thresholds,scoreMatrixDLD,smSolverDLD)
        #pragma omp parallel default(shared) private(a,d,l,scoreMatrixDLD,smSolverDLD)
        {
            #pragma omp master
            nThreads = omp_get_num_threads();

            // Enter the loop to compute the accuracy volume --- #pragma omp parallel for
            //#pragma omp for
            //#pragma omp parallel for private(scoreMatrixDLD,smSolverDLD)       
            #pragma omp for        
            for(a = 0; a < angleNum; ++a)
                for(d = 0; d < midPDisNum; ++d)
                    for(l = 0; l < helixLenNum; ++l)
                    {
                        std::vector<double> tempThresholds(7,0.0);
                        tempThresholds[0] = 0; // wf
                        tempThresholds[3] = midPDisList[i][0] + d*midPDisList[i][2];   // td
                        tempThresholds[4] = angleList[i][0] + a*angleList[i][2];       // ta
                        tempThresholds[5] = angleList[i][0] + a*angleList[i][2];       // ta
                        tempThresholds[6] = helixLenList[i][0] + l*helixLenList[i][2]; // tl

                        std::vector<int> tempMatchTruth = matchTruthList[i];

                        // Update the match truth based on different parameter setting. 
                        for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
                        {
                            if( 0 != tempMatchTruth[h] )
                            {
                                // The truth has the same format as in Mathematica --- index starts from 1
                                float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                                float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                                if( fabs(moduleLen - datalen) > tempThresholds[6] && tempThresholds[6] > 0 )
                                    tempMatchTruth[h] = 0;
                            }
                        }

                        //std::cout<<"Thresholds: "<<thresholds[0]<<", "<<thresholds[3]<<", "<<thresholds[4]<<", "<<thresholds[5]<<", "<<thresholds[6]<<"\n";
                        // Compute accuracy
                        ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, tempThresholds, scoreFunFlag);
                        // Solve for match result --- choose Simple Greedy Binarize method
                        SolveForResultDLD(smSolverDLD, scoreMatrixDLD, 1);
                        // Compute and record accuracy
                        accuracyVolume[a*angleOffset + d*midPDisOffset + l] = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
                    }// End for loop

            #pragma omp barrier
        } // End Parallel Block

        //// Export the result to the disk.
        // Make output name and choose output folder

        std::string exportFile = "";
        if(0 == scoreFunList[i])
        {
            exportFile = outPutDirGau + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + outPutExt;
        }
        else if(1 == scoreFunList[i])
        {
            exportFile = outPutDirEpan + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + outPutExt;
        }

        std::ofstream outputFile(exportFile, std::ios::trunc);

        // Use to convert float to string
        std::stringstream stream;
        // Initialize
        stream.str(std::string());

        if( !outputFile.is_open() )
            std::cout<<"Unable to open "<<exportFile<<"...\n";
        else
        {
            // Output volume
            for(int a = 0; a < angleNum; ++a)
            {
                // angle dimension
                stream<<"\"{";
                for(int d = 0; d < midPDisNum; ++d)
                {
                    // midpoint distance dimension
                    stream<<"{";
                    for(int l = 0; l < helixLenNum; ++l)
                    {
                        // helix length dimension
                        stream<<accuracyVolume[a*angleOffset + d*midPDisOffset + l];
                        if(l < helixLenNum - 1)
                            stream<<", ";
                    }

                    stream<<"}";
                    if(d < midPDisNum - 1)
                        stream<<", ";
                }

                stream<<"}\"";
                if(a < angleNum - 1)
                    stream<<",";
            }

            stream<<"\n";

            // Output Parameters --- a, d, l
            // To match for unit in Mathematica, scale by 0.01
            float tempScale = 0.01f;
            stream<<"\"{"<<angleList[i][0]<<", "<<angleList[i][1]<<", "<<angleList[i][2]<<"}\",";
            stream<<"\"{"<<midPDisList[i][0]*tempScale<<", "<<midPDisList[i][1]*tempScale<<", "<<midPDisList[i][2]*tempScale<<"}\",";
            stream<<"\"{"<<helixLenList[i][0]*tempScale<<", "<<helixLenList[i][1]*tempScale<<", "<<helixLenList[i][2]*tempScale<<"}\"";
        }

        // Actually write into the file
        outputFile<<stream.str();
        outputFile.close();

        // Free the memory
        delete[] accuracyVolume;

        // Compute elapsed time
        elapsedTime = cpuTimer.Tick();
        std::cout<<"Done with "<<i+1<<"th sample set. It takes "<<elapsedTime*(1e-3)<<"s. "<<nThreads<<" OpenMP threads were used.\n";

        // Move to the next sample set
        ++i;
    }

    std::cout<<"Finish processing the batch!\n";
}


// Compute the accuracy volume (angle, midpoint distance and helix length)
// --- wrap this script into a class. ToDo.
// Also, wrap run for one pair script into a class. ToDo.
void RunForAccuracyVolumeADL()
{
    // Wrap all this information into a class. ToDo. Flag Ha.
    
    std::cout<<"Begin to compute acuracy volume for ADL (IPB)...\n";

    std::string inputFile = "../Data/AccuracyBatch.txt";
    std::string dataDir   = "../Data/";
    std::string fileExt   = ".pdb";

    std::string outPutType    = "volume";
    std::string outPutDirGau  = "../Data/IterativePicking/GaussianAccuracy/";
    std::string outPutDirEpan = "../Data/IterativePicking/EpanAccuracy/";
    std::string outPutExt     = ".csv";

    std::vector< std::vector<std::string> > fileList;
    std::vector< std::vector<int> >         matchTruthList;
    std::vector< int >                      scoreFunList;
    std::vector< std::vector<double> >      angleList;
    std::vector< std::vector<double> >      midPDisList;
    std::vector< std::vector<double> >      helixLenList;

    // Read in all the batch data from AccuracyBatch.txt
    // The information in the file has to be line by line, in the order above.
    // The angles in the AccuracyBatch file are measured by degree
    std::ifstream file(inputFile);
    std::string line;

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ':');

            if("fileName" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<std::string> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( content[i] );
                }
                fileList.push_back(temp);
            }
            else if("matchTruth" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<int> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stoi(content[i]) );
                }
                matchTruthList.push_back(temp);
            }
            else if("scoreFun" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                scoreFunList.push_back(std::stoi(content[0]));
            }
            else if("angle" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( DToR( std::stod(content[i]) ) );
                }
                angleList.push_back(temp);
            }
            else if("midPointDistance" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                midPDisList.push_back(temp);
            }
            else if("helixLenght" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                helixLenList.push_back(temp);
            }            
        }

        std::cout<<"Finish reading in the file!\n";
    }

    file.close();        

    unsigned int i = 0;
    while(i < fileList.size())
    {
        std::cout<<"Begin to process the "<<i+1<<"th sample set...\n";

        //// Time it
        cpuTimer.Start();

        //// Read in file and retrieve geometry
        // Read in module set
        std::string fullPath = dataDir + fileList[i][0] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        char chainID = *(fileList[i][1].c_str());
        RetrieveHelixGeometry(moduleHelix, pdbReader, chainID, moduleOffset);

        // Read in data set
        fullPath = dataDir + fileList[i][2] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        chainID = *(fileList[i][3].c_str());
        RetrieveHelixGeometry(dataHelix, pdbReader, chainID, dataOffset);
        // Make sure module set has more helices
        SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

        // Initialize the accuracy volume --- at least 1x1x1 volume        
        int angleNum    = static_cast<int>( (angleList[i][1]    - angleList[i][0])    / angleList[i][2]    + indexEpsilon ) + 1;
        int midPDisNum  = static_cast<int>( (midPDisList[i][1]  - midPDisList[i][0])  / midPDisList[i][2]  + indexEpsilon ) + 1;
        int helixLenNum = static_cast<int>( (helixLenList[i][1] - helixLenList[i][0]) / helixLenList[i][2] + indexEpsilon ) + 1;

        scoreFunFlag = scoreFunList[i];
        //std::cout<<angleNum<<", "<<midPDisNum<<", "<<helixLenNum<<"\n";

        // Allocate memory to hold the volume
        float *accuracyVolume = new float[angleNum*midPDisNum*helixLenNum];

        int angleOffset   = midPDisNum*helixLenNum;
        int midPDisOffset = helixLenNum;

        int nThreads = 0;
        //int askedThreads = 2;
       //omp_set_num_threads(askedThreads);

        int a = 0; int d = 0; int l = 0;
        ScoreMatrixDoubleLineDir scoreMatrixDLD;
        SMSolverDLD                 smSolverDLD;
        //unsigned int h = 0;
        //#pragma omp parallel default(shared) private(a,d,l,thresholds,scoreMatrixDLD,smSolverDLD)
        #pragma omp parallel default(shared) private(a,d,l,scoreMatrixDLD,smSolverDLD)
        {
            #pragma omp master
            nThreads = omp_get_num_threads();

            // Enter the loop to compute the accuracy volume --- #pragma omp parallel for
            //#pragma omp for
            //#pragma omp parallel for private(scoreMatrixDLD,smSolverDLD)       
            #pragma omp for
            for(a = 0; a < angleNum; ++a)
                for(d = 0; d < midPDisNum; ++d)
                    for(l = 0; l < helixLenNum; ++l)
                    {
                        // Set up new parameters --- No distance weighting so far
                        //thresholds[0] = 0; // wf
                        //thresholds[3] = midPDisList[i][0] + d*midPDisList[i][2];   // td
                        //thresholds[4] = angleList[i][0] + a*angleList[i][2];       // ta
                        //thresholds[5] = thresholds[4]; // ta
                        //thresholds[6] = helixLenList[i][0] + l*helixLenList[i][2]; // tl

                        std::vector<double> tempThresholds(7,0.0);
                        tempThresholds[0] = 0; // wf
                        tempThresholds[3] = midPDisList[i][0] + d*midPDisList[i][2];   // td
                        tempThresholds[4] = angleList[i][0] + a*angleList[i][2];       // ta
                        tempThresholds[5] = angleList[i][0] + a*angleList[i][2];       // ta
                        tempThresholds[6] = helixLenList[i][0] + l*helixLenList[i][2]; // tl

                        std::vector<int> tempMatchTruth = matchTruthList[i];

                        // Update the match truth based on different parameter setting. 
                        for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
                        {
                            if( 0 != tempMatchTruth[h] )
                            {
                                // The truth has the same format as in Mathematica --- index starts from 1
                                float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                                float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                                if( fabs(moduleLen - datalen) > tempThresholds[6] && tempThresholds[6] > 0 )
                                    tempMatchTruth[h] = 0;
                            }
                        }

                        //std::cout<<"Thresholds: "<<thresholds[0]<<", "<<thresholds[3]<<", "<<thresholds[4]<<", "<<thresholds[5]<<", "<<thresholds[6]<<"\n";
                        // Compute accuracy
                        ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, tempThresholds, scoreFunFlag);
                        // Solve for match result
                        SolveForResultDLD(smSolverDLD, scoreMatrixDLD);
                        // Compute and record accuracy
                        accuracyVolume[a*angleOffset + d*midPDisOffset + l] = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
                        //if(a == 0 && d == 0 && l == 1)
                        //{
                        //    //scoreMatrixDLD.GetMatrix();
                        //    std::cout<<"Matrix Sum: "<<scoreMatrixDLD.GetMatrixSum()<<"\n";

                        //    for(auto it = tempMatchTruth.begin(); it != tempMatchTruth.end(); ++it)
                        //    {
                        //        std::cout<<*it<<", ";
                        //    }
                        //    std::cout<<"\n";

                        //    for(auto it = scoreMatrixDLD.GetParaList().begin(); it != scoreMatrixDLD.GetParaList().end(); ++it)
                        //    {
                        //        std::cout<<*it<<", ";
                        //    }
                        //    std::cout<<"\n";

                        //    std::cout<<"para: "<<tempThresholds[4]<<", "<<tempThresholds[3]<<", "<<tempThresholds[6]<<"\n";

                        //    for(auto it = smSolverDLD.GetMatchResult().begin(); it != smSolverDLD.GetMatchResult().end(); ++it)
                        //    {
                        //        std::cout<<*it<<", ";
                        //    }
                        //    std::cout<<"\n";

                        //    for(auto it = smSolverDLD.GetMatchFlag().begin(); it != smSolverDLD.GetMatchFlag().end(); ++it)
                        //    {
                        //        std::cout<<*it<<", ";
                        //    }
                        //    std::cout<<"\n";

                        //    std::cout<<"value: "<<accuracyVolume[a*angleOffset + d*midPDisOffset + l]<<"\n";
                        //}

                    }// End for loop

            #pragma omp barrier
        } // End Parallel Block

        //// Export the result to the disk.
        // Make output name and choose output folder
        std::string exportFile = "";
        if(0 == scoreFunList[i])
        {
            exportFile = outPutDirGau + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + outPutExt;
        }
        else if(1 == scoreFunList[i])
        {
            exportFile = outPutDirEpan + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + outPutExt;
        }
        
        std::ofstream outputFile(exportFile, std::ios::trunc);

        // Use to convert float to string
        std::stringstream stream;
        // Initialize
        stream.str(std::string());

        if( !outputFile.is_open() )
            std::cout<<"Unable to open "<<exportFile<<"...\n";
        else
        {
            // Output volume
            for(int a = 0; a < angleNum; ++a)
            {
                // angle dimension
                stream<<"\"{";
                for(int d = 0; d < midPDisNum; ++d)
                {
                    // midpoint distance dimension
                    stream<<"{";
                    for(int l = 0; l < helixLenNum; ++l)
                    {
                        // helix length dimension
                        stream<<accuracyVolume[a*angleOffset + d*midPDisOffset + l];
                        if(l < helixLenNum - 1)
                            stream<<", ";
                    }

                    stream<<"}";
                    if(d < midPDisNum - 1)
                        stream<<", ";
                }

                stream<<"}\"";
                if(a < angleNum - 1)
                    stream<<",";
            }

            stream<<"\n";

            // Output Parameters --- a, d, l
            // To match for unit in Mathematica, scale by 0.01
            float tempScale = 0.01f;
            stream<<"\"{"<<angleList[i][0]<<", "<<angleList[i][1]<<", "<<angleList[i][2]<<"}\",";
            stream<<"\"{"<<midPDisList[i][0]*tempScale<<", "<<midPDisList[i][1]*tempScale<<", "<<midPDisList[i][2]*tempScale<<"}\",";
            stream<<"\"{"<<helixLenList[i][0]*tempScale<<", "<<helixLenList[i][1]*tempScale<<", "<<helixLenList[i][2]*tempScale<<"}\"";
        }

        // Actually write into the file
        outputFile<<stream.str();
        outputFile.close();

        // Free the memory
        delete[] accuracyVolume;

        // Compute elapsed time
        elapsedTime = cpuTimer.Tick();
        std::cout<<"Done with "<<i+1<<"th sample set. It takes "<<elapsedTime*(1e-3)<<"s. "<<nThreads<<" OpenMP threads were used.\n";

        // Move to the next sample set
        ++i;
    }

    std::cout<<"Finish processing the batch!\n";
}


// Compute the accuracy volume (distance weight, midpoint distance and helix length)
// Fix angle. Absolute distance weight
void RunForAccuracyVolumeDWDL()
{
    // Wrap all this information into a class. ToDo.

    std::cout<<"Begin to compute acuracy volume for DWDL (IPB)...\n";

    std::string inputFile = "../Data/AccuracyBatchDW.txt";
    std::string dataDir   = "../Data/";
    std::string fileExt   = ".pdb";

    std::string outPutType    = "volume_dw";
    std::string outPutDirGau  = "../Data/IterativePicking/GaussianAccuracyDW/";
    std::string outPutDirEpan = "../Data/IterativePicking/EpanAccuracyDW/";
    std::string outPutExt     = ".csv";

    std::vector< std::vector<std::string> > fileList;
    std::vector< std::vector<int> >         matchTruthList;
    std::vector< int >                      scoreFunList;
    // Now the angle list serves as distance weight list
    std::vector< std::vector<double> >      angleList;
    std::vector< std::vector<double> >      midPDisList;
    std::vector< std::vector<double> >      helixLenList;

    // Read in all the batch data from AccuracyBatch.txt
    // The information in the file has to be line by line, in the order above.
    // The angles in the AccuracyBatch file are measured by degree
    std::ifstream file(inputFile);
    std::string line;

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ':');

            if("fileName" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<std::string> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( content[i] );
                }
                fileList.push_back(temp);
            }
            else if("matchTruth" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<int> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stoi(content[i]) );
                }
                matchTruthList.push_back(temp);
            }
            else if("scoreFun" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                scoreFunList.push_back(std::stoi(content[0]));
            }
            else if("angle" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    // Angle means distance weight threshold here
                    temp.push_back( std::stod(content[i]) );
                }
                angleList.push_back(temp);
            }
            else if("midPointDistance" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                midPDisList.push_back(temp);
            }
            else if("helixLenght" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                helixLenList.push_back(temp);
            }            
        }

        std::cout<<"Finish reading in the file!\n";
    }

    file.close();

    double fixAngle = 1.466078; // 1.466078 (84.0001 degree)
    unsigned int i = 0;
    while(i < fileList.size())
    {
        std::cout<<"Begin with the "<<i+1<<"th sample set...\n";

        //// Time it
        cpuTimer.Start();

        //// Read in file and retrieve geometry
        // Read in module set
        std::string fullPath = dataDir + fileList[i][0] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        char chainID = *(fileList[i][1].c_str());
        RetrieveHelixGeometry(moduleHelix, pdbReader, chainID, moduleOffset);

        // Read in data set
        fullPath = dataDir + fileList[i][2] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        chainID = *(fileList[i][3].c_str());
        RetrieveHelixGeometry(dataHelix, pdbReader, chainID, dataOffset);
        // Make sure module set has more helices
        SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

        // Compute the bounding box
        //ComputeBBD(inputModuleHelix, moduleBBD);
        //ComputeBBD(inputDataHelix, dataBBD);
        //std::cout<<"module and data BBD: "<<moduleBBD<<", "<<dataBBD<<"\n";

        // Initialize the accuracy volume --- at least 1x1x1 volume        
        int angleNum    = static_cast<int>( (angleList[i][1]    - angleList[i][0])    / angleList[i][2]    + indexEpsilon ) + 1;
        int midPDisNum  = static_cast<int>( (midPDisList[i][1]  - midPDisList[i][0])  / midPDisList[i][2]  + indexEpsilon ) + 1;
        int helixLenNum = static_cast<int>( (helixLenList[i][1] - helixLenList[i][0]) / helixLenList[i][2] + indexEpsilon ) + 1;

        scoreFunFlag = scoreFunList[i];
        //std::cout<<angleNum<<", "<<midPDisNum<<", "<<helixLenNum<<"\n";

        // Allocate memory to hold the volume
        float *accuracyVolume = new float[angleNum*midPDisNum*helixLenNum];

        int angleOffset   = midPDisNum*helixLenNum;
        int midPDisOffset = helixLenNum;

        //double angleStart = angleList[i][0];
        //double angleStep  = angleList[i][2];

        int nThreads = 0;
        //int askedThreads = 2;
        //omp_set_num_threads(askedThreads);

        int a = 0; int d = 0; int l = 0;
        //unsigned int h = 0;
        //#pragma omp parallel default(shared) private(a,d,l,thresholds,scoreMatrixDLD,smSolverDLD)
        #pragma omp parallel default(shared) private(a,d,l,scoreMatrixDLD,smSolverDLD)
        {
            #pragma omp master
            nThreads = omp_get_num_threads();

            // Enter the loop to compute the accuracy volume --- #pragma omp parallel for
            //#pragma omp for
            //#pragma omp parallel for private(scoreMatrixDLD,smSolverDLD)       
            #pragma omp for
            // Enter the loop to compute the accuracy volume
            for(a = 0; a < angleNum; ++a)
                for(d = 0; d < midPDisNum; ++d)
                    for(l = 0; l < helixLenNum; ++l)
                    {
                        std::vector<double> tempThresholds(7,0.0);
                        // Set up new parameters --- No distance weighting so far
                        tempThresholds[0] = 1.0; // wf    angleList[i][0] + a*angleList[i][2]
                        //thresholds[1] = moduleBBD*(angleStart + a*angleStep);
                        //thresholds[2] = dataBBD*(angleStart + a*angleStep);
                        tempThresholds[1] = angleList[i][0] + a*angleList[i][2];
                        tempThresholds[2] = angleList[i][0] + a*angleList[i][2];
                        tempThresholds[3] = midPDisList[i][0] + d*midPDisList[i][2];   // td
                        tempThresholds[4] = fixAngle; // ta
                        tempThresholds[5] = fixAngle; // ta
                        tempThresholds[6] = helixLenList[i][0] + l*helixLenList[i][2]; // tl

                        std::vector<int> tempMatchTruth = matchTruthList[i];

                        // Update the match truth based on different parameter setting.
                        for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
                        {
                            if( 0 != tempMatchTruth[h] )
                            {
                                // The truth has the same format as in Mathematica --- index starts from 1
                                float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                                float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                                if( fabs(moduleLen - datalen) > tempThresholds[6] && tempThresholds[6] > 0 )
                                    tempMatchTruth[h] = 0;
                            }
                        }

                        //std::cout<<"Thresholds: "<<thresholds[0]<<", "<<thresholds[3]<<", "<<thresholds[4]<<", "<<thresholds[5]<<", "<<thresholds[6]<<"\n";
                        // Compute accuracy
                        ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, tempThresholds, scoreFunFlag);
                        // Solve for match result
                        SolveForResultDLD(smSolverDLD, scoreMatrixDLD);
                        // Compute and record accuracy
                        accuracyVolume[a*angleOffset + d*midPDisOffset + l] = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
                    }// End for loop
            #pragma omp barrier
        } // End Parallel Block

        //// Export the result to the disk.
        // Make output name and choose output folder

        std::string exportFile = "";
        if(0 == scoreFunList[i])
        {
            exportFile = outPutDirGau + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + outPutExt;
        }
        else if(1 == scoreFunList[i])
        {
            exportFile = outPutDirEpan + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + outPutExt;
        }

        std::ofstream outputFile(exportFile, std::ios::trunc);

        // Use to convert float to string
        std::stringstream stream;
        // Initialize
        stream.str(std::string());

        if( !outputFile.is_open() )
            std::cout<<"Unable to open "<<exportFile<<"...\n";
        else
        {
            // Output volume
            for(int a = 0; a < angleNum; ++a)
            {
                // angle dimension
                stream<<"\"{";
                for(int d = 0; d < midPDisNum; ++d)
                {
                    // midpoint distance dimension
                    stream<<"{";
                    for(int l = 0; l < helixLenNum; ++l)
                    {
                        // helix length dimension
                        stream<<accuracyVolume[a*angleOffset + d*midPDisOffset + l];
                        if(l < helixLenNum - 1)
                            stream<<", ";
                    }

                    stream<<"}";
                    if(d < midPDisNum - 1)
                        stream<<", ";
                }

                stream<<"}\"";
                if(a < angleNum - 1)
                    stream<<",";
            }

            stream<<"\n";

            // Output Parameters --- a, d, l
            // To match for unit in Mathematica, scale by 0.01
            float tempScale = 0.01f;
            stream<<"\"{"<<angleList[i][0]<<", "<<angleList[i][1]<<", "<<angleList[i][2]<<"}\",";
            stream<<"\"{"<<midPDisList[i][0]*tempScale<<", "<<midPDisList[i][1]*tempScale<<", "<<midPDisList[i][2]*tempScale<<"}\",";
            stream<<"\"{"<<helixLenList[i][0]*tempScale<<", "<<helixLenList[i][1]*tempScale<<", "<<helixLenList[i][2]*tempScale<<"}\"";
        }

        // Actually write into the file
        outputFile<<stream.str();
        outputFile.close();

        // Free the memory
        delete[] accuracyVolume;

        // Compute elapsed time
        elapsedTime = cpuTimer.Tick();
        std::cout<<"Done with "<<i+1<<"th sample set. It takes "<<elapsedTime*(1e-3)<<"s. "<<nThreads<<" OpenMP threads were used.\n";


        // Move to the next sample set
        ++i;
    }

    std::cout<<"Finish processing the batch!\n";
}


// Compute the accuracy volume (distance weight, midpoint distance and helix length)
// Fix angle. BBD ratio distance weight
void RunForAccuracyVolumeDWDLBBD()
{
    // Wrap all this information into a class. ToDo.

    std::cout<<"Begin to compute acuracy volume for DWDL (IPB)...\n";

    std::string inputFile = "../Data/AccuracyBatchDW.txt";
    std::string dataDir   = "../Data/";
    std::string fileExt   = ".pdb";

    std::string outPutType    = "volume_dw";
    std::string outPutDirGau  = "../Data/IterativePicking/GaussianAccuracyDW/";
    std::string outPutDirEpan = "../Data/IterativePicking/EpanAccuracyDW/";
    std::string outPutExt     = ".csv";

    std::vector< std::vector<std::string> > fileList;
    std::vector< std::vector<int> >         matchTruthList;
    std::vector< int >                      scoreFunList;
    // Now the angle list serves as distance weight list
    std::vector< std::vector<double> >      angleList;
    std::vector< std::vector<double> >      midPDisList;
    std::vector< std::vector<double> >      helixLenList;

    // Read in all the batch data from AccuracyBatch.txt
    // The information in the file has to be line by line, in the order above.
    // The angles in the AccuracyBatch file are measured by degree
    std::ifstream file(inputFile);
    std::string line;

    if( !file.is_open() )
        std::cout<<"Unable to open "<<inputFile<<"...\n";
    else
    {
        std::cout<<"Read in the file...\n";

        while( std::getline(file, line) )
        {
            // Protect from empty lines
            if( "" == line )
                continue;

            std::vector<std::string> tokens;
            ParseByCharater(line, tokens, ':');

            if("fileName" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<std::string> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( content[i] );
                }
                fileList.push_back(temp);
            }
            else if("matchTruth" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<int> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stoi(content[i]) );
                }
                matchTruthList.push_back(temp);
            }
            else if("scoreFun" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                scoreFunList.push_back(std::stoi(content[0]));
            }
            else if("angle" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    // Angle means distance weight threshold here
                    temp.push_back( std::stod(content[i]) );
                }
                angleList.push_back(temp);
            }
            else if("midPointDistance" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                midPDisList.push_back(temp);
            }
            else if("helixLenght" == tokens[0])
            {
                // Parse each parameter out by space
                std::vector<std::string> content;
                ParseByCharater(tokens[1], content, ' ');

                std::vector<double> temp;
                for(int i = 0; i < (int)content.size(); ++i)
                {
                    temp.push_back( std::stod(content[i]) );
                }
                helixLenList.push_back(temp);
            }            
        }

        std::cout<<"Finish reading in the file!\n";
    }

    file.close();

    double fixAngle = 1.466078; // 1.466078 (84.0001 degree)
    unsigned int i = 0;
    while(i < fileList.size())
    {
        std::cout<<"Begin with the "<<i+1<<"th sample set...\n";

        //// Time it
        cpuTimer.Start();

        //// Read in file and retrieve geometry
        // Read in module set
        std::string fullPath = dataDir + fileList[i][0] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        char chainID = *(fileList[i][1].c_str());
        RetrieveHelixGeometry(moduleHelix, pdbReader, chainID, moduleOffset);

        // Read in data set
        fullPath = dataDir + fileList[i][2] + fileExt;
        ReadInTheFile(fullPath.c_str(), pdbReader);
        chainID = *(fileList[i][3].c_str());
        RetrieveHelixGeometry(dataHelix, pdbReader, chainID, dataOffset);
        // Make sure module set has more helices
        SwtichModuleData(inputModuleHelix, inputDataHelix, &moduleHelix, &dataHelix);

        // Compute the bounding box
        double moduleBBD = 0.0;
        double dataBBD   = 0.0;
        ComputeBBD(inputModuleHelix, moduleBBD);
        ComputeBBD(inputDataHelix, dataBBD);
        //std::cout<<"module and data BBD: "<<moduleBBD<<", "<<dataBBD<<"\n";

        // Initialize the accuracy volume --- at least 1x1x1 volume        
        int angleNum    = static_cast<int>( (angleList[i][1]    - angleList[i][0])    / angleList[i][2]    + indexEpsilon ) + 1;
        int midPDisNum  = static_cast<int>( (midPDisList[i][1]  - midPDisList[i][0])  / midPDisList[i][2]  + indexEpsilon ) + 1;
        int helixLenNum = static_cast<int>( (helixLenList[i][1] - helixLenList[i][0]) / helixLenList[i][2] + indexEpsilon ) + 1;

        scoreFunFlag = scoreFunList[i];
        //std::cout<<angleNum<<", "<<midPDisNum<<", "<<helixLenNum<<"\n";

        // Allocate memory to hold the volume
        float *accuracyVolume = new float[angleNum*midPDisNum*helixLenNum];

        int angleOffset   = midPDisNum*helixLenNum;
        int midPDisOffset = helixLenNum;

        double angleStart = angleList[i][0];
        double angleStep  = angleList[i][2];

        int nThreads = 0;
        //int askedThreads = 2;
        //omp_set_num_threads(askedThreads);

        int a = 0; int d = 0; int l = 0;
        //unsigned int h = 0;
        //#pragma omp parallel default(shared) private(a,d,l,thresholds,scoreMatrixDLD,smSolverDLD)
        #pragma omp parallel default(shared) private(a,d,l,scoreMatrixDLD,smSolverDLD)
        {
            #pragma omp master
            nThreads = omp_get_num_threads();

            // Enter the loop to compute the accuracy volume --- #pragma omp parallel for
            //#pragma omp for
            //#pragma omp parallel for private(scoreMatrixDLD,smSolverDLD)       
            #pragma omp for
            // Enter the loop to compute the accuracy volume
            for(a = 0; a < angleNum; ++a)
                for(d = 0; d < midPDisNum; ++d)
                    for(l = 0; l < helixLenNum; ++l)
                    {
                        std::vector<double> tempThresholds(7,0.0);
                        // Set up new parameters --- No distance weighting so far
                        tempThresholds[0] = 1.0; // wf
                        tempThresholds[1] = moduleBBD*(angleStart + a*angleStep);
                        tempThresholds[2] = dataBBD*(angleStart + a*angleStep);
                        tempThresholds[3] = midPDisList[i][0] + d*midPDisList[i][2];   // td
                        tempThresholds[4] = fixAngle; // ta
                        tempThresholds[5] = fixAngle; // ta
                        tempThresholds[6] = helixLenList[i][0] + l*helixLenList[i][2]; // tl

                        std::vector<int> tempMatchTruth = matchTruthList[i];

                        // Update the match truth based on different parameter setting.
                        for(unsigned int h = 0; h < tempMatchTruth.size(); ++h)
                        {
                            if( 0 != tempMatchTruth[h] )
                            {
                                // The truth has the same format as in Mathematica --- index starts from 1
                                float moduleLen = inputModuleHelix->GetLineSegment(tempMatchTruth[h]-1).GetLen();
                                float datalen   = inputDataHelix->GetLineSegment(h).GetLen();

                                if( fabs(moduleLen - datalen) > tempThresholds[6] && tempThresholds[6] > 0 )
                                    tempMatchTruth[h] = 0;
                            }
                        }

                        //std::cout<<"Thresholds: "<<thresholds[0]<<", "<<thresholds[3]<<", "<<thresholds[4]<<", "<<thresholds[5]<<", "<<thresholds[6]<<"\n";
                        // Compute accuracy
                        ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, tempThresholds, scoreFunFlag);
                        // Solve for match result
                        SolveForResultDLD(smSolverDLD, scoreMatrixDLD);
                        // Compute and record accuracy
                        accuracyVolume[a*angleOffset + d*midPDisOffset + l] = ComputeAccuracyDLD(tempMatchTruth, smSolverDLD.GetMatchResult(), smSolverDLD.GetMatchFlag());
                    }// End for loop
            #pragma omp barrier
        } // End Parallel Block

        //// Export the result to the disk.
        // Make output name and choose output folder

        std::string exportFile = "";
        if(0 == scoreFunList[i])
        {
            exportFile = outPutDirGau + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + outPutExt;
        }
        else if(1 == scoreFunList[i])
        {
            exportFile = outPutDirEpan + outPutType + "_" + fileList[i][0] + "_" + fileList[i][2] + outPutExt;
        }

        std::ofstream outputFile(exportFile, std::ios::trunc);

        // Use to convert float to string
        std::stringstream stream;
        // Initialize
        stream.str(std::string());

        if( !outputFile.is_open() )
            std::cout<<"Unable to open "<<exportFile<<"...\n";
        else
        {
            // Output volume
            for(int a = 0; a < angleNum; ++a)
            {
                // angle dimension
                stream<<"\"{";
                for(int d = 0; d < midPDisNum; ++d)
                {
                    // midpoint distance dimension
                    stream<<"{";
                    for(int l = 0; l < helixLenNum; ++l)
                    {
                        // helix length dimension
                        stream<<accuracyVolume[a*angleOffset + d*midPDisOffset + l];
                        if(l < helixLenNum - 1)
                            stream<<", ";
                    }

                    stream<<"}";
                    if(d < midPDisNum - 1)
                        stream<<", ";
                }

                stream<<"}\"";
                if(a < angleNum - 1)
                    stream<<",";
            }

            stream<<"\n";

            // Output Parameters --- a, d, l
            // To match for unit in Mathematica, scale by 0.01
            float tempScale = 0.01f;
            stream<<"\"{"<<angleList[i][0]<<", "<<angleList[i][1]<<", "<<angleList[i][2]<<"}\",";
            stream<<"\"{"<<midPDisList[i][0]*tempScale<<", "<<midPDisList[i][1]*tempScale<<", "<<midPDisList[i][2]*tempScale<<"}\",";
            stream<<"\"{"<<helixLenList[i][0]*tempScale<<", "<<helixLenList[i][1]*tempScale<<", "<<helixLenList[i][2]*tempScale<<"}\"";
        }

        // Actually write into the file
        outputFile<<stream.str();
        outputFile.close();

        // Free the memory
        delete[] accuracyVolume;

        // Compute elapsed time
        elapsedTime = cpuTimer.Tick();
        std::cout<<"Done with "<<i+1<<"th sample set. It takes "<<elapsedTime*(1e-3)<<"s. "<<nThreads<<" OpenMP threads were used.\n";


        // Move to the next sample set
        ++i;
    }

    std::cout<<"Finish processing the batch!\n";
}


/**************************************************************************************************************************/
/************************************************* Functions **************************************************************/
/**************************************************************************************************************************/

// Compute Bounding box diagonal
void ComputeBBD(HelixGeometry* helixGeo, double & diagonal)
{
    // Get the first point
    Vector3f tempPoint = (helixGeo->GetLineSegList()).at(0).GetP0();
    Vector3f maxPoint = tempPoint;
    Vector3f minPoint = tempPoint;

    // Find the AABB bounding box
    for(auto it = helixGeo->GetLineSegList().begin(); it != helixGeo->GetLineSegList().end(); ++it)
    {
        tempPoint = it->GetP0();
        maxPoint[0] = max(maxPoint[0], tempPoint[0]);
        maxPoint[1] = max(maxPoint[1], tempPoint[1]);
        maxPoint[2] = max(maxPoint[2], tempPoint[2]);
        minPoint[0] = min(minPoint[0], tempPoint[0]);
        minPoint[1] = min(minPoint[1], tempPoint[1]);
        minPoint[2] = min(minPoint[2], tempPoint[2]);
        tempPoint = it->GetP1();
        maxPoint[0] = max(maxPoint[0], tempPoint[0]);
        maxPoint[1] = max(maxPoint[1], tempPoint[1]);
        maxPoint[2] = max(maxPoint[2], tempPoint[2]);
        minPoint[0] = min(minPoint[0], tempPoint[0]);
        minPoint[1] = min(minPoint[1], tempPoint[1]);
        minPoint[2] = min(minPoint[2], tempPoint[2]);
    }

    // Compute the diagonal
    diagonal = EuDis(maxPoint, minPoint);
}

// Given match truth and match result, compute the correctness rate.
// For double line direction score matrix
float ComputeAccuracyDLD(const std::vector<int> & matchTruth, const std::vector<int> & matchResult, const std::vector<int> & m_matchFlag)
{
    int rightMatch = 0;
    int totalMatch = 0;

    for(unsigned int i = 0; i < matchResult.size(); ++i)
    {
        if( 0 != matchTruth[i] )
        {
            ++totalMatch;
            //if( matchTruth[i] == matchResult[i] && 0 == m_matchFlag[i] )
            //    ++rightMatch;
            // Try the accuracy not including flag --- debug
            if( matchTruth[i] == matchResult[i] )
                ++rightMatch;
        }
    }

    // The accuracy is the correctness ratio.
    return (float)rightMatch/(float)totalMatch;
}

// For single line direction score matrix
float ComputeAccuracySLD(const std::vector<int> & matchTruth, const std::vector<int> & matchResult)
{
    int rightMatch = 0;
    int totalMatch = 0;

    for(unsigned int i = 0; i < matchResult.size(); ++i)
    {
        if( 0 != matchTruth[i] )
        {
            ++totalMatch;
            if( matchTruth[i] == matchResult[i] )
                ++rightMatch;
        }
    }

    // The accuracy is the correctness ratio.
    return (float)rightMatch/(float)totalMatch;
}

void SetParameterForScoreMatrix(std::vector<double> & thresholds, std::string parameters, int & scoreFunFlag, int sfFlag)
{
    // Parse each parameter out by space
    std::vector<std::string> tokens;
    ParseByCharater(parameters, tokens, ' ');
    // Set up all the thresholds
    for(int i = 0; i < (int)tokens.size(); ++i)
    {
        thresholds[i] = std::stod(tokens[i]);
    }

    // Set up the score function
    scoreFunFlag = sfFlag;
}


// Dump match result (including flag) to disk for Mathematica
void ExportMatchResult()
{
    // Output .csv file
    std::ofstream outputFile(exportResultFileName, std::ios::trunc);

    if( !outputFile.is_open() )
        std::cout<<"Unable to open "<<exportResultFileName<<"...\n";
    else
    {
        // Use to convert float to string
        std::stringstream stream;
        // Initialize
        stream.str(std::string());

        /** Output match result **/
        for(auto it = smSolverDLD.GetMatchResult().begin(); it != smSolverDLD.GetMatchResult().end()-1; ++it)
        {
            stream<<*it<<",";
        }
        auto lastIt = smSolverDLD.GetMatchResult().end()-1;
        stream<<*lastIt<<"\n";

        /** Output result flag **/
        for(auto it = smSolverDLD.GetMatchFlag().begin(); it != smSolverDLD.GetMatchFlag().end()-1; ++it)
        {
            stream<<*it<<",";
        }
        lastIt = smSolverDLD.GetMatchFlag().end()-1;
        stream<<*lastIt<<"\n";

        /** Write into the file **/
        outputFile<<stream.str();
        outputFile.close();
    }
}


// Dump module and data geometry to disk for Mathematica
// As a file: ModuleDataHelix.csv
void ExportModuleAndDataGeometry(HelixGeometry *inputModuleHelix, HelixGeometry *inputDataHelix, bool switchFlag)
{
    // Output .csv file
    std::ofstream outputFile(exportModelFileName, std::ios::trunc);

    if( !outputFile.is_open() )
        std::cout<<"Unable to open "<<exportModelFileName<<"...\n";
    else
    {
        // Use to convert float to string
        std::stringstream stream;
        // Set precision of float to string
        stream.precision(10);
        // Initialize
        stream.str(std::string());

        /** Output Information whether module and data are switched **/
        if(true == switchFlag)
            stream<<"Module and data have been switched.\n";
        else
            stream<<"Module and data have NOT been switched.\n";

        std::cout<<"Here O\n";

        /** Output module geometry **/
        int tc = 0;
        std::cout<<inputModuleHelix->GetLineSegList().size()<<"\n";
        for( auto it = inputModuleHelix->GetLineSegList().begin(); it != inputModuleHelix->GetLineSegList().end()-1; ++it )
        {
            stream<<std::fixed<<"\"{{"<<it->GetP0()[0]<<","<<it->GetP0()[1]<<","<<it->GetP0()[2]<<"},{"<<it->GetP1()[0]<<","<<it->GetP1()[1]<<","<<it->GetP1()[2]<<"}}\",";
            //++tc;
            //std::cout<<tc<<"\n";
        }

        std::cout<<"Here A\n";

        // Add the last point
        auto lastIt = inputModuleHelix->GetLineSegList().end()-1;
        stream<<std::fixed<<"\"{{"<<lastIt->GetP0()[0]<<","<<lastIt->GetP0()[1]<<","<<lastIt->GetP0()[2]<<"},{"<<lastIt->GetP1()[0]<<","<<lastIt->GetP1()[1]<<","<<lastIt->GetP1()[2]<<"}}\"\n";

        std::cout<<"Here A\n";

        ///** Output module sequence **/
        //for( auto it = inputModuleHelix->GetSeqList().begin(); it != inputModuleHelix->GetSeqList().end()-1; ++it )
        //{
        //    for(auto it2 = it->begin(); it2 != it->end(); ++it2)
        //    {
        //        stream<<*it2;
        //    }
        //    stream<<",";
        //}

        std::cout<<"Here B\n";

        //// Add the sequence for last helix
        //auto lastIt2 = inputModuleHelix->GetSeqList().end()-1;
        //for(auto it2 = lastIt2->begin(); it2 != lastIt2->end(); ++it2)
        //{
        //    stream<<*it2;
        //}
        //stream<<"\n";

        std::cout<<"Here C\n";

        /** Output data geometry **/
        for( auto it = inputDataHelix->GetLineSegList().begin(); it != inputDataHelix->GetLineSegList().end()-1; ++it )
        {
            stream<<std::fixed<<"\"{{"<<it->GetP0()[0]<<","<<it->GetP0()[1]<<","<<it->GetP0()[2]<<"},{"<<it->GetP1()[0]<<","<<it->GetP1()[1]<<","<<it->GetP1()[2]<<"}}\",";
        }

        std::cout<<"Here D\n";

        // Add the last point
        lastIt = inputDataHelix->GetLineSegList().end()-1;
        stream<<std::fixed<<"\"{{"<<lastIt->GetP0()[0]<<","<<lastIt->GetP0()[1]<<","<<lastIt->GetP0()[2]<<"},{"<<lastIt->GetP1()[0]<<","<<lastIt->GetP1()[1]<<","<<lastIt->GetP1()[2]<<"}}\"\n";

        

        ///** Output data sequence **/
        //for( auto it = inputDataHelix->GetSeqList().begin(); it != inputDataHelix->GetSeqList().end()-1; ++it )
        //{
        //    for(auto it2 = it->begin(); it2 != it->end(); ++it2)
        //    {
        //        stream<<*it2;
        //    }
        //    stream<<",";
        //}
        //// Add the sequence for last helix
        //lastIt2 = inputDataHelix->GetSeqList().end()-1;
        //for(auto it2 = lastIt2->begin(); it2 != lastIt2->end(); ++it2)
        //{
        //    stream<<*it2;
        //}
        //stream<<"\n";

        std::cout<<"Here E\n";

        /** Write into the file **/
        outputFile<<stream.str();
        // Close the file
        outputFile.close();
    }

}

// Make sure that the module set contains more helices
bool SwtichModuleData(HelixGeometry* & inputModuleHelix, HelixGeometry* & inputDataHelix, HelixGeometry* moduleHelix, HelixGeometry* dataHelix)
{
    inputModuleHelix = moduleHelix;
    inputDataHelix   = dataHelix;
    if( moduleHelix->GetHelixNum() < dataHelix->GetHelixNum() )
    {
        std::cout<<"The input module and data are switched.\n";
        inputModuleHelix = dataHelix;
        inputDataHelix   = moduleHelix;

        //switchFlag = true;
        return true;
    }

    return false;
}

// Read in the PDB file
void ReadInTheFile(const char* fileName, PDBObject & pdbReader)
{
    // Begin to time
//    cpuTimer.Start();

    // Read the file into a PDB object
    pdbReader.ReadInPDBFile(fileName);
    
    // Output elapsed time
//    elapsedTime = cpuTimer.Tick();
//    std::cout<<"Time to read in PDB file: "<<elapsedTime*(1e-3)<<"s\n";
}


// Generate helix geometry (line segment) multiple chains
void RetrieveHelixGeometry(HelixGeometry& helixGeo, PDBObject& pdbReader, Vector3f modelOffset, std::string chainIDs)
{
    // Empty vector to hold the chains
    std::vector<ProteinChain> chains;

    // Read in all existing chains
    if("All" == chainIDs)
    {
        // Retrieve all existing chains
        //chains = pdbReader.GetChains();
        for(auto it = pdbReader.GetChains().begin(); it != pdbReader.GetChains().end(); ++it)
        {
            chains.push_back( it->second );
        }
    }
    // Read in specific chains --- case sensitive
    else
    {
        // Retrieve all specified chains
        for(auto chainID = chainIDs.begin(); chainID != chainIDs.end(); ++chainID)
        {
            // Chain ID protection.
            if( *chainID < 'A' || *chainID > 'Z' )
            {
                std::cout<<"Include invalid chain ID: "<<*chainID<<"...\n";
                continue;
            }
            else
            {
                // ATTENTION!!!! So far there is no protection on chain index
                // If the "pdbReader" does not contain that chain, we will get
                // an empty chain
                chains.push_back( pdbReader.GetChain(*chainID) );
            }
        }
    }

    // Allocate the memory to hold all the helices --- empty chains are allowed and will contribute nothing
    int helixNum = 0;
    for(auto chain = chains.begin(); chain != chains.end(); ++chain)
    {
        helixNum += chain->GetSSE().GetHelices().size();
    }
    helixGeo.Reset(helixNum);

    int i = 0;
    for(auto chain = chains.begin(); chain != chains.end(); ++chain)
    {
        auto & helixList    = chain->GetSSE().GetHelices();
        auto & residuePool  = chain->GetAllResidues();

        for(auto it = helixList.begin(); it != helixList.end(); ++it, ++i)
        {
            // Find all the atoms' residue indices
            auto & indices = it->GetResidueIndices();
            // Retrieve all the atoms
            pointList atomCoordList(indices[1]-indices[0]+1, Vector3f(0.0f,0.0f,0.0f));
            for(int h = indices[0], k = 0; h < indices[1]+1; ++h, ++k)
            {
                auto residue = residuePool.find(h);
                // If the residue is missing then continue
                if(residue != residuePool.end())
                {
                    // Record the atom. residue (map): first key, second value
                    atomCoordList[k] = residue->second.GetCAlpha().GetCoord();

                    // Record the atom's corresponding residue sequence in the helix list
                    helixGeo.AddSeq(i, residue->second.GetAminoAbbrev());
                }
            }
            // Fit the line from all the atoms and add the line into moduleHelix list
            lineSeg lineS = FitLine(atomCoordList);
            // Add the line segment into helix list
            helixGeo.AddLineSeg(i, LineSegment(lineS[0], lineS[1]));        

        }// End for of current chain

    } // End for of all chains

    // Update model center --- This is mainly for visualization in Mathematica
    helixGeo.ComputeModelCenter();
    helixGeo.SetModelOffset(modelOffset);
    helixGeo.TransformModel();
}

// Generate helix geometry (line segment) --- only one chain
void RetrieveHelixGeometry(HelixGeometry& helixGeo, PDBObject& pdbReader, char chainID, Vector3f modelOffset)
{
    // Begin to time
//    cpuTimer.Start();

    // Grab the helix list from the chain
    auto & chain        = pdbReader.GetChain(chainID);
    auto & helixList    = chain.GetSSE().GetHelices();
    auto & residuePool  = chain.GetAllResidues();
    // Initiate the helixGeometry --- Release old stuff and allocate new memory
    //std::cout<<"helixList.size(): "<<helixList.size()<<"\n";
    helixGeo.Reset(helixList.size());
    // Record all the helices, one by one
    int i = 0;
    for(auto it = helixList.begin(); it != helixList.end(); ++it, ++i)
    {
        // Find all the atoms' residue indices
        auto & indices = it->GetResidueIndices();
        // Retrieve all the atoms
        pointList atomCoordList(indices[1]-indices[0]+1, Vector3f());
        for(int h = indices[0], k = 0; h < indices[1]+1; ++h, ++k)
        {
            auto residue = residuePool.find(h);
            if(residue != residuePool.end())
            {
                // Record the atom. residue (map): first key, second value
                atomCoordList[k] = residue->second.GetCAlpha().GetCoord();

                // Record the atom's corresponding residue sequence in the helix list
                helixGeo.AddSeq(i, residue->second.GetAminoAbbrev());
            }

        }
        // Fit the line from all the atoms and add the line into moduleHelix list
        lineSeg lineS = FitLine(atomCoordList);
        // Add the line segment into helix list
        helixGeo.AddLineSeg(i, LineSegment(lineS[0], lineS[1]));        

    }

    // Update model center
    // This is mainly for visualization in Mathematica
    helixGeo.ComputeModelCenter();
    helixGeo.SetModelOffset(modelOffset);
    helixGeo.TransformModel();

    // Output elapsed time
//    elapsedTime = cpuTimer.Tick();
//    std::cout<<"Time to generate helix geometry: "<<elapsedTime*(1e-3)<<"s\n";
}


// Construct the score matrix for double line direction
void ConstructScoreMatrixDLD(ScoreMatrixDoubleLineDir & scoreMatrixDLD, HelixGeometry *moduleHelix, HelixGeometry *dataHelix, const std::vector<double> & thresholds, int & scoreFunFlag)
{
    // Begin to time
//    cpuTimer.Start();

    // Set up the module and data for the score matrix
    scoreMatrixDLD.SetModuleAndDataSet(moduleHelix, dataHelix);
    // Set up score funtion --- 0:Gaussian; 1:Epan
    scoreMatrixDLD.SetScoreFun(scoreFunFlag);
    // Set up parameters
    scoreMatrixDLD.SetParameters(thresholds[0],thresholds[1],thresholds[2],thresholds[3],thresholds[4],thresholds[5],thresholds[6]);
    // Generate the score matrix
    scoreMatrixDLD.GenScoreMatrix();

    // Output elapsed time
//    elapsedTime = cpuTimer.Tick();
//    std::cout<<"Time to construct score matrix: "<<elapsedTime*(1e-3)<<"s\n";

    //std::cout<< std::setprecision(15)<<scoreMatrixDLD.GetMatrix().col(0)[54]<<"\n";
    //for(int i = 0; i < scoreMatrixDLD.GetMatrix().row(0).size(); ++i)
    //    std::cout<<scoreMatrixDLD.GetMatrix().col(0)[i]<<", ";
    //std::cout<<"\n";
    //for(int i = 0; i < 10; ++i)
    //{
    //    std::cout<<scoreMatrixDLD.GetMatrix()(i,i)<<", ";
    //}
    //std::cout<<"\n";
}

// Solve the matrix for principle eigenvector and binarize it to obtain indicator vector
void SolveForResultDLD(SMSolverDLD & smSolverDLD, const ScoreMatrixDoubleLineDir & matrix, int methodFlag, double eigenRatio)
{
    // Begin to time
//    cpuTimer.Start();

    //smSolverDLD.SolveForEigenvectors(matrix.GetMatrix());
    switch (methodFlag)
    {
        case 0:
            smSolverDLD.IterativePicking(matrix);
    	    break;
        case 1:
            smSolverDLD.SimpleGreedyPick(matrix);
            break;
        case 2:
            smSolverDLD.IterativePickingNewStop(matrix, eigenRatio);
            break;
        default:
            std::cout<<"Invalid binarize method....\n";
    }
    

    // Output elapsed time
//    elapsedTime = cpuTimer.Tick();
//    std::cout<<"Time to solve for eigenvectors: "<<elapsedTime*(1e-3)<<"s\n";

    //for(int i = 10; i < 20; ++i)
    //    std::cout<<"eigenVec: "<<smSolverDLD.GetPrinEigenVec()[i]<<"\n";
}


// Construct the score matrix
void ConstructScoreMatrixSLD(ScoreMatrixSingleLineDir & scoreMatrixSLD, HelixGeometry *moduleHelix, HelixGeometry *dataHelix, const std::vector<double> & thresholds, int & scoreFunFlag)
{
    // Begin to time
    //    cpuTimer.Start();

    // Set up the module and data for the score matrix
    scoreMatrixSLD.SetModuleAndDataSet(moduleHelix, dataHelix);
    // Set up score funtion --- 0:Gaussian; 1:Epan
    scoreMatrixSLD.SetScoreFun(scoreFunFlag);
    // Set up parameters
    scoreMatrixSLD.SetParameters(thresholds[0],thresholds[1],thresholds[2],thresholds[3],thresholds[4],thresholds[5],thresholds[6]);
    // Generate the score matrix
    scoreMatrixSLD.GenScoreMatrix();

    // Output elapsed time
    //    elapsedTime = cpuTimer.Tick();
    //    std::cout<<"Time to construct score matrix: "<<elapsedTime*(1e-3)<<"s\n";
}

// Construct the score matrix for orient single line direction
void ConstructScoreMatrixSLDOrient(ScoreMatrixSingleLineDir & scoreMatrixSLD, HelixGeometry *moduleHelix, HelixGeometry *dataHelix, const std::vector<double> & thresholds, int & scoreFunFlag)
{
    // Begin to time
    //    cpuTimer.Start();

    // Set up the module and data for the score matrix
    scoreMatrixSLD.SetModuleAndDataSet(moduleHelix, dataHelix);
    // Set up score funtion --- 0:Gaussian; 1:Epan
    scoreMatrixSLD.SetScoreFun(scoreFunFlag);
    // Set up parameters
    scoreMatrixSLD.SetParameters(thresholds[0],thresholds[1],thresholds[2],thresholds[3],thresholds[4],thresholds[5],thresholds[6]);
    // Generate the score matrix
    scoreMatrixSLD.GenScoreMatrixOrientLine();

    // Output elapsed time
    //    elapsedTime = cpuTimer.Tick();
    //    std::cout<<"Time to construct score matrix: "<<elapsedTime*(1e-3)<<"s\n";
}


// Solve the matrix for principle eigenvector and binarize it to obtain indicator vector
void SolveForResultSLD(SMSolverSLD & smSolverSLD, const ScoreMatrixSingleLineDir & matrix, int methodFlag, double eigenRatio)
{
    // Begin to time
//    cpuTimer.Start();

    //smSolverDLD.SolveForEigenvectors(matrix.GetMatrix());
    switch (methodFlag)
    {
        case 0:
            smSolverSLD.IterativePicking(matrix);
    	    break;
        case 1:
            smSolverSLD.SimpleGreedyPick(matrix);
            break;
        case 2:
            smSolverSLD.HungarianPick(matrix);
            break;
        case 3:
            smSolverSLD.IterativePickingNewStop(matrix, eigenRatio);
            break;
        default:
            std::cout<<"Invalid binarize method....\n";
    }
    

    // Output elapsed time
//    elapsedTime = cpuTimer.Tick();
//    std::cout<<"Time to solve for eigenvectors: "<<elapsedTime*(1e-3)<<"s\n";

    //for(int i = 10; i < 20; ++i)
    //    std::cout<<"eigenVec: "<<smSolverDLD.GetPrinEigenVec()[i]<<"\n";
}


// Construct the score matrix for 2D point cloud
void ConstructScoreMatrix2DP(ScoreMatrix2DPointCloud & scoreMatrix2DP, PointCloud2D* modulePoints, PointCloud2D* dataPoints, const std::vector<double> & thresholds, int & scoreFunFlag)
{
    // Set up the module and data for the score matrix
    scoreMatrix2DP.SetModuleAndDataSet(modulePoints, dataPoints);
    // Set up score funtion --- 0:Gaussian; 1:Epan
    scoreMatrix2DP.SetScoreFun(scoreFunFlag);
    // Set up parameters
    scoreMatrix2DP.SetParameters(thresholds[0],thresholds[1],thresholds[2],thresholds[3],thresholds[4],thresholds[5],thresholds[6]);
    // Generate the score matrix
    scoreMatrix2DP.GenScoreMatrix();
}

// Solve the matrix for principle eigenvector and binarize it to obtain indicator vector
void SolveForResult2DP(SMSolver2DP & smSolver2DP, const ScoreMatrix2DPointCloud & matrix, int methodFlag)
{
    switch (methodFlag)
    {
    case 0:
        smSolver2DP.IterativePicking(matrix);
        break;
    case 1:
        smSolver2DP.SimpleGreedyPick(matrix);
        break;
    case 2:
        smSolver2DP.HungarianPick(matrix);
        break;
    default:
        std::cout<<"Invalid binarize method....\n";
    }
}

// Old test main function
void oldTestMain()
{
    // Set the file name
    char moduleFileName[100] = "../Data/1oel_A.pdb";
    char dataFileName[100]   = "../Data/2c7c_A.pdb";

    /************************* Read in the File and build geometry *************************/
    // Read the file into PDBObject. Without specific chainID,
    // the PDBReader will all the available chains

    ReadInTheFile(moduleFileName, pdbReader);
    // Retrieve Helix model 
    RetrieveHelixGeometry(moduleHelix, pdbReader, 'A', moduleOffset);

    ReadInTheFile(dataFileName, pdbReader);
    // Retrieve Helix model 
    RetrieveHelixGeometry(dataHelix, pdbReader, 'A', dataOffset);

    //auto lines  = moduleHelix.GetAllLineSegs();
    //int lineNUm = moduleHelix.GetlineNum();
    //Vector3f moduleCenter = moduleHelix.GetModelCenter();

    //for(int i = 0 ; i < lineNUm; ++i)
    //{
    //    Vector3f point1 = lines[i][0] - moduleCenter + moduleOffset;
    //    Vector3f point2 = lines[i][1] - moduleCenter + moduleOffset;
    //    std::cout<<point1(0)<<","<<point1(1)<<","<<","<<point1(2)<<"\n";
    //    std::cout<<point2(0)<<","<<point2(1)<<","<<","<<point2(2)<<"\n";
    //    std::cout<<"\n";
    //}

    /************************* Construct Score Matrix *************************/
    // Switch so that moduleHelix is always larger than dataHelix


    ConstructScoreMatrixDLD(scoreMatrixDLD, inputModuleHelix, inputDataHelix, thresholds, scoreFunFlag);

    //scoreMatrixDLD.GenScoreMatrix();
    //std::cout<<scoreMatrixDLD.GetMatrix();

    /************************* Solve For Eigenvector *************************/

    SolveForResultDLD(smSolverDLD, scoreMatrixDLD);

    //SelfAdjointEigenSolver<MatrixXd> solver;
    //MatrixXd m0(1040,1040);
    //for(int i=0; i<1040; ++i)
    //    for(int j=0; j<1040; ++j)
    //        m0(i,j) = 0;

    //m0(0,0)       = 1;
    //m0(0,1039)    = 4;
    //m0(1039,0)    = 4;
    //m0(1039,1039) = 1;

    //// Begin to time
    //cpuTimer.Start();

    //solver.compute(m0);

    //// Output elapsed time
    //elapsedTime = cpuTimer.Tick();
    //std::cout<<"Time to solve for eigenvalues: "<<elapsedTime*(1e-3)<<"s\n";

    //int eigenNum = solver.eigenvalues().size();
    //std::cout<<"eigen num: "<<eigenNum<<"\n";

    //std::cout<<"\neigenvalues:\n"<<solver.eigenvalues()[1039]<<"\n";
    //std::cout<<"\neigenvalues:\n"<<solver.eigenvalues()[1038]<<"\n";
    //std::cout<<"\neigenvalues:\n"<<solver.eigenvalues()[0]<<"\n";

    //std::vector<double> v(5);
    //v[0] = 1.0;
    //v[1] = 2.0;
    //v[2] = 3.0;
    //v[3] = 2.0;
    //v[4] = 1.0;


    //int largestEleIndex = std::distance(v.begin(), std::max_element(v.begin(), v.end()));
    //std::cout<<largestEleIndex<<"\n";
    //MatrixXd m1(4,4);

    //m1 << 1,2,3,4,
    //    5,6,7,8,
    //    9,10,11,12,
    //    13,14,15,16;


    //std::cout<<m1<<"\n";

    //m1.conservativeResize(2,2);

    //std::cout<<"\n matrix:\n "<<m1<<"\n";
    //std::cout<<"\n matrix size: "<<m1.cols()<<"\n";

    //std::cout<<"matrix size: "<<m1.cols()<<"\n";

    //VectorXd v1 = m1.col(0);
    //std::cout<<((ArrayXd)m1.col(0) > 2).all()<<"\n";    
    //std::cout<<"\n"<< m0 <<"\n";

    //std::cout<<"\neigenvectors:\n";

    //for(int i = 0; i < eigenNum; ++i)
    //    std::cout<<"\n"<<solver.eigenvectors().col(i)<<"\n";

    //std::cout<<"\neigenvalues:\n"<<solver.eigenvalues()<<"\n";

    //SYSTEM_INFO sysinfo;
    //GetSystemInfo( &sysinfo );
    //int numCPU = sysinfo.dwNumberOfProcessors;
    //std::cout<<"Number of processers: "<<numCPU<<std::endl;
    ////OMP_NUM_THREADS=n ./my_program;
    //omp_set_num_threads(12);
    //Eigen::setNbThreads(12);
    ////Eigen::setNbThreads(2);
    //int num_threads = Eigen::nbThreads();
    //std::cout<<"Eigen num_threads: "<<num_threads<<std::endl;

    //SelfAdjointEigenSolver<MatrixXd> solver;
    //MatrixXd m0 = MatrixXd::Random(1040, 1040);

    // Begin to time
    //cpuTimer.Start();

    //VectorXd tempV = m0.selfadjointView<Eigen::Lower>().eigenvalues();

    //solver.compute( m0.selfadjointView<Eigen::Lower>() );
    //solver.compute( m0 );

    //for(int i = 0; i < 19; ++i)
    //{
    //    m0 = MatrixXd::Random(2*(20-i)*(26-i), 2*(20-i)*(26-i));
    //    solver.compute(m0);
    //} 

    // Output elapsed time
    //elapsedTime = cpuTimer.Tick();
    //std::cout<<"Time to solve for eigenvalues: "<<elapsedTime*(1e-3)<<"s\n";



    /************************* Export result to disk for Mathematica *************************/


    // End timing
    //cpuTimer.GetTime();
    std::cout<<"Please press any key..."<<"\n";
    getchar();
}