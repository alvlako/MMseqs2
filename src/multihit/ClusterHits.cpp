#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "itoa.h"
#include "Util.h"


#ifdef OPENMP
#include <omp.h>
#endif

struct hit{
    unsigned int qPos;
    unsigned int tPos;
    bool qStrand;
    bool tStrand;
};

//std::sort by default sorts the hits by qPos in ascending order
bool operator<(const hit& a, const hit& b)
{
    return a.qPos < b.qPos;
}

//Log of binomial coefficient as defined in combinepvalperset
double LBinCoeff(double* lookup, int M, int k);


double hypergeoDensity(double* lookup, int k, int n, int N, int K){
    return exp(LBinCoeff(lookup, n, k) + LBinCoeff(lookup, N - n, K - k) - LBinCoeff(lookup, N, K));
}

double hypergeoDistribution(double* lookup, int k, int n, int N, int K){
    double sum = 0;
    for (int i = 0; i < k + 1; i++){
        sum += hypergeoDensity(lookup, i, n, N, K);
    }
    return sum;
}

double clusterPval(double* lookup, int k, int m, int K, int Nq, int Nt) {
    int minKm = (K < m) ? K : m;
    double sum = 0;
    for (int Kp = k; Kp < minKm; Kp++){
        double pHG = hypergeoDensity(lookup, (Kp - 1), (K - 1), (Nq - 1), (m - 1));
        double PHG = hypergeoDistribution(lookup, (k - 2), (Kp - 1), (Nt - 1), (m - 1));
        sum += pHG * (1 - pow(PHG, k));
    }
    return 1 - pow((1 - sum), (K - k + 1));
}


double clusterPval_fast(int m, int K, int Nq, int Nt) {
    double power = 2.0 * pow((m - 1), 2) * (K - 1);
    return 1 - pow((1 - 1.0 * K / (Nq * Nt)), power);
}


double orderingPval(int k, int m){
    return (1 - 1.0 * m / k)/(pow(2,m) * exp(lgamma(m+1)));
}


int findSpan(const std::vector<hit> &cluster){
    unsigned int iMax = 0;
    unsigned int iMin = INT_MAX;
    unsigned int jMax = 0;
    unsigned int jMin = INT_MAX;
    for (size_t l = 0; l < cluster.size(); l++){
        iMax = (cluster[l].qPos > iMax) ? cluster[l].qPos : iMax;
        iMin = (cluster[l].qPos < iMin) ? cluster[l].qPos : iMin;
        jMax = (cluster[l].tPos > jMax) ? cluster[l].tPos : jMax;
        jMin = (cluster[l].tPos < jMin) ? cluster[l].tPos : jMin;
    }
    //index 1 and 2 has a span of 2
    int spanI =iMax - iMin + 1;
    int spanJ =jMax - jMin + 1;
    return (spanI > spanJ) ? spanI : spanJ;
}

int findConservedPairs(std::vector<hit> &cluster){
    //re-index hits in clusters with pos in query set in ascending order
    std::sort(cluster.begin(), cluster.end());
    int m = 0;
    for (size_t l = 0; l < cluster.size()-1; l++){
        bool isSameOrder = (cluster[l+1].tPos > cluster[l].tPos);
        bool isSameStrand = (cluster[l].qStrand == cluster[l].tStrand);
        bool isSameStrand2 = (cluster[l+1].qStrand == cluster[l+1].tStrand);
        if((isSameStrand == isSameOrder) && (isSameStrand2 == isSameOrder)){
            m++;
        }
    }
    return m;
}

double orderingPval_fast(const std::vector<hit> &cluster){
    bool isSameOrder = (cluster[1].qPos > cluster[0].qPos) ? (cluster[1].tPos > cluster[0].tPos) : (cluster[1].tPos < cluster[0].tPos) ;
    bool isConservedPair = ((cluster[0].qStrand == cluster[0].tStrand) == isSameOrder ) && ((cluster[1].qStrand == cluster[1].tStrand) == isSameOrder );
    return isConservedPair ? 0.25 : 1.0;
}


double clusterMatchScore(double* lookup, std::vector<hit> &cluster, int K, int Nq, int Nt){
    int span = findSpan(cluster);
    int k = cluster.size();
    double pClu;
    double pOrd;
    if(k == 2){
        pClu = clusterPval_fast(span, K, Nq, Nt);
        pOrd = orderingPval_fast(cluster);
    }
    else{
        pClu = clusterPval(lookup, k, span, K, Nq, Nt);
        int m = findConservedPairs(cluster);
        pOrd = orderingPval(k, m);
    }

    //full score would be: -log(pClu) - log(pOrd) + log(1 -log(pClu) - log(pOrd))
    return -log(pClu) - log(pOrd); 
}

std::vector<hit> groupNodes(const std::vector<std::vector<int>> &nodeList, const std::vector<hit> &matchList, int i, int j){
    std::vector<hit> cluster;
    for(size_t m = 0; m < nodeList[i].size(); m++){
        cluster.push_back(matchList[nodeList[i][m]]);
    }
    for(size_t n = 0; n < nodeList[j].size(); n++){
        cluster.push_back(matchList[nodeList[j][n]]);
    }
    return cluster;
}


int clusterhits(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);


    std::string sizeDBName = std::string(par.db1) + "_set_size";
    std::string sizeDBIndex = std::string(par.db1) + "_set_size.index";
    DBReader<unsigned int> querySizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	querySizeReader.open(DBReader<unsigned int>::NOSORT);

    sizeDBName = std::string(par.db2) + "_set_size";
    sizeDBIndex = std::string(par.db2) + "_set_size.index";
    DBReader<unsigned int> targetSizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	targetSizeReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> qlookupReader(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_LOOKUP);
    qlookupReader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::LookupEntry* qlookup = qlookupReader.getLookup();

    DBReader<unsigned int> tlookupReader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_LOOKUP);
    tlookupReader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::LookupEntry* tlookup = tlookupReader.getLookup();

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> headerReader(par.hdr3.c_str(), par.hdr3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, resultReader.getDbtype());
    dbw.open();

    //find the max number of genes/ORFs of Nq&Nt
    //TODO: can we do that faster by scanning the header file?
    unsigned int maxOrfCount = 0;
    for (size_t i = 0; i < querySizeReader.getSize(); ++i) { 
        unsigned int currentCount = Util::fast_atoi<unsigned int>(querySizeReader.getData(i, 0));
        if (currentCount > maxOrfCount) {
            maxOrfCount = currentCount;
        };
    }
    for (size_t i = 0; i < targetSizeReader.getSize(); ++i) { 
        unsigned int currentCount = Util::fast_atoi<unsigned int>(targetSizeReader.getData(i, 0));
        if (currentCount > maxOrfCount) {
            maxOrfCount = currentCount;
        };
    }

    //create a lookup table for all possible log gamma values (n-1)!, n! will be lookup[n + 1]
    double* lGammaLookup = new double[maxOrfCount + 2];
    for (size_t i = 0; i < maxOrfCount + 2; ++i) { 
        lGammaLookup[i] = lgamma(i);
    }


Debug::Progress progress(resultReader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(10 * 1024);
        std::string header;
        header.reserve(1024);
        const char *entry[255];
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            progress.updateProgress();

            std::vector<hit> match;

            //from result header file read Nq & Nt
            header = headerReader.getData(i, thread_idx);
            std::vector<std::string>  hdrcolumns = Util::split(header, "\t");
            unsigned int Nq = Util::fast_atoi<size_t>(hdrcolumns[4].c_str()); //query set size
            unsigned int Nt = Util::fast_atoi<size_t>(hdrcolumns[5].c_str()); //target set size

            //read all hits from a match and extract pos and strand, append them to std::vector<hit> match
            char *data = resultReader.getData(i, thread_idx);
            while (*data != '\0'){
                //TODO: read first two column for qid and tid, get entryName from lookups
                size_t columns = Util::getWordsOfLine(data, entry, 255);
                data = Util::skipLine(data);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid alignment result record\n";
                    EXIT(EXIT_FAILURE);
                }
                unsigned int qid = Util::fast_atoi<unsigned int>(entry[0]);
                unsigned int tid = Util::fast_atoi<unsigned int>(entry[1]);

                //TODO: from seqid get entryName in lookup, split by "_" and retrieve info
                std::vector<std::string> qcolumns = Util::split(qlookup[qid].entryName, "_");
                std::vector<std::string> tcolumns = Util::split(tlookup[tid].entryName, "_");
                int qStart = Util::fast_atoi<size_t>(qcolumns[qcolumns.size()-2].c_str());
                int qEnd = Util::fast_atoi<size_t>(qcolumns.back().c_str());
                int tStart = Util::fast_atoi<size_t>(tcolumns[tcolumns.size()-2].c_str());
                int tEnd = Util::fast_atoi<size_t>(tcolumns.back().c_str());
                hit tmpBuffer;

                //pos is the protein index in the genome, strand is determined by start and end coordinates
                tmpBuffer.qPos = Util::fast_atoi<size_t>(qcolumns[qcolumns.size()-3].c_str());
                tmpBuffer.tPos = Util::fast_atoi<size_t>(tcolumns[tcolumns.size()-3].c_str());
                tmpBuffer.qStrand = (qStart < qEnd) ? 1 : 0;
                tmpBuffer.tStrand = (tStart < tEnd) ? 1 : 0;

                match.push_back(tmpBuffer);
            }

            //TODO: total number of hits can be read from the header
            //int K = Util::fast_atoi<int>(hdrcolumns[3].c_str());//total number of hits
            int K = match.size();

            //initiallize distance matrix to be [K][K]
            double** DistMat = new double*[K]; // Rows
            for (int i = 0; i < K; i++)
            {
                DistMat[i] = new double[K]; // Columns
            }

            std::vector<int> dmin(K); //index of closest cluster/highest score
            std::vector<std::vector<int>> nodes(K);
            //assign each node with the index of the singleton cluster
            for(int n = 0; n < K; n++){
                nodes[n].push_back(n);
            }

            for(int i = 0; i < K; i++){
                for(int j = 0; j < K; j++){
                    if(i == j){
                        DistMat[i][j] = 0.0; //set score = 0 to self similarities
                    }
                    else{
                        std::vector<hit> tmpCluster = groupNodes(nodes,match,i,j);
                        DistMat[i][j] = clusterMatchScore(lGammaLookup,tmpCluster,K,Nq,Nt); //score(i,j)
                    }
                    dmin[i] = (DistMat[i][j] > DistMat[i][dmin[i]]) ? j : dmin[i];
                }
            }

            double maxScore = DBL_MAX;
            while(maxScore > 1){ //par.sMin
                int i1 = 0;
                int i2;
                //find closest pair of clusters (i1,i2), i.e pair of clusters with highest score
                for(int s = 0; s < K-1; s++){
                    i1 = 0;
                    for (int i = 0; i < K; i++){
                        i1 = (DistMat[i][dmin[i]] > DistMat[i1][dmin[i1]]) ? i : i1;
                    }
                    i2 = dmin[i1];
                }

                maxScore = DistMat[i1][i2];

                //delete node i2 and add i2 to i1
                nodes[i2].clear();
                nodes[i1].push_back(i2);
                    
                //overwrite row and column i1 with dist[i1,i2][j]
                for(int j = 0; j < K; j++){
                    if(i1 == j){
                        DistMat[i1][j] = 0.0; 
                    }
                    else{
                        std::vector<hit> tmpCluster = groupNodes(nodes,match,i1,j);
                        DistMat[i1][j] = clusterMatchScore(lGammaLookup,tmpCluster,K,Nq,Nt);
                        DistMat[j][i1] = DistMat[i1][j];
                    }

                    //set scores in old row and column i2 to 0
                    DistMat[i2][j] = 0.0;
                    DistMat[j][i2] = 0.0;

                    //for row i1 find index dmin[i1] of closest cluster
                    dmin[i1] = (DistMat[i1][j] > DistMat[i1][dmin[i1]]) ? j : dmin[i1];
                    
                    //for every other row j != i1 && j != i2 compare the new score DistMat[j][i1] with DistMat[j][dmin[j]]
                    if(j != i1 && j != i2) {
                        dmin[j] = (DistMat[j][i1] > DistMat[j][dmin[j]]) ? i1 : dmin[j];
                    }

                }
                
            }


            
            delete[] DistMat;
            match.clear();

            dbw.writeData(buffer.c_str(), buffer.length(), resultReader.getDbKey(i), thread_idx);
            buffer.clear();
        }
    }
    dbw.close(true);
    resultReader.close();
    headerReader.close();


    qlookupReader.close();
    tlookupReader.close();
    querySizeReader.close();
    targetSizeReader.close();
    delete[] lGammaLookup;
    return EXIT_SUCCESS;
}

