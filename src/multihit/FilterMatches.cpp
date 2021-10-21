#include "Parameters.h"
//#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include <algorithm>
#include <climits>

#ifdef OPENMP
#include <omp.h>
#endif

int filtermatches(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string sizeDBName = std::string(par.db1) + "_set_size";
    std::string sizeDBIndex = std::string(par.db1) + "_set_size.index";
    DBReader<unsigned int> sizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    sizeReader.open(DBReader<unsigned int>::NOSORT);

    sizeDBName = std::string(par.db2) + "_set_size";
    sizeDBIndex = std::string(par.db2) + "_set_size.index";
    DBReader<unsigned int> tsizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    tsizeReader.open(DBReader<unsigned int>::NOSORT);

    std::string setDBName = std::string(par.db1) + "_member_to_set";
    std::string setDBIndex = std::string(par.db1) + "_member_to_set.index";
    DBReader<unsigned int> qsetReader(setDBName.c_str(), setDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qsetReader.open(DBReader<unsigned int>::NOSORT);

    setDBName = std::string(par.db2) + "_member_to_set";
    setDBIndex = std::string(par.db2) + "_member_to_set.index";
    DBReader<unsigned int> tsetReader(setDBName.c_str(), setDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    tsetReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> matchReader(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    matchReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbwriter(par.db5.c_str(), par.db5Index.c_str(), par.threads, par.compressed, resultReader.getDbtype());
    dbwriter.open();

    DBWriter headerWriter(par.hdr5.c_str(), par.hdr5Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

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

        unsigned int match_idx = 0;

        std::string cEval;
        cEval.reserve(255);
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < matchReader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int matchKey = matchReader.getDbKey(i);
            size_t resId = resultReader.getId(matchKey);
            if (resId == UINT_MAX) {
                Debug(Debug::WARNING) << "Missing alignment result for key " << matchKey << " \n";
                continue;
            }
            char *resultData = resultReader.getData(resId, thread_idx);
            char *matchData = matchReader.getData(i, thread_idx);
            // unsigned int resultKey = resultReader.getDbKey(i);
            // size_t qsetKey = Util::fast_atoi<size_t>(qsetReader.getDataByDBKey(resultKey, thread_idx));
            size_t qsetKey = resultReader.getDbKey(i);
            unsigned int setSize = Util::fast_atoi<unsigned int>(sizeReader.getDataByDBKey(qsetKey, thread_idx));
            double logPvalThr = log(par.alpha/(setSize + 1));
            while (*matchData != '\0'){
                size_t columns = Util::getWordsOfLine(matchData, entry, 255);
                matchData = Util::skipLine(matchData);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid alignment result record\n";
                    EXIT(EXIT_FAILURE);
                }
                unsigned int targetSetKey = Util::fast_atoi<unsigned int>(entry[0]);
                unsigned int tsetSize = Util::fast_atoi<unsigned int>(tsizeReader.getDataByDBKey(targetSetKey, thread_idx));
                //Decide if to filter out all hits if qsetkey == tsetkey
                //TODO: delete the following if statement when self matches are filtered out in search
                //TODO: disable if no self matches
                if(qsetKey != targetSetKey) {
                    //Add hit count, qset_size and tset_size to new header 
                    size_t hitCount = 0;
                    header.append(SSTR(qsetKey));
                    header.append("\t");
                    header.append(SSTR(targetSetKey));
                    header.append("\t");
                    header.append(entry[1], entry[2] - entry[1]);
                    char *data = resultData;//resultReader.getDataByDBKey(matchKey, thread_idx);
                    while (*data != '\0') {
                        //entry is rewritten here
                        columns = Util::getWordsOfLine(data, entry, 255);
                        unsigned int tid = Util::fast_atoi<unsigned int>(entry[1]);
                        unsigned int tsetKey = Util::fast_atoi<unsigned int>(tsetReader.getDataByDBKey(tid, thread_idx));
                        if(tsetKey != targetSetKey) {
                            data = Util::skipLine(data);
                            continue;
                        }
                        double logPval = strtod(entry[2], NULL);
                        if (logPval >= logPvalThr && logPval >= log(par.bhPvalThr)) {
                            data = Util::skipLine(data);
                            continue;
                        }
                        //(qsetid, tsetid,) qid, tid, p/e-value + rest of alignment info
                        data = Util::skipLine(data);
                        buffer.append(entry[0], entry[1] - entry[0] - 1);
                        buffer.append("\t");
                        buffer.append(entry[1], entry[2] - entry[1] - 1);
                        buffer.append("\t");
                        buffer.append(SSTR(exp(logPval)));
                        buffer.append("\t");
                        buffer.append(entry[3], data - entry[3]);
                        hitCount++;
                    }
                    if (buffer.back() != '\n'){
                        buffer.append("\n");
                    }
                    header.append("\t");
                    header.append(SSTR(hitCount));
                    header.append("\t");
                    header.append(SSTR(setSize)); //q set size
                    header.append("\t");
                    header.append(SSTR(tsetSize)); //t set size
                    header.append("\n");
                    dbwriter.writeData(buffer.c_str(), buffer.length(), match_idx, thread_idx);
                    headerWriter.writeData(header.c_str(), header.length(), match_idx, thread_idx);
                    match_idx++;
                    buffer.clear();
                    header.clear();
                }

            }
        }
    }
    dbwriter.close();
    headerWriter.close(true);
    matchReader.close();
    resultReader.close();
    qsetReader.close();
    tsetReader.close();
    sizeReader.close();
    return EXIT_SUCCESS;
}

