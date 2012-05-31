/**
 * Code used to test compressing difference positions using Bzet
 * Bzet available at https://github.com/achow/Bzet
 *
 * Define NODE_ELS to be the flavor of Bzet to use (4-64, powers of 2)
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>

#define BZET_IMPL_
#include "Bzet.h"

using namespace std;

enum DiffType { SNP = 0, DEL = 1, INS = 2 };
enum Base { B_A = 0, B_C = 1, B_T = 2, B_G = 3 };

struct GenDiff {
    DiffType difftype;
    string tag;
    size_t pos;
    string ref_alt;
};

void encode(std::vector<GenDiff>::iterator itstart, std::vector<GenDiff>::iterator itend) {
    BZET b;
    std::vector<GenDiff>::iterator it = itstart;
    for (; itstart != itend; itstart++) {
        b.seqset(itstart->pos);
    }
    cout << it->difftype << "," << it->tag << ", size: " << b.size() << endl;
}

int main() {
    /*
    cout << "Loading gwah shared library" << endl;
    void *hdl = dlopen(GWAHFILE, RTLD_LAZY | RTLD_GLOBAL);
    if (!hdl) {
        cerr << "Load failed" << endl;
        return 1;
    }
    dlerror();

    cout << "Loading symbols" << endl;
 */   
    cout << NODE_ELS << endl;

    std::vector<GenDiff> diffs;
    while (cin.good()) {
        string line;
        getline(cin, line);
        if (line.empty())
            break;
        char *cline = (char *) line.c_str();
        GenDiff cdiff;
        cdiff.difftype = (DiffType) atoi(strtok(cline, ","));
        cdiff.tag = strtok(NULL, ",");
        cdiff.pos = atoi(strtok(NULL, ","));
        cdiff.ref_alt = strtok(NULL, ",");
        diffs.push_back(cdiff);
    }

    cout << "Done reading input" << endl;

    clock_t TMS_CLK_TCK = sysconf(_SC_CLK_TCK);
    struct tms tmsstart, tmsend;

    string prevtag = diffs.begin()->tag;
    DiffType prevdiff = diffs.begin()->difftype;
    std::vector<GenDiff>::iterator chrstart = diffs.begin();
    times(&tmsstart);
    for (std::vector<GenDiff>::iterator it = diffs.begin(); it != diffs.end(); it++) {
        if (it->difftype != prevdiff || it->tag != prevtag || (it + 1) == diffs.end()) {
            encode(chrstart, it);
            chrstart = it;
            prevtag = it->tag;
            prevdiff = it->difftype;
        }
    }
    times(&tmsend);
    cout << "TIME: " << ((tmsend.tms_utime - tmsstart.tms_utime) / (double) TMS_CLK_TCK) << endl;
}
