/*
 * This implementation is pretty crude. It does the job, but don't look at the 
 * code since it's ugly. There aren't comments in here anyways. Just compile
 * and run.
 *
 * Input from stdin
 * Output to stdout (compressed form)
 * Diagnostic messges to stderr
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>

#define STEP_BYTES 2

#define BZET_IMPL_
#include "Bzet.h"

#define INIT_ALLOC 1024
#define OUTPUT
#define COUTPUT

using namespace std;

enum DiffType { SNP = 0, DEL = 1, INS = 2 };
enum Base { B_A = 0, B_C = 1, B_T = 2, B_G = 3 };

struct GenDiff {
    DiffType difftype;
    string tag;
    size_t pos;
    string ref_alt;
};

size_t encsnps(vector<GenDiff>& vgf,
             vector<GenDiff>::iterator itstart,
             vector<GenDiff>::iterator itend) {
    BZET snpdiffs;
    unsigned char *snps = (unsigned char *) malloc(INIT_ALLOC * sizeof(unsigned char));
    memset(snps, 0x00, INIT_ALLOC * sizeof(unsigned char));
    size_t nsnps = 0;
    size_t snpbufsize = INIT_ALLOC * sizeof(unsigned char);

    for (vector<GenDiff>::iterator it = itstart; it != itend; it++) {
        size_t pos = it->pos;
        char rep = it->ref_alt.at(2);
        Base b;
        if (rep == 'A')
            b = B_A;
        else if (rep == 'C')
            b = B_C;
        else if (rep == 'T')
            b = B_T;
        else if (rep == 'G')
            b = B_G;

        if (nsnps == snpbufsize * 4) {
            snps = (unsigned char *)  realloc(snps, snpbufsize * 2);
            if (!snps) {
                cerr << "SNP array realloc failed" << endl;
                exit(1);
            }

            memset(snps + snpbufsize, 0x00, snpbufsize);
            snpbufsize *= 2;
        }

        snpdiffs.seqset(pos);
        snps[nsnps / 4] |= (unsigned char) b << (3 - (nsnps % 4));
        nsnps++;
    }

#ifdef COUTPUT
    cout << itstart->difftype << ",";
    cout << itstart->tag << ",";
    cout << NODE_ELS << ",";
    cout << snpdiffs.size() << ",";
    cout << ((nsnps + 3) / 4) << endl;
    cout.flush();

    /*unsigned char *bzbuf = new unsigned char[snpdiffs.size()];
    snpdiffs.hex(bzbuf);
    fwrite(bzbuf, 1, snpdiffs.size(), stdout);
    delete[] bzbuf;*/
    snpdiffs.exportTo(stdout);
    fwrite(snps, 1, (nsnps + 3) / 4, stdout);
#endif

#ifdef OUTPUT_INDEXES
    string fname = itstart->tag + ".snp.indexes";
    FILE *fi = fopen(fname.c_str(), "w");
    for (vector<GenDiff>::iterator it = itstart; it != itend; it++) {
        char intbuf[20];
        memset(intbuf, 0x00, sizeof(intbuf));
        sprintf(intbuf, "%d", it->pos);
        string buf(intbuf);
        buf += "\n";
        fwrite(buf.c_str(), strlen(buf.c_str()), 1, fi);
    }
    fclose(fi);
#endif

#ifdef VERIFY
    FILE *f = fopen("tempfile.txt", "w+b");
    snpdiffs.exportTo(f);
    rewind(f);
    BZET v;
    v.importFrom(f, snpdiffs.size());
    cerr << "VERIFY: " << v.count() << " diffs" << endl;
    cerr << "v == snpdiffs: " << (v == snpdiffs) << endl;
    cerr << "v size: " << v.size() << endl;
    fseek(f, 0, SEEK_END);
    size_t fsize = ftell(f);
    cerr << "file size: " << fsize << endl;
    assert(v.count() == snpdiffs.count());
    fclose(f);
#endif

#ifdef OUTPUT
    cerr << "Size is " << (snpdiffs.size() + (nsnps + 3) / 4) << endl; 
    cerr << "Number of SNPs: " << snpdiffs.count() << " = " << nsnps << endl;
    cerr << "snpdiffs size: " << snpdiffs.size() << endl;
    cerr << "-> " << itstart->tag << " snpdiff Bzet" << NODE_ELS << " is " << snpdiffs.size();
    fprintf(stderr, " (%.2f KiB)\n", snpdiffs.size() / (double) 1024);
#endif
    return snpdiffs.size() + (nsnps + 3) / 4;
}

size_t encdels(vector<GenDiff>& vgf,
               vector<GenDiff>::iterator itstart,
               vector<GenDiff>::iterator itend) {
    BZET dels;
    for (vector<GenDiff>::iterator it = itstart; it != itend; it++) {
        size_t pos = it->pos;
        size_t ndels = it->ref_alt.find('/');
        for (int i = 0; i < ndels; i++) {
            dels.seqset(pos + i);
        }
    }

#ifdef COUTPUT
    cout << itstart->difftype << ",";
    cout << itstart->tag << ",";
    cout << NODE_ELS << ",";
    cout << dels.size() << ",";
    cout << 0 << endl;
    cout.flush();

    /*unsigned char *bzbuf = new unsigned char[dels.size()];
    dels.hex(bzbuf);
    fwrite(bzbuf, 1, dels.size(), stdout);
    delete[] bzbuf;*/
    dels.exportTo(stdout);
#endif

#ifdef OUTPUT_INDEXES
    string fname = itstart->tag + ".del.indexes";
    FILE *fi = fopen(fname.c_str(), "w");
    for (vector<GenDiff>::iterator it = itstart; it != itend; it++) {
        char intbuf[20];
        memset(intbuf, 0x00, sizeof(intbuf));
        sprintf(intbuf, "%d", it->pos);
        string buf(intbuf);
        buf += "\n";
        fwrite(buf.c_str(), strlen(buf.c_str()), 1, fi);
    }
    fclose(fi);
#endif

#ifdef OUTPUT
    cerr << "Dels size: " << dels.size() << endl;
    cerr << "-> " << itstart->tag << " del Bzet" << NODE_ELS << " is " << dels.size();
    fprintf(stderr, " (%.2f KiB)\n", dels.size() / (double) 1024);
#endif
    return dels.size();
}

size_t writevint(unsigned int n, unsigned char *buf) {
    if (n < 128) {
        *buf = (unsigned char) n + 128;
        return 1;
    }

    size_t bytesreqd = 0;
    unsigned int tempn = n;
    while (tempn > 0) {
        bytesreqd++;
        tempn /= 128;
    }
    
    size_t fac = 128;
    for (size_t i = bytesreqd - 1; i >= 0; i--) {
        buf[i] = (n % fac);
        if (i == bytesreqd - 1)
            buf[i] += 128;
        fac *= 128;
    }    

    return bytesreqd;
}

size_t encins(vector<GenDiff>& vgf,
               vector<GenDiff>::iterator itstart,
               vector<GenDiff>::iterator itend) {
    BZET ins;
    unsigned char *insbuf = (unsigned char *) malloc(INIT_ALLOC * sizeof(unsigned char));
    size_t insbufloc = 0;
    size_t insbufsize = INIT_ALLOC * sizeof(unsigned char);
    for (vector<GenDiff>::iterator it = itstart; it != itend; it++) {
        size_t pos = it->pos;
        size_t nins = it->ref_alt.find('/');
        size_t tmpn = nins;

        size_t bytesreqd = 0;
        while (tmpn > 0) {
            bytesreqd++;
            tmpn /= 128;
        }
        bytesreqd += (nins + 3) / 4;

        if (insbufloc + bytesreqd >= insbufsize) {
            insbuf = (unsigned char *) realloc(insbuf, insbufsize * 2);
            if (!insbuf) {
                cerr << "Unable to realloc in insenc" << endl;
                exit(1);
            }
            insbufsize *= 2;
        }

        size_t offs = writevint(nins, insbuf + insbufloc);
        insbufloc += offs;
        size_t baseswritten = 0;
        unsigned char basebuf = 0x00;
        for (size_t i = 0; i < nins; i++) {
            char insbase = it->ref_alt.at(i);
            if (insbase == 'A') {
                basebuf |= ((unsigned char) B_A << (3 - baseswritten));
            }
            else if (insbase == 'T') {
                basebuf |= ((unsigned char) B_T << (3 - baseswritten));
            }
            else if (insbase == 'C') {
                basebuf |= ((unsigned char) B_C << (3 - baseswritten));
            }
            else if (insbase == 'G') {
                basebuf |= ((unsigned char) B_G << (3 - baseswritten));
            }
            baseswritten = (baseswritten + 1) % 4;
            if (baseswritten == 0) {
                insbuf[insbufloc] = basebuf;
                insbufloc++;
                basebuf = 0x00;
            }
        }
        if (baseswritten != 0) {
            insbuf[insbufloc] = basebuf;
            insbufloc++;
        }

        ins.seqset(pos);
    }

#ifdef COUTPUT
    cout << itstart->difftype << ",";
    cout << itstart->tag << ",";
    cout << NODE_ELS << ",";
    cout << ins.size() << ",";
    cout << insbufloc << endl;
    cout.flush();

    /*unsigned char *bzbuf = new unsigned char[ins.size()];
    ins.hex(bzbuf);
    fwrite(bzbuf, 1, ins.size(), stdout);
    delete[] bzbuf;*/
    ins.exportTo(stdout);
    fwrite(insbuf, 1, insbufloc, stdout);
#endif

#ifdef OUTPUT_INDEXES
    string fname = itstart->tag + ".ins.indexes";
    FILE *fi = fopen(fname.c_str(), "w");
    for (vector<GenDiff>::iterator it = itstart; it != itend; it++) {
        char intbuf[20];
        memset(intbuf, 0x00, sizeof(intbuf));
        sprintf(intbuf, "%d", it->pos);
        string buf(intbuf);
        buf += "\n";
        fwrite(buf.c_str(), strlen(buf.c_str()), 1, fi);
    }
    fclose(fi);
#endif


#ifdef OUTPUT
    cerr << "Ins size: " << (ins.size() + insbufloc) << endl;
    cerr << "Insbufloc: " << insbufloc << endl;
    cerr << "Ins bzet size: " << ins.size() << endl;
    cerr << "-> " << itstart->tag << " ins Bzet" << NODE_ELS << " is " << ins.size();
    fprintf(stderr, " (%.2f KiB)\n", ins.size() / (double) 1024);
#endif
    return ins.size() + insbufloc;
}

bool compdiff(const GenDiff& left, const GenDiff& right) {
    if (left.tag == right.tag) {
        if (left.difftype == right.difftype)
            return (left.pos < right.pos);

        return (left.difftype < right.difftype);
    }

    return (left.tag < right.tag);
}

int main() {
    // Load diffs
    vector<GenDiff> diffs;
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

    // Sort
    //sort(diffs.begin(), diffs.end(), &compdiff);
    clock_t TMS_CLK_TCK = sysconf(_SC_CLK_TCK);
    struct tms tmsstart, tmsend;

    size_t totalsize = 0;
    size_t nsnps = 0, nins = 0, ndel = 0;
    size_t maxdel = 0, maxins = 0;
    string prevtag = diffs.begin()->tag;
    DiffType prevdiff = diffs.begin()->difftype;
    vector<GenDiff>::iterator chrstart = diffs.begin();
    times(&tmsstart);
    for (vector<GenDiff>::iterator it = diffs.begin(); it != diffs.end(); it++) {
        size_t count = it->ref_alt.find('/');
        if (it->difftype == INS) {
            nins += count;
            if (count > maxins)
                maxins = count;
        }
        else if (it->difftype == DEL) {
            ndel += count;
            if (count > maxdel)
                maxdel = count;
        }
        else 
            nsnps++;

        //if (it->difftype == INS)
        //    continue;

        if (it->difftype != prevdiff || it->tag != prevtag || (it + 1) == diffs.end()) {
            if ((it - 1)->difftype == SNP) {
#ifdef OUTPUT
                cerr << "encsnp for tag " << prevtag << endl;
#endif
                totalsize += encsnps(diffs, chrstart, it);
            }
            else if ((it - 1)->difftype == DEL) {
#ifdef OUTPUT
                cerr << "encdel for tag " << prevtag << endl;
#endif
                totalsize += encdels(diffs, chrstart, it);
            }
            else if ((it - 1)->difftype == INS) {
#ifdef OUTPUT
                cerr << "encins for tag " << prevtag << endl;
#endif
                totalsize += encins(diffs, chrstart, it);
            }
            chrstart = it;
            prevtag = it->tag;
            prevdiff = it->difftype;
        }
    }
    times(&tmsend);

    cerr << "Deletions: " << ndel << ", max " << maxdel << endl;
    cerr << "Insertions: " << nins << ", max " << maxins << endl;
    cerr << "SNPs: " << nsnps << endl;
    cerr << "Total: " << (ndel + nins + nsnps) << endl;
    cerr << "End size: " << totalsize << endl;
    cerr << "Bzet" << NODE_ELS << " time: " << ((tmsend.tms_utime - tmsstart.tms_utime) / (double) TMS_CLK_TCK) << endl;
}
