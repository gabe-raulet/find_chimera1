#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include "kppseq.h"

typedef std::tuple<int, int, int, int> seed_t;
void *Amem;

char comp(char c)
{
    switch (c)
    {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
    }
    return '\0';
}

std::string revcomp(const std::string& s)
{
    std::string s2(s);
    std::transform(s2.cbegin(), s2.cend(), s2.begin(), comp);
    std::reverse(s2.begin(), s2.end());
    return s2;
}

int main(int argc, char *argv[])
{
    int c, mat = 1, mis = -1, gap = -1, xdrop = 15, seedlen = 31, coverage_min = 0;

    while ((c = getopt(argc, argv, "x:A:B:O:l:c:")) >= 0)
    {
        if (c == 'x') xdrop = atoi(optarg);
        else if (c == 'A') mat = atoi(optarg);
        else if (c == 'B') mis = -atoi(optarg);
        else if (c == 'O') gap = -atoi(optarg);
        else if (c == 'l') seedlen = atoi(optarg);
        else if (c == 'c') coverage_min = atoi(optarg);
    }

    if (optind + 2 > argc)
    {
        std::cerr << "Usage: find_chimera [options] <in.fa|in.fq|-> <seeds.tsv>\n\n"
                  << "Options:\n"
                  << "    -x INT    x-drop value [15]\n"
                  << "    -l INT    seed length [31]\n"
                  << "    -A INT    matching score [1]\n"
                  << "    -B INT    mismatch penalty [1]\n"
                  << "    -O INT    gap penalty [1]\n"
                  << "    -c INT    if coverage is minus or equal to this create a gap [0]\n" << std::endl;
        return 1;
    }

    const char *seq_fname = argv[optind++];
    const char *seeds_fname = argv[optind];

    kseq ks = strcmp(seq_fname, "-")? kseq::open(seq_fname) : kseq::open(stdin);
    seqstore store(ks);

    std::vector<seed_t> seeds;
    std::ifstream seedstream(seeds_fname);

    std::string line;
    while (getline(seedstream, line))
    {
        seed_t seed;
        std::istringstream record(line);
        record >> std::get<0>(seed) >> std::get<1>(seed) >> std::get<2>(seed) >> std::get<3>(seed);
        --std::get<0>(seed), --std::get<1>(seed);
        seeds.push_back(seed);
    }

    Amem = static_cast<void*>(new char[3*(store.get_maxlen()+1)]);

    for (auto itr = seeds.begin(); itr != seeds.end(); ++itr)
    {
        std::string qs = store.query_seq(std::get<0>(*itr));
        std::string ts = store.query_seq(std::get<1>(*itr));
        std::string qn = store.query_name(std::get<0>(*itr));
        std::string tn = store.query_name(std::get<1>(*itr));
        int ql = qs.size();
        int tl = ts.size();
        int qb = std::get<2>(*itr);
        int tb = std::get<3>(*itr);
        int rc;
        int q_ext_l;
        int q_ext_r;
        int t_ext_l;
        int t_ext_r;

        std::string qseed = qs.substr(qb, seedlen);
        std::string tseed = ts.substr(tb, seedlen);

        if (qseed != tseed && qseed != revcomp(tseed)) std::cout << qseed << "\t" << tseed << std::endl;
    }

    return 0;
}
