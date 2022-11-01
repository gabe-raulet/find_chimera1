#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include "kppseq.h"

typedef std::tuple<int, int, int, int> seed_t;

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

    return 0;
}
