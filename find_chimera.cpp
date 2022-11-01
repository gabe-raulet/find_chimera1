#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <unistd.h>
#include <assert.h>
#include <limits>
#include "kppseq.h"

typedef std::tuple<int, int, int, int> seed_t;
typedef std::tuple<int, int> interval_t;
typedef std::tuple<uint32_t, std::vector<interval_t>> cov_t;

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

int extend_seed
(
    void *Amem,
    const char *qs, /* query sequence */
    const char *ts, /* target sequence */
    int ql, /* query sequence length */
    int tl, /* target sequence length */
    int qb, /* query seed position */
    int tb, /* target seed position */
    int seedlen, /* seed length */
    int mat, /* match reward */
    int mis, /* mismatch penalty */
    int gap, /* gap penalty */
    int xdrop, /* X-drop dropoff value */
    int& q_ext_l, /* coordinate of best left extension on query */
    int& q_ext_r, /* coordinate of best right extension on query */
    int& t_ext_l, /* coordinate of best left extension on target */
    int& t_ext_r  /* coordinate of best right extension on target */
);

int extend_seed_lr
(
    void *Amem, /* memory for anti-diagonals */
    const char *qs, /* query sequence */
    const char *ts, /* target sequence */
    int ql, /* query sequence length */
    int tl, /* target sequence length */
    int qb, /* query begin position */
    int tb, /* target begin position */
    int mat, /* match reward */
    int mis, /* mismatch penalty */
    int gap, /* gap penalty */
    int xdrop, /* X-drop dropoff value */
    bool extleft, /* true means extend left, false means extend right */
    int& q_ext, /* coordinate of best extension on query */
    int& t_ext  /* coordinate of best extension on target */
);

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

    std::vector<cov_t> coverage;
    for (int i = 0; i < store.get_numseqs(); ++i)
    {
        uint32_t len = store.query_seq(i).size();
        coverage.push_back({len, std::vector<interval_t>()});

    }

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

    void *Amem = static_cast<void*>(new int[3*(store.get_maxlen()+1)]);
    uint16_t *covered_bases = new uint16_t[store.get_maxlen()];


    for (auto itr = seeds.begin(); itr != seeds.end(); ++itr)
    {
        int qi = std::get<0>(*itr);
        int ti = std::get<1>(*itr);
        int qb = std::get<2>(*itr);
        int tb = std::get<3>(*itr);
        std::string qs = store.query_seq(qi);
        std::string ts = store.query_seq(ti);
        std::string qn = store.query_name(qi);
        std::string tn = store.query_name(ti);
        int ql = qs.size();
        int tl = ts.size();
        bool rc;
        int q_ext_l;
        int q_ext_r;
        int t_ext_l;
        int t_ext_r;

        std::string qseed = qs.substr(qb, seedlen);
        std::string tseed = ts.substr(tb, seedlen);

        if (qseed == tseed) rc = false;
        else if (qseed == revcomp(tseed)) rc = true;
        else
        {
            std::cerr << qseed << " and " << tseed << " do not match! skipping..." << std::endl;
            continue;
        }

        if (rc) { tseed = revcomp(tseed); tb = tl - tb - seedlen; }

        int score = extend_seed(Amem, qs.c_str(), ts.c_str(), ql, tl, qb, tb, seedlen, mat, mis, gap, xdrop, q_ext_l, q_ext_r, t_ext_l, t_ext_r);

        std::get<1>(coverage[qi]).push_back({q_ext_l, q_ext_r});
        std::get<1>(coverage[ti]).push_back({t_ext_l, t_ext_r});

        // std::cout << std::get<0>(coverage[qi]) << "\t" << std::get<0>(coverage[ti]) << std::endl;
        // std::cout << ql << "\t" << tl << std::endl;
        assert(std::get<0>(coverage[qi]) == static_cast<uint32_t>(ql));
        assert(std::get<0>(coverage[ti]) == static_cast<uint32_t>(tl));

        // char strand = rc? '-' : '+';
        // std::cout << qn << "\t" << ql << "\t" << q_ext_l << "\t" << q_ext_r << "\t" << strand << "\t" << tn << "\t" << tl << "\t" << t_ext_l << "\t" << t_ext_r << "\t" << score << std::endl;
    }

    return 0;
}

int extend_seed
(
    void *Amem,
    const char *qs, /* query sequence */
    const char *ts, /* target sequence */
    int ql, /* query sequence length */
    int tl, /* target sequence length */
    int qb, /* query seed position */
    int tb, /* target seed position */
    int seedlen, /* seed length */
    int mat, /* match reward */
    int mis, /* mismatch penalty */
    int gap, /* gap penalty */
    int xdrop, /* X-drop dropoff value */
    int& q_ext_l, /* coordinate of best left extension on query */
    int& q_ext_r, /* coordinate of best right extension on query */
    int& t_ext_l, /* coordinate of best left extension on target */
    int& t_ext_r  /* coordinate of best right extension on target */
)
{
    int ll = ql+tl+2;

    int rscore, lscore;

    if (qb + seedlen >= ql || tb + seedlen >= tl)
    {
        q_ext_r = qb + seedlen - 1;
        t_ext_r = tb + seedlen - 1;
        rscore = 0;
    }
    else
    {
        rscore = extend_seed_lr(Amem, qs, ts, ql, tl, qb+seedlen, tb+seedlen, mat, mis, gap, xdrop, false, q_ext_r, t_ext_r);
    }

    if (qb <= 0 || tb <= 0)
    {
        q_ext_l = qb;
        t_ext_l = tb;
        lscore = 0;
    }
    else
    {
        lscore = extend_seed_lr(Amem, qs, ts, ql, tl, qb-1, tb-1, mat, mis, gap, xdrop, true, q_ext_l, t_ext_l);
    }

    return rscore + lscore + mat*seedlen;
}

int extend_seed_lr
(
    void *Amem, /* memory for anti-diagonals */
    const char *qs, /* query sequence */
    const char *ts, /* target sequence */
    int ql, /* query sequence length */
    int tl, /* target sequence length */
    int qb, /* query begin position */
    int tb, /* target begin position */
    int mat, /* match reward */
    int mis, /* mismatch penalty */
    int gap, /* gap penalty */
    int xdrop, /* X-drop dropoff value */
    bool extleft, /* true means extend left, false means extend right */
    int& q_ext, /* coordinate of best extension on query */
    int& t_ext  /* coordinate of best extension on target */
)
{
    assert((qb >= 0 && tb >= 0 && qb < ql && tb < tl));

    int m = extleft? qb + 1 : ql - qb; /* number of rows */
    int n = extleft? tb + 1 : tl - tb; /* number of columns */
    int *A = static_cast<int*>(Amem); /* scoring anti-diagonals */
    int undef = std::min((std::numeric_limits<int>::min() / (m+n)) - gap, std::numeric_limits<int>::min() - gap); /* un-explored value (-inf) */
    int best = 0; /* best score seen so far */
    int cur_best = undef; /*best score seen on current anti-diagonal */
    int last_best = undef; /* best score seen on previous anti-diagonal */
    int L = 0; /* starting row index on current anti-diagonal */
    int U = 1; /* ending row index on current anti-diagonal */
    int Lprev = -undef; /* starting row index on previous anti-diagonal */
    int Uprev = undef; /* ending row index on previous anti-diagonal */
    int tr_i = 0; /* best query extension index */
    int tr_j = 0; /* best target extension index */
    int k = 0; /* current anti-diagonal index, i.e. anti-diagonal where i+j == k */
    int *ad1 = A; /* anti-diagonal 2 rounds ago */
    int *ad2 = ad1 + (m+1); /* anti-diagonal 1 round ago */
    int *ad3 = ad2 + (m+1); /* current anti-diagonal */

    ad2[0] = 0;

    do
    {
        k++;
        L = std::max(L, k-n);
        U = std::min(U, m);

        int Lnext = -undef;
        int Unext = undef;
        cur_best = undef;

        for (int i = L; i <= U; ++i)
        {
            int j = k-i;
            int sw = (i < U && ad2[i] > undef)? ad2[i] + gap : undef;
            int sn = (i-1 >= L && ad2[i-1] > undef)? ad2[i-1] + gap : undef;
            int align = extleft? (qs[qb-i+1]==ts[tb-j+1]) : (qs[qb+i-1]==ts[tb+j-1]);
            int snw = (Lprev <= i-1 && i-1 < Uprev && ad1[i-1] > undef)? ad1[i-1] + (align? mat : mis) : undef;

            ad3[i] = std::max({snw, sw, sn});

            if (ad3[i] < best - xdrop) ad3[i] = undef;
            if (ad3[i] > best) tr_i = i, tr_j = j;
            if (ad3[i] > undef && Lnext == -undef) Lnext = i;
            if (ad3[i] > undef && i+1 > Unext) Unext = i+1;

            cur_best = std::max(ad3[i], cur_best);
        }

        if (cur_best == undef && last_best == undef) break;

        best = std::max(best, cur_best);
        Lprev = L;
        Uprev = U;
        L = Lnext;
        U = Unext;
        last_best = cur_best;

        int *tmp = ad1;
        ad1 = ad2;
        ad2 = ad3;
        ad3 = tmp;

    } while (L <= U);

    q_ext = extleft? qb - tr_i + 1 : qb + tr_i - 1;
    t_ext = extleft? tb - tr_j + 1 : tb + tr_j - 1;

    return best;
}
