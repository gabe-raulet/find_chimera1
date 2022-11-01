#pragma once

#include <stdio.h>
#include <zlib.h>
#include <ctype.h>
#include <string.h>
#include <string>

template <int BUFLEN>
class kstream
{
private:
    unsigned char *buf;
    int beg, end;
    bool is_eof;
    gzFile f;

public:
    const static int SEP_SPACE = 0; // isspace(): \t, \n, \v, \f, \r
    const static int SEP_TAB = 1;   // isspace() && !' '
    const static int SEP_LINE = 2;  // '\n'
    const static int SEP_MAX = 2;

    kstream(gzFile f) : f(f), beg(0), end(0), is_eof(false), buf(new unsigned char[BUFLEN]) {}
    ~kstream() { delete[] buf; }

    int getc()
    {
        if (is_eof && beg >= end) return -1;
        if (beg >= end) {
            beg = 0;
            end = gzread(f, buf, BUFLEN);
            if (end < BUFLEN) is_eof = true;
            if (end == 0) return -1;
        }

        return (int)buf[beg++];
    }

    int getuntil(int delimiter, std::string& s, int* dret, bool append = false)
    {
        if (dret) *dret = 0;
        if (!append) s.clear();
        if (beg >= end && is_eof) return -1;

        for(;;)
        {
            int i;
            if (beg >= end) {
                if (!is_eof) {
                    beg = 0;
                    end = gzread(f, buf, BUFLEN);
                    if (end < BUFLEN) is_eof = true;
                    if (end == 0) break;
                } else break;
            }

            if (delimiter == SEP_LINE) {
                unsigned char *sep = static_cast<unsigned char*>(memchr(buf + beg, '\n', end - beg));
                i = sep? sep - buf : end;
            } else if (delimiter > SEP_MAX) {
                for (i = beg; i < end; ++i)
                    if (buf[i] == delimiter) break;
            } else if (delimiter == SEP_SPACE) {
                for (i = beg; i < end; ++i)
                    if (isspace(buf[i])) break;
            } else if (delimiter == SEP_TAB) {
                if (isspace(buf[i]) && buf[i] != ' ') break;
            } else i = 0; /* never come to here! */

            s.append(reinterpret_cast<char*>(buf + beg), i - beg);
            beg = i + 1;
            if (i < end) {
                if (dret) *dret = buf[i];
                break;
            }
        }

        return s.size();
    }
};

class kseq
{
private:
    std::string name, comment, seq, qual;
    int last_char;
    bool is_fastq;
    kstream<16384> ks;

public:
    kseq(gzFile fd) : ks(fd), is_fastq(false) {}

    const std::string& getname() const { return name; }
    const std::string& getcomment() const { return comment; }
    const std::string& getseq() const { return seq; }
    const std::string& getqual() const { return qual; }
    const bool isfastq() const { return is_fastq; }

    static kseq open(FILE *f) { return kseq(gzdopen(fileno(f), "r")); }
    static kseq open(std::string fname) { return kseq(gzopen(fname.c_str(), "r")); }

    int read()
    {
        int c, r;
        if (last_char == 0) {
            while ((c = ks.getc()) >= 0 && c != '>' && c != '@');
            if (c < 0) return c;
            last_char = c;
        }

        name.clear(), comment.clear(), seq.clear(), qual.clear();

        if ((r = ks.getuntil(ks.SEP_SPACE, name, &c)) < 0) return r;
        if (c != '\n') ks.getuntil(ks.SEP_LINE, comment, nullptr);
        while ((c = ks.getc()) >= 0 && c != '>' && c != '+' && c != '@') {
            if (c == '\n') continue;
            seq.push_back(c);
            ks.getuntil(ks.SEP_LINE, seq, nullptr, true);
        }
        if (c == '>' || c == '@') last_char = c;
        is_fastq = (c == '+');
        if (!is_fastq) return seq.size();
        while ((c = ks.getc()) >= 0 && c != '\n');
        if (c < 0) return c;
        while ((c = ks.getuntil(ks.SEP_LINE, qual, nullptr, true)) >= 0 && qual.size() < seq.size());
        if (c < 0) return c;
        last_char = 0;
        if (qual.size() != seq.size()) return -1;
        return seq.size();
    }
};

#include <unordered_map>
#include <vector>

class seqstore
{
private:
    std::vector<std::string> names;
    std::vector<std::string> comments;
    std::vector<std::string> seqs;
    std::vector<std::string> quals;
    std::unordered_map<std::string, int> idmap;
    int maxlen;

public:
    seqstore(kseq& ks) : maxlen(0)
    {
        int i = 0;
        while (ks.read() >= 0)
        {
            names.push_back(ks.getname());
            comments.push_back(ks.getcomment());
            seqs.push_back(ks.getseq());
            if (ks.isfastq()) quals.push_back(ks.getqual());
            idmap.insert({ks.getname(), i++});
            maxlen = ks.getseq().size() > maxlen? ks.getseq().size() : maxlen;
        }
    }

    int query_id(const std::string& name)
    {
        auto q = idmap.find(name);
        return q != idmap.end()? static_cast<int>(q->second) : -1;
    }

    int get_maxlen() const { return maxlen; }

    const std::string& query_name(int id) { return names[id]; }
    const std::string& query_comment(int id) { return comments[id]; }
    const std::string& query_seq(int id) { return seqs[id]; }
    const std::string& query_qual(int id) { return quals[id]; }

    const std::string& query_comment(const std::string& name) { return query_comment(query_id(name)); }
    const std::string& query_seq(const std::string& name) { return query_seq(query_id(name)); }
    const std::string& query_qual(const std::string& name) { return query_qual(query_id(name)); }
};
