#include <iostream>
using namespace std;

#include <sstream>

#include <cstdlib>

#include <cctype>
using std::toupper;

#include <cstring>

#include <string>
using std::string;

#include <sstream>

#include <getopt.h>

#include <fstream>
using std::ifstream;

#include <cmath>
using std::log;
using std::floor;


#include <vector>
using std::vector;

#include <algorithm>
using std::sort;
using std::reverse;


#include <google/sparse_hash_map> //memmory efficent but relatively slow
using google::sparse_hash_map;

#include <bitset>

#include <zlib.h>
#include <stdexcept>
#include <bitset>
#include <climits>


///2013.3.13, no calling external hash
using std::hash;

/* or
using ext::hash;reverse
using __gnu_cxx::hash;
*/



///Found these here http://mail-archives.apache.org/mod_mbox/trafficserver-dev/201110.mbox/%3CCACJPjhYf=+br1W39vyazP=ix
///eQZ-4Gh9-U6TtiEdReG3S4ZZng@mail.gmail.com%3E
#define MOD_GZIP_ZLIB_WINDOWSIZE 15
#define MOD_GZIP_ZLIB_CFACTOR    9
#define MOD_GZIP_ZLIB_BSIZE      8096

///Found this one here: http://panthema.net/2007/0328-ZLibString.html, author is Timo Bingmann
/** Compress a STL string using zlib with given compression level and return
  * the binary data. */
std::string compress_deflate(const std::string& str, int compressionlevel = Z_BEST_COMPRESSION)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();           // set the z_stream's input

    int ret=0;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw(std::runtime_error(oss.str()));
    }

    //cout << outstring.size() << endl;

    return outstring;
}


/* //to test and verifiy outputs
std::string decompress_deflate(const std::string& str)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (inflateInit(&zs) != Z_OK)
        throw(std::runtime_error("inflateInit failed while decompressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // get the decompressed bytes blockwise using repeated calls to inflate
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = inflate(&zs, 0);

        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer, zs.total_out - outstring.size());
        }

    } while (ret == Z_OK);

    inflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib decompression: (" << ret << ") "
            << zs.msg;
        throw(std::runtime_error(oss.str()));
    }

    return outstring;
}
*/


struct compare_char
{
    bool operator()(char c_1, char c_2) const
    {
        return (c_1==c_2);
        //return ((c_1==c_2) || (c_1 && c_2 && strcmp(c_1, c_2)== 0));

    }

};


struct compare_string //compare string
{
    bool operator()(string s1, string s2) const
    {
        return (s1 == s2); //|| (s1 && s2 && strcmp(s1, s2) == 0);

    }
};



string rev_comp_str_convert(string &forward_index_str)
{
    string rev_comp_str="";

    for (string::iterator it=forward_index_str.begin(); it!=forward_index_str.end(); ++it)
    {
        switch (*it)
        {
            case 'A':
                rev_comp_str+="T";
                break;

            case 'T':
                rev_comp_str+="A";
                break;

            case 'G':
                rev_comp_str+="C";
                break;

            case 'C':
                rev_comp_str+="G";
                break;

            case 'R':
                rev_comp_str+="Y";
                break;

            case 'Y':
                rev_comp_str+="R";

        }

    }

    reverse(rev_comp_str.begin(), rev_comp_str.end());

    return rev_comp_str;

}


char ry_convert(char read_char)
{
    if (read_char=='A' or read_char=='G')
    {
        return 'R'; ///purine

    } else if (read_char=='C' or read_char=='T')
    {
        return 'Y'; ///pyrimidine

    } else
    {
        return read_char;

    }
}


//define key_str
//nucleotide; AGCT
//amino acid; GASTCVLIMPFYWDENQHKR
void key_hash_register(string &key_str, sparse_hash_map<char, long long, hash<char>, compare_char> &key_reg_hash)
{
    long long key_index_cnt=0;

    key_reg_hash.resize(key_str.length());

    ///register key hash
    for (string::iterator it=key_str.begin(); it!=key_str.end(); ++it)
    {
        key_reg_hash[*it] = key_index_cnt;
        key_index_cnt++; //key index start from 0 ~
    }


}


///count and aribtrary create tokens to reduce the size of feature strings
string apply_word_token(sparse_hash_map<string, long long, hash<string>, compare_string> &compress_token_hash, string &str_block)
{
    bitset<8> char_replace=0;
    long long delimit_pos=0;
    string str_block_bk="";
    int str_compress_rate_length = int(str_block.length() * 0.25); //for speed up 1/4 size shrink as threshold

    for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator it=compress_token_hash.begin(); it!=compress_token_hash.end(); ++it)
    {
        str_block_bk="";

        while (str_block!=str_block_bk && str_compress_rate_length <= str_block.length()) //maximum compression, but speed decrease. Speed matters.
        {
            str_block_bk = str_block;

            delimit_pos = str_block.find(it->first);

            if (delimit_pos!=-1) //if found
            {
                char_replace = bitset<8>(it->second);
                char_replace[7]=1; //signature
                //replace
                str_block.replace(delimit_pos, (it->first).size(), reinterpret_cast<char *>(&char_replace));

            }

        }

    }


    return str_block;
}


void compress_word_token(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash,
                            sparse_hash_map<string, long long, hash<string>, compare_string> &compress_token_hash, bool ratio_flag,
                            double &feature_hit_cnt, int str_key_frag_size, stringstream &output_stream)
{
    char* stream_buf = NULL;//char stream_buf[70]; //give enough space
    size_t stream_size_t=0;

    string str_key_frag=""; //fragment of str_key(as a whole)
    string str_value="";
    long long min_frequency=0;

    for (sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>::iterator it=prim_index_hash.begin(); it!=prim_index_hash.end(); ++it)
    {
        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=(it->second).begin(); sub_it!=(it->second).end(); ++sub_it)
        {
            str_key_frag = (sub_it->first).substr(0, str_key_frag_size); //get less(or equal) than half segment

            if (ratio_flag==true)
            {
                stream_size_t=snprintf(NULL, 0, "%.4e", (double)(sub_it->second) / feature_hit_cnt);
                stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));
                snprintf(stream_buf, size_t(stream_buf), "%.4e", (double)(sub_it->second) / feature_hit_cnt);

            } else //convert longlong to string type
            {
                stream_size_t=snprintf(NULL, 0, "%lld", sub_it->second);
                stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));
                snprintf(stream_buf, size_t(stream_buf), "%lld", sub_it->second);

            }


            if (compress_token_hash.size()>compress_token_hash.max_size()-3)
            {
                break;

            }


            compress_token_hash[str_key_frag]++; //should check size limit
            compress_token_hash[string(stream_buf)]++; //should check size limit

        }

    }

    ///sort and remove low count tokens; maximum number of tokens = 128
    long long min_value=1;
    compress_token_hash.set_deleted_key(string()); //causing problem here

    while (compress_token_hash.size() > 128)
    {
        //find smallest
        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator it=compress_token_hash.begin(); it!=compress_token_hash.end(); it++)
        {
            if (it->second==min_value)
            {
                //cout << it->first << "\t" << it->second << endl;
                compress_token_hash.erase(it);
            }

            if (compress_token_hash.size()==128)
            {
                break;

            }

        }

        min_value++;

    }


    compress_token_hash.clear_deleted_key();
    compress_token_hash.resize(0); //compact the hashtable because erase doesn't remove(element) nor reduce hashtable

    //replace original value to (long long) index order and print?
    long long order_cnt=0;


    for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator it=compress_token_hash.begin(); it!=compress_token_hash.end(); it++)
    {
        it->second = order_cnt;

        stream_size_t=snprintf(NULL, 0, "#%s|%lld\n", (it->first).c_str(), order_cnt);
        stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char));
        //str_compress_key = (string)stream_buf;

        snprintf(stream_buf, size_t(stream_buf), "#%s|%lld\n", (it->first).c_str(), order_cnt);
        //output_stream << string(stream_buf);
        output_stream << stream_buf;

        order_cnt++;

    }

    free(stream_buf);
    stream_buf=NULL;

}

double feature_string_entropy(string &key_str, string &str_key) ///determine feature entropy
{

    vector<int> str_count_vector;
    str_count_vector.resize(key_str.size());

    for (string::iterator it=str_key.begin(); it!=str_key.end(); ++it)
    {
        str_count_vector[key_str.find(*it)]++;

    }

    double str_entropy=0.0;
    double each_frequency=0.0;
    int key_str_length = key_str.size();
    int str_key_length = str_key.size();

    for (vector<int>::iterator it=str_count_vector.begin(); it!=str_count_vector.end(); ++it)
    {
        if (*it!=0)
        {
            each_frequency=(double)(*it) / str_key_length;
            str_entropy-=each_frequency * (log(each_frequency) / log(key_str_length)); ///log base is the size of key_str. possible alphabets


        } else
        {
            each_frequency=0;
        }

    }

    return str_entropy;

}


///filter feature counts
void filter_feature_count(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash,
    long long bottom_count_limit, long long top_count_limit, double &feature_hit_cnt)
{
    for (sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>::iterator it=prim_index_hash.begin(); it!=prim_index_hash.end(); ++it)
    {
        (it->second).set_deleted_key(string());

        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=(it->second).begin(); sub_it!=(it->second).end(); ++sub_it)
        {
            if ((sub_it->second < bottom_count_limit) || (sub_it->second > top_count_limit && top_count_limit!=0))
            {
                feature_hit_cnt-=(sub_it->second); //substract selected feature count from total
                (it->second).erase(sub_it); //or (it->second).erase(sub_it->first); //erase by key

            }

        }

        (it->second).clear_deleted_key();
        (it->second).resize(0); //compact the hashtable because erase doesn't remove(element) nor reduce hashtable.

    }
    //prim_index_hash.resize(0); //necessary?

}


void feature_container_input(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string> &prim_index_hash,
    vector<string> &max_index_key_vector, string str_key)
{
    bool key_insert_flag=false;

    for (vector<string>::iterator it=max_index_key_vector.begin(); it!=max_index_key_vector.end(); ++it)
    {
        if ((str_key == *it) || (str_key < *it && prim_index_hash[*it].size() < INT_MAX)) ///alternatives?: prim_index_hash[*it].max_size(), INT_MAX, LONG_MAX
        {
            prim_index_hash[*it][str_key]++;
            //cout << str_key << "\t" << max_index_key_vector.size() << endl;
            key_insert_flag=true;
            break;

        } else if (str_key < *it)
        {
            prim_index_hash[str_key].resize(INT_MAX); ///reserve max capacity
            prim_index_hash[str_key][str_key]=1;

            prim_index_hash[*it].set_deleted_key(string()); //or string()

            for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=prim_index_hash[*it].begin(); sub_it!=prim_index_hash[*it].end(); ++sub_it)
            {
                if (sub_it->first <= str_key)
                {
                    prim_index_hash[str_key][sub_it->first]+=sub_it->second; //+= for self add
                    prim_index_hash[*it].erase(sub_it); //or prim_index_hash[*it].erase(sub_it->first)

                }

            }

            ///compact hashtable *it to a smllest valid size
            prim_index_hash[*it].clear_deleted_key(); //off delete key call
            prim_index_hash[*it].resize(0); //compact the hashtable because erase doesn't remove(element) nor reduce hashtable, invalidate iterator

            max_index_key_vector.push_back(str_key); ///add new represent key to front, and deque iteration become invalid? and invalid after all?

            key_insert_flag=true;
            break;

        }

    }


    if (key_insert_flag==false) ///new feature
    {
        prim_index_hash[str_key].resize(INT_MAX); ///reserve max capacity
        prim_index_hash[str_key][str_key]=1; //<string, long long>
        
        max_index_key_vector.push_back(str_key); //in vector, push_back does not invalidate iterator

    }

    sort(max_index_key_vector.begin(), max_index_key_vector.end()); //sorting invalidate iterator

}


void feature_container_output(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string> &prim_index_hash,
    vector<string> &max_index_key_vector, sparse_hash_map<string, long long, hash<string>, compare_string> &compress_token_hash,
    bool ratio_output_flag, int str_token_size, double &feature_hit_cnt, stringstream &output_stream)
{
    char* stream_buf=NULL;
    size_t stream_size_t=0;
    string str_block="";

    vector<string> sub_index_vector;

    prim_index_hash.set_deleted_key(string());

    for (vector<string>::iterator it=max_index_key_vector.begin(); it!=max_index_key_vector.end(); ++it)
    {
        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=prim_index_hash[*it].begin(); sub_it!=prim_index_hash[*it].end(); ++sub_it)
        {
            sub_index_vector.push_back(sub_it->first);

        }

        sort(sub_index_vector.begin(), sub_index_vector.end());
        prim_index_hash[*it].set_deleted_key(string());

        for (vector<string>::iterator vector_it=sub_index_vector.begin(); vector_it!=sub_index_vector.end(); ++vector_it)
        {
            if (ratio_output_flag==true)
            {
                //printf("org; %s|%.4e\n", (*vector_it).c_str(), (double)(prim_index_hash[*it])[*vector_it] / feature_hit_cnt);
                stream_size_t = snprintf(NULL, 0, "%s|%.4e\n", (*vector_it).c_str(), (double)(prim_index_hash[*it])[*vector_it] / feature_hit_cnt);
                stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

                snprintf(stream_buf, size_t(stream_buf), "%s|%.4e\n", (*vector_it).c_str(), (double)(prim_index_hash[*it])[*vector_it] / feature_hit_cnt);

            } else
            {
                //printf("org; %s|%lld\n", (*vector_it).c_str(), (prim_index_hash[*it])[*vector_it]);
                stream_size_t = snprintf(NULL, 0, "%s|%lld\n", (*vector_it).c_str(), (prim_index_hash[*it])[*vector_it]);
                stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

                snprintf(stream_buf, size_t(stream_buf), "%s|%lld\n", (*vector_it).c_str(), (prim_index_hash[*it])[*vector_it]);

            }

            str_block = string(stream_buf);

            output_stream << apply_word_token(compress_token_hash, str_block); ///compress a string using prepared tokens

            prim_index_hash[*it].erase(*vector_it); //erase element after output for decrease next hash search time

        }

        prim_index_hash[*it].clear_deleted_key();
        prim_index_hash[*it].resize(0); //compact the hashtable because erase doesn't remove(element) nor reduce hashtable.

        prim_index_hash.erase(*it); //by const key

        sub_index_vector.clear();

    }

    free(stream_buf);
    stream_buf=NULL;

    prim_index_hash.clear_deleted_key();
    prim_index_hash.resize(0); //compact the hashtable because erase doesn't remove(element) nor reduce hashtable.
}



///2012.6.26, compact version, separate with output function
int seq_read_window(stringstream &read_f, int feature_length, bool backward_flag, double &feature_hit_cnt, bool unmasked_treat_flag, bool ratio_output_flag, bool ry_code_flag,
    sparse_hash_map<char, long long, hash<char>, compare_char> &key_reg_hash, sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash,
    vector<string> &max_index_key_vector, string &key_str, double bottom_entropy_limit, double top_entropy_limit)
{

    string forward_index_str="";
    string backward_index_str="";
    string key_index_str=""; ///final feature key product after conversions

    char read_char='0';
    double str_entropy=0.0;

    bool fasta_head_flag=false;

    while (read_f.get(read_char))
    {
        if (read_char=='>' || read_f.eof())
        {
            fasta_head_flag=true;
            forward_index_str.clear();

        } else if (fasta_head_flag==true && read_char=='\n')
        {
            fasta_head_flag=false;

        } else if (fasta_head_flag==false && read_char!='\n')
        {

            if (unmasked_treat_flag==true) ///accept lower case letters; softmasked
            {
                read_char = toupper(read_char);

            }

            if (ry_code_flag==true) ///convert AGCT to RY code
            {
                read_char = ry_convert(read_char);

            }

            if(key_reg_hash.find(read_char)!=key_reg_hash.end())
            {
                forward_index_str+=read_char;

                if (forward_index_str.length() == feature_length)
                {

                    if (bottom_entropy_limit==0.0 && top_entropy_limit==1.0)
                    {
                        str_entropy = 0.5; ///force fix; continue

                    } else ///or calculate
                    {
                        str_entropy = feature_string_entropy(key_str, forward_index_str);

                    }

                    if (str_entropy >= bottom_entropy_limit && str_entropy <= top_entropy_limit)
                    {
                        if (backward_flag==true) ///consider complement strands
                        {
                            backward_index_str = rev_comp_str_convert(forward_index_str);

                            if (forward_index_str > backward_index_str) ///select based on feature lexicographical order (the one comes first)
                            {
                                key_index_str = backward_index_str;

                            } else
                            {
                                key_index_str = forward_index_str;

                            }

                        } else ///no complement strands
                        {
                            key_index_str = forward_index_str;

                        }

                        feature_container_input(prim_index_hash, max_index_key_vector, key_index_str);

                        feature_hit_cnt+=1;

                    }

                    forward_index_str.erase(0, 1);
                    //cout << forward_index_str << endl;
                }


            } else
            {
                forward_index_str.clear();
                backward_index_str.clear();

            }

        }

    }
    forward_index_str.clear();


    if (!prim_index_hash.empty())
    {
        return 1;

    }

    return 0; //if prim_index_hash is empty

}




long long vocab_size_measure(sparse_hash_map<string, sparse_hash_map<string,  long long, hash<string>, compare_string>, hash<string>, compare_string> &prim_index_hash)
{
    long long vocab_size=0; //long long or double?

    for (sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>::iterator it=prim_index_hash.begin(); it!=prim_index_hash.end(); ++it)
    {
        vocab_size+=(long long)(it->second).size();

    }

    return vocab_size;

}



unsigned long count_bits_per_alphabet(int n_key_type)
{
    unsigned long bits_per_alphabet=0;

    while(n_key_type > pow(2, bits_per_alphabet))
    {
        bits_per_alphabet++;

    }

    return (bits_per_alphabet);

}


void file_to_string_blocks(vector<string> &str_block_vector, ifstream &read_f)
{
    string read_buf="";

    while (getline(read_f, read_buf, '>')) ///'>' fasta header as a delimiter
    {
        if (!read_buf.empty())
        {
            str_block_vector.push_back(">"+read_buf);
            read_buf=string(); //clear

        }

    }

}


void show_help()
{
    cout << "Parameter usage [option][load_path][save_path]" << endl;
    cout << "-h, show_help()" << endl;
    cout << "-v, show_version()" << endl;
    cout << "-s [INT], feature size(dependency on memory size and bit size)" << endl;
    cout << "-e [INT], feature size end, limit feature_length++ while measuring vocabulary size" << endl;


    cout << endl << "Extend option;" << endl;
    cout << "-a, input is amino acid sequence(default is AGCT nucleotide)" << endl;
    cout << "-c, convert nucleotide (AGCT) code to RY code" << endl;
    cout << "-k [STR], manual define of key bases, ex; 'HJKL' -> recognize as (H, J, K, L), set of character without duplicates" << endl;
    cout << "-r, disable reverse complement counting(given -a will turn off reverse complement)" << endl;
    cout << "-n, output freqeuncy into ratio(of sum of all frequency)" << endl;
    cout << "-u, unmasked_treat, also accept lower letters(default; accept only upper letters)" << endl;

    cout << endl << "Vocabulary size measure;" << endl;
    cout << "-V, measure vocabular size; will stop when total size decreases; affect by 'feature count limit' and 'feature entropy limit' condition" << endl;

    cout << endl << "feature count limit (1 ~ 0 = maximum); remove a feature count below [-b] or above [-t]" << endl;
    cout << "-b [int], bottom_count_limit(default = 1)" << endl;
    cout << "-t [int], top_count_limit(default = 0 = maximum)" << endl;

    cout << endl << "feature entropy limit (0.0 ~ 1.0); remove feature string entropy below [-B] or above [-T]" << endl;
    cout << "-B [double], bottom_entropy_limit(default = 0.0)" << endl;
    cout << "-T [double], top_entropy_limit(default = 1.0)" << endl;
}


void show_profile()
{
    cout << "FFP string version 2v.2.1 (2018-9-21)\n";
    cout << "Value presentation: a point below 4 decimal places (%.4e)\n";

    cout << "Code by JaeJin Choi; https://github.com/jaejinchoi/FFP\n";
    cout << "compile; g++ -std=c++11 -o (output) (this script) -lz\n";
    cout << "Required; google-sparse-hash\n";
    cout << "Required; zlib 1.2.8+\n";

}



int main(int argc, char** argv)
{
    int feature_length=10; //default
    int feature_length_end=0;

    bool ratio_output_flag=false;
    bool backward_flag=true;
    bool max_vocab_find_flag=false;
    bool ry_code_flag=false;


    ///feature count range (1.0 ~ 0.0 = maximum)
    long long bottom_count_limit=1; ///in general, use bottom_count_limit=2 to find a point where vocabulary start to diverge
    long long top_count_limit=0; ///is a maximum

    ///feature entropy range (0.0~1.0)
    double top_entropy_limit=1.0;
    double bottom_entropy_limit=0.0;

    bool unmasked_treat_flag=false; ///default, the program does not accept lower case letters

    string key_str="AGCT"; ///default, nucleotides

    int opt;

    while ((opt = getopt(argc, argv, "hvs:e:ack:rnub:t:VB:T:")) !=EOF) // EOF = -1
    {
        switch (opt)
        {
            case 'h':
                show_help();
                exit(EXIT_SUCCESS);

            case 'v':
                show_profile();
                exit(1);

            case 's':
                feature_length = atoi(optarg);
                break;

            case 'e':
                feature_length_end = atoi(optarg);
                break;

            case 'a':
                key_str = "GASTCVLIMPFYWDENQHKR"; ///default, amino acids
                break;

            case 'c':
                ry_code_flag=true;
                break;

            case 'k': //user defined key_str
                key_str = optarg;
                break;

            case 'r':
                backward_flag=false;
                break;

            case 'n':
                ratio_output_flag=true;
                break;

            case 'u':
                unmasked_treat_flag=true;
                break;

            case 'b':
                bottom_count_limit=atoi(optarg);
                break;

            case 't':
                top_count_limit=atoi(optarg);
                break;

            case 'B':
                bottom_entropy_limit=atof(optarg);
                break;

            case 'T':
                top_entropy_limit=atof(optarg);
                break;

            case 'V':
                max_vocab_find_flag=true;
                break;

            default:
                show_help();
                exit(0);

        }
    }


    ///only nucleotide have reverse_complement option-double strand
    if (key_str!="AGCT") //including amino acides but not RY code
    {
        backward_flag=false;

    } else if (ry_code_flag==true)
    {
        key_str="RY";

    }


    if (bottom_count_limit > top_count_limit && top_count_limit!=0)
    {
        cout << "Wrong feature count filter boundry\n";
        cout << "bottom_count_limit <= top_count_limit\n";
        exit(1);

    } else if (bottom_entropy_limit > top_entropy_limit && top_entropy_limit==1.0)
    {
        cout << "Wrong feature entropy filter boundry\n";
        cout << "bottom_entropy_limit <= top_entropy_limit\n";
        exit(1);

    }

    ///what is the optimal token size to compress original feature strings?
    int str_token_size = min(feature_length, max(2, int(log(feature_length) / log(key_str.length()))));
    sparse_hash_map<string, long long, hash<string>, compare_string> compress_token_hash; ///str_token hash

    ///prepare necessary parameters / declare variables
    sparse_hash_map<char, long long, hash<char>, compare_char> key_reg_hash;

    sort(key_str.begin(), key_str.end()); ///sort given key_str by lexicographical order
    key_hash_register(key_str, key_reg_hash);

    sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string> prim_index_hash;

    vector<string> max_index_key_vector;

    ifstream read_f;
    vector<string> str_block_vector; ///store string blocks of a file

    ofstream write_f; //profile save

    stringstream read_str_f(ios::in|ios::out);
    stringstream output_stream(ios::out|ios::in|ios::ate);

    string output_total_string="";

    double feature_hit_cnt=0.0;
    long long past_vocab_size=0;
    long long recent_vocab_size=0;


    char* stream_buf = NULL;
    size_t stream_size_t=0;


    if (argc - optind == 2) ///[input file path][save file path]
    {
        read_f.open(argv[optind], ios::in); //read text file
        write_f.open(argv[optind+1], ios::out|ios::binary); //save binary file

        if (!read_f)
        {
            cout << "Unable to open file: " << argv[optind] << endl;
            return 1;

        } else
        {

            file_to_string_blocks(str_block_vector, read_f); ///dealing large files; fragmentate into each fasta
            read_f.close();

            while (1)
            {
                for (vector<string>::iterator it=str_block_vector.begin(); it!=str_block_vector.end(); ++it)
                {
                    read_str_f.rdbuf()->pubsetbuf(&(*it)[0], (*it).size()); //load content
                    read_str_f.seekg(0, ios_base::beg); //go to the beginning

                    seq_read_window(read_str_f, feature_length, backward_flag, feature_hit_cnt, unmasked_treat_flag,
                                        ratio_output_flag, ry_code_flag, key_reg_hash, prim_index_hash, max_index_key_vector,
                                        key_str, bottom_entropy_limit, top_entropy_limit);

                    read_str_f.str(string()); //clear content
                    read_str_f.clear(); //clear flag

                }

                ///feature count filter
                if (bottom_count_limit>1 || top_count_limit!=0)
                {
                    filter_feature_count(prim_index_hash, bottom_count_limit, top_count_limit, feature_hit_cnt);

                }


                if (max_vocab_find_flag==true && !prim_index_hash.empty())
                {
                    recent_vocab_size = vocab_size_measure(prim_index_hash);

                    stream_size_t = snprintf(NULL, 0, "l-mer; %d, vocab_size; %lld\n", feature_length, recent_vocab_size); //find string(c) size
                    stream_buf = (char*)realloc(stream_buf, (stream_size_t + 1)*sizeof(char *)); //+1 at the end(/0)
                    snprintf(stream_buf, size_t(stream_buf), "l-mer; %d, vocab_size; %lld\n", feature_length, recent_vocab_size); //format char*

                    output_stream << stream_buf;

                    if ((feature_length < feature_length_end || feature_length_end==0) && past_vocab_size <= recent_vocab_size)
                    {
                        ///prepare and clear variables and containers for next cycle
                        max_index_key_vector.clear();
                        prim_index_hash.clear();
                        prim_index_hash.resize(0);
                        feature_hit_cnt=0;

                        feature_length++;
                        
                        past_vocab_size = recent_vocab_size;

                    } else
                    {
                        write_f << output_stream.str(); //<< '\n'; //additional linebreak un necessary
                        break; //break while()

                    }

                } else if (max_vocab_find_flag==false && !prim_index_hash.empty()) ///write FFP to file
                {
                    compress_word_token(prim_index_hash, compress_token_hash, ratio_output_flag, feature_hit_cnt, str_token_size, output_stream); ///create tokens to reduce overall data size

                    feature_container_output(prim_index_hash, max_index_key_vector, compress_token_hash, ratio_output_flag, str_token_size, feature_hit_cnt, output_stream);

                    write_f << compress_deflate(output_stream.str(), Z_BEST_COMPRESSION);
                    //write_f << output_stream.str(); //<< '\n'; //without zlib compression
                    break; //break while(1)
                    
                } else
                {
                    cout << "Empty FFP ouput: " << argv[optind] << endl;
                    break;
                    
                }

            }

            ///clear containers/initialize
            max_index_key_vector.clear();
            prim_index_hash.clear();
            prim_index_hash.resize(0);
            feature_hit_cnt=0;

            output_stream.str(string());
            output_stream.clear();

            str_block_vector.clear();

        }

        write_f.close();


    } else
    {
        cout << "Please input load_path, and then save_path" << endl;
        exit(EXIT_SUCCESS);

					
    }

    free(stream_buf);
    stream_buf=NULL;
}



