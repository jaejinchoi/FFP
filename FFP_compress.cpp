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

//#include <deque>
//using std::deque;

#include <vector>
using std::vector;

#include <algorithm>
using std::sort;
using std::reverse;
//for using sort of deque

#include <google/sparse_hash_map> //memmory efficent but relatively slow
using google::sparse_hash_map;

#include <bitset>

#include <zlib.h>
#include <stdexcept>
#include <bitset>


///2013.3.13, no calling external hash
using std::hash;

/* or
using ext::hash;
using __gnu_cxx::hash;
*/


// Found these here http://mail-archives.apache.org/mod_mbox/trafficserver-dev/201110.mbox/%3CCACJPjhYf=+br1W39vyazP=ix
//eQZ-4Gh9-U6TtiEdReG3S4ZZng@mail.gmail.com%3E
#define MOD_GZIP_ZLIB_WINDOWSIZE 15
#define MOD_GZIP_ZLIB_CFACTOR    9
#define MOD_GZIP_ZLIB_BSIZE      8096


// Found this one here: http://panthema.net/2007/0328-ZLibString.html, author is Timo Bingmann
/** Compress a STL string using zlib with given compression level and return
  * the binary data. */
//std::string compress_deflate(const std::string& str, int compressionlevel = Z_BEST_COMPRESSION)
std::string compress_deflate(const std::string& str, int compressionlevel = Z_BEST_COMPRESSION)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();           // set the z_stream's input

    int ret;
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



struct compare_string //compare string
{
    bool operator()(string s1, string s2) const
    {
        return (s1 == s2); //|| (s1 && s2 && strcmp(s1, s2) == 0);
    }
}; //should put ; after struct declare


struct compare_char
{
    bool operator()(char c_1, char c_2) const
    {
        return (c_1==c_2);
        //return ((c_1==c_2) || (c_1 && c_2 && strcmp(c_1, c_2)== 0));

    }

};


//define key_str
//nucleotide; AGCT
//amino acid; GASTCVLIMPFYWDENQHKR
void register_key_hash(string &key_str, sparse_hash_map<char, long long, hash<char>, compare_char> &key_reg_hash)
{
    long long key_index_cnt=0;
    ///register key hash
    for (string::iterator it=key_str.begin(); it!=key_str.end(); ++it)
    {
        key_reg_hash[*it] = key_index_cnt;
        key_index_cnt++; //key index start from 0 ~
    }

    //determine length of key_index(by max size) = key_index_length, key_index_cnt>0
    //return (int)(floor(log(key_index_cnt)/log(10)+1));

}



long long key_hash_check(char &read_char, sparse_hash_map<char, long long, hash<char>, compare_char> &key_reg_hash)
{
    long long key_index_value=-1; //default; -1 = false

    for (sparse_hash_map<char, long long, hash<char>, compare_char>::iterator it=key_reg_hash.begin(); it!=key_reg_hash.end(); ++it)
    {
        if (it->first==read_char)
        {
            key_index_value = it->second;
            break;

        }

    }

    return key_index_value;

}


string rev_comp_str_convert(string &forward_index_str)
{
    string rev_comp_str;

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
        return 'R'; //purine

    } else if (read_char=='C' or read_char=='T')
    {
        return 'Y'; //pyrimidine

    } else
    {
        return read_char;

    }
}



void multi_hash_container_input(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash,
                                vector<string> &max_index_key_vector, string str_key)
{

    sort(max_index_key_vector.begin(), max_index_key_vector.end());

    bool str_key_input_flag=false;

    for (vector<string>::iterator it=max_index_key_vector.begin(); it!=max_index_key_vector.end(); ++it)
    {
        //if (str_key <= *it && prim_index_hash[*it].size() <= prim_index_hash[*it].max_size())
        if ((str_key == *it) || (str_key < *it && prim_index_hash[*it].size() < prim_index_hash[*it].max_size()))
        {
            prim_index_hash[*it][str_key]++;
            str_key_input_flag=true;
            //cout << str_key << "\t" << max_index_key_vector.size() << endl;
            break;

        } else if (str_key < *it) // && prim_index_hash[*it].size()==prim_index_hash[*it].max_size())
        {
            //cout << str_key << "\t" << "over" << endl;
            prim_index_hash[str_key][str_key]=1;
            max_index_key_vector.push_back(str_key); //add new represent key

            prim_index_hash[*it].set_deleted_key(string());

            for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=prim_index_hash[*it].begin(); sub_it!=prim_index_hash[*it].end(); ++sub_it)
            {
                if (sub_it->first <= str_key)
                {
                    prim_index_hash[str_key][sub_it->first]+=sub_it->second; //+= for self add
                    //prim_index_hash[*it].erase(sub_it->first); //erase by const key
                    prim_index_hash[*it].erase(sub_it); //erase by iterator position, also available

                }

            }
            ///compact hashtable *it to a smllest valid size
            prim_index_hash[*it].clear_deleted_key(); //off delete key call
            prim_index_hash[*it].resize(0); //not clear, but compact the hashtable because erase doesn't remove(element) and resize hashtable.
            //printf("done-resize\t%lld\n", prim_index_hash[*it].size());

            str_key_input_flag=true;
            break;

        }

    }


    if (str_key_input_flag==false) //sub_str > *it, add new
    {
        prim_index_hash[str_key][str_key]++; //<string, long long>
        max_index_key_vector.push_back(str_key);
    }


}




void text_compress_binary(sparse_hash_map<string, long long, hash<string>, compare_string> &compress_text_hash, string &str_block, stringstream &output_stream)
{
    bitset<8> char_replace=0;
    long long delimit_pos=0;
    string str_block_bk="";
    int str_compress_rate_length = int(str_block.length() * 0.25); //for speed up 1/4 size shrink as threshold

    for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator it=compress_text_hash.begin(); it!=compress_text_hash.end(); ++it)
    {
        str_block_bk="";

        //while (str_block!=str_block_bk)
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

    //write_f.write(str_block.c_str(), str_block.length());
    output_stream << str_block;

}


///first concate all item to single string and then print
void multi_hash_container_output(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash, vector<string> &max_index_key_vector,
                                 sparse_hash_map<string, long long, hash<string>, compare_string> &compress_text_hash, bool ratio_flag, int str_key_frag_size, long long &feature_hit_cnt, stringstream &output_stream)
{
    char stream_buf[70]; //give enough space
    int stream_size;
    string str_block;

    sort(max_index_key_vector.begin(), max_index_key_vector.end());

    vector<string> sub_index_vector;


    prim_index_hash.set_deleted_key(string());

    for (vector<string>::iterator it=max_index_key_vector.begin(); it!=max_index_key_vector.end(); ++it)
    {
        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=prim_index_hash[*it].begin(); sub_it!=prim_index_hash[*it].end(); ++sub_it)
        {
            sub_index_vector.push_back(sub_it->first);

        }

        sort(sub_index_vector.begin(), sub_index_vector.end());
        //prim_index_hash[*it].set_deleted_key(string());
        prim_index_hash[*it].set_deleted_key(string());

        for (vector<string>::iterator vector_it=sub_index_vector.begin(); vector_it!=sub_index_vector.end(); ++vector_it)
        {
            //2013.11.29, change delimiter tab to enter
            if (ratio_flag==true)
            {
                stream_size = snprintf(stream_buf, sizeof(stream_buf), "%s|%.4e\n", (*vector_it).c_str(), (double)(prim_index_hash[*it])[*vector_it] / (double)feature_hit_cnt);
                str_block = (string)stream_buf;

                text_compress_binary(compress_text_hash, str_block, output_stream);
                //output_stream << str_block;

            } else
            {
                stream_size = snprintf(stream_buf, sizeof(stream_buf), "%s|%lld\n", (*vector_it).c_str(), (prim_index_hash[*it])[*vector_it]);
                str_block = (string)stream_buf;

                text_compress_binary(compress_text_hash, str_block, output_stream);
                //output_stream << str_block;

            }

            memset(stream_buf, '\0', sizeof(stream_buf)); //release memory
            prim_index_hash[*it].erase(*vector_it); //erase element after output for decrease next hash search time
            prim_index_hash[*it].resize(0); //compact the hashtable because erase doesn't remove(element) and resize hashtable.

        }

        prim_index_hash[*it].clear_deleted_key();

        prim_index_hash.erase(*it); //by const key
        prim_index_hash.resize(0); //compact hash, save memory space, invalidate iterator

        sub_index_vector.clear();

    }

    prim_index_hash.clear_deleted_key();

}


///first concate all item to single string and then print
void compress_word_analyze(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash,
                            sparse_hash_map<string, long long, hash<string>, compare_string> &compress_text_hash, bool ratio_flag,
                            long long &feature_hit_cnt, int str_key_frag_size, stringstream &output_stream)
{
    char stream_buf[70]; //give enough space
    int stream_size;

    string str_key_frag; //fragment of str_key(as a whole)
    string str_value;
    long long min_frequency=0;

    for (sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>::iterator it=prim_index_hash.begin(); it!=prim_index_hash.end(); ++it)
    {
        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=(it->second).begin(); sub_it!=(it->second).end(); ++sub_it)
        {
            str_key_frag = (sub_it->first).substr(0, str_key_frag_size); //get less(or equal) than half segment

            if (ratio_flag==true)
            {
                stream_size=snprintf(stream_buf, sizeof(stream_buf), "%.4e", (double)(sub_it->second) / (double)feature_hit_cnt);
                str_value=(string)stream_buf;

            } else //convert longlong to string type
            {
                stream_size=snprintf(stream_buf, sizeof(stream_buf), "%lld", sub_it->second);
                str_value=(string)stream_buf;

            }


            if (compress_text_hash.size()>compress_text_hash.max_size()-3)
            {
                break;
            }


            compress_text_hash[str_key_frag]++; //should check size limit
            compress_text_hash[str_value]++; //should check size limit

            memset(stream_buf, '\0', sizeof(stream_buf)); //release memory

        }

    }

    //cout << compress_text_hash.size() << endl;

    ///sort and remove low frequency keys, maximum number of keys=128
    long long min_value=1;
    compress_text_hash.set_deleted_key(string()); //causing problem here

    while (compress_text_hash.size() > 128)
    {
        //find smallest
        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator it=compress_text_hash.begin(); it!=compress_text_hash.end(); it++)
        {
            if (it->second==min_value)
            {
                //cout << it->first << "\t" << it->second << endl;
                compress_text_hash.erase(it);
            }

            if (compress_text_hash.size()==128)
            {
                break;

            }

        }

        min_value++;

    }

    compress_text_hash.clear_deleted_key();
    //compress_text_hash.resize(0); //resize(0) invalidate hashtable

    //replace original value to (long long) index order and print?
    long long order_cnt=0;
    string str_compress_key;
    bitset<8> char_order_cnt=0;

    for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator it=compress_text_hash.begin(); it!=compress_text_hash.end(); it++)
    {
        it->second = order_cnt;
        char_order_cnt = bitset<8>(order_cnt);

        //stream_size=snprintf(stream_buf, 100, "#%s|%lld\n", (it->first).c_str(), order_cnt);
        char_order_cnt[7]=1;
        //stream_size=snprintf(stream_buf, 100, "#%s|%s\n", (it->first).c_str(), reinterpret_cast<char *>(&char_order_cnt));
        stream_size=snprintf(stream_buf, sizeof(stream_buf), "#%s|%lld\n", (it->first).c_str(), order_cnt);
        str_compress_key = (string)stream_buf;

        memset(stream_buf, '\0', sizeof(stream_buf)); //release memory

        //write_f.write(str_compress_key.c_str(), str_compress_key.length());
        output_stream << str_compress_key;

        order_cnt++;

    }


}



///2012.6.26, compact version, separate with output function
int read_sliding_process(stringstream &read_f, int feature_length, bool backward_flag, long long &feature_hit_cnt, bool unmasked_treat_flag, bool ratio_output_flag, bool ry_code_flag,
                            sparse_hash_map<char, long long, hash<char>, compare_char> &key_reg_hash, sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash,
                            vector<string> &max_index_key_vector)
{

    string forward_index_str="";
    string backward_index_str="";

    char read_char;

    //long long key_index_value=0;

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

            if (unmasked_treat_flag==true) //masked data provie masked region as lower case alphbet
            {
                read_char = toupper(read_char);

            }

            if (key_hash_check(read_char, key_reg_hash)!=-1) //char registered in key_reg_hash
            {

                if (ry_code_flag==true) //convert to ry code
                {
                    read_char = ry_convert(read_char);

                }

                forward_index_str+=read_char;

                if (forward_index_str.length() == feature_length)
                {

                    if (backward_flag==true)
                    {
                        backward_index_str = rev_comp_str_convert(forward_index_str);

                        if (forward_index_str > backward_index_str) //take smallest one
                        {
                            multi_hash_container_input(prim_index_hash, max_index_key_vector, backward_index_str);

                        } else //backward_index_str larger or equal
                        {
                            multi_hash_container_input(prim_index_hash, max_index_key_vector, forward_index_str);

                        }

                    } else //no backward consideration
                    {
                        multi_hash_container_input(prim_index_hash, max_index_key_vector, forward_index_str);
                    }

                    //2012.6.28, also summing +1
                    feature_hit_cnt+=1;

                    forward_index_str.erase(0, 1); //works as deque.pop()
                    //cout << forward_index_str << endl;
                }

            } else
            {
                forward_index_str.clear();
                backward_index_str.clear();

            }

        }

    }

    //refresh deque
    forward_index_str.clear();

    if (!prim_index_hash.empty())
    {
        return 1;

    }

    return 0; //in case prim_index_hash=empty

}



long long vocab_size_measure(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash, long long bottom_limit, long long top_limit)
{
    long long vocab_size=0;
    long long ret_value=0;

    for (sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>::iterator it=prim_index_hash.begin(); it!=prim_index_hash.end(); ++it)
    {
        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=(it->second).begin(); sub_it!=(it->second).end(); ++sub_it)
        {
            if ((sub_it->second >= bottom_limit || bottom_limit==0) && (sub_it->second <= top_limit || top_limit==0))
            {
                vocab_size++;

            }

        }
        //vocab_size+=(long long)(it->second).size();

    }

    return vocab_size;

}



void show_help()
{
    cout << "Parameter usage [option][load_path][save_path]" << endl;
    cout << "-h, show_help()" << endl;
    cout << "-s [INT], feature size(dependency on memory size and bit size)" << endl;

    cout << endl << "Extend option;" << endl;
    cout << "-a, input is amino acid sequence(default is AGCT nucleotide)" << endl;
    cout << "-c, convert nucleotide(AGCT) code to RY code" << endl;
    cout << "-k [STR], manual define of key_str, ex; 'HJKL' -> recognize as (H, J, K, L), set of character without repeat" << endl;
    cout << "-r, disable reverse complement counting(given -a will turn off reverse complement)" << endl;
    cout << "-n, output freqeuncy into ratio(of sum of all frequency)" << endl;
    cout << "-u, unmasked_treat, also accept lower letters(default; accept only upper letters)" << endl;

    cout << endl << "Vocabulary size measure;" << endl;
    cout << "-V, max vocab size measure" << endl;
    cout << "-b [long], bottom_limit(default = 2), minimum=0" << endl;
    cout << "-t [long], top_limit(default = 0 = maximum)" << endl;

}


void show_profile()
{
    cout << "2014.12.22, pfp_helen_rec\n";
    cout << "Code by JaeJin Choi\n";
    cout << "Support longer feature length using string key hash\n";
    cout << "Out-format (STR key | float(or long long) value) TAB\n";
    cout << "based on pfp_helen_rec, add RY code function\n";
    cout << "2014.12.22; Load whole file to memory(istringstream)\n";
    cout << "2015; add manual compression\n";
    cout << "2015.8.31; add zlib compression\n";
    cout << "compile; g++ -std=c++11 -o (output) (this script) -lz\n";
}



int main(int argc, char** argv)
{
    int feature_length=10; //default

    bool ratio_output_flag=false;
    bool part_output_flag=false;
    bool backward_flag=true;
    bool max_vocab_find_flag=false; //default;definition; count number of feature that (1<frequency)
    bool ry_code_flag=false;


    long long top_limit=0;
    long long bottom_limit=2; //default ==2

    //have trigger of upper case(using masked, unmasked exception)
    bool unmasked_treat_flag=false; //as default program don't accept lower case letter

    string key_str="AGCT"; //default
    int key_str_size=0;

    //parameter with option flag
    int opt;

    while ((opt = getopt(argc, argv, "hvs:ack:rnub:t:V")) !=EOF) // EOF = -1
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

            case 'a':
                key_str = "GASTCVLIMPFYWDENQHKR";
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
                bottom_limit=atoll(optarg);
                break;

            case 't':
                top_limit=atoll(optarg);
                break;

            case 'V':
                max_vocab_find_flag=true;
                break;

            default:
                show_help();
                exit(0);

        }
    }


    key_str_size = int(key_str.size());

    ///only nucleotide have reverse_complement option-double strand
    if (key_str!="AGCT") //including amino acides but not RY code
    {
        //str_key_frag_size = min(feature_length, max(2, int(feature_length * 0.2)+1));
        backward_flag=false;

    } else if (ry_code_flag==true)
    {
        key_str_size = 2; //RY code but can recognize 'ACGT' when slide reading

    }


    if (bottom_limit > top_limit && top_limit!=0)
    {
        cout << "Wrong arrangment of limit" << endl;
        cout << "bottom_limit <= top_limit or top_limit==0 or bottom_limit==0(minimum 0)" << endl;
        exit(1);

    }


    ///instead give specific fragment size per certain cases
    ///what is the optimal str_key_frag_size for a compression(regard speed and efficiency)?
    //int str_key_frag_size = min(feature_length, max(2, int(feature_length * 0.2)+1)); //30% of feature_length
    //int str_key_frag_size = 3; //as a default for AGCT coding
    int str_key_frag_size = min(feature_length, max(2, int(log(feature_length) / log(key_str_size) )));


    ///2015-5, optimize memory use; define hash contain size
    ///register key_str
    sparse_hash_map<char, long long, hash<char>, compare_char> key_reg_hash;
    key_reg_hash.resize(key_str.length());
    register_key_hash(key_str, key_reg_hash);

    ///hash table max size, let depend on the container size
    //long long sub_container_size_limit = (long long)pow(2, 23); //fixed/initial size size

    sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string> prim_index_hash;

    vector<string> max_index_key_vector;

    ///compression index
    sparse_hash_map<string, long long, hash<string>, compare_string> compress_text_hash;

    //cout << str_key_frag_size << endl;

    ifstream read_f;
    ofstream write_f; //profile save

    stringstream read_str_f(ios::in|ios::out);
    stringstream output_stream(ios::out|ios::in|ios::ate);

    string output_total_string="";

    long long feature_hit_cnt=0;
    long long past_vocab_size=0;
    long long recent_vocab_size=0;
    int original_feature_length = feature_length;
    int cnt=0;

    int ret_stream_size=0;
    char stream_buf[70];

    //while (argc - optind !=0)
    if (argc  - optind == 2) //direct file write version; load_path and save_path
    {
        read_f.open(argv[optind], ios::in); //read text file
        write_f.open(argv[optind+1], ios::out|ios::binary); //save binary file

        //while (!read_f.eof())
        if (!read_f)
        {
            cout << "Unable to open file: " << argv[optind] << endl;
            return 1;

        } else
        {
            read_str_f << read_f.rdbuf();
            read_str_f.seekg(0, ios::beg);
            read_f.close();

            while (!read_str_f.eof())
            {
                if (read_sliding_process(read_str_f, feature_length, backward_flag, feature_hit_cnt, unmasked_treat_flag, ratio_output_flag, ry_code_flag,
                    key_reg_hash, prim_index_hash, max_index_key_vector)==1)
                {
                    if (max_vocab_find_flag==true)
                    {
                        recent_vocab_size = vocab_size_measure(prim_index_hash, bottom_limit, top_limit);

                        if (past_vocab_size==0 || (past_vocab_size <= recent_vocab_size && recent_vocab_size > 0)) //not support part_output
                        {
                            //printf("l-mer; %d, vocab_size; %lld\n", feature_length, recent_vocab_size); //for vocabulary size trend notify
                            ret_stream_size = snprintf(stream_buf, sizeof(stream_buf), "l-mer; %d, vocab_size; %lld\n", feature_length, recent_vocab_size);
                            //ret_stream_size = snprintf(stream_buf, sizeof(stream_buf), "l-mer; %d, vocab_size; %lld, #amino; %lld\n", feature_length, recent_vocab_size, feature_hit_cnt + feature_length -1);
                            output_stream << stream_buf;
                            memset(stream_buf, '\0', sizeof(stream_buf));

                            past_vocab_size=recent_vocab_size;

                            ///clear and than seekg(clear procedure should be first)
                            read_str_f.clear(); //reset flag
                            read_str_f.seekg(0, ios_base::beg);

                            feature_length++;

                        } else
                        {

                            //printf("l-mer; %d, vocab_size; %lld\n", feature_length, recent_vocab_size); //for vocabulary size trend notify
                            ret_stream_size = snprintf(stream_buf, sizeof(stream_buf), "l-mer; %d, vocab_size; %lld\n", feature_length, recent_vocab_size);
                            //ret_stream_size = snprintf(stream_buf, sizeof(stream_buf), "l-mer; %d, vocab_size; %lld, #amino; %lld\n", feature_length, recent_vocab_size, feature_hit_cnt + feature_length -1);
                            output_stream << stream_buf;
                            memset(stream_buf, '\0', sizeof(stream_buf));

                            ret_stream_size = snprintf(stream_buf, sizeof(stream_buf),"#Point of MAX vocabulary size(%lld <= frequency <= %lld); %d\n", bottom_limit, top_limit, feature_length-1);
                            output_stream << stream_buf;
                            memset(stream_buf, '\0', sizeof(stream_buf));

                            max_index_key_vector.clear();
                            prim_index_hash.clear();
                            feature_hit_cnt=0;

                            write_f << output_stream.str(); //<< '\n'; //additional linebreak un necessary
                            write_f.close();

                            output_stream.str(string());
                            output_stream.clear();

                            break;

                        }

                    } else
                    {
                        ///release memory space
                        read_str_f.str(string());
                        read_str_f.clear();

                        //create compression index hash -> compress_text_hash
                        compress_word_analyze(prim_index_hash, compress_text_hash, ratio_output_flag, feature_hit_cnt, str_key_frag_size, output_stream);
                        //printf("input finish\n"); ///input check
                        multi_hash_container_output(prim_index_hash, max_index_key_vector, compress_text_hash, ratio_output_flag, str_key_frag_size, feature_hit_cnt, output_stream);

                        //ofstream write_f(argv[optind+1]);
                        write_f << compress_deflate(output_stream.str(), Z_BEST_COMPRESSION); //<< '\n'; //additional linebreak un necessary
                        write_f.close();

                        ///std output
                        //cout << compress_deflate(output_stream.str(), Z_BEST_COMPRESSION) << endl;
                        //cout << output_stream.str() << endl; //output without compress

                        output_stream.str(string());
                        output_stream.clear();
                    }

                    max_index_key_vector.clear();
                    prim_index_hash.clear();
                    compress_text_hash.clear();
                    feature_hit_cnt=0;

                }

            }

            read_str_f.str(string());
            read_str_f.clear();

        }

        //optind++; //increase optind

    } else
    {
        cout << "Please input load_path, and then save_path" << endl;
        //show_help();
        exit(EXIT_SUCCESS);

    }

}


