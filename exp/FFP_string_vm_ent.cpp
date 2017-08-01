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

//#include <deque>
//using std::deque;

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
#include <climits>


///2013.3.13, no calling external hash
using std::hash;

/* or
using ext::hash;reverse
using __gnu_cxx::hash;
*/


// Found these here http://mail-archives.apache.org/mod_mbox/trafficserver-dev/201110.mbox/%3CCACJPjhYf=+br1W39vyazP=ix
//eQZ-4Gh9-U6TtiEdReG3S4ZZng@mail.gmail.com%3E
#define MOD_GZIP_ZLIB_WINDOWSIZE 15
#define MOD_GZIP_ZLIB_CFACTOR    9reverse
#define MOD_GZIP_ZLIB_BSIZE      8096


// Found this one here: http://panthema.net/2007/0328-ZLibString.html, author is Timo Bingmann
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

/*
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
}; //should put ; after struct declare



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
        return 'R'; //purine

    } else if (read_char=='C' or read_char=='T')
    {
        return 'Y'; //pyrimidine

    } else
    {
        return read_char;

    }
}



string integer_to_bit_string(unsigned long &bits_per_alphabet, int int_size)
{
    string bit_string="";

    for (int cy1=bits_per_alphabet-1; cy1!=-1; cy1--)
    {
        //cout << cy1 << "\t" << pow(2, cy1) << endl;
        if (int_size >= pow(2, cy1))
        {
            int_size-=pow(2, cy1);
            bit_string+='1';

        } else
        {
            bit_string+='0';
        }


    }

    return bit_string;
}



void key_hash_register(string &key_str, sparse_hash_map<char, string, hash<char>, compare_char> &key_reg_hash,
    unsigned long &bits_per_alphabet, int additional_bit)
{
    int alphabet_cnt=additional_bit; //define base (eg 0 or 1)

    //sort(key_str.begin(), key_str.end()); //sort given key_str into lexical order
    //key_reg_hash.resize(key_str.length());

    for (string::iterator it=key_str.begin(); it!=key_str.end(); it++)
    {
        key_reg_hash[*it] = integer_to_bit_string(bits_per_alphabet, alphabet_cnt);
        alphabet_cnt++;

    }

}



string str_feature_to_compact_feature_str(string str_feature, unsigned long &bits_per_feature,
    sparse_hash_map<char, string, hash<char>, compare_char> &key_reg_hash)
{
    string key_bit_string="";

    for (string::iterator it=str_feature.begin(); it!=str_feature.end(); ++it)
    {
        key_bit_string+=key_reg_hash[*it];

    }

    string converted_key_string="";
    bitset<8> char_bit;

    ///no longer necessary
    /*
    ///2017-1 add one additional byte to prevent using deleted hash key (empty string())
    char_bit = bitset<8>("11111111"); //is not empty, to avoid deleted hash key ("")
    //converted_key_string+=reinterpret_cast<char *>(&char_bit);
    converted_key_string+=static_cast<unsigned char>(char_bit.to_ulong());
    */

    key_bit_string.resize(bits_per_feature, '1'); ///1 or 0?, use 1 to make consistent spacing for JSD
    //cout << key_bit_string.size() << "\t" <<  key_bit_string << endl;

    for (size_t it=0; it!=key_bit_string.length(); it+=8)
    {
        char_bit = bitset<8>(key_bit_string.substr(it, 8));
        //converted_key_string+=reinterpret_cast<char *>(&char_bit); //reinterpret_cast vs static_cast
        //converted_key_string+=static_cast<unsigned char>(char_bit.to_ulong()); //reinterpret_cast vs static_cast
        converted_key_string+=char(char_bit.to_ulong());
        //cout << it << "\t" << it+8 << "\t" << key_bit_string.substr(it, 8) << endl;

    }

    //cout << converted_key_string.size() << "\t" << converted_key_string << endl;
    //exit(0);

    return converted_key_string;

}



double base_string_entropy(string &key_str, string &str_key)
{
    vector<int> str_count_vector;
    str_count_vector.resize(key_str.size());

    for (string::iterator it=str_key.begin(); it!=str_key.end(); ++it)
    {
        str_count_vector[key_str.find(*it)]++;
        //cout << *it << ":" << key_str.find(*it) << "\t";

    }

    double str_entropy=0.0;
    double each_frequency=0.0;
    int key_str_length = key_str.size();
    int str_key_length = str_key.size();

    for (vector<int>::iterator it=str_count_vector.begin(); it!=str_count_vector.end(); ++it)
    {
        //each_frequency=(double)(*it) / key_str_length;str_entropy-=each_frequency * (log(each_frequency) / log(key_str_length)); ///log base is the size of key_str. possible alphabets
        if (*it!=0)
        {
            each_frequency=(double)(*it) / str_key_length;
            str_entropy-=each_frequency * (log(each_frequency) / log(key_str_length)); ///log base is the size of key_str. possible alphabets
            //str_entropy-=each_frequency * (log(each_frequency) / log(2)); ///log base is the size of key_str. possible alphabets
            //str_entropy-=each_frequency * (log(each_frequency) / log(str_key_length)); ///log base is the size of key_str. possible alphabets

        } else
        {
            each_frequency=0;
        }

        //cout << *it << ":" << each_frequency << "\t";

    }
    //cout << endl;
    //cout << key_str << "\t" << str_entropy << endl;

    return str_entropy;

}



void multi_hash_container_input(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string> &prim_index_hash,
    vector<string> &max_index_key_vector, string str_key, string &deleted_hash_key)
{
    bool key_insert_flag=false;

    for (vector<string>::iterator it=max_index_key_vector.begin(); it!=max_index_key_vector.end(); ++it)
    {
        if ((str_key == *it) || (str_key < *it && prim_index_hash[*it].size() < INT_MAX))
        {
            prim_index_hash[*it][str_key]++;
            //cout << str_key << "\t" << max_index_key_vector.size() << endl;
            key_insert_flag=true;
            break;

        } else if (str_key < *it) // && prim_index_hash[*it].size()==prim_index_hash[*it].max_size())
        {
            prim_index_hash[str_key][str_key]=1;
            prim_index_hash[*it].set_deleted_key(deleted_hash_key); //or string()

            //for (sparse_hash_map< vector<bool>, long long, hash< <vector<bool> >, compare_bool_vector>::iterator sub_it=prim_index_hash[*it].begin(); sub_it!=prim_index_hash[*it].end(); ++sub_it)
            for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=prim_index_hash[*it].begin(); sub_it!=prim_index_hash[*it].end(); ++sub_it)
            {
                if (sub_it->first <= str_key)
                {
                    prim_index_hash[str_key][sub_it->first]+=sub_it->second; //+= for self add
                    //prim_index_hash[*it].erase(sub_it->first); //erase by key
                    prim_index_hash[*it].erase(sub_it); //erase by iterator position, also available and not invalidate iteration

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


    if (key_insert_flag==false)
    {
        prim_index_hash[str_key][str_key]=1; //<string, long long>
        max_index_key_vector.push_back(str_key); //in vector, push at back does not invalidate iterator

    }

    sort(max_index_key_vector.begin(), max_index_key_vector.end()); //invalidate iterator

}


void multi_hash_container_output(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string> &prim_index_hash,
    vector<string> &max_index_key_vector, unsigned long &bits_per_feature, double &feature_hit_cnt, bool ratio_flag, stringstream &output_stream)
{
    ///first, header
    ///first byte define the length of feature(as bits)

    unsigned long bytes_per_feature = bits_per_feature / 8; //convert unit bits to bytes
    unsigned long bytes_per_value = 0;

    //record device defined variable size, at the point of write
    if (ratio_flag==true) //double
    {
        bytes_per_value = sizeof(double);
        ///double = 8 bytes (xe-308 ~ xe+308, 15 digit precision), float = 4 bytes (xe-38 ~ xe+38, 6 digit precision)

    } else //long long
    {
        bytes_per_value = sizeof(long long);

    }

    if (!prim_index_hash.empty())
    {
        output_stream << static_cast<unsigned char>(bytes_per_feature); //feature bytes size, 1 bytes
        output_stream << static_cast<unsigned char>(bytes_per_value); //value bytes size, 1 bytes

    }

    ///second, actual data, focused on each soild feature (of given length)
    ///in formatted blocks of [feature (string)][value (long long or double)]

    double value_is_double=0.0;
    long long value_is_lld=0;

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
            ///put feature as string
            output_stream << (*vector_it);
            //output_stream.write((*vector_it).c_str(), (*vector_it).size());
            //cout << *vector_it << "\t" << (*vector_it).length() << "\t";

            ///no delimiter but by block size
            ///put value
            if (ratio_flag==true)
            {
                value_is_double = (double)(prim_index_hash[*it])[*vector_it] / feature_hit_cnt;
                output_stream.write(reinterpret_cast<char*>(&value_is_double), sizeof(double));
                //output_stream.write(reinterpret_cast<char*>(&value_is_double), sizeof(value_is_double));
                //cout << value_is_double << endl;
                //cout << (*vector_it) << "\t" << value_is_double << endl;

            } else
            {
                value_is_lld = (prim_index_hash[*it])[*vector_it];
                output_stream.write(reinterpret_cast<char*>(&value_is_lld), sizeof(long long));
                //output_stream.write(reinterpret_cast<char*>(&value_is_lld), sizeof(value_is_lld));
                //cout << (*vector_it) << "\t" << value_is_lld << endl;

            }

            prim_index_hash[*it].erase(*vector_it); //erase element after output for decrease next hash search time

        }

        prim_index_hash[*it].clear_deleted_key();
        prim_index_hash[*it].resize(0); //compact the hashtable because erase doesn't remove(element) nor reduce hashtable.

        prim_index_hash.erase(*it); //by const key
        sub_index_vector.clear();

    }

    prim_index_hash.clear_deleted_key();
    prim_index_hash.resize(0); //compact the hashtable because erase doesn't remove(element) nor reduce hashtable.


    /*
    string tmp_compressed_str = compress_deflate(output_stream.str(), Z_BEST_COMPRESSION); //<< '\n'; //additional linebreak unnecessary
    output_stream.str(string());
    output_stream.clear();

    output_stream << decompress_deflate(tmp_compressed_str);
    tmp_compressed_str.clear();


    ///reverse check (file read)
    output_stream.seekg(ios::beg);
    output_stream.clear();

    unsigned long key_length = output_stream.get() + 0;
    unsigned long value_length = output_stream.get() + 0;

    //output_stream.read(reinterpret_cast<char*>(&key_feature), key_length);
    //output_stream.read(reinterpret_cast<double>(&key_value), value_length);
    //char key_char[key_length];
    char key_buf[256]; //to add null at the end
    string key_str;

    double key_value;

    cout << "----------------------" << endl;

    while (!output_stream.eof())
    {
        while (output_stream.read(key_buf, key_length))
        //while (!output_stream.eof())
        {
            //output_stream.read(key_buf, key_length);
            key_str = string(key_buf);
            output_stream.read(reinterpret_cast<char*>(&key_value), value_length);

            //cout << string(key_buf) << "\t" << key_value << endl; //"\t";
            cout << key_str << "\t" << key_value << endl;
        }


    }
    //cout << sizeof(key_buf) << "\t" << size_t(key_buf) << endl;
    */


}


///2012.6.26, compact version, separate with output function
int read_sliding_process(stringstream &read_f, int feature_length, bool backward_flag, double &feature_hit_cnt, bool unmasked_treat_flag, bool ratio_output_flag, bool ry_code_flag,
    sparse_hash_map<char, string, hash<char>, compare_char> &key_reg_hash, sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash,
    vector<string> &max_index_key_vector, unsigned long &bits_per_feature, string &key_str, string &deleted_hash_key, double bottom_limit, double top_limit)
{

    string forward_index_str="";
    string backward_index_str="";
    string key_index_str=""; //final hash key product of convertion

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

            if (unmasked_treat_flag==true) //masked data provie masked region as lower case alphbet
            {
                read_char = toupper(read_char);

            }

            if (ry_code_flag==true) //convert to ry code
            {
                read_char = ry_convert(read_char);

            }

            //if (key_hash_check(read_char, key_reg_hash)==true) //char registered in key_reg_hash
            if(key_reg_hash.find(read_char)!=key_reg_hash.end())
            {
                forward_index_str+=read_char;

                if (forward_index_str.length() == feature_length)
                {
                    str_entropy = base_string_entropy(key_str, forward_index_str);
                    //if (base_string_entropy(key_str, forward_index_str) >= bottom_limit && base_string_entropy(key_str, forward_index_str) <= top_limit)
                    if (str_entropy >= bottom_limit && str_entropy <= top_limit)
                    {
                        if (backward_flag==true)
                        {
                            backward_index_str = rev_comp_str_convert(forward_index_str);

                            if (forward_index_str > backward_index_str) //take smallest one
                            {
                                key_index_str = str_feature_to_compact_feature_str(backward_index_str, bits_per_feature, key_reg_hash);
                                //multi_hash_container_input(prim_index_hash, max_index_key_vector, strdup(backward_index_str.c_str()));
                                //multi_hash_container_input(prim_index_hash, max_index_key_vector, backward_index_str.data());

                            } else //backward_index_str larger or equal
                            {
                                key_index_str = str_feature_to_compact_feature_str(forward_index_str, bits_per_feature, key_reg_hash);
                                //multi_hash_container_input(prim_index_hash, max_index_key_vector, strdup(forward_index_str.c_str()));
                                //multi_hash_container_input(prim_index_hash, max_index_key_vector, forward_index_str.data());

                            }

                        } else //no backward consideration
                        {
                            key_index_str = str_feature_to_compact_feature_str(forward_index_str, bits_per_feature, key_reg_hash);

                        }

                        //multi_hash_container_input(prim_index_hash, max_index_key_vector, strdup(key_index_str.c_str()), deleted_hash_key);
                        multi_hash_container_input(prim_index_hash, max_index_key_vector, key_index_str, deleted_hash_key);

                        //multi_hash_container_input(prim_index_hash, max_index_key_vector, strdup(forward_index_str.c_str()));
                        //multi_hash_container_input(prim_index_hash, max_index_key_vector, key_index_str);
                        feature_hit_cnt+=1;

                    } /*else if (base_string_entropy(key_str, forward_index_str) >= top_limit)
                    {
                        cout << "exc; " << forward_index_str << "\t" << base_string_entropy(key_str, forward_index_str) << endl;


                    }
                    */
                    //2012.6.28, also summing +1
                    //feature_hit_cnt+=1;

                    forward_index_str.erase(0, 1); //works as deque.pop()
                    //cout << forward_index_str << endl;
                }


            } else
            {
                forward_index_str.clear();
                backward_index_str.clear();
                //key_index_str.clear();

            }

        }

    }
    forward_index_str.clear();
    //key_index_str.clear();

    if (!prim_index_hash.empty())
    {
        return 1;

    }

    return 0; //in case prim_index_hash=empty

}


/* ///filter by frequency is unfair and is arbitary. Frequency greatly relies on genome size
void filter_frequency(sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>  &prim_index_hash,
    long long bottom_limit, long long top_limit, double &feature_hit_cnt, string &deleted_hash_key)
{
    for (sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>::iterator it=prim_index_hash.begin(); it!=prim_index_hash.end(); ++it)
    {
        (it->second).set_deleted_key(deleted_hash_key); //or string()

        for (sparse_hash_map<string, long long, hash<string>, compare_string>::iterator sub_it=(it->second).begin(); sub_it!=(it->second).end(); ++sub_it)
        {
            //if (sub_it->second < bottom_limit || sub_it->second > top_limit)
            if ((sub_it->second < bottom_limit && bottom_limit!=1) || (sub_it->second > top_limit && top_limit!=0))
            {
                feature_hit_cnt-=(sub_it->second); //substract selected feature count
                (it->second).erase(sub_it); //erase by iterator position?, also available and not invalidate iteration
                //(it->second).erase(sub_it->first); //erase by key

            }

        }

        (it->second).clear_deleted_key(); //off delete key call
        (it->second).resize(0); //compact the hashtable because erase doesn't remove(element) nor reduce hashtable.

    }
    //prim_index_hash.resize(0); //necessary?

}
*/



long long vocab_size_measure(sparse_hash_map<string, sparse_hash_map<string,  long long, hash<string>, compare_string>, hash<string>, compare_string> &prim_index_hash)
{
    long long vocab_size=0; //long long or double?

    for (sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string>::iterator it=prim_index_hash.begin(); it!=prim_index_hash.end(); ++it)
    {
        vocab_size+=(long long)(it->second).size();

    }

    return vocab_size;

}


///testing function
/*
void test_func()
{
    ///test swap
    sparse_hash_map<char, int, hash<char>, compare_char> temp_hash_a;
    //temp_hash_a = {{'a',1}, {'b',1}, {'c',1}, {'d',1}};
    temp_hash_a['a']=1;
    temp_hash_a['b']=2;
    temp_hash_a['c']=3;
    temp_hash_a['d']=4;

    sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string> prim_index_hash;
    cout << prim_index_hash.max_size() << endl;

    deque<string> index_deque;
    cout << index_deque.max_size() << endl;

    sparse_hash_map<string, sparse_hash_map< vector<bool>, long long, hash< vector<bool> >, compare_bool_vector>, hash< vector<bool> >, compare_bool_vector> prim_index_hash_vb;
    cout << prim_index_hash_vb.max_size() << endl;

    deque< vector<bool> > index_deque_vb;
    cout << index_deque_vb.max_size() << endl;

    sparse_hash_map<char, int, hash<char>, compare_char> temp_hash_b;
    temp_hash_b.resize(4);

    cout << "Before" << "\t" << temp_hash_a.size() <<"\t"<< temp_hash_b.size() << endl;

    //swap(temp_hash_a, temp_hash_b);
    //swap_ranges(temp_hash_a.begin(), temp_hash_a.begin(), temp_hash_b.begin());

    move(temp_hash_a.begin(), temp_hash_a.begin(), inserter(temp_hash_b, temp_hash_b.begin()));


    cout << "After" << "\t" << temp_hash_a.size() <<"\t"<< temp_hash_b.size() << endl;


    vector<bool> vec_bol;
    cout << vec_bol.max_size() << endl;

    vector<int> vec_int;
    cout << vec_int.max_size() << endl;

    vector<char> vec_char;

    cout << vec_char.max_size() << endl;

    exit(0);


}
*/

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

    while (getline(read_f, read_buf, '>')) //sizeof(read_buf)))
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
    cout << "-s [INT], feature size(dependency on memory size and bit size)" << endl;
    cout << "-e [INT], feature size end, limit feature_length++ while measuring vocabular_size" << endl;


    cout << endl << "Extend option;" << endl;
    cout << "-a, input is amino acid sequence(default is AGCT nucleotide)" << endl;
    cout << "-c, convert nucleotide (AGCT) code to RY code" << endl;
    cout << "-k [STR], manual define of key bases, ex; 'HJKL' -> recognize as (H, J, K, L), set of character without duplicates" << endl;
    cout << "-r, disable reverse complement counting(given -a will turn off reverse complement)" << endl;
    cout << "-n, output freqeuncy into ratio(of sum of all frequency)" << endl;
    cout << "-u, unmasked_treat, also accept lower letters(default; accept only upper letters)" << endl;

    cout << endl << "Vocabulary size measure;" << endl;
    cout << "-V, max vocab size measure" << endl;
    //cout << "-b [long], bottom_limit(default = 2), minimum=0" << endl;
    //cout << "-t [long], top_limit(default = 0 = maximum)" << endl;

    cout << endl << "string entropy limit;" << endl;
    cout << "-b [double], bottom_limit(default = 0.0), minimum=0" << endl;
    cout << "-t [double], top_limit(default = 1.0), maximum=1" << endl;

}


void show_profile()
{
    cout << "Code by JaeJin Choi\n";
    cout << "compile; g++ -std=c++11 -o (output) (this script) -lz\n";
    cout << "Google sparsehash\n";
    cout << "zlib 1.2.8+";

}



int main(int argc, char** argv)
{
    int feature_length=10; //default
    int feature_length_end=0;

    bool ratio_output_flag=false;
    bool backward_flag=true;
    bool max_vocab_find_flag=false; //default;definition; count number of feature that (1<frequency)
    bool ry_code_flag=false;


    double top_limit=1.0; //default==1.0==maximum
    double bottom_limit=0.0; //default==0.0=minimum

    //have trigger of upper case(using masked, unmasked exception)
    bool unmasked_treat_flag=false; //as default program don't accept lower case letter

    string key_str="AGCT"; //default is nucleotide

    //parameter with option flag
    int opt;

    while ((opt = getopt(argc, argv, "hvs:e:ack:rnub:t:V")) !=EOF) // EOF = -1
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
                bottom_limit=atof(optarg);
                break;

            case 't':
                top_limit=atof(optarg);
                break;

            case 'V':
                max_vocab_find_flag=true;
                break;

            default:
                show_help();
                exit(0);

        }
    }


    int additional_bit=1; //avoid 00000000 (empty char) so start from 1

    ///only nucleotide have reverse_complement option-double strand
    if (key_str!="AGCT") //including amino acides but not RY code
    {
        //str_key_frag_size = min(feature_length, max(2, int(feature_length * 0.2)+1));
        backward_flag=false;

    } else if (ry_code_flag==true)
    {
        key_str="RY";

    }


    if (bottom_limit > top_limit && top_limit!=0)
    {
        cout << "Wrong arrangment of limit\n";
        cout << "bottom_limit <= top_limit or top_limit==0(maximum) or bottom_limit==0(minimum 0)\n";
        exit(1);

    }


    ///instead give specific fragment size per certain cases
    ///what is the optimal str_key_frag_size for a compression(regard speed and efficiency)?
    //int str_key_frag_size = min(feature_length, max(2, int(feature_length * 0.2)+1)); //30% of feature_length
    //int str_key_frag_size = 3; //as a default for AGCT coding
    //int str_key_frag_size = min(feature_length, max(2, int(log(feature_length) / log(key_str_size) )));


    ///2015-5, optimize memory use; define hash contain size
    ///register key_str
    //unsigned long bits_per_alphabet = (unsigned long)(ceil((double)key_str.length()+1 / 2)); //key_str.length() + 1, add 1 for empty string '0000~' for deleted hash key reserved
    unsigned long bits_per_alphabet = count_bits_per_alphabet(key_str.length() + additional_bit);
    unsigned long bits_per_feature = (unsigned long)(ceil(((double)(bits_per_alphabet*feature_length) / 8)) * 8);

    sparse_hash_map<char, string, hash<char>, compare_char> key_reg_hash;
    key_hash_register(key_str, key_reg_hash, bits_per_alphabet, additional_bit);

    sparse_hash_map<string, sparse_hash_map<string, long long, hash<string>, compare_string>, hash<string>, compare_string> prim_index_hash;

    //string deleted_hash_key="ÿ"; //is "11111111" == char(255)
    string deleted_hash_key=""; //is empty string

    vector<string> max_index_key_vector;

    ifstream read_f;
    vector<string> str_block_vector; ///store string blocks of a file

    ofstream write_f; //profile save


    stringstream read_str_f(ios::in|ios::out);
    //stringstream output_stream(ios::out|ios::in|ios::binary|ios::ate);
    stringstream output_stream(ios::out|ios::in|ios::ate);

    string output_total_string="";

    double feature_hit_cnt=0.0;
    long long past_vocab_size=0;
    long long recent_vocab_size=0;

    //char stream_buf[70]; //replaced by  char* malloc
    char* stream_buf = NULL;
    size_t stream_size_t=0;

    //while (argc - optind !=0)
    if (argc  - optind == 2) //direct file write version; load_path and save_path
    {
        read_f.open(argv[optind], ios::in); //read text file
        write_f.open(argv[optind+1], ios::out|ios::binary); //save binary file

        if (!read_f)
        {
            cout << "Unable to open file: " << argv[optind] << endl;
            return 1;

        } else
        {

            file_to_string_blocks(str_block_vector, read_f);
            read_f.close();

            while (1)
            {

                for (vector<string>::iterator it=str_block_vector.begin(); it!=str_block_vector.end(); ++it)
                {
                    read_str_f.rdbuf()->pubsetbuf(&(*it)[0], (*it).size());

                    read_str_f.clear(); //clear flag
                    read_str_f.seekg(0, ios_base::beg);
                    //cout << (*it).size() << endl;
                    read_sliding_process(read_str_f, feature_length, backward_flag, feature_hit_cnt, unmasked_treat_flag,
                                        ratio_output_flag, ry_code_flag, key_reg_hash, prim_index_hash, max_index_key_vector,
                                        bits_per_feature, key_str, deleted_hash_key, bottom_limit, top_limit);

                }


                /*
                if (bottom_limit > 1 || top_limit!=0)
                {
                    filter_frequency(prim_index_hash, bottom_limit, top_limit, feature_hit_cnt, deleted_hash_key);

                }
                */


                if (max_vocab_find_flag==true && feature_length <= feature_length_end && feature_length_end!=0)
                {
                    recent_vocab_size = vocab_size_measure(prim_index_hash);

                    stream_size_t = snprintf(NULL, 0, "l-mer; %d, vocab_size; %lld\n", feature_length, recent_vocab_size); //find string(c) size
                    ///formated size larger than given stream_buff?
                    stream_buf = (char*)realloc(stream_buf, (stream_size_t + 1)*sizeof(char *)); //+1 at the end(/0)
                    snprintf(stream_buf, size_t(stream_buf), "l-mer; %d, vocab_size; %lld\n", feature_length, recent_vocab_size); //format char*

                    output_stream << stream_buf;

                    if (feature_length < feature_length_end)
                    {
                        max_index_key_vector.clear();
                        prim_index_hash.clear(); // == hash.erase(hash.begin(), hash.end()). if so, actually elements remain and take memory space
                        prim_index_hash.resize(0); //free memory
                        feature_hit_cnt=0;

                        feature_length++;
                        bits_per_feature = (unsigned long)(ceil(((double)(bits_per_alphabet*feature_length) / 8)) * 8); //recalculate


                    } else
                    {
                        write_f << output_stream.str(); //<< '\n'; //additional linebreak un necessary
                        break; //break while()

                    }

                } else if (max_vocab_find_flag==false) ///write FFP to file
                {
                    ///no output for now 2016-12-19
                    //printf("Bit-wise output not supported yet, 2016-12-19\n");
                    //printf("index_deque-size: %lld", (long long)max_index_key_vector.size());
                    multi_hash_container_output(prim_index_hash, max_index_key_vector, bits_per_feature, feature_hit_cnt, ratio_output_flag, output_stream);

                    write_f << compress_deflate(output_stream.str(), Z_BEST_COMPRESSION); //<< '\n'; //additional linebreak unnecessary
                    //write_f << output_stream.str(); //<< '\n'; //additional linebreak unnecessary

                    break; //break while()

                }


            }


            max_index_key_vector.clear();
            prim_index_hash.clear(); // == hash.erase(hash.begin(), hash.end()). if so, actually elements remain and take memory space
            prim_index_hash.resize(0); //free memory
            feature_hit_cnt=0;

            output_stream.str(string());
            output_stream.clear();

            str_block_vector.clear();

        }

        write_f.close();
        //optind++; //increase optind

    } else
    {
        cout << "Please input load_path, and then save_path" << endl;
        //show_help();
        exit(EXIT_SUCCESS);

    }

    free(stream_buf);
    stream_buf=NULL;
}



