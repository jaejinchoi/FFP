#include <iostream>
using namespace std;

#include <cstdlib>

#include <string>
using std::string;

#include <cstring>
//to use memset

#include <sstream>

#include <getopt.h>

#include <fstream>

#include <cmath>

#include <vector>
using std::vector;

#include <future> //retrieve return value from thread(template), require environment option check, -std=c++11 -lpthread
using std::future;
using std::future_status;
//using std::future_errc;
//using std::future_error;
#include <chrono>
//using std::chrono::milliseconds;

#include <algorithm>
using std::sort;

#include <bitset>


#include <zlib.h>
#include <stdexcept>

///Found zlib function here http://mail-archives.apache.org/mod_mbox/trafficserver-dev/201110.mbox/%3CCACJPjhYf=+br1W39vyazP=ix
///eQZ-4Gh9-U6TtiEdReG3S4ZZng@mail.gmail.com%3E
#define MOD_GZIP_ZLIB_WINDOWSIZE 15
#define MOD_GZIP_ZLIB_CFACTOR    9
#define MOD_GZIP_ZLIB_BSIZE      8096

std::string decompress_deflate(const std::string& str, int &zs_ret)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (inflateInit(&zs) != Z_OK)
        throw(std::runtime_error("inflateInit failed while decompressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();

    //int ret;
    char outbuffer[32768];
    std::string outstring;

    // get the decompressed bytes blockwise using repeated calls to inflate
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        zs_ret = inflate(&zs, 0);

        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer, zs.total_out - outstring.size());
        }

    } while (zs_ret == Z_OK);

    inflateEnd(&zs);

    if (zs_ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib decompression: (" << zs_ret << ") "
            << zs.msg;
        throw(std::runtime_error(oss.str()));
        exit(0); //exit when error occur during decompression

    }


    return outstring;
}



void read_binary_block(stringstream &t_stream
    , unsigned long &bytes_per_feature
    , unsigned long &bytes_per_value
    , int &feature_length
    , string &n_key
    , double &n_value
    , unsigned long &vocab_size
    )
{
	char t_chars[64];

    if (bytes_per_feature==0 && bytes_per_value==0) ///read header information
	{
        // determine vocab size can be done using file size (in bytes)
        // in bytes. (file_size - (header: sizeof(bytes_per_feature) + sizeof(bytes_per_value) + sizeof(feature_length)) ) / (bytes_per_feature + bytes_per_value)
        t_stream.seekg(0, ios::end);
        vocab_size = (unsigned long)(t_stream.tellg()) - sizeof(bytes_per_feature) - sizeof(bytes_per_value) - sizeof(feature_length);
        t_stream.seekg(0, ios::beg);


		t_stream.read(t_chars, sizeof(bytes_per_feature));
		bytes_per_feature = *reinterpret_cast<const unsigned long*>(t_chars);
		memset(&(t_chars[0]), 0, sizeof(t_chars)); //is this necessary?

		t_stream.read(t_chars, sizeof(bytes_per_value));
		bytes_per_value = *reinterpret_cast<const unsigned long*>(t_chars); //( t_stream.get(sizeof(bytes_per_value)).data() ); //expand 127 limit to 255 limit (4 letters -> 800, 20 letters -> 320)
		memset(&(t_chars[0]), 0, sizeof(t_chars));

		t_stream.read(t_chars, sizeof(feature_length));
		feature_length = *reinterpret_cast<const int*>(t_chars); //( t_stream.get(sizeof(feature_length)).data() ); //expand 127 limit to 255 limit (4 letters -> 800, 20 letters -> 320)
		memset(&(t_chars[0]), 0, sizeof(t_chars));

		n_key.resize(bytes_per_feature);
        
        // calculate vocab_size
        vocab_size = vocab_size / (bytes_per_feature + bytes_per_value);

	}

    t_stream.read(reinterpret_cast<char*>(&n_key[0]), bytes_per_feature); ///read feature
    t_stream.read(reinterpret_cast<char*>(&n_value), bytes_per_value);

}


void JSD_divergence(string &p_key
    , string &q_key
    , double &p_value
    , double &q_value
    , double &Hp
    , double &Hq
    , double &Hm
    )
{

    if ((p_key > q_key && q_value!=0) || p_value==0) ///add q_value
    {
        Hq-=q_value * log2(q_value);
        Hm-=q_value * log2(0.5 * q_value); ///p_value==0

        q_value=0;

    } else if ((p_key < q_key && p_value!=0) || q_value==0) ///add p_value
    {
        Hp-=p_value * log2(p_value);
        Hm-=p_value * log2(0.5 * p_value); ///q_value==0

        p_value=0;

    } else if (p_key==q_key && p_value!=0 && q_value!=0) ///p_key==q_key neither values are 0
    {
        Hm-=(p_value + q_value) * log2(0.5 * (p_value + q_value));
        Hp-=p_value * log2(p_value);
        Hq-=q_value * log2(q_value);

        p_value=0;
        q_value=0;

    }

}


void JSD_divergence_size_adjustment(
    string &p_key
    , string &q_key
    , double &p_value
    , double &q_value
    , double &Hp
    , double &Hq
    , double &Hm
    //, double &Hm_p
    //, double &Hm_q
    , double &p_vratio
    , double &q_vratio
    )
{

    if ((p_key > q_key && q_value!=0) || p_value==0) ///add q_value
    {
        //Hq-=q_value * log(q_value) / log(2);
        //Hm-=q_value * log(q_vratio * q_value) / log(2); ///p_value==0

        Hq-=q_value * log2(q_value);
        Hm-=q_vratio * q_value * log2(q_vratio * q_value); ///p_value==0

        q_value=0;
        
        //q_cnt+=1.0;

    } else if ((p_key < q_key && p_value!=0) || q_value==0) ///add p_value
    {
        //Hp-=p_value * log(p_value) / log(2);
        //Hm-=p_value * log( p_vratio * p_value) / log(2); ///q_value==0

        Hp-=p_value * log2(p_value);
        Hm-=p_vratio * p_value * log2(p_vratio * p_value);

        p_value=0;
        
        //p_cnt+=1.0;

    } else if (p_key==q_key && p_value!=0 && q_value!=0) ///p_key==q_key neither values are 0
    {
        //Hm-=(p_value + q_value) * log(0.5 * (p_value + q_value)) / log(2);
        //Hm-=( (p_vratio * p_value) + (q_vratio * q_value) ) * log( (p_vratio * p_value) + (q_vratio * q_value) ) / log(2);

        Hm-=( (p_vratio * p_value) + (q_vratio * q_value) ) * log2( (p_vratio * p_value) + (q_vratio * q_value) );

        Hp-=p_value * log2(p_value);
        Hq-=q_value * log2(q_value);

        p_value=0;
        q_value=0;
        
        //p_cnt+=1.0;
        //q_cnt+=1.0;

    }

}


///Kullback-Leibler distance. symmetrized relative entropy: ( KL(A||B) + KL(B||A) ) / 2
// by definition KL cannot handle if P(i) or Q(i) is 0 inside log, for instance, log2(P(i)) or log2(Q(i))
void KL_distance(string &p_key
    , string &q_key
    , double &p_value
    , double &q_value
    , double &Hp
    , double &Hq
    , double &Hm
    )
{

    //cout << std::numeric_limits<double>::min() << endl; //hte smallest double number defined
    //to use KL, set 0 probability to the smallest double number (laplace's pseudo number)
    ///*
    if ((p_key > q_key && q_value!=0) || p_value==0) ///add q_value, p_value
    {
        Hq-=q_value * log2(q_value);

        Hp-=std::numeric_limits<double>::min() * log2(q_value);
        Hm-=q_value * log2(std::numeric_limits<double>::min()) + std::numeric_limits<double>::min()*log2(q_value);

        q_value=0;

    } else if ((p_key < q_key && p_value!=0) || q_value==0) ///add p_value
    {
        Hp-=p_value * log2(p_value);

        Hq-=std::numeric_limits<double>::min() * log2(p_value);
        Hm-=p_value * log2(std::numeric_limits<double>::min()) + std::numeric_limits<double>::min()*log2(p_value);

        p_value=0;

    } else if (p_key==q_key) // && p_value!=0 && q_value!=0) ///p_key==q_key neither values are 0
    {
        Hp-=p_value * log2(p_value);
        Hq-=q_value * log2(q_value);
        Hm-=(p_value * log2(q_value) + q_value * log2(p_value));

        p_value=0;
        q_value=0;

    }
    //*/

}


///Jaccard distance: 1 - ( n(intersection set) / n(union set) )
void JACCARD_distance(string &p_key
    , string &q_key
    , double &p_value
    , double &q_value
    , double &Hp
    , double &Hq
    , double &Hm
    )
{

    if ((p_key > q_key && q_value!=0) || p_value==0) ///only in set(q)
    {
        Hq+=1.0;

        q_value=0;

    } else if ((p_key < q_key && p_value!=0) || q_value==0) ///only in set(p)
    {
        Hp+=1.0;

        p_value=0;

    } else if (p_key==q_key && p_value!=0 && q_value!=0) ///set(p) and set(q) share
    {
        Hm+=1.0;

        p_value=0;
        q_value=0;

    }

}



double calculate_distance(string p_path, string q_f_buf
    , int delimiter_int
    , string distance_type
    )
{
    ///p -> istringstream, q -> cref, istringstream
    double p_value=0.0;
    double q_value=0.0;

    string p_key="";
    string q_key="";

    double Hm=0.0;
    //size adjustment variatns
    double Hm_p=0.0;
    double Hm_q=0.0;

    double Hp=0.0;
    double Hq=0.0;
    double r_value=0.0; //return distance value

    int p_zs_ret=0;
    unsigned long p_bytes_per_feature=0;
    unsigned long p_bytes_per_value=0;
	int p_feature_length=0;
    // double p_cnt=0.0; //p(x) vocabulary size
    unsigned long p_vocab_size = 0;

	//int q_zs_ret=0;
    unsigned long q_bytes_per_feature=0;
    unsigned long q_bytes_per_value=0;
	int q_feature_length=0;
    // double q_cnt=0.0; //q(x) vocabulary size
    unsigned long q_vocab_size = 0;

    ///set read_p_f
    ifstream read_f;
    string p_f_buf;

    read_f.open(p_path.c_str(), ios::in|ios::binary);
    read_f.seekg(0, ios::end);
    p_f_buf.resize(read_f.tellg()); //reach max size.
    read_f.seekg(0, ios::beg);
    read_f.read(&p_f_buf[0], p_f_buf.size());
    read_f.close();

    stringstream read_p_f(decompress_deflate(p_f_buf, p_zs_ret), ios::in|ios::out|ios::binary);
    read_p_f.seekg(0, ios::beg);
    p_f_buf.clear(); // possible?


    ///set read_q_f
    stringstream read_q_f(q_f_buf, ios::in|ios::out|ios::binary);
    read_q_f.seekg(0, ios::beg);

    ///read header info and check feature length
    read_binary_block(read_p_f, p_bytes_per_feature, p_bytes_per_value, p_feature_length, p_key, p_value, p_vocab_size);
    read_binary_block(read_q_f, q_bytes_per_feature, q_bytes_per_value, q_feature_length, q_key, q_value, q_vocab_size);

    // should be converted to double before division
    double p_vratio = double(p_vocab_size) / double(p_vocab_size + q_vocab_size);
    double q_vratio = double(q_vocab_size)/ double(p_vocab_size + q_vocab_size);

    ///anything wrong during inflation/decompression step, or comparing with different feature lengths, should output an error value (-1)
    if (p_feature_length!=q_feature_length) //fool proof
    {
        cerr << "Different feature_lengths compared: " << p_feature_length << " vs. " << q_feature_length << endl;

        read_p_f.str(string());
        read_p_f.clear(); //clear flag
        p_f_buf.clear();

        read_q_f.str(string());
        read_q_f.clear(); //clear flag

        //r_value=-1.0;
        exit(0);

    }

    while (!read_p_f.eof() || !read_q_f.eof()) ///run until both files reach EOF()
    {
        if (p_value==0 && !read_p_f.eof())
        {
            read_binary_block(read_p_f, p_bytes_per_feature, p_bytes_per_value, p_feature_length, p_key, p_value, p_vocab_size);

        }

        if (q_value==0 && !read_q_f.eof())
        {
            read_binary_block(read_q_f, q_bytes_per_feature, q_bytes_per_value, q_feature_length, q_key, q_value, q_vocab_size);

        }

        if (p_value!=0 || q_value!=0) //requires to avoid nan
        {
            ///support various distances (or dissimilarity); can add more in future
            if (distance_type=="jsdiv" || distance_type=="jsdist") //JS-divergence
            {
                JSD_divergence(p_key, q_key, p_value, q_value, Hp, Hq, Hm);

            } else if (distance_type=="jacc") //JS-distance
            {
                JACCARD_distance(p_key, q_key, p_value, q_value, Hp, Hq, Hm);

            } else if (distance_type=="jsda") //JS-distance
            {
                JSD_divergence_size_adjustment(p_key, q_key, p_value, q_value, Hp, Hq, Hm, p_vratio, q_vratio);

            } else if (distance_type=="kls") //KL relative entropy symmetrized
            {
                KL_distance(p_key, q_key, p_value, q_value, Hp, Hq, Hm);

            }

        }
    }

    ///clear variables and containers
    read_p_f.str(string());
    read_p_f.clear(); //clear flag
    p_f_buf.clear();

    read_q_f.str(string());
    read_q_f.clear(); //clear flag
    ///q_f_buf is a constant reference


	if (r_value==-1)
	{
    	return r_value;
	}


    if (distance_type=="jsdiv")
    {
        r_value = 0.5*(Hm - Hp - Hq); //[0,1] bound

    } else if (distance_type=="jsdist") //JS-divergence
    {
        r_value = sqrt(0.5*(Hm - Hp - Hq)); //[0,1] bound

    } else if (distance_type=="jacc") //JS-distance
    {
        r_value = (Hp + Hq) / (Hp + Hq + Hm);

    } else if (distance_type=="jsda") //size-weighted JS-distance
    {
        r_value = Hm - (p_vratio*Hp) - (q_vratio*Hq);

    } else if (distance_type=="kls") //KL relative entropy symmetrized
    {
        r_value = 0.5*(Hm - Hp - Hq); //[0, x) lower bounded

    }

	return r_value;

}



struct future_handle
{
    int n_row;
    int n_col;
    future<double> n_fut;
    bool in_act=false;

};


void to_square_matrix_output(stringstream &output_stream, bool item_tab_flag) ///convert a low triangular matrix to a square matrix
{

    output_stream.seekg(0, ios::beg);
    string read_line="";

    vector<string> line_vector;
    vector<vector<string>> matrix_vector;

    int item_count=0;
    string item_count_str="";

    std:getline(output_stream, item_count_str); ///first line is a number of items
    item_count = atoi(item_count_str.c_str());

    string item_name="";

    string value_string=""; ///values + tabs string chunck
    string each_value_string=""; ///individual value string


    while (getline(output_stream, read_line)) ///read each line; line break as a delimiter
    {
        if (item_tab_flag==true)
        {
            item_name = read_line.substr(0, read_line.find("\t"));
            value_string = read_line.substr(read_line.find("\t")+1);

        } else //fixed phylip format
        {
            item_name = read_line.substr(0, 10); ///front 10 letters are item_name
            value_string = read_line.substr(10);
        }

        line_vector.push_back(item_name);

        for (string::iterator s_it=value_string.begin(); s_it!=value_string.end(); s_it++) ///tab tokenize
        {
            if (*s_it=='\t' && each_value_string!="") ///tab as a delimiter
            {
                line_vector.push_back(each_value_string);
                each_value_string=""; ///each_value_string.clear(); clear

            } else
            {
                each_value_string+=(*s_it); ///sum

            }

        }

        if (each_value_string!="")
        {
            line_vector.push_back(each_value_string);
            each_value_string=""; ///each_value_string.clear(); clear

        }

        line_vector.resize(item_count+1, ""); ///resize container (item_name + values); extra spaces remain empty ""
        matrix_vector.push_back(line_vector);

        line_vector.clear();

    }


    output_stream.clear();
    output_stream.str(""); //clear previous contents
    output_stream.seekg(0, ios::beg);

    output_stream << item_count_str << endl; ///a number of items

    ///check and reform to a sqaure matrix (symmetric)
    for (int r_cnt=0; r_cnt!=item_count; r_cnt++) ///row
    {
        output_stream << matrix_vector[r_cnt][0]; ///first column is item_name (not value)

        if (item_tab_flag==true) output_stream << "\t"; //add tab after the item_name

        for (int c_cnt=1; c_cnt!=item_count+1; c_cnt++) ///1 based
        {
            if (matrix_vector[r_cnt][c_cnt]=="" && matrix_vector[c_cnt-1][r_cnt+1]!="") ///if empty find a corresponding value; the corresponding value shouldn't be empty (either low or high triangular matrix)
            {
                output_stream << matrix_vector[c_cnt-1][r_cnt+1];

            } else
            {
                output_stream << matrix_vector[r_cnt][c_cnt];

            }

            if (c_cnt!=item_count) ///if not the last
            {
                output_stream << "\t"; ///add delimiters between values

            }

        }

        output_stream << "\n"; ///add line breaks

    }

}



void print_value_vector_str(stringstream &output_stream
    , vector< vector<double> > &fut_value_vector
    , vector<string> &load_path_vector
    , int start_item_n
    , bool symmetric_flag
    , bool item_tab_flag
    )
{
    size_t stream_size_t=0;
    char* stream_buf=NULL;

    string base_name_str; //basename of path(file name)

    if (start_item_n==0)
    {
        stream_size_t = snprintf(NULL, 0, "%d\n", (int)load_path_vector.size());
        stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

        snprintf(stream_buf, size_t(stream_buf), "%d\n", (int)load_path_vector.size());
        output_stream << stream_buf; //item_size
		//free(stream_buf);
    }


    for (vector< vector<double> >::size_type r_it=start_item_n; r_it<load_path_vector.size(); r_it++)
    {
        base_name_str = load_path_vector[r_it].substr(load_path_vector[r_it].rfind('/')+1);

        ///9 characters maximum; for PHYLIP format
        if (item_tab_flag==false) //go for phylip format
        {
            if (base_name_str.length()>9)
            {
                base_name_str.erase(base_name_str.begin()+9, base_name_str.end()); //in c++, string variable end with 'string/' so require -1

            }

            stream_size_t = snprintf(NULL, 0, "%-10s", base_name_str.c_str()); ///item names, with a character limit 10(or 9)
            stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

            snprintf(stream_buf, size_t(stream_buf), "%-10s", base_name_str.c_str());

        } else //use tab
        {
            stream_size_t = snprintf(NULL, 0, "%s\t", base_name_str.c_str()); ///item names, with no length limit and use tab as a separator
            stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

            snprintf(stream_buf, size_t(stream_buf), "%s\t", base_name_str.c_str());
        }

        output_stream << stream_buf;


        for (vector<double>::size_type c_it=0; c_it<fut_value_vector[r_it].size(); c_it++)
        {

            ///print, a point below 8 decimal places
            stream_size_t = snprintf(NULL, 0, "%.8g\t", fut_value_vector[r_it][c_it]);
            stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));
            snprintf(stream_buf, size_t(stream_buf), "%.8g\t", fut_value_vector[r_it][c_it]);

            output_stream << stream_buf;

        }

        //diagonal self distance, is 0 and not actually calculated
        stream_size_t = snprintf(NULL, 0, "%.8g\n", 0.0);
        stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

        snprintf(stream_buf, size_t(stream_buf), "%.8g\n", 0.0);
        output_stream << stream_buf;

        fut_value_vector[r_it].clear(); //empty vector
    }

    free(stream_buf);

    //fu`t_value_vector.clear();

    if (symmetric_flag==true)
    {
        to_square_matrix_output(output_stream, item_tab_flag); ///convert a low triangular matrix to a square matrix

    }

    cout << output_stream.str(); ///output a matrix

}


///calculate JS Divergence using multiple threads
void multi_thread_manage(vector< vector<double> > &fut_value_vector
    , vector<string> &load_path_vector
    , int thread_n_limit
    , int delimiter_int
    , int start_item_n
    , string distance_type
    )
{
    ///prepare multi-threading
    future_handle fut_struct[thread_n_limit];
    fut_value_vector.resize(load_path_vector.size());

    int n_row=0;
    int n_col=0;

    int rev_pos=0;
    bool input_occur=false;

    ///prepare a constant reference
    ifstream read_f;
    string q_f_buf;
    string q_decompress_buf;
    int q_zs_ret=0;


    for (vector<vector<double>>::size_type r_it=start_item_n; r_it<load_path_vector.size(); ++r_it)
    {
        fut_value_vector[r_it].resize(r_it);

        ///shared across threads
        if (r_it!=0)
        {
            ///efficient method to copy a whole file content to memory
            ///http://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
            read_f.open(load_path_vector[r_it].c_str(), ios::in|ios::binary);
            read_f.seekg(0, ios::end);
            q_f_buf.resize(read_f.tellg()); //reach max size.
            read_f.seekg(0, ios::beg);
            read_f.read(&q_f_buf[0], q_f_buf.size());
            read_f.close();

        }

        ///independent per thread
        for (int cy1=0; cy1!=thread_n_limit; ++cy1)
        {
            if (fut_struct[cy1].in_act==true)
            {
                n_row=fut_struct[cy1].n_row;
                n_col=fut_struct[cy1].n_col;
                fut_value_vector[n_row][n_col]=fut_struct[cy1].n_fut.get();

                fut_struct[cy1].in_act=false;

            }

        }

        if (!q_f_buf.empty()) ///clear a constant reference
        {
            q_decompress_buf = decompress_deflate(q_f_buf, q_zs_ret);

            if (q_zs_ret!=Z_STREAM_END) //q_file decompress sanity check
            {
                cerr << "Z_STREAM | decompression unsuccessful: " << load_path_vector[r_it] << endl;
                exit(0);
            }

            q_f_buf.clear();

        }


        for (vector<vector<double>>::size_type c_it=0; c_it<r_it; ++c_it)
        {
            input_occur=false;

            while (input_occur==false)
            {
                if (rev_pos>=thread_n_limit)
                {
                    rev_pos=0;
                }

                for (int cy1=rev_pos; cy1!=thread_n_limit; ++cy1)
                {
                    ///check running-finished task
                    if (fut_struct[cy1].in_act==true && fut_struct[cy1].n_fut.wait_for(std::chrono::milliseconds(0))==future_status::ready)
                    {
                        n_row=fut_struct[cy1].n_row;
                        n_col=fut_struct[cy1].n_col;

                        fut_value_vector[n_row][n_col]=fut_struct[cy1].n_fut.get(); //now future status vecome valid()

                        fut_struct[cy1].in_act=false;

                    }

                    ///input new task
                    if (fut_struct[cy1].in_act==false)
                    {
                        fut_struct[cy1].n_row=r_it;
                        fut_struct[cy1].n_col=c_it;

                        ///input is [feature][feature value] format
                        fut_struct[cy1].n_fut=async(launch::async, calculate_distance, load_path_vector[c_it], std::cref(q_decompress_buf), delimiter_int, distance_type);

                        fut_struct[cy1].in_act=true;
                        input_occur=true;

                        rev_pos++; ///reserve next thread number

                        break;

                    }

                }

            }

        }

		///final check to finsih remaining running thread
		for (int cy1=0; cy1!=thread_n_limit; ++cy1)
		{
			if (fut_struct[cy1].in_act==true)
			{
				n_row=fut_struct[cy1].n_row;
				n_col=fut_struct[cy1].n_col;
				fut_value_vector[n_row][n_col]=fut_struct[cy1].n_fut.get();

				fut_struct[cy1].in_act=false;
			}

		}

    }

    q_decompress_buf.clear();
    //print_value_vector(fut_value_vector, load_path_vector);

}

int read_reserved_matrix(stringstream &output_stream
    , vector<string> &load_path_vector
    , string &reserve_matrix_path
    , bool item_tab_flag
    )
{
    ifstream read_f;
    string read_line;
    int start_item_n=0;

    read_f.open(reserve_matrix_path.c_str(), ios::in);

    getline(read_f, read_line, '\n'); ///first line of matrix indicate the number of items
    start_item_n = stoi(read_line);

    string item_name_str;
    vector<string> resv_load_path_vector; /// in general item name is file name as taxonID

    output_stream << to_string(load_path_vector.size()) << '\n'; //item number

	std::size_t str_delim = 0;


    while (getline(read_f, read_line, '\n')) //2021-1 needed to change
    {
		str_delim = (item_tab_flag==false) ? 10 : read_line.find("\t"); ///phylip, the first 10 letters, or no length limit, separated by tab

		item_name_str = read_line.substr(0, str_delim); //first 10 characters
        item_name_str.erase(remove(item_name_str.begin(), item_name_str.end(), ' '), item_name_str.end()); //trim spaces

        resv_load_path_vector.push_back(item_name_str);

        output_stream << read_line << '\n';

    }

    read_f.close();


    vector<string>::iterator it_pos;
    int replaced_path_cnt=0;
    int reserve_matrix_cnt=resv_load_path_vector.size();

    for (vector<string>::iterator it=load_path_vector.begin(); it!=load_path_vector.end(); ++it)
    {
        item_name_str = (*it).substr((*it).rfind('/')+1);

        if (item_name_str.length()>9 && item_tab_flag==false)
        {
            item_name_str.erase(item_name_str.begin()+9, item_name_str.end()); ///in c++, string variable end with 'string/' so require -1

        }

        it_pos = find(resv_load_path_vector.begin(), resv_load_path_vector.end(), item_name_str);

        if (it_pos!=resv_load_path_vector.end()) //exist
        {
            resv_load_path_vector[int(distance(resv_load_path_vector.begin(), it_pos))] = *it; //replace item name to file path
            replaced_path_cnt++;

        } else //new
        {
            resv_load_path_vector.push_back(*it); //add new file path

        }

    }

    if (replaced_path_cnt!=reserve_matrix_cnt) //means not fully covering reserved matrix items
    {
        start_item_n=-1;
    }


    load_path_vector = resv_load_path_vector; ///replace original load_path_vector to order reserved(including new) load_path_vector

    return start_item_n; ///return starting point, 0 base
}



void show_help()
{
    cout << "Parameters: [option][load_paths(requires full path)]\n";
    cout << "--help | -h, show help\n";
    cout << "--version | -v, show help\n";
    cout << "--thread | -t [int], set number of threads (default = 5)\n";
    cout << "--reserved | -r [path], input reserved distance matrix\n";
    cout << "--symmetric | -s, output a symmetric matrix; default is a low triangular matrix\n";
    cout << "-T, Use TAB as a separator between row names and distances. Default is PHYLIP format that limit row names up to 9 characters\n";
    
    // various distances
    cout << "--distance | -d [str], type of distance or dissimilarity metric\n";

    cout << "\tjsdiv : Jensen-Shannon Divergence (JS divergence)\n";
    
    cout << "\tjsdist: Jensen-Shannon Distance (JS distance = square_root(JS divergence))\n";
    
    cout << "\tjacc : Jaccard Distance\n";
    
    // https://stats.stackexchange.com/questions/97938/calculate-the-kullback-leibler-divergence-in-practice
    cout << "\tkls : Symmetrized (Kullback-Leibler) relative entropy; [0, inf)\n";
    cout << "\t\tKL assumes (P(i) or Q(i) !=0 or 1) and cannot define 0 probability\n";
    cout << "\t\tThus, use Laplace's pseudonumber to improvise (std::numeric_limits<double>::min() !=0)\n";

    //experimental feature
    cout << "\t(Experimental distances)\n";
    cout << "\tjsda : Size-weighted Jensen-Shannon Divergence\n";
    cout << "\t\tWeighted by vocabulary size; Dist(P||Q) = a*KL(P||a*P + b*Q) + b*KL(P||a*P + b*Q); a+b = 1.0, a <> b, a and b are vocabulary size ratio\n";

}


void show_profile()
{
    cout << "JSD distance calculate; 2v.4.0\n";
    cout << "Value presentation: a poin below 8 decimal places (%.8g)\n";
    cout << "Code by JaeJin Choi; https://github.com/jaejinchoi/FFP\n";
    // cout << "Compile; g++ -std=c++11 -pthread -o (output) (this script) -lz\n";
    cout << "Require; zlib 1.2.8+\n";

}


int main(int argc, char** argv)
{

    int opt;
    //int opt_index;
    int delimiter_int = 10; //char(10)=="\n"
    int thread_n_limit=5; ///default, 5 threads

    bool symmetric_flag=false; ///default, output a low triangular matrix
    bool js_distance_flag=false; ///default, output JS divergence
    string reserve_matrix_path="";

    bool item_tab_flag=false; //default is false, limit item name length by 9 characters, in PHYLIP format. true accept full item names and use tab as a separator

    ///check supported distance_type
    string distance_type="jsdiv"; //default is jsd : js-divergence. using string form
    vector<string> able_distance_type = {
        "jsdiv" //JS-divergence
        , "jsdist" //JS-distance
        , "jacc" //Jaccard
        , "jsda" //JS-divergence weighted (experimental)
        , "kls" //KL relative entropy symmetrized
        };

    //using getopt_long (requires struct defining long form options)
    static struct option long_options[] = {
        {"thread", required_argument, NULL, 't'}, //NULL or nullptr?
        {"reserved", required_argument, NULL, 'r'}, //using long-form only
        {"distance", required_argument, NULL, 'd'},
        {"symmetric", no_argument, NULL, 's'},
        {"tab", no_argument, NULL, 'T'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // while ((opt = getopt(argc, argv, "ht:r:d:vsT")) !=EOF) #simple argument take
    while ((opt = getopt_long(argc, argv, "ht:r:d:vsT", long_options, NULL)) !=EOF) // EOF = -1
    {
        switch (opt)
        {
            case 'h':
                show_help();
                exit(EXIT_SUCCESS);

            case 'v':
                show_profile();
                exit(EXIT_SUCCESS);

            case 't':
                thread_n_limit=max(1, atoi(optarg)); ///least 1 thread or not run
                break;

            ///*
            case 'r':
                reserve_matrix_path=string(optarg);
                break;
            //*/

            case 'd':
                distance_type = string(optarg);
                break;

            case 's':
                symmetric_flag=true;
                break;

            case 'T':
                item_tab_flag=true; //enable TAB as an item name separator; disable PHYLIP format
                break;

            //*
            case '?':
            default:
                //show_help();
                exit(0);

            //*/
        }
    }

    if (std::find(able_distance_type.begin(), able_distance_type.end(), distance_type)==able_distance_type.end()) //not in
    {
        cerr << "Unsupported distance type requested: " << distance_type << endl;
        exit(0);
    }


    vector<string> load_path_vector;
    vector< vector<double> > fut_value_vector;

    int start_item_n=0; ///a reserved_matrix parameter

    stringstream output_stream(ios::out|ios::in|ios::ate);
    output_stream.str(string());
    int input_load_path_cnt=0;

    if (argc - optind > 1)
    {
        while (argc!=optind)
        {
            load_path_vector.push_back(argv[optind]);
            optind++;

        }

        input_load_path_cnt=load_path_vector.size();

        sort(load_path_vector.begin(), load_path_vector.end());

        if (reserve_matrix_path!="")
        {
            ///print reserved matrix(or given), replace original load_path_vector order and element
            start_item_n = read_reserved_matrix(output_stream, load_path_vector, reserve_matrix_path, item_tab_flag);

        }

        if (start_item_n==-1)
        {
            printf("Check if items in reserved matrix are included in input files\n");

        } else
        {
            multi_thread_manage(fut_value_vector, load_path_vector, thread_n_limit, delimiter_int, start_item_n, distance_type); ///multi-threading; calculate JS Divergence

            print_value_vector_str(output_stream, fut_value_vector, load_path_vector, start_item_n, symmetric_flag, item_tab_flag); ///final output (standard output)

        }

        output_stream.str(string());
        output_stream.clear();


    } else
    {
        cerr << "A number of input items should be least a pair" << endl;

    }

}

