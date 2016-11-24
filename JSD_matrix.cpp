#include <iostream>
using namespace std;

#include <cstdlib>

#include <string>
using std::string;

#include <cstring>
//to use memset

#include <sstream>
#include <getopt.h>
//#include <unordered_map>

#include <fstream>
//using ifstream
#include <cmath>
//#include <list>
#include <vector>
using std::vector;
//#include <thread> //for multithreading-however require c++11,
//using std::thread;
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

//#include <pthread.h> //basic thread option, thread is advanced, c++11, thread  is more convenient

#include <zlib.h>
#include <stdexcept>

// Found these here http://mail-archives.apache.org/mod_mbox/trafficserver-dev/201110.mbox/%3CCACJPjhYf=+br1W39vyazP=ix
//eQZ-4Gh9-U6TtiEdReG3S4ZZng@mail.gmail.com%3E
#define MOD_GZIP_ZLIB_WINDOWSIZE 15
#define MOD_GZIP_ZLIB_CFACTOR    9
#define MOD_GZIP_ZLIB_BSIZE      8096


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


/// string key, double value
void split_str(string &str_buffer, string &str_key, double &re_value)
{
    //str_buffer="5c3c|23612";
    //str_buffer="22016|5.9207e-04";
    long long delimit_position=0; //fix delimitor as '|'
    delimit_position = str_buffer.find('|');

    str_key = str_buffer.substr(0, delimit_position);
    re_value = atof(str_buffer.substr(delimit_position+1).data());

}


string decompress_key_text(stringstream &t_stream, vector<string> &t_vector, int delimiter_int)
{
    ///t_vector -> independent, use index as a key; and value as replacement
    string decompressed_str;
    string str_key;
    double re_value;
    char get_char;

    bitset<8> current_bitset=0;
    bool key_sign_flag=false;

    while (t_stream.get(get_char))
    {
        current_bitset = bitset<8>(get_char);

        if (get_char=='#') //'#' as key signature
        {
            key_sign_flag=true;

        } else if (get_char==char(delimiter_int) && decompressed_str!="") //as linedelimiter
        {
            if (key_sign_flag==true) //if key line, save to t_vector
            {
                split_str(decompressed_str, str_key, re_value);
                t_vector[(int)(re_value)]=str_key;
                decompressed_str.clear();
                key_sign_flag=false;

            } else
            {
                //cout << decompressed_str << endl;
                return decompressed_str;

            }

        } else if (current_bitset[7]==1) //compress marked
        {
            current_bitset[7]=0;
            decompressed_str+=t_vector[(int)(current_bitset.to_ulong())];

        } else //normal char
        {
            decompressed_str+=get_char;

        }

    }

}

double jsd_distance_key_sort(string p_path, string q_decompress_buf, int delimiter_int)
{
    ///p -> istringstream, q -> cref, istringstream
    string p_buffer="";
    string q_buffer="";

    double p_value=0;
    double q_value=0;

    string p_key="";

    string index_q_key="";
    string q_key="";

    double Hm=0.0;
    double Hp=0.0;
    double Hq=0.0;

    //decompress key vector
    vector<string> p_vector;
    p_vector.resize(128);

    vector<string> q_vector;
    q_vector.resize(128);

    ///set read_p_f
    ifstream read_f;
    string p_f_buf;

    read_f.open(p_path.c_str(), ios::in|ios::binary);
    read_f.seekg(0, ios::end);
    p_f_buf.resize(read_f.tellg()); //reach max size.
    read_f.seekg(0, ios::beg);
    read_f.read(&p_f_buf[0], p_f_buf.size());
    read_f.close();

    stringstream read_p_f(decompress_deflate(p_f_buf), ios::in|ios::out|ios::binary);
    //istringstream read_p_f(p_f_buf, ios::in|ios::binary);
    read_p_f.seekg(0, ios::beg);
    p_f_buf.clear(); // possible?

    ///set read_q_f
    stringstream read_q_f(q_decompress_buf, ios::in|ios::out|ios::binary);
    //istringstream read_q_f(q_f_buf, ios::in|ios::binary);
    read_q_f.seekg(0, ios::beg);

    while (!read_p_f.eof() || !read_q_f.eof())
    {

        if (p_value==0 && !read_p_f.eof())
        {
            p_buffer = decompress_key_text(read_p_f, p_vector, delimiter_int);

            if (p_buffer!="\n" and p_buffer!="")
            {
                split_str(p_buffer, p_key, p_value);

            } else
            {
                p_value=0;

            }

        }

        if (q_value==0 && !read_q_f.eof())
        {
            q_buffer = decompress_key_text(read_q_f, q_vector, delimiter_int);

            if (q_buffer!="\n" and q_buffer!="")
            {
                split_str(q_buffer, q_key, q_value);

            } else
            {
                q_value=0;
            }

        }


        //printf("p: %s:%f \t q: %s:%f\n", p_key.c_str(), p_value, q_key.c_str(), q_value);
        if (p_value==0 && q_value==0)
        {
            continue;

        } else if (p_key==q_key)
        {
            Hm-=(p_value + q_value) * log(0.5 * (p_value + q_value)) / log(2);
            Hp-=p_value * log(p_value) / log(2);
            Hq-=q_value * log(q_value) / log(2);

            p_value=0;
            q_value=0;

        } else if ((p_key > q_key && q_value!=0) || p_value==0) //add q_value
        {
            Hq-=q_value * log(q_value) / log(2);
            Hm-=q_value * log(0.5 * q_value) / log(2); //p_value==0

            q_value=0;

        } else if ((p_key < q_key && p_value!=0) || q_value==0) //add p_value
        {
            Hp-=p_value * log(p_value) / log(2);
            Hm-=p_value * log(0.5 * p_value) / log(2); //q_value==0

            p_value=0;
        }

        //printf("Hp: %f \t Hq: %f \t Hm: %f\n", Hp, Hq, Hm);

    }

    ///2015-3 empty container first and then clear the flag is proper?
    read_p_f.str(string()); //empty container/buffer?
    read_p_f.clear(); //clear flag
    p_f_buf.clear(); //

    read_q_f.str(string()); //empty container/buffer?
    read_q_f.clear(); //clear flag


    return 0.5*(Hm - Hp - Hq);

}



struct future_handle
{
    int n_row;
    int n_col;
    future<double> n_fut;
    bool in_act=false;

};



void print_value_vector_str(stringstream &output_stream, vector< vector<double> > &fut_value_vector, vector<string> &load_path_vector, int start_item_n, bool js_distance_flag)
{
    size_t stream_size_t=0;
    char* stream_buf=NULL;

    //stringstream output_stream(ios::out|ios::in|ios::ate);
    //output_stream.str(string());

    string base_name_str; //basename of path(file name)
    //output matrix-in neighbor joining input format
    //cout << load_path_vector.size() << endl; //number of items

    if (start_item_n==0)
    {
        stream_size_t = snprintf(NULL, 0, "%d\n", (int)load_path_vector.size());
        stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

        snprintf(stream_buf, size_t(stream_buf), "%d\n", (int)load_path_vector.size());
        output_stream << stream_buf; //item_size


    }


    for (vector< vector<double> >::size_type r_it=start_item_n; r_it<load_path_vector.size(); r_it++)
    {
        base_name_str = load_path_vector[r_it].substr(load_path_vector[r_it].rfind('/')+1);
        ///require trim at the end if base_name_str length over 9, to prevent value error
        if (base_name_str.length()>9)
        {
            base_name_str.erase(base_name_str.begin()+9, base_name_str.end()); //in c++, string variable end with 'string/' so require -1

        }

        stream_size_t = snprintf(NULL, 0, "%-10s", base_name_str.c_str()); ///item names. character limits 10(or 9)
        stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

        snprintf(stream_buf, size_t(stream_buf), "%-10s", base_name_str.c_str());
        output_stream << stream_buf;


        for (vector<double>::size_type c_it=0; c_it<fut_value_vector[r_it].size(); c_it++)
        {

            if (js_distance_flag==true) //sqrt(JS divergence)=JS distance
            {
                stream_size_t = snprintf(NULL, 0, "%.4e\t", sqrt(fut_value_vector[r_it][c_it]));
                stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

                snprintf(stream_buf, size_t(stream_buf), "%.4e\t", sqrt(fut_value_vector[r_it][c_it]));

            } else
            {
                stream_size_t = snprintf(NULL, 0, "%.4e\t", fut_value_vector[r_it][c_it]);
                stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

                snprintf(stream_buf, size_t(stream_buf), "%.4e\t", fut_value_vector[r_it][c_it]);

            }

            output_stream << stream_buf;

        }

        stream_size_t = snprintf(NULL, 0, "%.4e\n", 0.0);
        stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));

        snprintf(stream_buf, size_t(stream_buf), "%.4e\n", 0.0);
        output_stream << stream_buf;

        fut_value_vector[r_it].clear(); //empty vector

    }

    free(stream_buf);

    //refresh vector
    fut_value_vector.clear();

    cout << output_stream.str(); //output in one

    //output_stream.str(string());
    //output_stream.clear();

}


///2014-12-10, memory hybrid calculation
//void multi_thread_manage_mem(vector<string> &load_path_vector, int thread_n_limit, int delimiter_int, bool input_have_key)
void multi_thread_manage(vector< vector<double> > &fut_value_vector, vector<string> &load_path_vector, int thread_n_limit, int delimiter_int, int start_item_n)
{
    //create future handle
    future_handle fut_struct[thread_n_limit];
    //create value container
    fut_value_vector.resize(load_path_vector.size());

    int n_row=0;
    int n_col=0;

    int rev_pos=0;
    bool input_occur=false;

    //shared string define
    ifstream read_f;
    string q_f_buf;
    string q_decompress_buf;

    //create value matrix
    for (vector<vector<double>>::size_type r_it=start_item_n; r_it<load_path_vector.size(); ++r_it)
    {
        fut_value_vector[r_it].resize(r_it);

        if (r_it!=0) // and input_have_key==true)
        {
            ///http://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html <-show benchmark
            read_f.open(load_path_vector[r_it].c_str(), ios::in|ios::binary);
            read_f.seekg(0, ios::end);
            q_f_buf.resize(read_f.tellg()); //reach max size.
            read_f.seekg(0, ios::beg);
            read_f.read(&q_f_buf[0], q_f_buf.size());
            read_f.close();

            //q_decompress_buf = decompress_deflate(q_f_buf);
            //q_f_buf.clear();
        } ///give more time to finish any task below, meanwhile do independent work


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

        //q_decompress_buf.clear(); //necessary?

        if (!q_f_buf.empty())
        //if (r_it!=0)
        {
            q_decompress_buf = decompress_deflate(q_f_buf);
            q_f_buf.clear();

        }


        for (vector<vector<double>>::size_type c_it=0; c_it<r_it; ++c_it)
        {
            ///2007-7-10 balance calling multi-threading
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

                        ///input is KEY|VALUE, as a default
                        fut_struct[cy1].n_fut=async(launch::async, jsd_distance_key_sort, load_path_vector[c_it], std::cref(q_decompress_buf), delimiter_int);

                        fut_struct[cy1].in_act=true;
                        input_occur=true;

                        rev_pos++; //reserve next thread number

                        break;

                    }

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

    q_decompress_buf.clear();

    ///print vector value
    //print_value_vector(fut_value_vector, load_path_vector);

}


int reserve_matrix_item_analyze(stringstream &output_stream, vector<string> &load_path_vector, string &reserve_matrix_path)
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

    while (getline(read_f, read_line, '\n'))
    {
        item_name_str = read_line.substr(0, 10); //first 10 characters
        item_name_str.erase(remove(item_name_str.begin(), item_name_str.end(), ' '), item_name_str.end());

        resv_load_path_vector.push_back(item_name_str);

        //reserve matrix input
        output_stream << read_line << '\n';

    }

    read_f.close();


    vector<string>::iterator it_pos;
    int replaced_path_cnt=0;
    int reserve_matrix_cnt=resv_load_path_vector.size();

    for (vector<string>::iterator it=load_path_vector.begin(); it!=load_path_vector.end(); ++it)
    {
        item_name_str = (*it).substr((*it).rfind('/')+1);

        if (item_name_str.length()>9)
        {
            item_name_str.erase(item_name_str.begin()+9, item_name_str.end()); //in c++, string variable end with 'string/' so require -1

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
    cout << "Parameter, [option][load_paths(require full path)]\n";
    cout << "-h, show help\n"; //<< endl;
    cout << "-c [int], specific integer code of single or escape character (default = TAB = char(9))\n";
    //cout << "-k, disable input file which have key and value pair, format='key|value(ratio), but key should be in-order'\n";
    cout << " KEY | VALUE, align output\n";
    //cout << "-l, enable full matrix, default is lower trianglar matrix\n";
    cout << "-t [int], set number of threads(default=5)\n";
    cout << "-r [path], input reserved distance matrix\n";
    cout << "-d, output JS-distance matrix, square-root(JS-divergence) = JS distance\n";

}


void show_profile()
{
    cout << "JSD distance calculation-in matrix(thread support), 2014.11.10" << endl;
    cout << "Code by JaeJin Choi\n";
    cout << "compile; g++ -std=c++11 -pthread -o (output) (this script) -lz\n";
    cout << "string key, double value\n";
    cout << "lower trianglar matrix-without redundant calculation" << endl;
    cout << "need caution in item name should be less than 10 characters(maximum 9), to prevent output error-trim at the end of string length 10" << endl;
    cout << "require '-std=c++11 -lpthread' for correct compile, !may not work if g++ version not support c++11" << endl;
	cout << "2014-11-10, diagonal cell calculation to avoid HDD over-read" << endl;
	cout << "2015-4-20, avoid data race by independently through-out each calculation(independent memory load)\n";
    cout << "2016-1-25; latest zlib acquired; 1.2.8.(May, 2013)\n";
    cout << "2016-1-27, 28; optimize thread managing function and distance calculation function\n";
    cout << "2016-2-2; accept reserved distance matrix input to avoid redo any finished pair-wise calculation\n";
}


int main(int argc, char** argv) //in case of output JSD distance(incomplete)
{

    //parameter with option flag
    int opt;
    int delimiter_int = 10;
    int thread_n_limit=5; //8 as a default number of thread, to prevent overload
    bool input_have_key=true;
    bool js_distance_flag=false; //default is JS divergence, sqr(JS divergence)=JS distance
    string reserve_matrix_path="";

    while ((opt = getopt(argc, argv, "hc:t:r:dv")) !=EOF) // EOF = -1
    {
        switch (opt)
        {
            case 'h':
                show_help();
                exit(EXIT_SUCCESS);

            case 'c':
                delimiter_int = atoi(optarg);
                /*
                if (atoi(optarg)!=0)
                {
                    delimiter = char(atoi(optarg));

                } else
                {
                    cout << "-c [int], input integer of single or escape character(numberic)" << endl;
                    exit(0);

                }
                */
                break;

            case 't':
                thread_n_limit=max(1, atoi(optarg)); //minimum 1
                break;

            case 'r':
                reserve_matrix_path=optarg; //string
                break;

            case 'd':
                js_distance_flag=true;
                break;

            case 'v':
                show_profile();
                exit(EXIT_SUCCESS);

            default:
                show_help();
                exit(0);
        }
    }

    vector<string> load_path_vector;
    //vector<string> resv_load_path_vector;
    vector< vector<double> > fut_value_vector;
    int start_item_n=0; //default = 0

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
        //sort vector
        sort(load_path_vector.begin(), load_path_vector.end());

        if (reserve_matrix_path!="")
        {
            ///print reserved matrix(or given), replace original load_path_vector order and element
            start_item_n = reserve_matrix_item_analyze(output_stream, load_path_vector, reserve_matrix_path);

        }

        if (start_item_n==-1)
        {
            printf("Please check if load_path files also include items in reserved matrix?\n");

        } else
        {
            ///manage all thread works and printif
            multi_thread_manage(fut_value_vector, load_path_vector, thread_n_limit, delimiter_int, start_item_n);
            //cout << "finish_distance" << "\t" << load_path_vector.size() << "\t" << start_item_n << "\t" << fut_value_vector.size() << endl;
            print_value_vector_str(output_stream, fut_value_vector, load_path_vector, start_item_n, js_distance_flag); //print out and refresh all variable

        }

        output_stream.str(string());
        output_stream.clear();

    } else
    {
        cout << "Inputs should be more or equal to two" << endl;

    }

}
