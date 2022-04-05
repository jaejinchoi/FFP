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



void read_binary_block(stringstream &t_stream, unsigned long &bytes_per_feature, unsigned long &bytes_per_value, string &n_key, double &n_value)
{
    //int feature_length=0;

    if (bytes_per_feature==0 && bytes_per_value==0) ///read header information
    {
        bytes_per_feature = t_stream.get() + 0; ///read one byte; feature size
        bytes_per_value = t_stream.get() + 0; ///read one byte; feature value size (long long or double)
        t_stream.get() + 0; ///feature_length; not used but reserved for 1 byte

        n_key.resize(bytes_per_feature);

    }

    t_stream.read(reinterpret_cast<char*>(&n_key[0]), bytes_per_feature); ///read feature
    t_stream.read(reinterpret_cast<char*>(&n_value), bytes_per_value);

    /*
    if (!t_stream.read(reinterpret_cast<char*>(&n_key[0]), bytes_per_feature) || !t_stream.read(reinterpret_cast<char*>(&n_value), bytes_per_value))
    {
        //n_key="";
        n_value=0.0;
    }
    */

}


double jsd_distance(string p_path, string q_f_buf, int delimiter_int)
{
    ///p -> istringstream, q -> cref, istringstream
    double p_value=0.0;
    double q_value=0.0;

    string p_key="";
    string q_key="";

    double Hm=0.0;
    double Hp=0.0;
    double Hq=0.0;


    unsigned long p_bytes_per_feature=0;
    unsigned long p_bytes_per_value=0;

    unsigned long q_bytes_per_feature=0;
    unsigned long q_bytes_per_value=0;


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

    read_p_f.seekg(0, ios::beg);
    p_f_buf.clear(); // possible?


    ///set read_q_f
    stringstream read_q_f(q_f_buf, ios::in|ios::out|ios::binary);
    read_q_f.seekg(0, ios::beg);


    while (!read_p_f.eof() || !read_q_f.eof()) ///run until both files reach EOF()
    {
        if (p_value==0 && !read_p_f.eof())
        {
            read_binary_block(read_p_f, p_bytes_per_feature, p_bytes_per_value, p_key, p_value);

        }

        if (q_value==0 && !read_q_f.eof())
        {
            read_binary_block(read_q_f, q_bytes_per_feature, q_bytes_per_value, q_key, q_value);

        }


        if (p_value==0 && q_value==0)
        {
            continue;

        } else if ((p_key > q_key && q_value!=0) || p_value==0) ///add q_value
        {
            Hq-=q_value * log(q_value) / log(2);
            Hm-=q_value * log(0.5 * q_value) / log(2); ///p_value==0

            q_value=0;

        } else if ((p_key < q_key && p_value!=0) || q_value==0) ///add p_value
        {
            Hp-=p_value * log(p_value) / log(2);
            Hm-=p_value * log(0.5 * p_value) / log(2); ///q_value==0

            p_value=0;

        } else if (p_key==q_key && p_value!=0 && q_value!=0) ///p_key==q_key neither values are 0
        {
            Hm-=(p_value + q_value) * log(0.5 * (p_value + q_value)) / log(2);
            Hp-=p_value * log(p_value) / log(2);
            Hq-=q_value * log(q_value) / log(2);

            p_value=0;
            q_value=0;

        }

    }

    ///clear variables and containers
    read_p_f.str(string());
    read_p_f.clear(); //clear flag
    p_f_buf.clear();

    read_q_f.str(string());
    read_q_f.clear(); //clear flag
    ///q_f_buf is a constant reference

    return 0.5*(Hm - Hp - Hq);

}



struct future_handle
{
    int n_row;
    int n_col;
    future<double> n_fut;
    bool in_act=false;

};



void to_square_matrix_output(stringstream &output_stream) ///convert a low triangular matrix to a square matrix
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
        item_name = read_line.substr(0, 10); ///front 10 letters are item_name
        line_vector.push_back(item_name);

        value_string = read_line.substr(10);

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
    , bool js_distance_flag
    , bool symmetric_flag
    , bool item_tab_flag
    )
{
    size_t stream_size_t=0;
    char* stream_buf=NULL;

    string base_name_str; //basename of path(file name)

    double n_value=0.0;

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

            if (js_distance_flag==true) ///sqrt(JS divergence)=JS distance; output JS distance instead of JS divergence
            {
                n_value = sqrt(fut_value_vector[r_it][c_it]);

            } else
            {
                n_value = fut_value_vector[r_it][c_it];

            }

            ///print, a point below 8 decimal places
            stream_size_t = snprintf(NULL, 0, "%.8g\t", n_value);
            stream_buf = (char*)realloc(stream_buf, (stream_size_t+1)*sizeof(char *));
            snprintf(stream_buf, size_t(stream_buf), "%.8g\t", n_value);

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
        to_square_matrix_output(output_stream); ///convert a low triangular matrix to a square matrix

    }

    cout << output_stream.str(); ///output a matrix

}


///calculate JS Divergence using multiple threads
void multi_thread_manage(vector< vector<double> > &fut_value_vector, vector<string> &load_path_vector, int thread_n_limit, int delimiter_int, int start_item_n)
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


        if (!q_f_buf.empty()) ///clear a constant reference
        {
            q_decompress_buf = decompress_deflate(q_f_buf);
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
                        fut_struct[cy1].n_fut=async(launch::async, jsd_distance, load_path_vector[c_it], std::cref(q_decompress_buf), delimiter_int);

                        fut_struct[cy1].in_act=true;
                        input_occur=true;

                        rev_pos++; ///reserve next thread number

                        break;

                    }

                }

            }

        }

        //vary q_f_buf valid-time can cause unwanted consequences, so the position moved to the current position.
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
    cout << "Parameter usage, [option][load_paths]\n";
    cout << "-h, show help\n";
    cout << "-t [int], set number of threads (default=5)\n";
    cout << "-r [path], input reserved distance matrix\n";

    cout << "-s, output a symmetric matrix; default is a low triangular matrix\n";
    cout << "-T, use TAB as a row name separator, instead of PHYLIP format that limit row names up to the first 9 characters\n";

    cout << "-d, distance metric convert: JS-Divergence to JS-Distance; square-root(JS-divergence) = JS distance\n";
}


void show_profile()
{
    cout << "FFP distance calculate; 2v.3.2\n";
    cout << "Code by JaeJin Choi; https://github.com/jaejinchoi/FFP\n";
    cout << "Value presentation: a point below 8 decimal places (%.8g)\n";

    //cout << "Compile; g++ -std=c++11 -pthread -o (output) (this script) -lz\n";
    //cout << "Required; zlib 1.2.8+\n";
}


int main(int argc, char** argv)
{

    int opt;
    int delimiter_int = 10;
    int thread_n_limit=5; ///default, 5 threads

    bool symmetric_flag=false; ///default, output a low triangular matrix
    bool js_distance_flag=false; ///default, output JS divergence
    bool item_tab_flag=false; //default is false: PHYLIP format that limit row names by 9 characters; true: use tab as a separate to use any lengths of row names

    string reserve_matrix_path="";

    while ((opt = getopt(argc, argv, "ht:r:dvsT")) !=EOF) // EOF = -1
    {
        switch (opt)
        {
            case 'h':
                show_help();
                exit(EXIT_SUCCESS);

            case 't':
                thread_n_limit=max(1, atoi(optarg)); ///least 1 thread or not run
                break;

            case 'r':
                reserve_matrix_path=optarg;
                break;

            case 'd':
                js_distance_flag=true;
                break;

            case 's':
                symmetric_flag=true;
                break;

            case 'T':
                item_tab_flag=true; //enable TAB as an item name separator; disable PHYLIP format
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
            multi_thread_manage(fut_value_vector, load_path_vector, thread_n_limit, delimiter_int, start_item_n); ///multi-threading; calculate JS Divergence

            print_value_vector_str(output_stream, fut_value_vector, load_path_vector, start_item_n, js_distance_flag, symmetric_flag, item_tab_flag); ///final output (standard output)

        }

        output_stream.str(string());
        output_stream.clear();

    } else
    {
        cout << "Input items should more than one" << endl;

    }

}

