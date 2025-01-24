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
#include <assert.h>

///Found zlib function here http://mail-archives.apache.org/mod_mbox/trafficserver-dev/201110.mbox/%3CCACJPjhYf=+br1W39vyazP=ix
///eQZ-4Gh9-U6TtiEdReG3S4ZZng@mail.gmail.com%3E
#define MOD_GZIP_ZLIB_WINDOWSIZE 15
#define MOD_GZIP_ZLIB_CFACTOR    9
#define MOD_GZIP_ZLIB_BSIZE      8096


#define BLOCK_SIZE 32768 //1024*32 = 32KB (random access block size = 32KB)
///larger BLOCK_SIZE takes more time (e.g., 64KB, 8MB, 16MB)


void decompress_fraction(z_stream &zs, int &ret, int &fflush, stringstream &feed_stream, string &outstring)
{
    char inbuffer[BLOCK_SIZE];
	char outbuffer[BLOCK_SIZE];

    zs.avail_in = feed_stream.read(inbuffer, BLOCK_SIZE).gcount();
    zs.next_in = reinterpret_cast<Bytef*>(inbuffer);

    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = BLOCK_SIZE;

        ret = inflate(&zs, fflush);
        assert(ret!=Z_STREAM_ERROR);

        outstring.append(outbuffer, BLOCK_SIZE - zs.avail_out);

    } while(zs.avail_out==0);

}


void decompress_record_block(stringstream &feed_stream
    , int &feed_status
    , z_stream &zs
    , string &outstring
    , unsigned long &bytes_per_feature
    , unsigned long &bytes_per_value
    , int &feature_length
    , string &n_key
    , double &n_value
)
{
	int ret, fflush;
    fflush = Z_NO_FLUSH; //Z_NO_FLUSH only for inflation/decompression process

	if (feed_status==0)
	{
    	memset(&zs, 0, sizeof(zs)); //initialize
		ret = inflateInit(&zs);

		feed_status=1; //initiation failed

	    if (ret != Z_OK)
		{
			cerr << "inflateInit failed while decompressing" << endl;
			feed_status=-1; //initiation failed
		}

	}



    if (outstring.size() <= bytes_per_feature + bytes_per_value) // should exclude following condition to run correctly (for zlib end flag)!, && !feed_stream.eof())
	{
		decompress_fraction(zs, ret, fflush, feed_stream, outstring);

	}


	if (bytes_per_feature==0 && bytes_per_value==0) ///read header information
    {
        ///convert string to defined numeric variables, using *reinterpret_cast
		bytes_per_feature = *reinterpret_cast<const unsigned long*>( outstring.substr(0, sizeof(bytes_per_feature)).data() ); //expand 127 limit to 255 limit (4 letters -> 800, 20 letters -> 320)
		bytes_per_value = *reinterpret_cast<const unsigned long*>( outstring.substr(sizeof(bytes_per_feature), sizeof(bytes_per_value)).data() ); //expand 127 limit to 255 limit (4 letters -> 800, 20 letters -> 320)
		feature_length = *reinterpret_cast<const int*>( outstring.substr(sizeof(bytes_per_feature) + sizeof(bytes_per_feature), sizeof(feature_length)).data() ); //expand 127 limit to 255 limit (4 letters -> 800, 20 letters -> 320)

        n_key.resize(bytes_per_feature);
		//cout << "bytes_per_feature | bytes_per_value | feature_length : " << bytes_per_feature << "\t" << bytes_per_value << "\t" << feature_length << endl;
		outstring.erase(0, sizeof(bytes_per_feature) + sizeof(bytes_per_feature) + sizeof(feature_length));
    }

    ///read key and value blocks and then return
    n_key = string(outstring.substr(0, bytes_per_feature));
    n_value = *reinterpret_cast<const double*>( outstring.substr(bytes_per_feature, bytes_per_value).data() );
    outstring.erase(0, bytes_per_feature + bytes_per_value);

	///*
	if (feed_stream.eof() && outstring.size()==0) //finished decompressing and reading
	{
		//cout << "zs.total_out: " << zs.total_out << endl;
		inflateEnd(&zs); //flush

		if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
			cerr << "Exception during zlib decompression: (" << ret << ") " << zs.msg << endl;
			outstring.clear(); //should be 0 at this point
			feed_status=-1; //inflation finished incorrectly

		} else
		{
			feed_status=2; //inflation finished correctly and fully consumed decompressed outstring

		}

	}
	//*/
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


double calculate_distance(string p_path, const string &q_f_buf, bool &in_read, int delimiter_int, int distance_type)
{
    ///p -> istringstream, q -> cref = const referencing, istringstream

    ///Q_case
    string q_key="";
    double q_value=0.0;
    unsigned long q_bytes_per_feature=0;
    unsigned long q_bytes_per_value=0;
    int q_feature_length=0;


    ///setup Q_case stringstream
    stringstream read_q_f(q_f_buf, ios::in|ios::out|ios::binary); //does this copy q_f_buf to stringstream?
    read_q_f.seekg(0, ios::beg);

    z_stream q_zs;
    int q_feed_status=0;
	string q_outstring; //decompressed string

    in_read = false; ///instructions after this point can be done independently or simultaneously


    ///P_case
    string p_key="";
    double p_value=0.0;
    unsigned long p_bytes_per_feature=0;
    unsigned long p_bytes_per_value=0;
    int p_feature_length=0;

    ///setup P_case stringstream
    ifstream read_f;
    read_f.open(p_path.c_str(), ios::in|ios::binary);
	stringstream read_p_f(ios::in|ios::out|ios::binary);
    read_p_f << read_f.rdbuf(); //takes much time but no copy and use less memory (and less complex)
	read_f.close();
    read_p_f.seekg(0, ios::beg);

	z_stream p_zs;
	int p_feed_status=0;
	string p_outstring; //decompressed strin


    ///JSD value related variables
    double Hm=0.0;
    double Hp=0.0;
    double Hq=0.0;

	double r_value=0.0; //default; return error value = -1.0, else 0 <= X <= 1

    do
    {
        if (p_value==0 && p_feed_status!=2) ///conditions to read P_case
        {
			decompress_record_block(read_p_f, p_feed_status, p_zs, p_outstring, p_bytes_per_feature, p_bytes_per_value, p_feature_length, p_key, p_value);

        }

        if (q_value==0 && q_feed_status!=2) ///conditions to read P_case
        {
			decompress_record_block(read_q_f, q_feed_status, q_zs, q_outstring, q_bytes_per_feature, q_bytes_per_value, q_feature_length, q_key, q_value);

        }


        ///anything wrong during inflation/decompression step, or comparing with different feature lengths, should output an error value (-1)
		if (p_feed_status==-1 || q_feed_status==-1 || (p_feature_length!=q_feature_length)) //fool proof
		{
            cerr << "zlib process error or different feature_lengths compared" << endl;
			r_value=-1.0; ///note, ant distance value is >=0, value = -1 indicates error/fail
			break;

		}


		///various distances (or dissimilarity)
		switch(distance_type)
		{
            case 0: //case 0 and 1 are using JSD_divergence metric (case 1 for JSD_distance)
            case 1:
                JSD_divergence(p_key, q_key, p_value, q_value, Hp, Hq, Hm);
                break;

            case 2:
                JACCARD_distance(p_key, q_key, p_value, q_value, Hp, Hq, Hm);
                break;
		}

    } while (p_outstring.size()!=0 || q_outstring.size()!=0 || p_value!=0 || q_value!=0);

    ///clear variables and containers
    read_p_f.clear(); //clear flag
    read_p_f.str(string());
	p_outstring.clear();

    read_q_f.clear(); //clear flag
    read_q_f.str(string());
	q_outstring.clear();
    ///q_f_buf is constant referencing and should not mutate or clear

	if (r_value==-1)
	{
    	return r_value;
	}


    switch(distance_type)
    {
        case 0: //JSD divergence
            r_value = 0.5*(Hm - Hp - Hq);
            break;

        case 1: //JSD distance  = sqrt(JSD_divergence)
            r_value = sqrt(0.5*(Hm - Hp - Hq));
            break;

        case 2: //Jaccard distance
            r_value = (Hp + Hq) / (Hp + Hq + Hm);
            break;

        default: //fool proof
            r_value = -1;
            break;
    }

	return r_value;

}



struct future_handle
{
    int n_row;
    int n_col;
    future<double> n_fut;
    bool in_act=false;
	bool in_read=false;

};

/// considering to remove the option [-s] symmetric output option
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
    output_stream.str(string()); //flush previous contents
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


    if (symmetric_flag==true)
    {
        to_square_matrix_output(output_stream); ///convert a low triangular matrix to a square matrix

    }

    cout << output_stream.str(); ///output a matrix

}


///calculate JS Divergence using multiple threads
void multi_thread_manage(vector< vector<double> > &fut_value_vector, vector<string> &load_path_vector, int thread_n_limit, int delimiter_int, int start_item_n, int distance_type)
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
    //string q_decompress_buf;
	int in_read_cnt=0;

    for (vector<vector<double>>::size_type r_it=start_item_n; r_it<load_path_vector.size(); ++r_it)
    {
        fut_value_vector[r_it].resize(r_it);

        ///shared across threads
        if (r_it!=0)
        {
            ///efficient method to copy a whole file content to memory (as string); http://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
            read_f.open(load_path_vector[r_it].c_str(), ios::in|ios::out|ios::binary);
            read_f.seekg(0, ios::end);
            q_f_buf.resize(read_f.tellg()); //reach max size.
            read_f.seekg(0, ios::beg);
            read_f.read(&q_f_buf[0], q_f_buf.size()); //overwrite
            read_f.close();

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

                for (int cy1=rev_pos; cy1!=thread_n_limit; ++cy1) //input work in order
                {
                    ///check running-finished task
                    if (fut_struct[cy1].in_act==true && fut_struct[cy1].n_fut.wait_for(std::chrono::milliseconds(0))==future_status::ready)
                    {
                        n_row=fut_struct[cy1].n_row;
                        n_col=fut_struct[cy1].n_col;
                        fut_value_vector[n_row][n_col]=fut_struct[cy1].n_fut.get(); //now future status become valid()

                        fut_struct[cy1].in_act=false;

                    }

                    //*
                    ///feed new task
                    if (fut_struct[cy1].in_act==false)
                    {
                        fut_struct[cy1].n_row=r_it;
                        fut_struct[cy1].n_col=c_it;

						fut_struct[cy1].in_read=true;

                        ///input is [feature][feature value] format
                        fut_struct[cy1].n_fut=async(launch::async, calculate_distance, load_path_vector[c_it], std::cref(q_f_buf), std::ref(fut_struct[cy1].in_read), delimiter_int, distance_type); //solved 2020-2-22

						fut_struct[cy1].in_act=true;
                        input_occur=true;

                        rev_pos++; ///reserve next thread number

                        break;

                    }
                    //*/
                }

            }

        }

        ////// check and retrieve any tasks in action, to go next row
		for (int cy1=0; cy1!=thread_n_limit; ++cy1)
		{
			if (fut_struct[cy1].in_act==true)
			{
				n_row=fut_struct[cy1].n_row;
				n_col=fut_struct[cy1].n_col;
				fut_value_vector[n_row][n_col]=fut_struct[cy1].n_fut.get();

				//cout << n_row << "\t" << n_col << "\t" << fut_value_vector[n_row][n_col] << endl;

				fut_struct[cy1].in_act=false;
				fut_struct[cy1].in_read=false;
			}

		}

    }

	q_f_buf.clear();

}


int read_reserved_matrix(stringstream &output_stream, vector<string> &load_path_vector, string &reserve_matrix_path, bool item_tab_flag)
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

        if (item_name_str.length()>9)
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
    cout << "-h, show help\n";

    cout << "-t [int], set number of threads (default = 5)\n";
    cout << "-r [path], input reserved distance matrix\n";

    cout << "-s, output a symmetric matrix; default is a low triangular matrix\n";
    cout << "-T, Disable PHYLIP matrix format limiting item names to the first 9 characters, and accept full item names using tab as a separator\n";

    cout << "-d [int], type of distance metric\n";
    cout << "\t0 : Jensen-Shannon Divergence (JS divergence), is default\n";
    cout << "\t1 : Jensen-Shannon Distance (JS distance = square_root(JS divergence))\n";
    cout << "\t2 : Jaccard Distance\n";
    //cout << "distance metric convert: JS-Divergence to JS-Distance; square-root(JS-divergence) = JS distance\n";

}


void show_profile()
{
    cout << "FFP_distance_calculate; 2v.4.1\n";
    cout << "Code by JaeJin Choi; https://github.com/jaejinchoi/FFP\n";
    //cout << "Value presentation: a poin below 8 decimal places (%.8g)\n";
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
    string reserve_matrix_path="";

    int distance_type=0; //default is 0:JSD divergence
    bool item_tab_flag=false; //default is false, limit item name length by 9 characters, in PHYLIP format. true accept full item names and use tab as a separator

    while ((opt = getopt(argc, argv, "ht:r:d:vsT")) !=EOF) // EOF = -1
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

            case 'r':
                reserve_matrix_path=optarg;
                break;

            case 'd':
                distance_type = atoi(optarg);
                break;

            case 's':
                symmetric_flag=true;
                break;

            case 'T':
                item_tab_flag=true; //enable TAB as an item name separator; disable PHYLIP format
                break;

            default:
                show_help();
                exit(0);
        }
    }

    ///check supported distance_type
    if (distance_type < 0 || distance_type > 2)
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

