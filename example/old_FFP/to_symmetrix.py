import os, sys, getopt
import StringIO

        
def read_assymetric_matrix(load_path=''):
    taxon_name_dict={}
    taxon_value_list=[] #double list

    read_f = open(load_path, 'r')

    taxon_cnt=0 #start from 1
    taxon_name=''

    #first row of a matrix is the number of items, so skip
    read_f.readline()

    for read_line in read_f:
        read_line = read_line.strip()

        if (read_line==''):
            continue

        else:
            #tab as a delimiter?, some tree programs limit item name to 10 characters (e.g, phylip neighbor)
            #read_str_list = read_line.split('\t') #[item name, value 1, value 2, -, -, ]
            #taxon_value_list.append(read_str_list[1:])
            
            taxon_name = read_line[:10].strip() #vast item used to be from JGI

            taxon_name_dict[taxon_cnt] = taxon_name
            taxon_cnt+=1
            
            taxon_value_list.append(read_line[10:].split('\t'))

                
    read_f.close()

    return taxon_name_dict, taxon_value_list



#convert asymmetric -> symmetric matrix
def reform_to_symmetric(taxon_value_list=[]):

    #list container read left->right(index 0->n) in order
    for cy1 in range(0, len(taxon_value_list)):
        
        for cy2 in range(0, len(taxon_value_list[cy1])):
            
            if (cy1==cy2):
                continue

            elif (cy1>cy2 and cy1 > len(taxon_value_list[cy2])-1):
                taxon_value_list[cy2].append(taxon_value_list[cy1][cy2])


def output_matrix(taxon_name_dict={}, taxon_value_list=[], save_path=''):

    str_file = StringIO.StringIO()

    str_file.write("%d\n" % (len(taxon_name_dict))) #is the number of items
    #print len(taxon_name_dict) #is the number of items

    for taxon_index, taxon_name in taxon_name_dict.iteritems():
        #print "%-10s\t%s" % (taxon_name, "\t".join(taxon_value_list[taxon_index]))
        str_file.write("%-10s\t%s\n" % (taxon_name, "\t".join(taxon_value_list[taxon_index])))

    if (save_path!=''):
        write_f = open(save_path, 'w')
        write_f.write(str_file.getvalue())
        write_f.close()

    else: #standard output
        print str_file.getvalue()

    str_file.close()
    

def show_help():
    print '''
    Convert asymmetric to symmetric distance matrix
    [program][load_path][save_path (optional or standard output)]
    -h, show_help()
    '''
    sys.exit()


if __name__=="__main__":

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h')

    except:
        sys.exit()

    for opt, arg in opts:

        if (opt=='-h'):
            show_help()
            
        else:
            show_help()

    taxon_name_dict={}
    taxon_value_list=[]

    load_path=''
    save_path=''

    if (len(args) > 0):

        load_path = os.path.abspath(args[0])

        if (len(args)==2):
            save_path = os.path.abspath(args[1])
            print "output will be saved to %s" % (save_path)
            
        taxon_name_dict, taxon_value_list = read_assymetric_matrix(load_path)

        reform_to_symmetric(taxon_value_list)

        output_matrix(taxon_name_dict, taxon_value_list, save_path)

    else:

        print "No input"

        

