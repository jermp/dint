import sys, os

collection_basename = sys.argv[1]
output_prefix_name = sys.argv[2]
query_log1 = sys.argv[3]
query_log2 = sys.argv[4]

strategies = ["opt",
              "block_interpolative",
              "block_qmx",
              "block_simple16",
              "block_optpfor",
              "block_vbyte",
              "block_varintgb",
              "block_varintg8iu",
              "streamvbyte_block",
              "maskedvbyte_block"]

results_directory = "./results." + output_prefix_name
if not os.path.exists(results_directory):
    os.makedirs(results_directory)

bins_directory = "./bins." + output_prefix_name
if not os.path.exists(bins_directory):
    os.makedirs(bins_directory)

for type in strategies:
    # build index
    output = output_prefix_name + "." + type
    cmd = "./create_freq_index " + type + " " + collection_basename + " " + bins_directory + "/" + output + ".bin"
    cmd += " 2> " + results_directory + "/" + output + ".log"
    # print cmd
    os.system(cmd)

    # perform queries
    cmd = "./queries " + type + " and " + bins_directory + "/" + output + ".bin < " + query_log1 + " >> " + results_directory + "/" + output + ".querytime"
    for i in xrange(0, 4):
        os.system(cmd)
    cmd = "./queries " + type + " and " + bins_directory + "/" + output + ".bin < " + query_log2 + " >> " + results_directory + "/" + output + ".querytime"
    for i in xrange(0, 4):
        os.system(cmd)
