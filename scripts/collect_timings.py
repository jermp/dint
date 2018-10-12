import sys, os

path_to_basename = sys.argv[1] # e.g., '/data2/inverted_indexes/gov2/gov2.sorted-text'
path_to_binaries = sys.argv[2] # e.g., './bin'
path_to_results = sys.argv[3]  # e.g., './results'
prefix_name = sys.argv[4]      # e.g., 'gov2'
query_log = sys.argv[5]

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

for type in strategies:
    cmd = "./queries " + type + " and " + path_to_binaries + "/" + prefix_name + "." + type + ".bin < " + query_log + " >> " + path_to_results + "/" + prefix_name + "." + type + ".querytime"
    # print cmd
    for i in xrange(0, 3):
        os.system(cmd)
            