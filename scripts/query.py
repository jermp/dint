import sys, os

bins_directory = sys.argv[1]
output_prefix_name = sys.argv[2]
query_log = sys.argv[3]
output = sys.argv[4]

strategies = [
              "opt",
              "block_interpolative",
              "block_qmx",
              # "block_simple16",
              "block_optpfor",
              # "block_vbyte",
              "block_varintgb",
              "block_varintg8iu",
              "block_streamvbyte",
              "block_maskedvbyte"
              ]

for type in strategies:
    # perform queries
    index = output_prefix_name + "." + type + ".bin"
    cmd = "./queries " + type + " and " + bins_directory + "/" + index + " < " + query_log + " >> " + output
    for i in xrange(0, 1):
        os.system(cmd)