import sys, os

collection_basename = sys.argv[1]
output_prefix_name = sys.argv[2]
query_logs_basename = sys.argv[3]

strategies = [
              "opt",
              "block_interpolative",
              "block_qmx",
              "block_simple16",
              "block_optpfor",
              "block_vbyte",
              "block_varintgb",
              "block_varintg8iu",
              "block_streamvbyte",
              "block_maskedvbyte"
              ]

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
    os.system(cmd)

    # perform queries
    for suffix in [".0.mapped.1k", ".0.mapped.selective", ".1.mapped.1k", ".1.mapped.selective"]:
      cmd = "./queries " + type + " and " + bins_directory + "/" + output + ".bin < " + query_logs_basename + suffix + " >> " + results_directory + "/" + output + ".querytime"
      for i in xrange(0, 4):
          os.system(cmd)

    os.system("rm " + bins_directory + "/" + output + ".bin")
