import sys, os

collection_basename = sys.argv[1]
ouput_prefix_name = sys.argv[2]

strategies = ["block_varintgb",
              "block_varintg8iu",
              "block_vbyte",
              "block_qmx",
              "block_simple16",
              "block_optpfor",
              "opt",
              "block_interpolative"]

for strategy in strategies:
    output = ouput_prefix_name + "." + strategy
    cmd = "./create_freq_index " + strategy + " " + collection_basename
    cmd += " 2> ./competitors/" + output + ".log"
    # print cmd
    os.system(cmd)
