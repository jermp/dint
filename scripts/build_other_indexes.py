import sys, os

collection_basename = sys.argv[1]
ouput_prefix_name = sys.argv[2]

strategies = ["block_vbyte", "block_varint", "block_optpfor", "block_interpolative", "block_qmx", "block_simple16", "opt"]

for strategy in strategies:
    output = ouput_prefix_name + "." + strategy
    cmd = "./create_freq_index " + strategy + " " + collection_basename
    cmd += " 2> ./competitors/" + output + ".log"
    # print cmd
    os.system(cmd)
