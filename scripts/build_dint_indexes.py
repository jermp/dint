import sys, os

collection_basename = sys.argv[1]
ouput_prefix_name = sys.argv[2]

strategies = ["DSF", "PDF", "LSS"]

for strategy in strategies:
    output = ouput_prefix_name + "." + strategy
    cmd = "./create_freq_index " + strategy + "_block_dint " + collection_basename + " " + output + ".bin"
    cmd += " 2> " + output + ".log"
    os.system(cmd)
