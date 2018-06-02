import sys, os

collection_basename = sys.argv[1]
ouput_prefix_name = sys.argv[2]

strategies = ["DSF", "PDF", "LSS"]

for strategy in strategies:
    cmd = "./create_freq_index " + strategy + "_block_dint " + collection_basename + " " + ouput_prefix_name + "." + strategy + ".out"
    os.system("cmd")
