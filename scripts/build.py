import sys, os

collection_basename = sys.argv[1]
output_prefix_name = sys.argv[2]
bins_directory = sys.argv[3]

strategies = [
              # "opt",
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

for type in strategies:
    # build index
    output = output_prefix_name + "." + type
    cmd = "./create_freq_index " + type + " " + collection_basename + " " + bins_directory + "/" + output + ".bin"
    os.system(cmd)