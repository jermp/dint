import sys, os

collection_basename = sys.argv[1]
ouput_prefix_name = sys.argv[2]

strategies = ["vbyte_block",
              "varint_G8IU_block",
              "optpfor_block",
              "interpolative_block",
              "qmx_block",
              "simple16_block",
              "opt",
              ]

for strategy in strategies:
    output = ouput_prefix_name + "." + strategy
    cmd = "./create_freq_index " + strategy + " " + collection_basename
    cmd += " 2> " + output + ".log"
    # print cmd
    os.system(cmd)
