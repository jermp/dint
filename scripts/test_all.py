import sys, os

dataset_filename = sys.argv[1]
output_prefix = sys.argv[2]
dictionary_filename = sys.argv[3]

codecs = ["interpolative", "optpfor", "varintg8iu", "qmx", "vbyte", "u32",
          "simple16", "streamvbyte", "maskedvbyte", "varintgb", "dint"]

for codec in codecs:
    output = output_prefix + "." + codec + ".out"
    enc_cmd = "./encode " + codec + " " + dataset_filename + " --out " + output
    # dec_cmd = "./decode " + codec + " " + output
    if codec == "dint":
        enc_cmd += " " + "--dict " + dictionary_filename
        dec_cmd += " " + "--dict " + dictionary_filename
    os.system(enc_cmd)
    # for i in xrange(0, 3):
    #     os.system(dec_cmd)
    # os.system("rm " + output)
