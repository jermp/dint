import sys, os

dataset_filename = sys.argv[1]
output_prefix = sys.argv[2]

codecs = ["interpolative",
          "qmx",
          "optpfor",
          "simple16",
          "varintgb",
          "varintg8iu",
          "vbyte",
          "maskedvbyte",
          "streamvbyte"]

for codec in codecs:
    output = output_prefix + "." + codec + ".out"
    enc_cmd = "./encode " + codec + " " + dataset_filename + " --out " + output
    dec_cmd = "./decode " + codec + " " + output
    os.system(enc_cmd)
    for i in xrange(0, 3):
        os.system(dec_cmd)
    os.system("rm " + output)
