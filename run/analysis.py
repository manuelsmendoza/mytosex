""" Analysis
"""
import os
import sys
import tensorflow as tf
from src.func_analysis import *
from src.func_setup import load_settings
from src.func_util import tnow, pass_file
from tensorflow import keras

# Load settings
settings = load_settings(os.getenv("MYTOSEX_SETTINGS"))
tmp_dir = os.path.join(settings["output_dir"], "tmp")
data_dir = os.path.join(settings["output_dir"], "data")
figs_dir = os.path.join(settings["output_dir"], "figures")

# Building the references of
print(tnow() + " INFO: Merging both mitogenomes of reference", file=sys.stdout)
for ext in [".fasta", ".gff", ".bed"]:
    open(os.path.join(tmp_dir, settings["reference"]["alias"] + ext), "w").close()
    for mt_type in ["mtf", "mtm"]:
        if settings["reference"][mt_type + "_ncbi"]:
            pass_file(
                os.path.join(tmp_dir, settings["reference"]["alias"] + "_" + mt_type + ext),
                os.path.join(tmp_dir, settings["reference"]["alias"] + ext)
            )
        else:
            pass_file(
                os.path.join(settings["reference"][mt_type]),
                os.path.join(tmp_dir, settings["reference"]["alias"] + ext)
            )

# Build the sequence index
# print(tnow() + " INFO: Building the reference index", file=sys.stdout)
# build_index(
#     os.path.join(tmp_dir, settings["reference"]["alias"] + ".fasta"),
#     os.path.join(tmp_dir, settings["reference"]["alias"]),
#     settings["numb_threads"]
# )

# Align the samples
# for sample in list(settings["samples"].keys()):
#     print(tnow() + " INFO: Aligning the reads of " + settings["samples"][sample]["alias"], file=sys.stdout)
#     if settings["samples"][sample]["layout"] == "paired" and settings["samples"][sample]["ncbi"]:
#         reads_align(
#             index=os.path.join(tmp_dir, settings["reference"]["alias"]),
#             alias=settings["samples"][sample]["alias"],
#             outdir=tmp_dir,
#             threads=settings["numb_threads"],
#             freads=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + "_1.fastq.gz"),
#             rreads=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + "_2.fastq.gz")
#         )
#     elif settings["samples"][sample]["layout"] == "paired" and not settings["samples"][sample]["ncbi"]:
#         reads_align(
#             index=os.path.join(tmp_dir, settings["reference"]["alias"]),
#             alias=settings["samples"][sample]["alias"],
#             outdir=tmp_dir,
#             threads=settings["numb_threads"],
#             freads=settings["samples"][sample]["forward"],
#             rreads=settings["samples"][sample]["reverse"]
#         )
#     elif settings["samples"][sample]["layout"] == "single" and settings["samples"][sample]["ncbi"]:
#         reads_align(
#             index=os.path.join(tmp_dir, settings["reference"]["alias"]),
#             alias=settings["samples"][sample]["alias"],
#             outdir=tmp_dir,
#             threads=settings["numb_threads"],
#             sreads=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".fastq.gz")
#         )
#     elif settings["samples"][sample]["layout"] == "single" and not settings["samples"][sample]["ncbi"]:
#         reads_align(
#             index=os.path.join(tmp_dir, settings["reference"]["alias"]),
#             alias=settings["samples"][sample]["alias"],
#             outdir=tmp_dir,
#             threads=settings["numb_threads"],
#             sreads=settings["samples"][sample]["single"]
#         )

# Filter the alignments and perform sex prediction (also extract more stats)
metrics_list = []
results_list = []
for sample in list(settings["samples"].keys()):
    # print(tnow() + " INFO: Filtering the alignments of " + settings["samples"][sample]["alias"], file=sys.stdout)
    # if settings["samples"][sample]["layout"] == "paired":
    #     filter_alignment(
    #         alignment=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".sam"),
    #         threads=settings["numb_threads"],
    #         layout="paired"
    #     )
    # elif settings["samples"][sample]["layout"] == "single":
    #     filter_alignment(
    #         alignment=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".sam"),
    #         threads=settings["numb_threads"],
    #         layout="single"
    #     )

#     print(tnow() + " INFO: Extracting alignment statistics of " + settings["samples"][sample]["alias"], file=sys.stdout)
#     extract_stats(
#         alignment=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".markdup.bam"),
#         features=os.path.join(tmp_dir, settings["reference"]["alias"] + ".bed"),
#         threads=settings["numb_threads"]
#     )
# 
#     sample_stats = alignment_stats(
#             coverage=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".cov.tsv"),
#             feat_coverage=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".bedcov.tsv"),
#             depth=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".depth.tsv"),
#             salias=settings["samples"][sample]["alias"],
#             ralias=settings["reference"]["alias"]
#         )
#     metrics_list.append(sample_stats)
# 
# print(tnow() + " INFO: Exporting alignment statistics", file=sys.stdout)
# align_metrics = pd.concat(metrics_list)
# align_metrics.to_csv(
#      os.path.join(data_dir, "align_stats.tsv"),
#      sep="\t",
#      index=False
# )
print(os.path.abspath(__file__))
print(tnow() + " INFO: Inferring the sex of the samples", file=sys.stdout)
model = tf.keras.models.load_model("model/med")
sex_prediction = model.predict(align_metrics.loc[:, ["mtfcov", "mtmcov", "mtfmd", "mtmmd", "mtfgi", "mtmgi"]])
sex_prediction = np.array([x[0] for x in np.array(sex_prediction).round()], dtype="int")

sex_result = pd.DataFrame.from_dict(
    {
        "sample": np.array(align_metrics.loc[:, "sample"], dtype="str"),
        "sex": sex_prediction
    }
)
sex_result["sex"].replace({0: "Female", 1: "Male"}, inplace=True)
sex_result.to_csv(
    os.path.join(settings["output_dir"], "sex_prediction.tsv"),
    sep="\t",
    index=False
)

print(sex_result)
