import numpy as np
import math
import cython
import pandas as pd
import multiprocessing as mp
import argparse
import deepdish
import matplotlib.pyplot as plt

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
PWM_SCHEMA = ["MOTIF_NAME", "MOTIF_PWM_FWD", "1", "2", "3", "4", "5", "6", "7", "8"]

def fetch_footprinting_args():
    parser=argparse.ArgumentParser(description="get marginal footprint scores")
    parser.add_argument("-tt", "--title", type=str, required=True, help="title for plot")
    parser.add_argument("-pwm_f", "--motifs_to_pwm", type=str, required=True)
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")
    args = parser.parse_args()
    return args

def get_footprint_characteristics(motif_width, corrected_footprint_1, nomotif_corrected_footprint_1, title, output_prefix, motif_name):


    fc_change = corrected_footprint_1 -  nomotif_corrected_footprint_1
    smoothed_fc = np.convolve(fc_change, np.ones(5)/5, mode='same')
    grad = np.gradient(smoothed_fc)

    rs = 500-motif_width//2-50
    ls = 500+motif_width//2+50

    midpoint=motif_width//2+50
    ldx = midpoint+np.argmax(grad[500:ls])
    rdx = np.argmin(grad[rs:500])
    minv = np.min(smoothed_fc[rs:ls][rdx:ldx])

    offset_right = np.min((grad[rs:ls][0:rdx][::-1]>0).nonzero())
    offset_left = np.min((grad[rs:ls][ldx:]<0).nonzero())

    right_max = rdx - offset_right
    left_max = ldx + offset_left

    maxv = np.max([smoothed_fc[rs:ls][right_max], smoothed_fc[rs:ls][left_max]])

    scores_path = "{}_{}_footprint_score.txt".format(output_prefix, motif_name)

    #f =  open(scores_path, "w")
    text = ",".join([motif_name, str(np.round(maxv-minv,2)), str(left_max-right_max), str(left_max-midpoint), str(midpoint-right_max), str(motif_width)])
    #f.write(text)
    #f.close()

    plt.figure()   
    plt.plot(smoothed_fc[rs:ls])
    vertical_lin_1 = [np.min(smoothed_fc[rs:ls])]*(ls-rs)
    vertical_lin_1[right_max] =np.max(smoothed_fc[rs:ls])
    plt.plot(vertical_lin_1, c="red", alpha=0.5)
    vertical_lin_2 = [np.min(smoothed_fc[rs:ls])]*(ls-rs)
    vertical_lin_2[left_max] =np.max(smoothed_fc[rs:ls])
    plt.plot(vertical_lin_2, c="red", alpha=0.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel('Counts',fontsize=18)
    plt.xlabel('motif_wdith+100bp bases around motif insertion',fontsize=18)
    plt.title(title+"_"+motif_name+"_"+text, fontsize=18)
    #plt.savefig(output_prefix+".{}.corrected.footprint.png".format(motif_name))
    plt.savefig(output_prefix+".{}.uncorrected.footprint.png".format(motif_name))

    return text


def get_footprint_characteristics_new(motif_width, corrected_footprint_1, nomotif_corrected_footprint_1, title, output_prefix, motif_name):


    fc_change = corrected_footprint_1 -  nomotif_corrected_footprint_1
    smoothed_fc = np.convolve(fc_change, np.ones(5)/5, mode='same')


    rs = 500-motif_width//2-50
    ls = 500+motif_width//2+50

    midpoint=motif_width//2+50
    ldx = midpoint+np.argmax(smoothed_fc[500:ls])
    rdx = np.argmax(smoothed_fc[rs:500])
    minv = np.min(smoothed_fc[rs:ls][rdx:ldx])

    right_max = rdx 
    left_max = ldx 

    maxv = np.max([smoothed_fc[rs:ls][right_max], smoothed_fc[rs:ls][left_max]])

    scores_path = "{}_{}_footprint_score.txt".format(output_prefix, motif_name)

    #f =  open(scores_path, "w")
    text = ",".join([motif_name, str(np.round(maxv-minv,2)), str(left_max-right_max), str(left_max-midpoint), str(midpoint-right_max), str(motif_width)])
    #f.write(text)
    #f.close()

    plt.figure()   
    plt.plot(smoothed_fc[rs:ls])
    vertical_lin_1 = [np.min(smoothed_fc[rs:ls])]*(ls-rs)
    vertical_lin_1[right_max] =np.max(smoothed_fc[rs:ls])
    plt.plot(vertical_lin_1, c="red", alpha=0.5)
    vertical_lin_2 = [np.min(smoothed_fc[rs:ls])]*(ls-rs)
    vertical_lin_2[left_max] =np.max(smoothed_fc[rs:ls])
    plt.plot(vertical_lin_2, c="red", alpha=0.5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel('Counts',fontsize=18)
    plt.xlabel('motif_wdith+100bp bases around motif insertion',fontsize=18)
    plt.title(title+"_"+motif_name+"_"+text, fontsize=18)
    plt.savefig(output_prefix+".{}.corrected.new1.footprint.png".format(motif_name))
    #plt.savefig(output_prefix+".{}.uncorrected.footprint.png".format(motif_name))

    return text


def main():
    args = fetch_footprinting_args()

    pwm_df = pd.read_csv(args.motifs_to_pwm, sep='\t',names=PWM_SCHEMA)

    motif_widths={}
    for i,r in pwm_df.iterrows():
        motif_widths[r["MOTIF_NAME"]] = len(r["MOTIF_PWM_FWD"])

    #scores_path = "{}_uncorrected_footprints.h5".format(args.output_prefix)
    scores_path = "{}_footprints.h5".format(args.output_prefix)
    scores = deepdish.io.load(scores_path)
    control_fotprint = scores["control"][0]

    values = []
    for motif in  scores.keys():
        if motif=="control":
            continue
        motif_footprint = scores[motif][0]
        width = motif_widths[motif]
        scores_text = get_footprint_characteristics_new(width, motif_footprint, control_fotprint, args.title, args.output_prefix, motif)
        values.append(scores_text.split(","))
    
    df = pd.DataFrame(values, columns = ['MOTIF_NAME', 'footprint_score', "footprint_width", "distance_from_center_left",  "distance_from_center_right", "motif_width"])
    df.to_csv(args.output_prefix+"_all_new_scores.csv", index=False)
    #df.to_csv(args.output_prefix+"_uncorrected_all_scores.csv", index=False)



if __name__ == '__main__':
    main()

