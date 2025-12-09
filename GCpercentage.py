import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns
import matplotlib.lines as mlines

"""
Please note that this code was created for a specific context and 
cannot be easily generalized. Read the description of Main() to use it.
"""

def one_line_fasta():
	file = "xxx.fasta"
	file2 = "xxx_oneline.fa"
	with open(file, "r") as infile, open(file2, "w") as outfile:
		infile = infile.readlines()
		for index, line in enumerate(infile):
			if line.startswith(">") and index == 0:
				outfile.write(line)
			elif line.startswith(">") and index > 0:
				new_line = "\n"+line
				outfile.write(new_line)
			else:
				line = line.strip()
				outfile.write(line)

def GC_calculator(line):
	len_line = len(line)
	line = line.strip()
	A_T = 0
	G_C = 0
	for char in line:
		if char == "A" or char == "T":
			A_T += 1
		else:
			G_C += 1
	A_Tper = A_T*100/len_line
	G_Cper = G_C*100/len_line

	return str(A_Tper), str(G_Cper)



def GC_percent():
	file = "xxx_oneline.fa"
	file2 = "xxx_oneline_GCpercent.tsv"
	with open(file, "r") as infile, open(file2, "w") as outfile:
		infile = infile.readlines()
		for index, line in enumerate(infile):
			new_line = ""
			if index%2 == 0:
				new_line += line.strip()
				A_T, G_C = GC_calculator(infile[index+1])
				new_line += "\t"+G_C+"\n"
			outfile.write(new_line)

def mean():
	input_file = "xxx_oneline_GCpercent.tsv"
	output_file = "xxx_oneline_GCpercent.txt"

	df = pd.read_csv(input_file, sep="\t", header=None)

	df.columns = ["Genome", "GC", "Species", "Annotation"]

	df_stats = df.groupby(["Genome", "Species", "Annotation"], as_index=False)["GC"].agg(["mean", "std"])

	df_stats = df_stats.rename(columns={"mean": "GC_Mean", "std": "GC_std"})

	df_stats = df_stats.reset_index(drop=True)

	df_stats = df_stats.sort_values(by="GC_Mean", ascending=False)

	df_stats.to_csv(output_file, sep="\t", index=False)
	return(output_file)

def plot_mean(file):
	input_file = file
	df = pd.read_csv(input_file, sep="\t")
	
	df = df.sort_values(by='GC_Mean', ascending=False).reset_index(drop=True)
	
	shape_dict = {
		"EV-A clade 2": "^",
		"EV-A clade 1": "s",
		"EV-A clade 3": "v",
		"EV-A clade 4": "D",
		"EV-B Hs": "s",
		"EV-B NHP": "D",
		"EV-J": "o",
		"EV-N": "o",
		"EV-H": "o",
		"EV-L": "o",
		"EV-C NC": "^",
		"EV-C Cult": "s",
		"EV-C resp": "v",
		"EV-D": "o",
		"EV-M": "o"
	}
	
	color_dict = {
		"EV-A clade 2": "#3dff32",
		"EV-A clade 1": "#3dff32",
		"EV-A clade 3": "#3dff32",
		"EV-A clade 4": "#3dff32",
		"EV-B Hs": "#0000ff",
		"EV-B NHP": "#0000ff",
		"EV-J": "#fe6fcf",
		"EV-N": "#ff8001",
		"EV-H": "#68cb6e",
		"EV-L": "#830cff",
		"EV-C NC": "#ff211b",
		"EV-C Cult": "#ff211b",
		"EV-C resp": "#ff211b",
		"EV-D": "#66ccff",
		"EV-M": "#820a07"
	}
	
	sort_legend = [
		"EV-A clade 1",
		"EV-A clade 2",
		"EV-A clade 3",
		"EV-A clade 4",
		"EV-B Hs",
		"EV-B NHP",
		"EV-C NC",
		"EV-C Cult",
		"EV-C resp",
		"EV-D",
		"EV-H",
		"EV-J",
		"EV-L",
		"EV-N",
		"EV-M"
	]
	

	plt.figure(figsize=(12, 6))
	
	for i, row in df.iterrows():
		plt.errorbar(
			x=row["Genome"], y=row["GC_Mean"], yerr=row["GC_std"], 
			fmt=shape_dict.get(row["Annotation"], "o"),  # Default shape: circle
			color=color_dict.get(row["Annotation"], "black"),  # Default color: black
			markersize=8, capsize=5
		)
	
	
	plt.xticks(rotation=90, fontsize=10)
	plt.yticks(fontsize=12)
	plt.ylabel("GC Mean (%)", fontsize=14)
	plt.title("%GC distribution", fontsize=16)
	plt.grid(True, linestyle="--", alpha=0.5)
	
	legend_elements = []
	
	present_annotations = df['Annotation'].unique()
	
	
	for annotation in sort_legend:
		if annotation in present_annotations:
			legend_elements.append(
				mlines.Line2D([], [], 
					color=color_dict.get(annotation, "black"), 
					marker=shape_dict.get(annotation, "o"), 
					linestyle='None', 
					markersize=10, 
					label=annotation)
			)
	
	
	plt.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1, 0.5), 
	          fontsize=10, frameon=True)
	
	plt.tight_layout()
	plt.show()


def main():
	"""
	README:
	This code allows to calculate the GC% of sequences from a fasta file.
	To use it, first launch the one_line_fasta() and GC_percent() functions. 
	Indicate the name of your file by replacing the xxx in the functions one_line_fasta(), GC_percent() and mean().
	Change the output file xxx_oneline_GCpercent.tsv to create 4 columns without header (Genome GC Species Annotation).
	Then you can launch line 184 and 185 to plot the results.

	usage: python3 GCpercentage.py
	"""
	#one_line_fasta()
	#GC_percent()
	file = mean()
	plot_average = plot_mean(file)


if __name__ == '__main__':
	main()