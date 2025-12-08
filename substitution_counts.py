import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import matplotlib.font_manager as fm

# Step 0: Define the Arial font path with sizes
arial_font = fm.FontProperties(fname="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf", size=16)
legend_font = fm.FontProperties(fname="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf", size=8)

# Set global font properties except for legend
plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16
})

# Step 1: Load the Excel file
df = pd.read_excel("substitution_counts.xlsx")  

# Step 2: Remove commas and convert numeric columns to integers
df["MR297"] = df["MR297"].astype(str).str.replace(",", "").astype(int)
df["ML-1"] = df["ML-1"].astype(str).str.replace(",", "").astype(int)

# Step 3: Plot the grouped bar chart
x = np.arange(len(df["Substitution"]))  # label locations
width = 0.20  # width of the bars

plt.figure(figsize=(12, 6))
ax = plt.gca()

# Force y-axis to use plain (non-scientific) format
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
ax.ticklabel_format(style='plain', axis='y')

# Plot bars with black outlines (edgecolor='black', linewidth=1)
plt.bar(x - width/2, df["MR297"], width, label="MR297", color="#d62728",
        edgecolor="black", linewidth=1)
plt.bar(x + width/2, df["ML-1"], width, label="ML-1", color="#66c2a5",
        edgecolor="black", linewidth=1)

# Step 4: Customize axes and appearance
plt.xticks(x, df["Substitution"], rotation=0, fontproperties=arial_font)
plt.yticks([0, 500000, 1000000, 1500000, 2000000, 2500000],fontproperties=arial_font)        
plt.ylabel("No. of SNPs", fontproperties=arial_font)
plt.title("Transition and Transversion Substitution Comparison", fontproperties=arial_font)

plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

# Step 5: Save and show the plot
plt.savefig("/mnt/d/Ubuntu_WGRS/ts_tv_comparison.png", dpi=300)
