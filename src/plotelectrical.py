import numpy as np
import csv
import configparser
import matplotlib.pyplot as plt
import progressbar

N = 0
size = 100
mode = "chaos"

"""
    Reads in necessary parameters from config.ini
"""
def read_config(mode="chaos",file="config.ini"):
    config = configparser.ConfigParser()
    config.read(file)

    global N, size
    size = int(config[mode]['size'])
    N_max = float(config[mode]['N_max'])
    N_output = float(config[mode]['N_output'])
    sample_rate = float(config[mode]['sample_rate'])
    N = int((N_max - N_output)/sample_rate)

print("Preparing data...")
read_config(mode=mode)


data = np.empty((1, N, size, size))
el_reader = csv.reader(open("electrical.csv", "r"), delimiter=",")

bar = progressbar.ProgressBar(max_value=N, redirect_stdout=True)

i = 0
k = 0
for el_line in el_reader:
    if not el_line:
        i += 1
        k = 0
        bar.update(i)
    else:
        arr_el = np.array(list(el_line))
        data[0,i,k] = arr_el
        k += 1
        
bar.finish()


################  plotting   ####################

import numpy as np
 
print("Start Plotting...")

for i in range(N):
        plt.imshow(data[0,i], cmap=plt.cm.Reds_r, interpolation="bilinear", origin="lower", vmin=0, vmax=1)
        plt.axis('off')  
        if i == 0: plt.colorbar(orientation="horizontal", shrink=0.52, pad=0.05)
        plt.savefig(f'../Results/electrical/electrical_{i}.png') 
        if (i%10) == 0: print(f"{i} pictures ready")
        if ((i+1)%N) == 0 and (i != 0): print("Yeeah, all pictures ready")


