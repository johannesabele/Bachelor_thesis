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
def read_config(mode='chaos',file='config.ini'):
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
mech_reader = csv.reader(open("mechanical.csv", "r"), delimiter=",")

bar = progressbar.ProgressBar(max_value=N, redirect_stdout=True)

i = 0
k = 0
for mech_line in mech_reader:
    if not mech_line:
        i += 1
        k = 0
        bar.update(i)
    else:
        arr_mech = np.array(list(mech_line))
        data[0,i,k] = arr_mech
        k += 1
        
bar.finish()


#####################  plotting  ########################

import matplotlib.pyplot as plt
 

print("Start Plotting...")

for i in range(N):
    plt.imshow(data[0,i], origin="lower", interpolation="bilinear", cmap=plt.cm.RdBu, vmin=-0.2, vmax=0.2)  
    plt.axis('off')  
    if i == 0: plt.colorbar(orientation="horizontal", shrink=0.52, pad=0.05)
    plt.savefig(f'../Results/mechanical/mechanical_{i}.png')  
    if (i%10) == 0: print(f"{i} pictures ready")
    if ((i+1)%N) == 0 and (i != 0): print("Yeeah, all pictures ready")


