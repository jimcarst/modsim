import matplotlib.pyplot as pt
import sys

ftype = "conf"
data = []
dir = "configurations/"
speed = 100000000

##Initialisation
file = open(dir + ftype + "_step0.output", "r")
for line in file:
    data.append([])
    line_list = list(line.split(" "))
    for l in line_list:
        data[-1].append(float(l))
file.close()

#Plotting of the whole model
def plotModel(input_data):
    xlen, ylen = len(input_data), len(input_data[0])
    spin_up_x, spin_up_y, spin_down_x, spin_down_y = [], [], [], []
    for x in range(xlen):
        for y in range(ylen):
            if input_data[x][y] == 1:
                spin_up_x.append(x)
                spin_up_y.append(y)
            elif input_data[x][y] == -1:
                spin_down_x.append(x)
                spin_down_y.append(y)
            else:
                print(input_data[x][y])


def readoutput(filename):
    file = open(filename, "r")
    ld = 0
    for line in file:
        dt = list(line.split(" "))
        for l in range(len(dt)):
            data[ld][l] = float(dt[l])
        ld += 1
    file.close()
    
    
pt.hot()
fig = pt.imshow(data,interpolation='none',vmin=-1.0,vmax=1.0)

tot_steps = 300000
step = 100
for i in range(1010, tot_steps, step):
    try:
        readoutput(dir + ftype + "_step{:d}.output".format(i))
    except:
        print("Error: could not open data file: " + dir + ftype + "_step{:d}.output".format(i))
        exit()
    plotModel(data)
    fig.set_data(data)
    pt.pause(1.E-4)
    pt.draw()
input()