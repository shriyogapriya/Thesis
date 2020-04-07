
# coding: utf-8

# In[78]:


import re
import csv
import numpy
import matplotlib.pyplot as plt


def get_relative_values(array):
    temp = []
    for i in range(len(array)):
        x = i % 2
        t = float(array[i]) - float(originalIdx[x]) 
        if float(originalIdx[x])*t >0:
            t -= 1
        t /= (float(seqLen[x])-1)
        temp.append(t)
    return temp


input_file = csv.DictReader(open("dc.csv"))

pred = []
orig = []
relativePositions = []

count=0
for row in input_file:
    count+=1
    query_target = re.split(r':', row['query:target'])
    pred = re.split(r'[:&]|\s+', row['prediction'])
    orig = re.split(r'[:&]|\s+', row['original'])
    seqLen = re.split(r':', row['seqLen'])
    originalIdx = re.split(r':', row['originalIdxPos0'])
    
    if (len(query_target)<2):
        continue

    for i in range(len(orig)):
        if orig[i].upper() == 'NA':
            orig[i] =float("nan")

    
    print('-' * 100)
    print('query_target:', query_target)
    print('pred:', pred)
    print('orig:', orig)
    print('seqLen:', seqLen)
    print('originalIdx:', originalIdx)
    print('-' * 100)
    map_index = {}

    relativePositionsPred = get_relative_values(pred)
    relativePositionsOrig = get_relative_values(orig)

    print('relativePositionsPred:', relativePositionsPred)
    print('relativePositionsOrig:', relativePositionsOrig)

    max1 = int(originalIdx[0]) + int(seqLen[0])
    max2 = int(originalIdx[1]) + int(seqLen[1])
    if (int(originalIdx[0])* max1) > 0 :
        max1 -= 1
    if (int(originalIdx[1])* max2) > 0 :
        max2 -= 1

    fig, ax1 = plt.subplots(1, 1, figsize=(6, 2), dpi=150)

    ax1.set_ylim(1, 2)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels([str(max2), originalIdx[1]])

    ax2 = ax1.twiny()
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels([originalIdx[0], str(max1)])

    for i in range(0, len(relativePositionsPred), 4):
        boxCorners = numpy.matrix([
            [relativePositionsPred[i], 1 - relativePositionsPred[i + 1],
             1 - relativePositionsPred[i + 3], relativePositionsPred[i + 2]],
            [2, 1, 1, 2]
        ])
        box = plt.Polygon(boxCorners.transpose(), closed=True,
                          fill=True, color="orange", alpha=0.5, label="prediction")
        ax1.add_patch(box)

    for i in range(0, len(relativePositionsOrig), 4):
        boxCorners = numpy.matrix([
            [relativePositionsOrig[i], 1 - relativePositionsOrig[i + 1],
             1 - relativePositionsOrig[i + 3], relativePositionsOrig[i + 2]],
            [2, 1, 1, 2]
        ])
        box = plt.Polygon(boxCorners.transpose(), closed=True,
                          fill=True, color="blue", alpha=0.5, label="original")
        ax1.add_patch(box)

    ax1.set_xlabel(query_target[1])
    ax2.set_xlabel(query_target[0])
    ax1.yaxis.set_visible(False)
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(
        handles=[handles[0], handles[2]],
        labels=[labels[0], labels[2]]
    )
    #plt.show()
    plt.savefig('rri-comparison-'+str(count)+'.pdf',bbox_inches='tight')

