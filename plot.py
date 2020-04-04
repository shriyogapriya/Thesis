
# coding: utf-8

# In[2]:


import re
import csv
import numpy
import matplotlib.pyplot as plt


def get_relative_values(array):
    temp = []
    for i in range(len(array)):
        if i > 0:
            if i % 2 == 1:
                temp.append(
                    (float(array[i]) - float(originalIdx[0])) / float(seqLen[0]))
            else:
                temp.append(
                    (float(array[i]) - float(originalIdx[1])) / float(seqLen[1]))
    return temp


input_file = csv.DictReader(open("dc.csv"))

pred = []
orig = []
relativePositions = []

for row in input_file:
    query_target = re.split(r'[: & ]\s*', row['query:target'])
    _pred = re.split(r'[: & ]\s*', row['prediction'])
    _orig = re.split(r'[: & ]\s*', row['original'])
    seqLen = re.split(r'[: & ]\s*', row['seqLen'])
    originalIdx = re.split(r'[: & ]\s*', row['originalIdxPos0'])

    pred = []
    orig = []

    for x, y in zip(_pred, _orig):
        if x.upper == 'NA' or y.upper() == 'NA':
            pred.append(0)
            orig.append(0)
        else:
            pred.append(x)
            orig.append(y)

    pred.insert(0, 0)
    orig.insert(0, 0)
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

    max1 = int(originalIdx[0]) + int(seqLen[0]) -            int(0 if int(originalIdx[0]) < 0 else 1)
    max2 = int(originalIdx[1]) + int(seqLen[1]) -            int(0 if int(originalIdx[1]) < 0 else 1)

    fig, ax1 = plt.subplots(1, 1, figsize=(6, 3), dpi=50)

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
    plt.show()

