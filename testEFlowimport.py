# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Connor Pecar

import tensorflow
from tensorflow import keras
import energyflow
from energyflow.archs import PFN
import numpy as np

modelname = "pfn_testEpic_000-2_vecQele_nHFS2_500_bs10k_bestValLoss"
model = keras.models.load_model(modelname)

def eflowPredict(feat, globalfeat):
    feat = np.asarray(feat)
    feat = np.reshape(feat, (1,len(feat),7))
    globalfeat = np.asarray(globalfeat)
    globalfeat = np.reshape(globalfeat, (1,10))    
    pred = model([feat,globalfeat])
    pred = np.reshape(pred,(4))
    return pred
