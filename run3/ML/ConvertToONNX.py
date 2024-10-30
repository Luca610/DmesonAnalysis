import os
import sys
import argparse
import pickle
import yaml
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml_converter.h4ml_converter import H4MLConverter

def main(): #pylint: disable=too-many-statements
    # read config file
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileNameML.yml', help='config file name for ml')
    args = parser.parse_args()

    print('Loading analysis configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading analysis configuration: Done!')

    PtBins = [[a, b] for a, b in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'])]

    output_dir = os.path.expanduser(inputCfg['conversion']['output_folder'])
    if os.path.isdir(output_dir):
            print((f'\033[93mWARNING: Output directory \'{output_dir}\' already exists,'
                   ' overwrites possibly ongoing!\033[0m'))
    else:
        os.makedirs(output_dir)

    for iBin, PtBin in enumerate(PtBins):
        print(f'\n\033[94mConverting model --- {PtBin[0]} < pT < {PtBin[1]} GeV/c\033[0m')
        ModelList = inputCfg['ml']['saved_models']
        ModelPath = ModelList[iBin]
        if not isinstance(ModelPath, str):
            print('\033[91mERROR: path to model not correctly defined!\033[0m')
            sys.exit()
        ModelPath = os.path.expanduser(ModelPath)
        print(f'Loaded saved model: {ModelPath}')
        ModelHandl = ModelHandler()
        ModelHandl.load_model_handler(ModelPath)

        model_converter = H4MLConverter(ModelHandl) # create the converter object
        model_onnx = model_converter.convert_model_onnx(1, len(inputCfg['ml']['training_columns']))
        model_converter.dump_model_onnx(f"{output_dir}/ModelHandler_pT_{PtBin[0]}_{PtBin[1]}.onnx") # dump the model in ONNX format

main()