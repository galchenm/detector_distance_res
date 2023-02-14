#!/usr/bin/env python3
# coding: utf8

"""

"""

import os
import sys

import numpy as np 
import glob

import re
import argparse
import logging

os.nice(0)

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_cmdline_args():
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter)
    parser.add_argument('-path','--path', type=str, help="The path of folder/s that contain/s GXPARM file")
    parser.add_argument('-f','--f', type=str, help="File with path to folders of interest")
    parser.add_argument('-tr',type=float, help="Treshold for resolution")
    return parser.parse_args()

def get_resolution(CORRECTLP):
    resolution = 0
    index1 = None
    index2 = None
    
    with open(CORRECTLP, 'r') as stream:
        lines = stream.readlines()
        for line in lines:
            
            if line.startswith(' RESOLUTION RANGE  I/Sigma  Chi^2  R-FACTOR  R-FACTOR  NUMBER ACCEPTED REJECTED\n'):
                index1 = lines.index(line)
            if line.startswith('   --------------------------------------------------------------------------\n'):
                index2 = lines.index(line)
                break
        
        if index1 is None or index2 is None:
            return -200000
        else:
            i_start = index1+3
            l = re.sub(r"\s+", "**", lines[i_start].strip('\n')).split('**')[1:]
            
            if len(l) == 9:
                res, I_over_sigma = float(l[0]), float(l[2])
                prev_res, prev_I_over_sigma = 0., 0.
                
                if I_over_sigma < 1:
                    return -6000
                
                while I_over_sigma >= 1. and i_start < index2:
                    if I_over_sigma == 1.:
                        return res

                    prev_res, prev_I_over_sigma = res, I_over_sigma
                    i_start += 1
                    l = re.sub(r"\s+", "**", lines[i_start].strip('\n')).split('**')[1:]
                    try:
                        res, I_over_sigma = float(l[0]), float(l[2])
                    except ValueError:
                        return prev_res
                
                
                
                k = round((I_over_sigma - prev_I_over_sigma)/(res - prev_res),3)
                b = round((res*prev_I_over_sigma - prev_res*I_over_sigma)/(res-prev_res),3)
                
                try:
                    resolution = round((1-b)/k,3)
                    
                except ZeroDivisionError:
                    return -1000
            else:
                print(f'Something wrong with data in {CORRECTLP}')
                return -5000

    return resolution

def DDR_main(GXPARMXDS):
    with open(GXPARMXDS, 'r') as fileflow1:
                            
        line = fileflow1.readline()
        line = fileflow1.readline()
        line = fileflow1.readline()
        k_vect = line.split()
        lambd = float(k_vect[0])
        kx,ky,kz = [float(i)*lambd for i in k_vect[1:]]
        
        line = fileflow1.readline()
        
        line = fileflow1.readline()
        line = fileflow1.readline()
        line = fileflow1.readline()
        line = fileflow1.readline()
        
        data_array = line.split()
        res_x = float(data_array[-2])
        res_y = float(data_array[-1])
        
        line = fileflow1.readline()
        X, Y, detector_distance = line.split()
        X = float(X)
        Y = float(Y)
        detector_distance = float(detector_distance)
        
        dX = kx * abs(detector_distance) / res_x
        dY = ky * abs(detector_distance) / res_y
    return X, Y, dX, dY, detector_distance, kz


def calculations(path, distances, delta_x, delta_y, dic_old):
    global FILE_NAME1
    global FILE_NAME
    
    logger = logging.getLogger('app')
    if os.path.exists(os.path.join(path,FILE_NAME1)):

        X, Y, dX, dY, detector_distance, kz = DDR_main(os.path.join(path,FILE_NAME1))

        distances.append(detector_distance)
        delta_x.append(dX)
        delta_y.append(dY)

        dic_old[path] = [np.array([X+dX, Y+dY]),np.array([dX,dY,kz])]
        logger.info(f'For {path} the distance is {detector_distance} and detector centre is ({X+dX}, {Y+dY})')
    else:
        print('In current path {} file with name {} does not exist'.format(path, FILE_NAME1))
    
    return distances, delta_x, delta_y, dic_old


def processing(all_paths, treshold):
    
    global FILE_NAME1
    global FILE_NAME
    

    dic_old = {}
    delta_x = []
    delta_y = []
    data_array = []

    distances = []
    processed_files = 0
    for path in all_paths:
        if os.path.exists(os.path.join(path,FILE_NAME)):
        
            resolution_data = get_resolution(os.path.join(path,FILE_NAME))
            print(f'{path} - {resolution_data}')
            if resolution_data < 0:
                print('With CORRECT.LP file in {} is something wrong'.format(path))
            else:
                if treshold and resolution_data < treshold:
                    distances, delta_x, delta_y, dic_old = calculations(path, distances, delta_x, delta_y, dic_old)
                    processed_files += 1                     
                elif treshold is None:
                    distances, delta_x, delta_y, dic_old = calculations(path, distances, delta_x, delta_y, dic_old)
        else:
            print('In current path {} file with name {} does not exist'.format(path, FILE_NAME))
        
    if len(delta_x)>0:
        dx_median = round(np.median(np.array(delta_x)),4)
        dy_median = round(np.median(np.array(delta_y)),4)

        sigma_dx = round(np.std(np.array(delta_x)),6)
        sigma_dy = round(np.std(np.array(delta_y)), 6)

        center_x = round(np.median([dic_old[i][0][0] for i in dic_old]),4)
        center_y = round(np.median([dic_old[i][0][1] for i in dic_old]),4)

        sigma_x = round(np.std([dic_old[i][0][0] for i in dic_old]),6)
        sigma_y = round(np.std([dic_old[i][0][1] for i in dic_old]),6)
        if treshold is not None:
            print('The number of processed files is {}'.format(processed_files))
        print('Detector center is ({}, {}),'.format(center_x,center_y))
        print('Standart deviation for center is ({}, {})'.format(sigma_x, sigma_y))
        print('Detector distance (mean) is {},'.format(round(np.mean(np.array(distances)),4)))
        print('Detector distance (median) is {},'.format(round(np.median(np.array(distances)),4)))
        print('Detector distance deviation is {},'.format(round(np.std(np.array(distances)),4)))
    else:
        if treshold is not None:
            print(f'No samples with res < {treshold}')
        else:
            print('Something wrong with data')
    
    return 0

if __name__ == "__main__":
    args = parse_cmdline_args()
    
    treshold = args.tr
    
    level = logging.INFO
    logger = logging.getLogger('app')
    logger.setLevel(level)
    log_file = 'file.log'
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    f_handler = logging.FileHandler(os.path.join(os.getcwd(), 'file.log'))
    logger.addHandler(f_handler)
    print("Log file is {}".format(os.path.join(os.getcwd(), 'file.log')))

    # !!!!Specify the full name of the correct file of interest to us!!!!
    FILE_NAME1 = 'GXPARM.XDS'
    FILE_NAME = 'CORRECT.LP'

    if args.path is not None:
        main_path = args.path
        
        # Initialize an empty list, where later we add all the paths to our files
        all_paths = []
        for path, dirs, all_files in os.walk(main_path):
            for file in all_files:
                if file == FILE_NAME1:
                    all_paths.append(path)
    elif args.f is not None:
        all_paths = [os.path.dirname(line.strip()) for line in open(args.f, 'r').readlines() if len(line.strip()) > 0]
        

    for path in all_paths:
        if os.access(path, os.R_OK) == 'false':
            exit()
    processing(all_paths, treshold)

    print('FINISH')