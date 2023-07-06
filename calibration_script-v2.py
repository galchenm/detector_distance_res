import os
import sys

import numpy as np 
import glob

import re
import argparse
import itertools

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_cmdline_args():
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter)
    parser.add_argument('-p','--p', type=str, help="The path of folder/s that contain/s GXPARM file")
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
                res1, I1 = float(l[0]), float(l[2])
                p_res1, p_I1 = 0., 0.

                while I1 >= 1. and i_start < index2:
                    if I1 == 1.:
                        return res1

                    p_res1, p_I1 = res1, I1
                    i_start += 1
                    l = re.sub(r"\s+", "**", lines[i_start].strip('\n')).split('**')[1:]
                    try:
                        res1, I1 = float(l[0]), float(l[2])
                    except ValueError:
                        return p_res1
                
                k = round((I1 - p_I1)/(res1 - p_res1),3)
                b = round((res1*p_I1-p_res1*I1)/(res1-p_res1),3)
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
        num, a, b, c, alpha, betta, gamma = line.split()
        a = float(a)
        b = float(b)
        c = float(c)
        
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
    return X, Y, dX, dY, detector_distance, a, b, c

def reading_XDS(XDS_INP):
    with open(XDS_INP, 'r') as file:
        distance = 0
        wavelength = 0   
        res_x = 0
        res_y = 0        
        for line in file:
            if not(line.startswith('!')) and line.strip().startswith('DETECTOR_DISTANCE='):
                distance = float(line.split('=')[1].strip())
                
            if not(line.startswith('!')) and line.strip().startswith('X-RAY_WAVELENGTH='):
                wavelength = float(line.split('=')[1].strip())
                
            if not(line.startswith('!')) and line.strip().startswith('QX='):
                res_x = float(line.split('=')[1].strip())
                
            if not(line.startswith('!')) and line.strip().startswith('QY='):
                res_y = float(line.split('=')[1].strip())
                               
            if wavelength != 0 and distance != 0 and res_x != 0 and res_y != 0:
                return distance, wavelength, res_x, res_y

def calculation(file1, z1, lambda1, file2, z2, res_x, res_y, X_centre, Y_centre, k_vect): #, L_vect_a, L_vect_b, L_vect_c): #file1, file2 - GXPARM; z1, z2 - detector distance from XDS.INP
    
    if z2 < z1:
        file1, z1, file2, z2 = file2, z2, file1, z1 
    delta_L = z2 - z1
    
    X1, Y1, dX1, dY1, L1, a1, b1, c1 = DDR_main(file1)
    X2, Y2, dX2, dY2, L2, a2, b2, c2 = DDR_main(file2)
    
    X_detector_center1 = X1 + dX1
    Y_detector_center1 = Y1 + dY1
    
    X_detector_center2 = X2 + dX2
    Y_detector_center2 = Y2 + dY2    

    
    prefix = os.path.commonprefix([file1,file2])
    print(f'For {file1.replace(prefix,"")} the distance is {L1} and detector centre is ({X_detector_center1}, {Y_detector_center1}), cell parameters - ({a1}, {b1}, {c1})')
    print(f'For {file2.replace(prefix,"")} the distance is {L2} and detector centre is ({X_detector_center2}, {Y_detector_center2}), cell parameters - ({a2}, {b2}, {c2})')
    
    La = a2 * delta_L / (a1 * L2 / L1 - a2)
    Lb = b2 * delta_L / (b1 * L2 / L1 - b2)
    Lc = c2 * delta_L / (c1 * L2 / L1 - c2)
    
    L_arr = [La, Lb, Lc]
    print('L for all three cell parameters', L_arr)
    
    k0 = [(X_detector_center2 - X_detector_center1)*res_x, (Y_detector_center2 - Y_detector_center1)*res_y, delta_L]
    
    k0 /= np.abs(np.sqrt(k0[0] ** 2 +k0[1] ** 2 + k0[2] **2)) #(np.linalg.norm(k0)*np.sign(k0))
    print(f'k0 vector is ({k0[0]}, {k0[1]}, {k0[2]})')
    
    cxo = X_detector_center1 - L1 * k0[0] / res_x
    cyo = Y_detector_center1 - L1* k0[1] / res_y
    
    print(f'Corrected detector centre is ({cxo}, {cyo})')
    
    k0 /= lambda1
    
    print(f'k0 vector is ({k0[0]}, {k0[1]}, {k0[2]})')
    
    X_centre.append(cxo) 
    Y_centre.append(cyo)
    k_vect.append(np.array(k0))

    #L_vect_a.append(La) 
    #L_vect_b.append(Lb)
    #L_vect_c.append(Lc)
    
    print('--------------------------------------------------------------------------------\n')
    return file1, z1, file2, z2, X_centre, Y_centre, k_vect, La, Lb, Lc #L_vect_a, L_vect_b, L_vect_c
    

def processing(all_paths, treshold):
    
    global GXPARM_XDS
    global CORRECT_LP
    
    X_centre = []
    Y_centre = []
    k_vect = []
    #L_vect_a = [] 
    #L_vect_b = []
    #L_vect_c = []
    
    processing_res = {}
    
    files_with_detector_distance = {} # filename: distance

    for path in all_paths:
        if os.path.exists(os.path.join(path,CORRECT_LP)):
            #print(os.path.join(path,CORRECT_LP))
            resolution_data = get_resolution(os.path.join(path,CORRECT_LP))
            if resolution_data < 0:
                print('With CORRECT.LP file in {} is something wrong'.format(path))
            else:
                if treshold and resolution_data < treshold:
                    XDS_INP = os.path.join(path, 'XDS.INP')
                    if os.path.exists(XDS_INP) and os.path.exists(os.path.join(path,GXPARM_XDS)):
                        distance, wavelength, res_x, res_y =  reading_XDS(XDS_INP)
                        files_with_detector_distance[os.path.join(path,GXPARM_XDS)] = {'distance': distance, 'wavelength': wavelength, 'res_x': res_x, 'res_y': res_y}
                    else:
                        print('In current path {} one of the file XDS.INP/GXPARM.XDS does not exist'.format(path))                    
                elif treshold is None:
                    XDS_INP = os.path.join(path, 'XDS.INP')
                    if os.path.exists(XDS_INP) and os.path.exists(os.path.join(path,GXPARM_XDS)):
                        distance, wavelength, res_x, res_y = reading_XDS(XDS_INP)
                        files_with_detector_distance[os.path.join(path,GXPARM_XDS)] = {'distance': distance, 'wavelength': wavelength, 'res_x': res_x, 'res_y': res_y}
                    else:
                        print('In current path {} one of the file XDS.INP/GXPARM.XDS does not exist'.format(path))            
        else:
            print('In current path {} file with name {} does not exist'.format(path, CORRECT_LP))
    
    list_of_files = files_with_detector_distance.keys()
    if len(list_of_files) != 0:
        for pair in list(itertools.combinations(list_of_files,2)):
            file1, file2 = pair
            z1 = files_with_detector_distance[file1]['distance']
            z2 = files_with_detector_distance[file2]['distance']
            lambda1 = files_with_detector_distance[file1]['wavelength']
            res_x = files_with_detector_distance[file1]['res_x'] 
            res_y = files_with_detector_distance[file1]['res_y']
            
            if abs(z2-z1)/max(z1,z2) >= 0.1:
                #X_centre, Y_centre, k_vect, L_vect_a, L_vect_b, L_vect_c = calculation(file1, z1, lambda1, file2, z2, res_x, res_y, X_centre, Y_centre, k_vect,  L_vect_a, L_vect_b, L_vect_c)
                file1, z1, file2, z2, X_centre, Y_centre, k_vect, La, Lb, Lc = calculation(file1, z1, lambda1, file2, z2, res_x, res_y, X_centre, Y_centre, k_vect)
                #processing_res[z1] = {'z2': z2, 'La': La,'Lb': Lb,'Lc': Lc}
                processing_res[min(z1,z2)] = {'z2': max(z1,z2), 'La': La,'Lb': Lb,'Lc': Lc}
    else:
        print(f'Check all parameters you used. Maybe there is no dataset with treshold = {treshold}')
    

    if len(X_centre)>0:
        x_median = np.round(np.median(np.array(X_centre)),6)
        y_median = np.round(np.median(np.array(Y_centre)),6)

        x_mean = np.round(np.mean(np.array(X_centre)),6)
        y_mean = np.round(np.mean(np.array(Y_centre)),6)

        sigma_x = np.round(np.std(np.array(X_centre)),6)
        sigma_y = np.round(np.std(np.array(Y_centre)),6)
        
        
        #key_distances = processing_res.keys()
        #min_z1 = min(key_distances)
        #max_z2 = max(key_distances)
        
        #delta = max_z2 - min_z1
        
        #La_max = processing_res[max_z2]['La'] - delta
        #Lb_max = processing_res[max_z2]['Lb'] - delta
        #Lc_max = processing_res[max_z2]['Lc'] - delta
        
        #mean_La = np.mean((np.array([La_max, processing_res[min_z1]['La']])))
        #std_La = np.std((np.array([La_max, processing_res[min_z1]['La']])))
        
        #mean_Lb = np.mean((np.array([Lb_max, processing_res[min_z1]['Lb']])))
        #std_Lb = np.std((np.array([Lb_max, processing_res[min_z1]['Lb']])))
        
        #mean_Lc = np.mean((np.array([Lc_max, processing_res[min_z1]['Lc']])))
        #std_Lc = np.std((np.array([Lc_max, processing_res[min_z1]['Lc']])))
        
        
        k_mean = np.round(np.mean(np.array(k_vect),axis=0),6)
        k_median = np.round(np.median(np.array(k_vect),axis=0),6)
        sigma_k = np.round(np.std(np.array(k_vect),axis=0),6)
        
        
        print('Final result:\n')
        
        message = f'\nMedian detector centre = ({x_median}, {y_median})\nMean detector centre = ({x_mean}, {y_mean})\nSTD for the detector centre = ({sigma_x}, {sigma_y})'
        print(message)
        
        #message = f'\nMean detector distance (a) = {mean_La}\nSTD for the detector distance (a) = {std_La}'
        #print(message)

        #message = f'\nMean detector distance (b) = {mean_Lb}\nSTD for the detector distance (b) = {std_Lb}'
        #print(message)        

        #message = f'\nMean detector distance (c) = {mean_Lc}\nSTD for the detector distance (c) = {std_Lc}'
        #print(message)
        
        message = f'\nMedian beam direction = {k_median}\nMean beam direction = {k_mean}\nSTD = {sigma_k}'
        print(message)

        
    
if __name__ == "__main__":
    args = parse_cmdline_args()
    
    treshold = args.tr
    
    GXPARM_XDS = 'GXPARM.XDS'
    CORRECT_LP = 'CORRECT.LP'

    if args.p is not None:
        main_path = args.p
        # Initialize an empty list, where later we add all the paths to our files
        all_paths = []
        for path, dirs, all_files in os.walk(main_path):
            for file in all_files:
                if file == GXPARM_XDS:
                    all_paths.append(path)
    elif args.f is not None:
        all_paths = [os.path.dirname(line.strip()) for line in open(args.f, 'r').readlines() if len(line.strip()) > 0]

    for path in all_paths:
        if os.access(path, os.R_OK) == 'false':
            exit()
    processing(all_paths, treshold)
