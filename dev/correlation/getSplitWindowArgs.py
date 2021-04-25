# --scripts for get the ROI coordinate of heatmap --
import argparse as ap
import cv2
import numpy as np
import os
import csv


def findroi(img, ret=False):
    """[summary]

    Args:
        img ([numpy.ndarray]): Heatmap image for operation
        ret (bool, optional): Define whether to return the image drawn contours. Defaults to False.

    Returns:
        [coordinate_list]: [the contours value of ROI]
    """
    # RGB image convert to HSV image
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    # set thresholds to filter the heatmap region
    minval = np.array([1, 0, 46])
    maxval = np.array([254, 255, 255])
    mask = cv2.inRange(hsv, minval, maxval)
    # get the color region
    color = cv2.bitwise_and(img, img, mask=mask)
    # image gray and binary
    gray = cv2.cvtColor(color, cv2.COLOR_BGR2GRAY)
    gray[:, int(gray.shape[1] * 0.9):] = 0
    ret, binary = cv2.threshold(gray, 1, 255, cv2.THRESH_BINARY)
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT,
                                       (3, 3))  # kernel for image open computation
    # close computation for fill holes
    image_close = cv2.morphologyEx(
        binary, cv2.MORPH_CLOSE,
        kernel)
    cnts, _ = cv2.findContours(image_close, cv2.RETR_EXTERNAL,
                               cv2.CHAIN_APPROX_SIMPLE)
    img_show = img.copy()
    coordinate_list = []
    for c in cnts:
        (x, y, w, h) = cv2.boundingRect(c)
        img = cv2.rectangle(img_show,  (x, y), (x + w, y + h),
                            (0, 255, 0), 2) if ret else None
        coordinate_list.append([x, y, w, h])
    print('find roi completed...')
    return coordinate_list, img_show


def getcoordinate(coordinate):
    """[summary]

    Args:
        coordinate ([list]): list with raw coordinate of the ROI

    Returns:
        coor_info : contains the init coordinate and widths and heights of every gap
    """
    coor_arr = np.array(coordinate)
    #print(f'{len(coor_arr)} ROI')
    #print(coor_arr)
    coor_arr= np.flip(coor_arr,axis=0)
    _, idx = np.unique(coor_arr[:, 0], return_index=True)
    x_coor = coor_arr[:, 0][np.sort(idx)].tolist()
    _, idx = np.unique(coor_arr[:, 1], return_index=True)
    y_coor = coor_arr[:, 1][np.sort(idx)].tolist()
    columns = len(x_coor)
    rows = len(y_coor)
    #print(f'{rows}rows,{columns}columns')
    #print(coor_arr)
    w_coor  = [coor_arr[i][2] for i in range(columns)]
    h_coor = [coor_arr[i][3] for i in range(0, columns*rows, columns)]
    init = np.array(coor_arr[0][:2]).tolist()
    #coor_info = np.array([init, w_coor, h_coor], dtype=object)
    print(f'get the coordinate...\nthere are {len(w_coor)} colums, {len(h_coor)} rows')
    return init, w_coor, h_coor


def get_args():
    parser = ap.ArgumentParser()
    parser.add_argument('-i', "--input", help="Input image path",
                        required=True)
    parser.add_argument('-o', "--output", help="Output path to save the info_file, csv",
                        required=True)

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    img = cv2.imread(os.path.join(args.input,'correlationPlot.png'))
    coor_list, img_cnts = findroi(img)
    init, w_coor, h_coor = getcoordinate(coor_list)
    with open(os.path.join(args.output,'splitWinArgs.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(init)
        writer.writerow(w_coor)
        writer.writerow(h_coor)
        csvfile.close()
    print(
        f'''mission completed, file save in {os.path.join(args.output,'splitWinArgs.csv')}''')
