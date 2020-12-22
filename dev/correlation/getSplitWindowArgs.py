import argparse as ap
import numpy as np
## using command: pip install opencv-python
## may use this code to install module 'skbuild': pip install scikit-build
import cv2
import pandas as pd

if __name__ == "__main__":
    # Argument Parser
    parser = ap.ArgumentParser()
    parser.add_argument('-p', "--path", help="Path to image",
            required=True)
    parser.add_argument('-l', "--row", type = int, help="Rows of splited regions",
            required=True)
    parser.add_argument('-g', "--columns", type = int, help="Columns of splited regions",
            required=True)
    parser.add_argument('-o', "--outpath", help="Output path",
            required=True)
    args = vars(parser.parse_args())
    
    
    #把热图抠出来
    path = args["path"]
    row = args['row']
    columns = args['columns']
    output = args['outpath']

    hotmap = cv2.imread(path + "correlationPlot.png")
    hsv = cv2.cvtColor(hotmap,cv2.COLOR_BGR2HSV)
    minval = np.array([1,0,46])
    maxval = np.array([124,255,255])
    mask = cv2.inRange(hsv,minval,maxval)
    color = cv2.bitwise_and(hotmap,hotmap,mask=mask)

    #扣出区块
    color[color==0]=255
    gray = cv2.cvtColor(color,cv2.COLOR_BGR2GRAY)
    gray[:,int(gray.shape[1]*0.9):]=255
    ret,binary = cv2.threshold(gray,254,255,cv2.THRESH_BINARY)
    contours, hierarchy = cv2.findContours(binary,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    img = cv2.drawContours(hotmap.copy(), contours, -1, (0,255,0), 1)

    #按序号取出不同方块
    inf = []
    for i in range(len(contours)):
        x,y,w,h = cv2.boundingRect(contours[i])
        inf.append([x,y,w,h])
    inf.sort(key=lambda x:x[1])
    inf.sort(key=lambda y:y[0])
    inf.pop(0)
    r = [i for i in range(1,row+1)]*columns
    c = [[j]*row for j in range(1,columns+1)]
    c = [j for m in c for j in m]
    info = np.array(inf)
    
    #保存数据
    df = pd.DataFrame(info,columns=['x','y','w','h'])
    df['row'] = r
    df['columns'] = c
    df.to_csv(output+'splitWinArgs.csv')

    
