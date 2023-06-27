#!/usr/bin/env python
from amescap.FV3_utils import sol2ls,ls2sol
from amescap.Script_utils import prYellow,prRed
import argparse #parsing arguments
import numpy as np


parser = argparse.ArgumentParser(description='Gives the solar longitude from a SOL or a SOL array (start stop, step), adapted from areols.py',
                                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('sol', nargs='+',type=float,
                             help='''Input is sol number, return solar longitude \n'''
                                  '''Usage: ./MarsCalendar.py 750. \n'''
                                  '''       ./MarsCalendar.py start stop step''')

parser.add_argument('-ls','--ls', action='store_true',
                    help="""Reverse operation. Inpout is Ls, output is sol \n"""
                    """> Usage: ./MarsCalendar.py start stop step' -ls \n"""
                    """ \n""")


parser.add_argument('-my', nargs='+',type=float,default=0.,
                             help=''' Mars Year For ls>sol add 668 if my=1,1336 if my=2 etc.. \n'''
                                  '''Usage: ./MarsCalendar.py  350 -ls -my 2 \n''')

parser.add_argument('-cum', action='store_true',
                             help='''For sol>ls return ls as cummulative 0>360>720... instead of [0-360] \n'''
                                  '''Usage: ./MarsCalendar.py  670 -cum \n''')


if __name__ == '__main__':


    #Load in Mars YEAR (if any, default is zero) and cummulative Ls
    my=np.asarray(parser.parse_args().my).astype(float)
    cum=False
    if parser.parse_args().cum:cum=True

    data_input=np.asarray(parser.parse_args().sol).astype(float)

    if len(data_input)==1:
        in_array=data_input
    elif len(data_input)==3:
        in_array=np.arange(data_input[0],data_input[1],data_input[2]) #start stop step
    else:
        prRed('Wrong number of arguments: enter [sol/ls] or [start stop step]')
        exit()

    #Requesting ls instead
    if parser.parse_args().ls:
        txt_multi='   Ls    |    Sol    '
        result=ls2sol(in_array)
    else:
        txt_multi='    SOL    |    Ls    '
        result=sol2ls(in_array,cummulative=cum)

    #Is scalar, turn as float
    if len(np.atleast_1d(result))==1:result=[result]
    #Display data
    print(txt_multi)
    print('-----------------------')
    for i in range(0,len(in_array)):
        print(' %7.2f   |    %7.3f  '%(in_array[i],result[i]+my*668.))


