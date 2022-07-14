#!/bin/bash

demo/GAN4BER.py "Crohns disease"
demo/GAN4BER.py "CRC"
demo/GAN4BER.py "Obesity"
demo/GAN4BER.py "Overweight"
demo/GAN4BER.py "T2D"
demo/Merge.py "demo/output/Non-healthy-BER.csv"
