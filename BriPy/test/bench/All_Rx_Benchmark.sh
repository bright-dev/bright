#!/bin/sh

./Rx_Benchmark.py -c LWR_NEA
cp LWR_NEA_Compare.txt All_Compare.txt

echo "-----------------" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt

./Rx_Benchmark.py -c LWR_VIS51
cat LWR_VIS51_Compare.txt >> All_Compare.txt

echo "-----------------" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt

./Rx_Benchmark.py -c FR_NEA
cat FR_NEA_Compare.txt >> All_Compare.txt

echo "-----------------" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt

./Rx_Benchmark.py -c FR_VISp1
cat FR_VISp1_Compare.txt >> All_Compare.txt

echo "-----------------" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt
echo "" >> All_Compare.txt

./Rx_Benchmark.py -c FR_VISp5
cat FR_VISp5_Compare.txt >> All_Compare.txt
